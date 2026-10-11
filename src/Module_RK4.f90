Module Runge_Kutta_4
  ! This is the standard RK4 explicit time marching module
  use SF_Constant,     only: ik,rk,weno5_order,fifth_order,first_order,Unsteady_NS_Analysis
  use SF_CFD_Global,   only: Iteration,DT0,Select_Scheme,ModelType

  ! Arrays and matrix
  use SF_CFD_Global,   only: UconsOld,ShockH_old,ShockV_old,&
                             dUconsRK1,dUconsRK2,dUconsRK3,dUconsRK4,dUcons0,&
                             Ucons_RK1,Ucons_RK2,Ucons_RK3,&
                             ShockAc_RK1,ShockAc_RK2,ShockAc_RK3,ShockAc_RK4,&
                             ShockV_RK1,ShockV_RK2,ShockV_RK3,&
                             Ucons0,ShockH,ShockV,ShockAc,&
                             dShockH,dShockV,dUcons
  use SFitting,        only: Calculate_ShockAC3D,Calculate_ShockAC3D_unsteady
  use OutputParaView,    only:output_results
  implicit none
  !
  contains

  ! This is the RK4 time marching subroutine
  SUBROUTINE IncreamentRK4
    implicit none
    real( kind = rk ):: RK_dt
    ! Calculate the increament  
    
     UconsOld = Ucons0;
    ShockH_old= ShockH;
    ShockV_old= ShockV;
   
   ! ================== 1st Step ==================  

    CALL Calculate_Flux

    CALL Calculate_ShockAC3D

  !call test_complete_flux

           RK_dt = 0.5_rk * DT0;
                                    ! K1
       Ucons_RK1 = UconsOld + RK_dt * dUcons0;
       
       dUconsRK1 = dUcons0;

     ShockAc_RK1 = ShockAc;

      ShockV_RK1 = ShockV_old + RK_dt * ShockAc;

    ! Update Variables, Jaco and Boundary Conditions 
          Ucons0 = Ucons_RK1;
          ShockH = ShockH_old + RK_dt * ShockV;
          ShockV = ShockV_RK1;
          

      if(ModelType == 1) then
          call Calculate_Jaco
      else
          call Singular_Calculate_Jaco
      endif

      CALL Update_Variables

   ! ==============================================
    
   ! ================== 2nd Step ==================
    CALL Calculate_Flux

    CALL Calculate_ShockAC3D

           RK_dt = 0.5_rk * DT0;
                                      ! K2
        Ucons_RK2 = UconsOld + RK_dt * dUcons0;

       dUconsRK2 = dUcons0;

     ShockAc_RK2 = ShockAc;

      ShockV_RK2 = ShockV_old + RK_dt * ShockAc;


    ! Update Variables, Jaco and Boundary Conditions 
          Ucons0 = Ucons_RK2;
          ShockH = ShockH_old + RK_dt * ShockV;
          ShockV = ShockV_RK2;
      
      if(ModelType == 1) then
          call Calculate_Jaco
      else
          call Singular_Calculate_Jaco
      endif
      CALL Update_Variables
      
   ! ==============================================

   ! ================== 3rd Step ==================
      CALL Calculate_Flux     ! Calculate the dUcons0

      CALL Calculate_ShockAC3D

          RK_dt = DT0;
                                      ! K3
      Ucons_RK3 = UconsOld + RK_dt * dUcons0;

      dUconsRK3 = dUcons0;

       ShockAc_RK3 = ShockAc;

        ShockV_RK3 = ShockV_old + RK_dt * ShockAc;

      ! Update Variables, Jaco and Boundary Conditions
         Ucons0 = Ucons_RK3;
         ShockH = ShockH_old + RK_dt * ShockV;
         ShockV = ShockV_RK3;
 
      if(ModelType == 1) then
          call Calculate_Jaco
      else
          call Singular_Calculate_Jaco
      endif
      CALL Update_Variables
      
   ! ==============================================

   ! ================== 4th Step ==================
      CALL Calculate_Flux

      CALL Calculate_ShockAC3D
   
        dUconsRK4 = dUcons0;
      ShockAc_RK4 = ShockAc;

      dUcons = DT0 * (  dUconsRK1+ 2.0_rk *   dUconsRK2+ 2.0_rk *   dUconsRK3+   dUconsRK4)/ 6.0_rk;
     dShockH = DT0 * ( ShockV_old+ 2.0_rk *  ShockV_RK1+ 2.0_rk *  ShockV_RK2+  ShockV_RK3)/ 6.0_rk;
     dShockV = DT0 * (ShockAc_RK1+ 2.0_rk * ShockAc_RK2+ 2.0_rk * ShockAc_RK3+ ShockAc_RK4)/ 6.0_rk;
 
     ! Reset the original variables to keep the same form of all time schemes
      ShockH = ShockH_old + dShockH;
      ShockV = ShockV_old + dShockV;
      Ucons0 = UconsOld   + dUcons;
      
   ! ==============================================
      

  END SUBROUTINE IncreamentRK4


  SUBROUTINE IncreamentRK4_unsteady
    use SF_CFD_Global,   only: Pert_Type,&
                               Time_old,Time_new,Time_RK0,Time_RK1,Time_RK2,Time_RK3,&
                               Rho_inf,U_inf,V_inf,W_inf,P_inf,&
                               Rho_pert_tau,U_pert_tau,V_pert_tau,W_pert_tau,P_pert_tau,&
                               Rho_free,U_free,V_free,W_free,P_free,T_free,&
                               shockXtau,C_inf,k_infty,mach_ref,epsilon,gamma,&
                               x_grid,nz_local,Ny,nx_local

    implicit none
    integer( kind=ik)::ic,jc,kc
    real( kind = rk ):: RK_dt
    
    ! Calculate the increament  
     !if(MyId == 0 ) write(*,*)"Reading Iteration Number"
     !open(113,file='RESU/IterationNumber.dat',form='unformatted',status='old')
     !read(113)Iteration
     !close(113)
     !write(*,*)Iteration
     Time_old = Time_new;
     
     UconsOld = Ucons0;
     ShockH_old= ShockH;
     ShockV_old= ShockV;
   
   ! ================== 1st Step ==================
    CALL Calculate_Flux

    CALL Update_Unsteady_Freestream(Time_old)
      
    CALL Calculate_ShockAC3D_unsteady
    !call test_complete_flux(Iteration-1)
           RK_dt = 0.5_rk * DT0;
                                    ! K1
       Ucons_RK1 = UconsOld + RK_dt * dUcons0;

       dUconsRK1 = dUcons0;

     ShockAc_RK1 = ShockAc;

      ShockV_RK1 = ShockV_old + RK_dt * ShockAc;

      TIME_RK1   = Time_old + RK_dt;
      
    ! Update Variables, Jaco and Boundary Conditions 
          Ucons0 = Ucons_RK1;
          ShockH = ShockH_old + RK_dt * ShockV;
          ShockV = ShockV_RK1;
    !
      if(ModelType == 1) then
        CALL Calculate_Jaco
      else
        CALL Singular_Calculate_Jaco
      endif
      CALL Update_Variables
!write(*,*)Ucons0
   ! ==============================================
    
   ! ================== 2nd Step ==================
    CALL Calculate_Flux

    CALL Update_Unsteady_Freestream(Time_RK1)
    
    CALL Calculate_ShockAC3D_unsteady

           RK_dt = 0.5_rk * DT0;
                                      ! K2
        Ucons_RK2 = UconsOld + RK_dt * dUcons0;

       dUconsRK2 = dUcons0;

     ShockAc_RK2 = ShockAc;

      ShockV_RK2 = ShockV_old + RK_dt * ShockAc;
      
      TIME_RK2   = Time_old + RK_dt;

    ! Update Variables, Jaco and Boundary Conditions 
          Ucons0 = Ucons_RK2;
          ShockH = ShockH_old + RK_dt * ShockV;
          ShockV = ShockV_RK2;
    !
      if(ModelType == 1) then
        CALL Calculate_Jaco
      else
        CALL Singular_Calculate_Jaco
      endif
      CALL Update_Variables
   ! ==============================================

   ! ================== 3rd Step ==================
      CALL Calculate_Flux     ! Calculate the dUcons0
!write(*,*)Time_RK2
      CALL Update_Unsteady_Freestream(Time_RK2)
      
      CALL Calculate_ShockAC3D_unsteady

          RK_dt = DT0;
                                      ! K3
      Ucons_RK3 = UconsOld + RK_dt * dUcons0;

      dUconsRK3 = dUcons0;

       ShockAc_RK3 = ShockAc;

        ShockV_RK3 = ShockV_old + RK_dt * ShockAc;
        
        TIME_RK3   = Time_old + RK_dt;

      ! Update Variables, Jaco and Boundary Conditions
         Ucons0 = Ucons_RK3;
         ShockH = ShockH_old + RK_dt * ShockV;
         ShockV = ShockV_RK3;

      if(ModelType == 1) then
        CALL Calculate_Jaco
      else
        CALL Singular_Calculate_Jaco
      endif
      CALL Update_Variables
   ! ==============================================

   ! ================== 4th Step ==================
      CALL Calculate_Flux

      CALL Update_Unsteady_Freestream(Time_RK3)
      
      CALL Calculate_ShockAC3D_unsteady
   
        dUconsRK4 = dUcons0;
      ShockAc_RK4 = ShockAc;

      dUcons = DT0 * (  dUconsRK1+ 2.0_rk *   dUconsRK2+ 2.0_rk *   dUconsRK3+   dUconsRK4)/ 6.0_rk;
     dShockH = DT0 * ( ShockV_old+ 2.0_rk *  ShockV_RK1+ 2.0_rk *  ShockV_RK2+  ShockV_RK3)/ 6.0_rk;
     dShockV = DT0 * (ShockAc_RK1+ 2.0_rk * ShockAc_RK2+ 2.0_rk * ShockAc_RK3+ ShockAc_RK4)/ 6.0_rk;
 
     ! Reset the original variables to keep the same form of all time schemes
      ShockH = ShockH_old + dShockH;
      ShockV = ShockV_old + dShockV;
      Ucons0 = UconsOld   + dUcons;
    
      
   ! ==============================================
  
  END SUBROUTINE IncreamentRK4_unsteady
  SUBROUTINE Update_Unsteady_Freestream(time)
    use SF_CFD_Global, only: Pert_Type,C_inf,k_infty,epsilon,Mach_Ref,Gamma, &
      Rho_inf,U_inf,V_inf,W_inf,P_inf,x_grid,shockXtau,nx_local,Ny,nz_local, &
      Rho_free,U_free,V_free,W_free,P_free,T_free, &
      Rho_pert_tau,U_pert_tau,V_pert_tau,W_pert_tau,P_pert_tau
    use MPI_GLOBAL, only: MPI_COMM_WORLD,ierr
    use mpi, only: MPI_Abort
    implicit none
    real(kind=rk), intent(in) :: time
    real(kind=rk) :: phase,wave,wave_tau,rho_amp,u_amp,v_amp,p_amp
    integer(kind=ik) :: ic,kc
    rho_amp = epsilon*Mach_Ref
    u_amp = 0.d0
    v_amp = 0.d0
    p_amp = 0.d0
    select case(Pert_Type)
    case(1,2)
      u_amp = epsilon
      if(Pert_Type == 2) u_amp = -epsilon
      p_amp = epsilon/Mach_Ref
    case(3)
      ! Entropy forcing has density but no pressure or velocity fluctuation.
    case(4)
      rho_amp = 0.d0
      v_amp = epsilon
    case default
      write(*,*)'Invalid unsteady perturbation type:',Pert_Type
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      stop
    end select
    do kc=1,nz_local
      do ic=1,nx_local
        phase = k_infty*((x_grid(ic,Ny,kc)-1.d0)-C_inf*time)
        wave = cos(phase)
        wave_tau = -k_infty*(shockXtau(ic,kc)-C_inf)*sin(phase)
        Rho_free(ic,Ny,kc) = Rho_inf+rho_amp*wave
        U_free(ic,Ny,kc) = U_inf+u_amp*wave
        V_free(ic,Ny,kc) = V_inf+v_amp*wave
        W_free(ic,Ny,kc) = W_inf
        P_free(ic,Ny,kc) = P_inf+p_amp*wave
        T_free(ic,Ny,kc) = Gamma*Mach_Ref**2*P_free(ic,Ny,kc)/Rho_free(ic,Ny,kc)
        Rho_pert_tau(ic,Ny,kc) = rho_amp*wave_tau
        U_pert_tau(ic,Ny,kc) = u_amp*wave_tau
        V_pert_tau(ic,Ny,kc) = v_amp*wave_tau
        W_pert_tau(ic,Ny,kc) = 0.d0
        P_pert_tau(ic,Ny,kc) = p_amp*wave_tau
      enddo
    enddo
  END SUBROUTINE Update_Unsteady_Freestream
End Module Runge_Kutta_4
