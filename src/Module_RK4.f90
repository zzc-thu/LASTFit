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
    real( kind = rk ), dimension(:,:,:),   allocatable :: Rho_per,U_per,V_per,W_per,P_per
    real( kind = rk ), dimension(:,:,:),   allocatable :: i_t,i_t_tau
    
    allocate(Rho_per(1:nx_local,1:Ny,1:nz_local))
    allocate(U_per(1:nx_local,1:Ny,1:nz_local))
    allocate(V_per(1:nx_local,1:Ny,1:nz_local))
    allocate(W_per(1:nx_local,1:Ny,1:nz_local))
    allocate(P_per(1:nx_local,1:Ny,1:nz_local))
    allocate(i_t(1:nx_local,1:Ny,1:nz_local))
    allocate(i_t_tau(1:nx_local,1:Ny,1:nz_local))
    ! Calculate the increament  
     !if(MyId == 0 ) write(*,*)"Reading Iteration Number"
     !open(113,file='RESU/IterationNumber.dat',form='unformatted',status='old')
     !read(113)Iteration
     !close(113)
     !write(*,*)Iteration
     Time_old = (Iteration-1) * DT0;
     
     UconsOld = Ucons0;
     ShockH_old= ShockH;
     ShockV_old= ShockV;
   
   ! ================== 1st Step ==================
    CALL Calculate_Flux

    jc = Ny;
      do kc = 1,nz_local
        do ic = 1,nx_local

    if (Pert_Type .EQ. 1)then !fast 
          i_t(ic,jc,kc) = cos(k_infty*((x_grid(ic,jc,kc)-1.d0) - C_inf * Time_old));
      i_t_tau(ic,jc,kc) = -k_infty*(shockXtau(ic,kc)-C_inf)*sin(k_infty*((x_grid(ic,jc,kc)-1.d0)-C_inf * Time_old));
       
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
           
        Rho_pert_tau(ic,jc,kc) = epsilon * mach_ref * i_t_tau(ic,jc,kc);
          U_pert_tau(ic,jc,kc) = epsilon * i_t_tau(ic,jc,kc);
          V_pert_tau(ic,jc,kc) = 0.d0;
          W_pert_tau(ic,jc,kc) = 0.d0;
          P_pert_tau(ic,jc,kc) = epsilon * mach_ref * gamma * i_t_tau(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
            
           Rho_free(ic,jc,kc)  = Rho_inf + Rho_per(ic,jc,kc);
             U_free(ic,jc,kc)  = U_inf   + U_per(ic,jc,kc);
             V_free(ic,jc,kc)  = V_inf   + V_per(ic,jc,kc);
             W_free(ic,jc,kc)  = W_inf   + W_per(ic,jc,kc);
             P_free(ic,jc,kc)  = P_inf   + P_per(ic,jc,kc);
             T_free(ic,jc,kc)  = Gamma * Mach_Ref * Mach_Ref * P_free(ic,jc,kc) / Rho_free(ic,jc,kc);
          elseif (Pert_Type .EQ. 2)then !slow
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = - epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
          elseif (Pert_Type .EQ. 3)then
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;   
            P_per(ic,jc,kc)    = 0.d0;
          !if (Pert_Type .EQ. 4)then
          else
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = 0.d0;                
          endif    
        enddo
      enddo
      
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
      CALL Calculate_Jaco
      CALL Update_Variables
!write(*,*)Ucons0
   ! ==============================================
    
   ! ================== 2nd Step ==================
    CALL Calculate_Flux

    jc = Ny;
      do kc = 1,nz_local
        do ic = 1,nx_local
     
    if (Pert_Type .EQ. 1)then !fast 
          i_t(ic,jc,kc) = cos(k_infty*((x_grid(ic,jc,kc)-1.d0) - C_inf * Time_RK1));
      i_t_tau(ic,jc,kc) = -k_infty*(shockXtau(ic,kc)-C_inf)*sin(k_infty*((x_grid(ic,jc,kc)-1.d0)-C_inf * Time_RK1));
         
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
           
          Rho_pert_tau(ic,jc,kc) = epsilon * mach_ref * i_t_tau(ic,jc,kc);
            U_pert_tau(ic,jc,kc) = epsilon * i_t_tau(ic,jc,kc);
            V_pert_tau(ic,jc,kc) = 0.d0;
            W_pert_tau(ic,jc,kc) = 0.d0;
            P_pert_tau(ic,jc,kc) = epsilon * mach_ref * gamma * i_t_tau(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
            
          Rho_free(ic,jc,kc)   = Rho_inf + Rho_per(ic,jc,kc);
            U_free(ic,jc,kc)   = U_inf   + U_per(ic,jc,kc);
            V_free(ic,jc,kc)   = V_inf   + V_per(ic,jc,kc);
            W_free(ic,jc,kc)   = W_inf   + W_per(ic,jc,kc);
            P_free(ic,jc,kc)   = P_inf   + P_per(ic,jc,kc);
          
          elseif (Pert_Type .EQ. 2)then !slow
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = - epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
          elseif (Pert_Type .EQ. 3)then
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;   
            P_per(ic,jc,kc)    = 0.d0;
          !if (Pert_Type .EQ. 4)then
          else
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = 0.d0;                
          endif    
        enddo
      enddo
    
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
      CALL Calculate_Jaco
      CALL Update_Variables
   ! ==============================================

   ! ================== 3rd Step ==================
      CALL Calculate_Flux     ! Calculate the dUcons0
!write(*,*)Time_RK2
      jc = Ny;
      do kc = 1,nz_local
        do ic = 1,nx_local
     
    if (Pert_Type .EQ. 1)then !fast 
          i_t(ic,jc,kc) = cos(k_infty*((x_grid(ic,jc,kc)-1.d0) - C_inf * Time_RK2));
      i_t_tau(ic,jc,kc) = -k_infty*(shockXtau(ic,kc)-C_inf)*sin(k_infty*((x_grid(ic,jc,kc)-1.d0)-C_inf * Time_RK2));
         
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
           
          Rho_pert_tau(ic,jc,kc)= epsilon * mach_ref * i_t_tau(ic,jc,kc);
            U_pert_tau(ic,jc,kc) = epsilon * i_t_tau(ic,jc,kc);
            V_pert_tau(ic,jc,kc) = 0.d0;
            W_pert_tau(ic,jc,kc) = 0.d0;
            P_pert_tau(ic,jc,kc) = epsilon * mach_ref * gamma * i_t_tau(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
            
          Rho_free(ic,jc,kc)   = Rho_inf + Rho_per(ic,jc,kc);
            U_free(ic,jc,kc)   = U_inf   + U_per(ic,jc,kc);
            V_free(ic,jc,kc)   = V_inf   + V_per(ic,jc,kc);
            W_free(ic,jc,kc)   = W_inf   + W_per(ic,jc,kc);
            P_free(ic,jc,kc)   = P_inf   + P_per(ic,jc,kc);
          
          elseif (Pert_Type .EQ. 2)then !slow
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = - epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
          elseif (Pert_Type .EQ. 3)then
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;   
            P_per(ic,jc,kc)    = 0.d0;
          !if (Pert_Type .EQ. 4)then
          else
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = 0.d0;                
          endif    
        enddo
      enddo
      
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

      CALL Calculate_Jaco
      CALL Update_Variables
   ! ==============================================

   ! ================== 4th Step ==================
      CALL Calculate_Flux

      jc = Ny;
      do kc = 1,nz_local
        do ic = 1,nx_local
     
    if (Pert_Type .EQ. 1)then !fast 
          i_t(ic,jc,kc) = cos(k_infty*((x_grid(ic,jc,kc)-1.d0) - C_inf * Time_RK3));
      i_t_tau(ic,jc,kc) = -k_infty*(shockXtau(ic,kc)-C_inf)*sin(k_infty*((x_grid(ic,jc,kc)-1.d0)-C_inf * Time_RK3));
         
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
           
          Rho_pert_tau(ic,jc,kc)= epsilon * mach_ref * i_t_tau(ic,jc,kc);
            U_pert_tau(ic,jc,kc) = epsilon * i_t_tau(ic,jc,kc);
            V_pert_tau(ic,jc,kc) = 0.d0;
            W_pert_tau(ic,jc,kc) = 0.d0;
            P_pert_tau(ic,jc,kc) = epsilon * mach_ref * gamma* i_t_tau(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
            
          Rho_free(ic,jc,kc)   = Rho_inf + Rho_per(ic,jc,kc);
            U_free(ic,jc,kc)   = U_inf   + U_per(ic,jc,kc);
            V_free(ic,jc,kc)   = V_inf   + V_per(ic,jc,kc);
            W_free(ic,jc,kc)   = W_inf   + W_per(ic,jc,kc);
            P_free(ic,jc,kc)   = P_inf   + P_per(ic,jc,kc);
          
          elseif (Pert_Type .EQ. 2)then !slow
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = - epsilon * i_t(ic,jc,kc);
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = epsilon * mach_ref * gamma* i_t(ic,jc,kc)/ ( gamma * mach_ref * mach_ref );
          elseif (Pert_Type .EQ. 3)then
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = 0.d0;
            W_per(ic,jc,kc)    = 0.d0;   
            P_per(ic,jc,kc)    = 0.d0;
          !if (Pert_Type .EQ. 4)then
          else
          Rho_per(ic,jc,kc)    = epsilon * mach_ref * i_t(ic,jc,kc);
            U_per(ic,jc,kc)    = 0.d0;
            V_per(ic,jc,kc)    = epsilon * i_t(ic,jc,kc);
            W_per(ic,jc,kc)    = 0.d0;            
            P_per(ic,jc,kc)    = 0.d0;                
          endif    
        enddo
      enddo
      
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
  Deallocate(Rho_per)
  Deallocate(U_per)
  Deallocate(V_per)
  Deallocate(W_per)
  Deallocate(P_per)
  Deallocate(i_t)
  Deallocate(i_t_tau)
  
  END SUBROUTINE IncreamentRK4_unsteady
End Module Runge_Kutta_4