Module Runge_Kutta_4
  ! This is the standard RK4 explicit time marching module
  use SF_Constant,     only: rk,weno5_order,fifth_order,first_order,Unsteady_NS_Analysis
  use SF_CFD_Global,   only: DT0,Select_Scheme

  ! Arrays and matrix
  use SF_CFD_Global,   only: UconsOld,ShockH_old,ShockV_old,&
                             dUconsRK1,dUconsRK2,dUconsRK3,dUconsRK4,dUcons0,&
                             Ucons_RK1,Ucons_RK2,Ucons_RK3,&
                             ShockAc_RK1,ShockAc_RK2,ShockAc_RK3,ShockAc_RK4,&
                             ShockV_RK1,ShockV_RK2,ShockV_RK3,&
                             Ucons0,ShockH,ShockV,ShockAc,&
                             dShockH,dShockV,dUcons,AnalysisType&
                             TotalTime

  use SFitting,        only: Calculate_ShockAC3D,Calculate_ShockAC3D_unsteady
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

    if(AnalysisType == Unsteady_NS_Analysis) then
    CALL Calculate_ShockAC3D_unsteady
    else
    CALL Calculate_ShockAC3D
    end if
    write(*,*)shockac
    
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
    !
      CALL Calculate_Jaco
      CALL Update_Variables
!write(*,*)Ucons0
   ! ==============================================
    
   ! ================== 2nd Step ==================
    CALL Calculate_Flux
    if(AnalysisType == Unsteady_NS_Analysis) then
    CALL Calculate_ShockAC3D_unsteady
    else
    CALL Calculate_ShockAC3D
    end if

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
    !
      CALL Calculate_Jaco
      CALL Update_Variables
   ! ==============================================

   ! ================== 3rd Step ==================
      CALL Calculate_Flux     ! Calculate the dUcons0
      if(AnalysisType == Unsteady_NS_Analysis) then
      CALL Calculate_ShockAC3D_unsteady
      else
      CALL Calculate_ShockAC3D
      end if

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

      CALL Calculate_Jaco
      CALL Update_Variables
   ! ==============================================

   ! ================== 4th Step ==================
      CALL Calculate_Flux
      if(AnalysisType == Unsteady_NS_Analysis) then
      CALL Calculate_ShockAC3D_unsteady
      else
      CALL Calculate_ShockAC3D
      end if
   
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


End Module Runge_Kutta_4