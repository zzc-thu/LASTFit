! This module is used to give the variables behind the shock surface
! by using the Rankine-Hugoniot relations, and the shock acceleration
! by using the Xiaolin Zhong's method (JCP 1999).
Module SFitting
    implicit none
    ! Local variables

    contains

SUBROUTINE ShockRelation3D
  ! This subroutine calculate the primitive variables after the shock wave	
  use SF_Constant,  only:ik,rk,OverLap
  use SF_CFD_Global,only:nx_local,Ny,nz_local,&
                         Gamma,Mach_Ref,&
                         U_inf,V_inf,W_inf,P_inf,Rho_inf
  ! arrays
  use SF_CFD_Global,only:ShockNormalX,ShockNormalY,ShockNormalZ,ShockV,&
                         WallNormalX,WallNormalY,WallNormalZ,&
                         U,V,W,P,Rho,T
  ! subs and functions

  implicit none
  real( kind = rk )::Delta_Inf,Masn2Inf,DeltaShkUinf
  
  integer( kind=ik)::ic,jc,kc
  
  jc = Ny ! Standing for the Shock Surface
  do kc = 1,nz_local
   do ic = 1,nx_local
    ! Velocity ahead of the shocks
    Delta_Inf = U_inf*ShockNormalX(ic,kc) + V_inf*ShockNormalY(ic,kc) + W_inf*ShockNormalZ(ic,kc) &
             &- ShockV(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX(ic,kc)&
                             &+WallNormalY(ic,kc)*ShockNormalY(ic,kc) &
                             &+WallNormalZ(ic,kc)*ShockNormalZ(ic,kc));
    ! Shock normal Mach number^2
      Masn2Inf = (Mach_Ref * Delta_Inf)**2.d0;
    ! Velocity normal to the shock after the shock
      DeltaShkUinf = 2.d0 * Delta_Inf / (Gamma + 1.d0) * (1.d0 / Masn2Inf - 1.d0);
      
      U(ic,jc,kc) = U_inf + DeltaShkUinf * ShockNormalX(ic,kc);

      V(ic,jc,kc) = V_inf + DeltaShkUinf * ShockNormalY(ic,kc);

      W(ic,jc,kc) = W_inf + DeltaShkUinf * ShockNormalZ(ic,kc);

      P(ic,jc,kc) = P_inf * (1.d0 + 2.d0 * Gamma / (Gamma + 1.d0) * (Masn2Inf - 1.d0));

    Rho(ic,jc,kc) = Rho_inf * (Gamma + 1.d0) * Masn2Inf / ((Gamma - 1.d0) * Masn2Inf + 2.d0);

      T(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / Rho(ic,jc,kc);
   enddo
  enddo
  

END SUBROUTINE ShockRelation3D

SUBROUTINE Calculate_ShockAC3D
  use SF_Constant,    only:ik,rk,NumVar
  use SF_CFD_Global,  only:nx_local,Ny,nz_local,Gamma

  ! arrays
  use SF_CFD_Global,  only:CV_INF,FINV_INF,GINV_INF,HINV_INF,invF,invG,invH,&
                           Rho,U,V,W,P,ShockNormalX,ShockNormalY,ShockNormalZ,ShockVN,&
                           ax_tau,ay_tau,az_tau,dUcons0,Ucons0,ShockAc,ShockXtau,ShockYtau,ShockZtau,&
                           WallNormalX,WallNormalY,WallNormalZ

  implicit none
  real( kind = rk )::locRho,locU,locV,locW,locP,locC,locUns
  real( kind = rk )::LAM5,L5(5_ik)
  real( kind = rk )::SUM_L5_U,SUM_L5_FG,FG
  integer( kind = ik )::ic,jc,kc,iVar
   
  jc = Ny;
  do kc = 1,nz_local
   do ic = 1,nx_local
    locRho = Rho(ic,jc,kc);
    locU   =   U(ic,jc,kc);
    locV   =   V(ic,jc,kc);
    locW   =   W(ic,jc,kc);
    locP   =   P(ic,jc,kc);
    locC   = sqrt(Gamma * locP / locRho);
    ! Calculate the local normal
    locUns = locU * ShockNormalX(ic,kc) + locV * ShockNormalY(ic,kc) + locW * ShockNormalZ(ic,kc);
    LAM5   = locUns - ShockVN(ic,kc) + locC;

    L5(1) = 0.5d0 * (Gamma - 1.d0)*( locU * locU + locV * locV + locW * locW ) - locUns * locC; 
    L5(2) =        -(Gamma - 1.d0)* locU + ShockNormalX(ic,kc) * locC;
    L5(3) =        -(Gamma - 1.d0)* locV + ShockNormalY(ic,kc) * locC;
    L5(4) =        -(Gamma - 1.d0)* locW + ShockNormalZ(ic,kc) * locC;
    L5(5) =         (Gamma - 1.d0);

    SUM_L5_U  = 0.d0;
    SUM_L5_FG = 0.d0;

    DO iVar = 1,NumVar
      FG = LAM5 * dUcons0(iVar,ic,jc,kc) + ax_tau(ic,kc) * (invF(iVar,ic,jc,kc) - FINV_INF(iVar)) &
                                         + ay_tau(ic,kc) * (invG(iVar,ic,jc,kc) - GINV_INF(iVar)) &
                                         + az_tau(ic,kc) * (invH(iVar,ic,jc,kc) - HINV_INF(iVar));

      SUM_L5_U  = SUM_L5_U  + L5(iVar) * (Ucons0(iVar,ic,jc,kc) - CV_INF(iVar));
      SUM_L5_FG = SUM_L5_FG + L5(iVar) * FG;                       
    ENDDO
    
    ! Calculate the ShockAc
    ShockAc(ic,kc) = SUM_L5_FG/SUM_L5_U - (ax_tau(ic,kc) * ShockXtau(ic,kc) &
                                        &+ ay_tau(ic,kc) * ShockYtau(ic,kc) &
                                        &+ az_tau(ic,kc) * ShockZtau(ic,kc));

    ShockAc(ic,kc) = ShockAc(ic,kc)/( WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                                   &+ WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                                   &+ WallNormalZ(ic,kc) * ShockNormalZ(ic,kc));
   enddo
  enddo

END SUBROUTINE Calculate_ShockAC3D

SUBROUTINE Calculate_ShockAC3D_Implicit
  use SF_Constant,    only:ik,rk,NumVar,overLAP
  use SF_CFD_Global,  only:nx_local,Ny,nz_local,Gamma

  ! arrays
  use SF_CFD_Global,  only:CV_INF,FINV_INF,GINV_INF,HINV_INF,invF,invG,invH,&
                           Rho,U,V,W,P,ShockNormalX,ShockNormalY,ShockNormalZ,ShockVN,&
                           ax_tau,ay_tau,az_tau,dUcons0,Ucons0,ShockAc,ShockXtau,ShockYtau,ShockZtau,&
                           WallNormalX,WallNormalY,WallNormalZ,CO_AC_UT,AsUs_P1,AsUs_P3,Ucons0

  ! subs and functions
  use MPI_GLOBAL,     only: Parallel_Exchange_Surface
  
  implicit none
  real( kind = rk )::locRho,locU,locV,locW,locP,locC,locUns
  real( kind = rk )::LAM5,L5(5_ik)
  real( kind = rk )::SUM_L5_U,SUM_L5_FG,FG,CO_TEMP
  real( kind = rk )::SUM_L5_A(NumVar),SUM_L5_B(NumVar),SUM_L5_C(NumVar)
  real( kind = rk )::Amat(NumVar,NumVar),Bmat(NumVar,NumVar),Cmat(NumVar,NumVar)
  real( kind = rk )::Kc_Gamma,localQ2,localHan
  integer( kind = ik )::ic,jc,kc,iVar,jVar
   
  Kc_Gamma = Gamma - 1.0_rk;

  jc = Ny;
  do kc = 1,nz_local
    do ic = 1,nx_local
      locRho = Rho(ic,jc,kc);
      locU   =   U(ic,jc,kc);
      locV   =   V(ic,jc,kc);
      locW   =   W(ic,jc,kc);
      locP   =   P(ic,jc,kc);
      locC   = sqrt(Gamma * locP / locRho);
      ! local Q2: velocity magnitude square
      localQ2 = locU * locU + locV * locV + locW * locW;
      ! local enthalpy
      localHan = (Ucons0(5,ic,jc,kc) + P(ic,jc,kc)) / Rho(ic,jc,kc);
  
      ! Calculate the local normal
      locUns = locU * ShockNormalX(ic,kc) + locV * ShockNormalY(ic,kc) + locW * ShockNormalZ(ic,kc);
      LAM5   = locUns - ShockVN(ic,kc) + locC;
  
      L5(1) = 0.5d0 * (Gamma - 1.d0)* localQ2 - locUns * locC; 
      L5(2) =       - (Gamma - 1.d0)* locU + ShockNormalX(ic,kc) * locC;
      L5(3) =       - (Gamma - 1.d0)* locV + ShockNormalY(ic,kc) * locC;
      L5(4) =       - (Gamma - 1.d0)* locW + ShockNormalZ(ic,kc) * locC;
      L5(5) =         (Gamma - 1.d0);
  
      SUM_L5_U  = 0.d0;
      SUM_L5_FG = 0.d0;

    DO iVar = 1,NumVar
      FG = LAM5 * dUcons0(iVar,ic,jc,kc) + ax_tau(ic,kc) * (invF(iVar,ic,jc,kc) - FINV_INF(iVar)) &
                                       & + ay_tau(ic,kc) * (invG(iVar,ic,jc,kc) - GINV_INF(iVar)) &
                                       & + az_tau(ic,kc) * (invH(iVar,ic,jc,kc) - HINV_INF(iVar));
      ! L5^T * (U0-UINF)
      SUM_L5_U  = SUM_L5_U  + L5(iVar) * (Ucons0(iVar,ic,jc,kc) - CV_INF(iVar));
      SUM_L5_FG = SUM_L5_FG + L5(iVar) * FG;                       
    ENDDO

    ! AsUs can be divided into 3 parts
    DO iVar = 1, NumVar
      AsUs_P1(iVar,ic,kc) = - SUM_L5_FG * L5(iVar)/( (WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
        &+ WallNormalY(ic,kc) * ShockNormalY(ic,kc) + WallNormalZ(ic,kc) * ShockNormalZ(ic,kc)) &
        &* SUM_L5_U * SUM_L5_U);
    ENDDO

    ! Part 2 Coefficients
    CO_TEMP = LAM5 /( (WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                   & + WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                   & + WallNormalZ(ic,kc) * ShockNormalZ(ic,kc)) * SUM_L5_U );

    DO iVar = 1, NumVar
      CO_AC_UT(iVar,ic,kc) = CO_TEMP * L5(iVar);
    ENDDO
    
    ! Part 3
        Amat = 0.d0;
        Bmat = 0.d0;
        Cmat = 0.d0;
    
    ! Linearized Jacobian matrix of the inviscid flux along x-direction 
        Amat(1, 2)   =  1.d0;
        
        Amat(2, 1)   = 0.5d0 * Kc_Gamma * localQ2 - locU * locU;
        Amat(2, 2)   = (3.d0 - Gamma) * locU;
        Amat(2, 3)   = (1.d0 - Gamma) * locV;
        Amat(2, 4)   = (1.d0 - Gamma) * locW;
        Amat(2, 5)   = Gamma - 1.d0;

        Amat(3, 1)   = -locU * locV;
        Amat(3, 2)   = locV;
        Amat(3, 3)   = locU;

        Amat(4, 1)   = -locU * locW;
        Amat(4, 2)   = locW;
        Amat(4, 4)   = locU;

        Amat(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * locU;
        Amat(5, 2)   = localHan - Kc_Gamma * locU * locU;
        Amat(5, 3)   =          - Kc_Gamma * locU * locV;
        Amat(5, 4)   =          - Kc_Gamma * locU * locW;
        Amat(5, 5)   = Kc_Gamma * locU;

        ! Linearized Jacobian matrix of the inviscid flux along y-direction
        Bmat(1, 3)   = 1.d0;

        Bmat(2, 1)   = -locU * locV;
        Bmat(2, 2)   =  locV;
        Bmat(2, 3)   =  locU;

        Bmat(3, 1)   = 0.5d0 * Kc_Gamma * localQ2 - locV * locV;
        Bmat(3, 2)   = (1.d0 - Gamma) * locU;
        Bmat(3, 3)   = (3.d0 - Gamma) * locV;
        Bmat(3, 4)   = (1.d0 - Gamma) * locW;
        Bmat(3, 5)   = Gamma - 1.d0;

        Bmat(4, 1)   = -locV * locW;
        Bmat(4, 3)   =  locW;
        Bmat(4, 4)   =  locV;

        Bmat(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * locV;
        Bmat(5, 2)   =          - Kc_Gamma * locV * locU;
        Bmat(5, 3)   = localHan - Kc_Gamma * locV * locV;
        Bmat(5, 4)   =          - Kc_Gamma * locV * locW;
        Bmat(5, 5)   = Kc_Gamma * locV;

        ! Linearized Jacobian matrix of the inviscid flux along z-direction
        Cmat(1, 4)   = 1.d0;

        Cmat(2, 1)   = -locU * locW;
        Cmat(2, 2)   =  locW;
        Cmat(2, 4)   =  locU;

        Cmat(3, 1)   = -locV * locW;
        Cmat(3, 3)   =  locW;
        Cmat(3, 4)   =  locV;

        Cmat(4, 1)   = 0.5d0 * Kc_Gamma * localQ2 - locW * locW;
        Cmat(4, 2)   = (1.d0 - Gamma) * locU;
        Cmat(4, 3)   = (1.d0 - Gamma) * locV;
        Cmat(4, 4)   = (3.d0 - Gamma) * locW;
        Cmat(4, 5)   = Gamma - 1.d0;

        Cmat(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * locW;
        Cmat(5, 2)   =          - Kc_Gamma * locW * locU;
        Cmat(5, 3)   =          - Kc_Gamma * locW * locV;
        Cmat(5, 4)   = localHan - Kc_Gamma * locW * locW;
        Cmat(5, 5)   = Kc_Gamma * locW;

        SUM_L5_A = 0.d0;
        SUM_L5_B = 0.d0; 
        SUM_L5_C = 0.d0;

    DO jVar = 1, NumVar
      DO iVar = 1, NumVar
        SUM_L5_A(jVar) = SUM_L5_A(jVar) + L5(iVar) * ax_tau(ic,kc) * Amat(iVar,jVar);
        SUM_L5_B(jVar) = SUM_L5_B(jVar) + L5(iVar) * ay_tau(ic,kc) * Bmat(iVar,jVar);
        SUM_L5_C(jVar) = SUM_L5_C(jVar) + L5(iVar) * az_tau(ic,kc) * Cmat(iVar,jVar);
      ENDDO
        AsUs_P3(jVar,ic,kc) = (SUM_L5_A(jVar) + SUM_L5_B(jVar) + SUM_L5_C(jVar))/&
                             &((WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                             &+ WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                             &+ WallNormalZ(ic,kc) * ShockNormalZ(ic,kc)) * SUM_L5_U);
    ENDDO
    
    ! Calculate the ShockAc
    ShockAc(ic,kc) = SUM_L5_FG/SUM_L5_U - (ax_tau(ic,kc) * ShockXtau(ic,kc) &
                                        &+ ay_tau(ic,kc) * ShockYtau(ic,kc) &
                                        &+ az_tau(ic,kc) * ShockZtau(ic,kc));

    ShockAc(ic,kc) = ShockAc(ic,kc)/( WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                                   &+ WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                                   &+ WallNormalZ(ic,kc) * ShockNormalZ(ic,kc));
   enddo
  enddo

  ! Exchange the shock related variables
  call Parallel_Exchange_Surface(ShockAc) 

  call Parallel_Exchange_Surface(AsUs_P1(1,:,:))
  call Parallel_Exchange_Surface(AsUs_P1(2,:,:))
  call Parallel_Exchange_Surface(AsUs_P1(3,:,:))
  call Parallel_Exchange_Surface(AsUs_P1(4,:,:))
  call Parallel_Exchange_Surface(AsUs_P1(5,:,:))

  call Parallel_Exchange_Surface(CO_AC_UT(1,:,:))
  call Parallel_Exchange_Surface(CO_AC_UT(2,:,:))
  call Parallel_Exchange_Surface(CO_AC_UT(3,:,:))
  call Parallel_Exchange_Surface(CO_AC_UT(4,:,:))
  call Parallel_Exchange_Surface(CO_AC_UT(5,:,:))

  call Parallel_Exchange_Surface(AsUs_P3(1,:,:))
  call Parallel_Exchange_Surface(AsUs_P3(2,:,:))
  call Parallel_Exchange_Surface(AsUs_P3(3,:,:))
  call Parallel_Exchange_Surface(AsUs_P3(4,:,:))
  call Parallel_Exchange_Surface(AsUs_P3(5,:,:))

END SUBROUTINE Calculate_ShockAC3D_Implicit


SUBROUTINE ShockRelation3D_per(Rho_free,U_free,V_free,W_free,P_free,sRho_pert,sU_pert,sV_pert,sW_pert,sP_pert,sT_pert)
  ! This subroutine calculate the primitive variables after the shock wave	
  use SF_Constant,  only:ik,rk,OverLap
  use SF_CFD_Global,only:nx_local,Ny,nz_local,&
                         Gamma,Mach_Ref,&
                         U_inf,V_inf,W_inf,P_inf,Rho_inf
  ! arrays
  use SF_CFD_Global,only:ShockNormalX,ShockNormalY,ShockNormalZ,ShockV,&
                         WallNormalX,WallNormalY,WallNormalZ,&
                         shockNormalX_steady,shockNormalY_steady,shockNormalZ_steady,shockV_steady
  ! subs and functions
  
  implicit none
  real( kind = rk ), dimension(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),intent(in)  :: Rho_free,U_free,V_free,W_free,P_free
  real( kind = rk ), dimension(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),intent(out) :: sRho_pert,sU_pert,sV_pert,sW_pert,sP_pert,sT_pert

  ! Local variables
  real( kind = rk ), dimension(:,:,:),allocatable:: sRho_un,sU_un,sV_un,sW_un,sP_un,sT_un
  real( kind = rk ), dimension(:,:,:),allocatable:: sRho_st,sU_st,sV_st,sW_st,sP_st,sT_st
  real( kind = rk )::Delta_Inf,Masn2Inf,DeltaShkUinf,Delta2_free,Delta2
  real( kind = rk )::Delta_Inf_free,Masn2Inf_free,DeltaShkUinf_free
  integer( kind = ik )::ic,jc,kc
  
  allocate(sRho_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sU_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sV_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),&
           sW_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sP_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sT_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
  allocate(sRho_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sU_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sV_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),&
           sW_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sP_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sT_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
  
  
  !write(*,*),'shockNormalX(1,1)',shockNormalX(41,1)
  !write(*,*),'shockNormalX_steady(1,1)',shockNormalX_steady(41,1)
  
  jc = Ny ! Standing for the Shock Surface
  do kc = 1-overLap,nz_local+overLap
   do ic = 1-overLap,nx_local+overLap
    ! Velocity ahead of the shocks
    Delta_Inf_free = U_free(ic,jc,kc)*ShockNormalX(ic,kc) + V_free(ic,jc,kc)*ShockNormalY(ic,kc) + W_free(ic,jc,kc)*ShockNormalZ(ic,kc) &
             &- ShockV(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX(ic,kc) &
                             &+WallNormalY(ic,kc)*ShockNormalY(ic,kc) &
                             &+WallNormalZ(ic,kc)*ShockNormalZ(ic,kc));

    Delta2_free = (Gamma - 1.d0) * Delta_inf_free/(Gamma + 1.d0) + 2.d0 * Gamma * P_free(ic,jc,kc)/((Gamma + 1.d0)*Rho_free(ic,jc,kc) * delta_inf_free);
   
      sU_un(ic,jc,kc) = U_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalX(ic,kc);

      sV_un(ic,jc,kc) = V_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalY(ic,kc);

      sW_un(ic,jc,kc) = W_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalZ(ic,kc);

      sP_un(ic,jc,kc) = 2.0 * Rho_free(ic,jc,kc) * delta_inf_free * delta_inf_free/(gamma + 1.0) - (gamma - 1.0) * P_free(ic,jc,kc) / (gamma + 1.0);

    sRho_un(ic,jc,kc) = delta_inf_free * Rho_free(ic,jc,kc) / delta2_free;

      sT_un(ic,jc,kc) = sP_un(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / sRho_un(ic,jc,kc);
   enddo
  enddo
 
  
  jc = Ny ! Standing for the Shock Surface
  do kc = 1-overLap,nz_local+overLap
   do ic = 1-overLap,nx_local+overLap
    ! Velocity ahead of the shocks
    Delta_Inf = U_inf * ShockNormalX(ic,kc) + V_inf * ShockNormalY(ic,kc) + W_inf * ShockNormalZ(ic,kc) &
              &- ShockV_steady(ic,kc) *(WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                                    & + WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                                    & + WallNormalZ(ic,kc) * ShockNormalZ(ic,kc));
    
    Delta2 = (Gamma - 1.d0) * Delta_inf/(Gamma + 1.d0) + 2.d0 * Gamma * P_inf/((Gamma + 1.d0)*Rho_inf * delta_inf);
    
      sU_st(ic,jc,kc) = U_inf + (Delta2 - Delta_Inf) * ShockNormalX(ic,kc);

      sV_st(ic,jc,kc) = V_inf + (Delta2 - Delta_Inf) * ShockNormalY(ic,kc);

      sW_st(ic,jc,kc) = W_inf + (Delta2 - Delta_Inf) * ShockNormalZ(ic,kc);

      sP_st(ic,jc,kc) = 2.0 * Rho_inf * delta_inf * delta_inf/(gamma + 1.0) - (gamma - 1.0) * P_inf / (gamma + 1.0);

    sRho_st(ic,jc,kc) = delta_inf * Rho_inf / delta2;

      sT_st(ic,jc,kc) = sP_st(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / sRho_st(ic,jc,kc);
   enddo
  enddo
  ! Calculate the basic varaibles of the perturbed state
  jc = Ny ! Standing for the Shock Surface
    do kc = 1-overLap,nz_local+overLap
      do ic = 1-overLap,nx_local+overLap
        sRho_pert(ic,jc,kc)   = sRho_un(ic,jc,kc) - sRho_st(ic,jc,kc)
          sU_pert(ic,jc,kc)   = sU_un(ic,jc,kc)   - sU_st(ic,jc,kc)
          sV_pert(ic,jc,kc)   = sV_un(ic,jc,kc)   - sV_st(ic,jc,kc)
          sW_pert(ic,jc,kc)   = sW_un(ic,jc,kc)   - sW_st(ic,jc,kc)
          sP_pert(ic,jc,kc)   = sP_un(ic,jc,kc)   - sP_st(ic,jc,kc)
          sT_pert(ic,jc,kc)   = sT_un(ic,jc,kc)   - sT_st(ic,jc,kc)
      enddo
    enddo

    deallocate(sRho_un,sU_un,sV_un,sW_un,sP_un,sT_un)
    deallocate(sRho_st,sU_st,sV_st,sW_st,sP_st,sT_st)
END SUBROUTINE ShockRelation3D_per

SUBROUTINE ShockRelation3D_per0(Rho_free,U_free,V_free,W_free,P_free,sRho_pert,sU_pert,sV_pert,sW_pert,sP_pert,sT_pert)
  ! This subroutine calculate the primitive variables after the shock wave	
  use SF_Constant,  only:ik,rk,OverLap
  use SF_CFD_Global,only:nx_local,Ny,nz_local,&
                         Gamma,Mach_Ref,&
                         U_inf,V_inf,W_inf,P_inf,Rho_inf
  ! arrays
  use SF_CFD_Global,only:ShockNormalX,ShockNormalY,ShockNormalZ,ShockV,&
                         WallNormalX,WallNormalY,WallNormalZ,&
                         shockNormalX_steady,shockNormalY_steady,shockNormalZ_steady,shockV_steady
  ! subs and functions
  
  implicit none
  real( kind = rk ), dimension(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),intent(in)  :: Rho_free,U_free,V_free,W_free,P_free
  real( kind = rk ), dimension(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),intent(out) :: sRho_pert,sU_pert,sV_pert,sW_pert,sP_pert,sT_pert

  ! Local variables
  real( kind = rk ), dimension(:,:,:),allocatable:: sRho_un,sU_un,sV_un,sW_un,sP_un,sT_un
  real( kind = rk ), dimension(:,:,:),allocatable:: sRho_st,sU_st,sV_st,sW_st,sP_st,sT_st
  real( kind = rk )::Delta_Inf,Masn2Inf,DeltaShkUinf,Delta2_free,Delta2
  real( kind = rk )::Delta_Inf_free,Masn2Inf_free,DeltaShkUinf_free
  integer( kind = ik )::ic,jc,kc
  
  allocate(sRho_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sU_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sV_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),&
           sW_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sP_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sT_un(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
  allocate(sRho_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sU_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sV_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),&
           sW_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sP_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap),sT_st(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
  jc = Ny ! Standing for the Shock Surface
  do kc = 1-overLap,nz_local+overLap
   do ic = 1-overLap,nx_local+overLap
    ! Velocity ahead of the shocks
    Delta_Inf_free = U_free(ic,jc,kc)*ShockNormalX_steady(ic,kc) + V_free(ic,jc,kc)*ShockNormalY_steady(ic,kc) + W_free(ic,jc,kc)*ShockNormalZ_steady(ic,kc) &
             &- ShockV_steady(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX_steady(ic,kc) &
                                    &+WallNormalY(ic,kc)*ShockNormalY_steady(ic,kc) &
                                    &+WallNormalZ(ic,kc)*ShockNormalZ_steady(ic,kc));

    Delta2_free = (Gamma - 1.d0) * Delta_inf_free/(Gamma + 1.d0) + 2.d0 * Gamma * P_free(ic,jc,kc)/((Gamma + 1.d0)*Rho_free(ic,jc,kc) * delta_inf_free);
   
      sU_un(ic,jc,kc) = U_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalX_steady(ic,kc);

      sV_un(ic,jc,kc) = V_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalY_steady(ic,kc);

      sW_un(ic,jc,kc) = W_free(ic,jc,kc) + (Delta2_free - Delta_Inf_free) * ShockNormalZ_steady(ic,kc);

      sP_un(ic,jc,kc) = 2.0 * Rho_free(ic,jc,kc) * delta_inf_free * delta_inf_free/(gamma + 1.0) - (gamma - 1.0) * P_free(ic,jc,kc) / (gamma + 1.0);

    sRho_un(ic,jc,kc) = delta_inf_free * Rho_free(ic,jc,kc) / delta2_free;

      sT_un(ic,jc,kc) = sP_un(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / sRho_un(ic,jc,kc);
   enddo
  enddo

  jc = Ny ! Standing for the Shock Surface
  do kc = 1-overLap,nz_local+overLap
   do ic = 1-overLap,nx_local+overLap
    ! Velocity ahead of the shocks
    Delta_Inf = U_inf * ShockNormalX_steady(ic,kc) + V_inf * ShockNormalY_steady(ic,kc) + W_inf * ShockNormalZ_steady(ic,kc) &
              &- ShockV_steady(ic,kc) *(WallNormalX(ic,kc) * ShockNormalX_steady(ic,kc) &
                                    & + WallNormalY(ic,kc) * ShockNormalY_steady(ic,kc) &
                                    & + WallNormalZ(ic,kc) * ShockNormalZ_steady(ic,kc));
    
    Delta2 = (Gamma - 1.d0) * Delta_inf/(Gamma + 1.d0) + 2.d0 * Gamma * P_inf/((Gamma + 1.d0)*Rho_inf * delta_inf);
    
      sU_st(ic,jc,kc) = U_inf + (Delta2 - Delta_Inf) * ShockNormalX_steady(ic,kc);

      sV_st(ic,jc,kc) = V_inf + (Delta2 - Delta_Inf) * ShockNormalY_steady(ic,kc);

      sW_st(ic,jc,kc) = W_inf + (Delta2 - Delta_Inf) * ShockNormalZ_steady(ic,kc);

      sP_st(ic,jc,kc) = 2.0 * Rho_inf * delta_inf * delta_inf/(gamma + 1.0) - (gamma - 1.0) * P_inf / (gamma + 1.0);

    sRho_st(ic,jc,kc) = delta_inf * Rho_inf / delta2;

      sT_st(ic,jc,kc) = sP_st(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / sRho_st(ic,jc,kc);
   enddo
  enddo
  ! Calculate the basic varaibles of the perturbed state
  jc = Ny ! Standing for the Shock Surface
    do kc = 1-overLap,nz_local+overLap
      do ic = 1-overLap,nx_local+overLap
        sRho_pert(ic,jc,kc)   = sRho_un(ic,jc,kc) - sRho_st(ic,jc,kc)
          sU_pert(ic,jc,kc)   = sU_un(ic,jc,kc)   - sU_st(ic,jc,kc)
          sV_pert(ic,jc,kc)   = sV_un(ic,jc,kc)   - sV_st(ic,jc,kc)
          sW_pert(ic,jc,kc)   = sW_un(ic,jc,kc)   - sW_st(ic,jc,kc)
          sP_pert(ic,jc,kc)   = sP_un(ic,jc,kc)   - sP_st(ic,jc,kc)
          sT_pert(ic,jc,kc)   = sT_un(ic,jc,kc)   - sT_st(ic,jc,kc)
      enddo
    enddo

    deallocate(sRho_un,sU_un,sV_un,sW_un,sP_un,sT_un)
    deallocate(sRho_st,sU_st,sV_st,sW_st,sP_st,sT_st)
END SUBROUTINE ShockRelation3D_per0

SUBROUTINE Calculate_ShockAC3D_unsteady
  use SF_Constant,    only:ik,rk,NumVar
  use SF_CFD_Global,  only:nx_local,Ny,nz_local,Gamma

  ! arrays
  use SF_CFD_Global,  only:CV_INF,FINV_INF,GINV_INF,HINV_INF,invF,invG,invH,&
                           Rho,U,V,W,P,ShockNormalX,ShockNormalY,ShockNormalZ,ShockVN,&
                           ax_tau,ay_tau,az_tau,dUcons0,Ucons0,ShockAc,ShockXtau,ShockYtau,ShockZtau,&
                           WallNormalX,WallNormalY,WallNormalZ,x_grid,&
                           Rho_inf,U_inf,V_inf,W_inf,P_inf,C_inf,CV_inf,&
                           Pert_Type,k_infty,epsilon,dt0,Gamma,Mach_Ref,&
                           Rho_free,U_free,V_free,W_free,P_free,&
                           Rho_pert_tau,U_pert_tau,V_pert_tau,W_pert_tau,P_pert_tau
  !subs and functions
  implicit none
  real( kind = rk )::locRho,locU,locV,locW,locP,locC,locUns
  real( kind = rk )::LAM5,L5(5_ik)
  real( kind = rk )::SUM_L5_U,SUM_L5_FG,FG
  real( kind = rk )::FG_f,SUM_L5_FG_f
  integer( kind = ik )::ic,jc,kc,iVar

  real( kind = rk ), dimension(:,:,:,:), allocatable :: CV_free,FCV_free,GCV_free,HCV_free
  real( kind = rk ), dimension(:,:,:,:), allocatable :: CV_free_tau,FCV_free_tau,GCV_free_tau,HCV_free_tau
  real( kind = rk ), dimension(:,:,:,:), allocatable :: FGCV_free_tau
  
    allocate(CV_free(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(FCV_free(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(GCV_free(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(HCV_free(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(CV_free_tau(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(FCV_free_tau(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(GCV_free_tau(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(HCV_free_tau(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    allocate(FGCV_free_tau(1:NumVar,1:nx_local,1:Ny,1:nz_local))
    
    ! Calculate the inf variables
  jc = Ny;
  do kc = 1,nz_local
   do ic = 1,nx_local
    locRho = Rho(ic,jc,kc);
    locU   =   U(ic,jc,kc);
    locV   =   V(ic,jc,kc);
    locW   =   W(ic,jc,kc);
    locP   =   P(ic,jc,kc);
    locC   =   sqrt(Gamma * locP / locRho);
    ! Calculate the local normal
    locUns = locU * ShockNormalX(ic,kc) + locV * ShockNormalY(ic,kc) + locW * ShockNormalZ(ic,kc);
    LAM5   = locUns - ShockVN(ic,kc) + locC;

    L5(1) = 0.5d0 * ( Gamma - 1.d0 )*( locU * locU + locV * locV + locW * locW ) - locUns * locC; 
    L5(2) =        -( Gamma - 1.d0 )*  locU + ShockNormalX(ic,kc) * locC;
    L5(3) =        -( Gamma - 1.d0 )*  locV + ShockNormalY(ic,kc) * locC;
    L5(4) =        -( Gamma - 1.d0 )*  locW + ShockNormalZ(ic,kc) * locC;
    L5(5) =         ( Gamma - 1.d0 );

    SUM_L5_U  = 0.d0;
    SUM_L5_FG = 0.d0;
    SUM_L5_FG_f = 0.d0;
    
   
       call UtoCV_CPG(Rho_free(ic,jc,kc),U_free(ic,jc,kc),V_free(ic,jc,kc),W_free(ic,jc,kc),P_free(ic,jc,kc),CV_free(:,ic,jc,kc))
       
       call UtoFlux_CPG3D(Rho_free(ic,jc,kc),U_free(ic,jc,kc),V_free(ic,jc,kc),W_free(ic,jc,kc),P_free(ic,jc,kc),FCV_FREE(:,ic,jc,kc),&
           &GCV_FREE(:,ic,jc,kc),HCV_FREE(:,ic,jc,kc))
       
       call UtoCV_tau(Rho_free(ic,jc,kc),U_free(ic,jc,kc),V_free(ic,jc,kc),W_free(ic,jc,kc),&
           &Rho_pert_tau(ic,jc,kc),U_pert_tau(ic,jc,kc),V_pert_tau(ic,jc,kc),W_pert_tau(ic,jc,kc),P_pert_tau(ic,jc,kc),CV_free_tau(:,ic,jc,kc))
       
       call UtoFlux_tau(Rho_free(ic,jc,kc),U_free(ic,jc,kc),V_free(ic,jc,kc),W_free(ic,jc,kc),P_free(ic,jc,kc),&
           &Rho_pert_tau(ic,jc,kc),U_pert_tau(ic,jc,kc),V_pert_tau(ic,jc,kc),W_pert_tau(ic,jc,kc),P_pert_tau(ic,jc,kc),&
           &FCV_FREE_tau(:,ic,jc,kc),GCV_FREE_tau(:,ic,jc,kc),HCV_FREE_tau(:,ic,jc,kc))
      
    DO iVar = 1,NumVar
      FG = LAM5 * dUcons0(iVar,ic,jc,kc) + ax_tau(ic,kc) * (invF(iVar,ic,jc,kc) - FCV_FREE(iVar,ic,jc,kc)) &
                                         + ay_tau(ic,kc) * (invG(iVar,ic,jc,kc) - GCV_FREE(iVar,ic,jc,kc)) &
                                         + az_tau(ic,kc) * (invH(iVar,ic,jc,kc) - HCV_FREE(iVar,ic,jc,kc));

      FGCV_Free_tau(iVar,ic,jc,kc) = (ShockNormalX(ic,kc)* FCV_free_tau(iVar,ic,jc,kc)) + &
                                   & (ShockNormalY(ic,kc)* GCV_free_tau(iVar,ic,jc,kc)) + &
                                   & (ShockNormalZ(ic,kc)* HCV_free_tau(iVar,ic,jc,kc));
       
      FG_f = FGCV_Free_tau(iVar,ic,jc,kc) - ShockVN(ic,kc) * CV_free_tau(iVar,ic,jc,kc);
      
      SUM_L5_U    = SUM_L5_U  + L5(iVar) * (Ucons0(iVar,ic,jc,kc) - CV_free(iVar,ic,jc,kc));
      SUM_L5_FG   = SUM_L5_FG + L5(iVar) * FG;   
      SUM_L5_FG_f = SUM_L5_FG_f + L5(iVar) * FG_f;
    ENDDO
    ! Calculate the ShockAc
    ShockAc(ic,kc) = SUM_L5_FG/SUM_L5_U - (ax_tau(ic,kc) * ShockXtau(ic,kc) &
                                        &+ ay_tau(ic,kc) * ShockYtau(ic,kc) &
                                        &+ az_tau(ic,kc) * ShockZtau(ic,kc))&
                                        &- (SUM_L5_FG_f/SUM_L5_U);

    ShockAc(ic,kc) = ShockAc(ic,kc)/( WallNormalX(ic,kc) * ShockNormalX(ic,kc) &
                                   &+ WallNormalY(ic,kc) * ShockNormalY(ic,kc) &
                                   &+ WallNormalZ(ic,kc) * ShockNormalZ(ic,kc));
   enddo
  enddo
  
!write(*,*),'shockac(:,1)',shockac(:,1) 

     deallocate(CV_free)
     deallocate(FCV_free)
     deallocate(GCV_free)
     deallocate(HCV_free)
     deallocate(CV_free_tau)
     deallocate(FCV_free_tau)
     deallocate(GCV_free_tau)
     deallocate(HCV_free_tau)
     deallocate(FGCV_free_tau)
  
END SUBROUTINE Calculate_ShockAC3D_unsteady

SUBROUTINE ShockRelation3D_unsteady_old
  ! This subroutine calculate the primitive variables after the shock wave	
  use SF_Constant,  only:ik,rk
  use SF_CFD_Global,only:nx_local,Ny,nz_local,&
                         Gamma,Mach_Ref
  ! arrays
  use SF_CFD_Global,only:ShockNormalX,ShockNormalY,ShockNormalZ,ShockV,&
                         WallNormalX,WallNormalY,WallNormalZ,&
                         U,V,W,P,Rho,T,Rho_free,U_free,V_free,W_free,P_free
  ! subs and functions

  implicit none
  real( kind = rk )::Delta_Inf,Masn2Inf,DeltaShkUinf
  integer( kind=ik)::ic,jc,kc
  
  jc = Ny ! Standing for the Shock Surface
  do kc = 1,nz_local
   do ic = 1,nx_local
    ! Velocity ahead of the shocks
    Delta_Inf = U_free(ic,jc,kc) * ShockNormalX(ic,kc) + V_free(ic,jc,kc) * ShockNormalY(ic,kc) + W_free(ic,jc,kc) * ShockNormalZ(ic,kc) &
             &- ShockV(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX(ic,kc) &
                             &+WallNormalY(ic,kc)*ShockNormalY(ic,kc) &
                             &+WallNormalZ(ic,kc)*ShockNormalZ(ic,kc));
    ! Shock normal Mach number^2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
      Masn2Inf = (Mach_Ref * Delta_Inf)**2.d0;
    ! Velocity normal to the shock after the shock
      DeltaShkUinf = 2.d0 * Delta_Inf / (Gamma + 1.d0) * (1.d0 / Masn2Inf - 1.d0);
      
      U(ic,jc,kc) = U_free(ic,jc,kc) + DeltaShkUinf * ShockNormalX(ic,kc);

      V(ic,jc,kc) = V_free(ic,jc,kc) + DeltaShkUinf * ShockNormalY(ic,kc);

      W(ic,jc,kc) = W_free(ic,jc,kc) + DeltaShkUinf * ShockNormalZ(ic,kc);

      P(ic,jc,kc) = P_free(ic,jc,kc) * (1.d0 + 2.d0 * Gamma / (Gamma + 1.d0) * (Masn2Inf - 1.d0));

    Rho(ic,jc,kc) = Rho_free(ic,jc,kc) * (Gamma + 1.d0) * Masn2Inf / ((Gamma - 1.d0) * Masn2Inf + 2.d0);

      T(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / Rho(ic,jc,kc);
   enddo
  enddo
!write(*,*)rho(1,Ny,1)
END SUBROUTINE ShockRelation3D_unsteady_old

SUBROUTINE ShockRelation3D_unsteady
  ! This subroutine calculate the primitive variables after the shock wave	
  use SF_Constant,  only:ik,rk
  use SF_CFD_Global,only:nx_local,Ny,nz_local,&
                         Gamma,Mach_Ref
  ! arrays
  use SF_CFD_Global,only:ShockNormalX,ShockNormalY,ShockNormalZ,ShockV,&
                         WallNormalX,WallNormalY,WallNormalZ,&
                         U,V,W,P,Rho,T,Rho_free,U_free,V_free,W_free,P_free
  ! subs and functions

  implicit none
  real( kind = rk )::Delta_Inf,Delta2
  integer( kind=ik)::ic,jc,kc
  
  jc = Ny ! Standing for the Shock Surface
  do kc = 1,nz_local
   do ic = 1,nx_local
    ! Velocity ahead of the shocks
    Delta_Inf = U_free(ic,jc,kc) * ShockNormalX(ic,kc) + V_free(ic,jc,kc) * ShockNormalY(ic,kc) + W_free(ic,jc,kc) * ShockNormalZ(ic,kc) &
             &- ShockV(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX(ic,kc) &
                             &+WallNormalY(ic,kc)*ShockNormalY(ic,kc) &
                             &+WallNormalZ(ic,kc)*ShockNormalZ(ic,kc));
    Delta2 = (Gamma - 1.d0) * Delta_inf/(Gamma + 1.d0) + 2.d0 * Gamma * P_free(ic,jc,kc)/((Gamma + 1.d0)*Rho_free(ic,jc,kc) * delta_inf);
    ! Velocity normal to the shock after the shock
      
      U(ic,jc,kc) = U_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalX(ic,kc);

      V(ic,jc,kc) = V_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalY(ic,kc);

      W(ic,jc,kc) = W_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalZ(ic,kc);

      P(ic,jc,kc) = 2.0 * Rho_free(ic,jc,kc) * delta_inf * delta_inf/(gamma + 1.0) - (gamma - 1.0) * P_free(ic,jc,kc) / (gamma + 1.0);

    Rho(ic,jc,kc) = delta_inf * Rho_free(ic,jc,kc) / delta2;

      T(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / Rho(ic,jc,kc);
   enddo
  enddo
!write(*,*)rho(1,Ny,1)
END SUBROUTINE ShockRelation3D_unsteady

END Module SFitting