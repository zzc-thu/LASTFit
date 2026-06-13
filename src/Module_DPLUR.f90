Module DPLUR
  ! This module is used to calculate the DPLUR method
  ! In this methods, we store the diagonal terms of the Jacobian matrix
  use SF_Constant,                only: ik,rk,NumVar,dplur_sub_iter,overLAP
  use SF_CFD_Global,              only: nx_local,Ny,nz_local,Gamma,&
                                      & Mach_Ref,Pr_Ref,Re_Ref,Csthlnd_Ref

  use MPI_GLOBAL,                 only: npx,npx0,npz,npz0

  ! some arrays and matrices
  use SF_CFD_Global,              only: Rho,U,V,W,P,T,Ucons0,Cs,&
                                      & dxidx,dxidy,dxidz,dxidt,&
                                      & detadx,detady,detadz,detadt,&
                                      & dzetadx,dzetady,dzetadz,dzetadt,&
                                      & DFhatDU,DHhatDU,DGhatDU,Jaco,&
                                      & Udx,Udy,Udz,Vdx,Vdy,Vdz,Wdx,Wdy,Wdz,&
                                      & Tdx,Tdy,Tdz,dUcons0,&
                                      & DvisFhatDp,DvisGhatDp,DvisHhatDp,&
                                      & DvisFhatDpxi,DvisGhatDpxi,DvisHhatDpxi,&
                                      & DvisFhatDpet,DvisGhatDpet,DvisHhatDpet,&
                                      & DvisFhatDpzt,DvisGhatDpzt,DvisHhatDpzt,&
                                      & Rds_Xi,Rds_Eta,Rds_Zeta,DRDH,DRDHxi,DRDHzeta,&
                                      & DRDV,DUsDH,DUsDHdxi,DUsDHdzeta,DUsDV,&
                                      & L_LUall,D_LUall,U_LUall,L_LUacUm,L_LUacUs,L_LUacH,&
                                      & L_LUacD,IPIV,DIAG,ADU,CDU,DU,DH,DV,DH_xi,DH_zeta,DV_xi,DV_zeta,&
                                      & ShockAcdHdxi,ShockAcdHdzeta,ShockAcdVdxi,ShockAcdVdzeta,ShockH,&
                                      & ShockV,SHockAc,DPLUR_RHS,Rds_Xi,Rds_Eta,Rds_Zeta,CO_AC_UT,AsUs_P1,AsUs_P3,&
                                      & ShockAcdH,ShockAcdV,DT0,dUcons,dShockH,dShockV,UconsOld,ShockH_Old,ShockV_Old,ModelType
  
  implicit none

contains

subroutine JacobiMatInv_CPG
  !===================================================================================!
  !  invFhat = (F * xix  + G * xiy  + H * xiz  + U0 * xit) 
  !  invGhat = (F * etax + G * etay + H * etaz + U0 * etat)
  !  invHhat = (F * zetax + G * zetay + H * zetaz + U0 * zetat)
  !
  !  DFhatDU,DHhatDU,DGhatDU are the derivatives of invFhat,invGhat,invHhat
  !  with respect to U0
  !
  !  DFhatDU = (mat_A * xix + mat_B * xiy + mat_C * xiz + mat_I * xit)
  !  DGhatDU = (mat_A * etax + mat_B * etay + mat_C * etaz + mat_I * etat)
  !  DHhatDU = (mat_A * zetax + mat_B * zetay + mat_C * zetaz + mat_I * zetat)
  !===================================================================================!
  implicit none
  integer(kind = ik) :: ic,jc,kc

  real( kind = rk ):: Kc_Gamma,local_Cs,localHan,localQ2

  ! The local Jacobian matrix of the inviscid flux along x-,y-,z-direction
  ! with respect to the conservative variables
  real(kind = rk), dimension(NumVar,NumVar) :: mat_A,mat_B,mat_C,mat_I
  

  Kc_Gamma = Gamma - 1.0_rk;
  local_Cs = 0.0_rk;
  localHan = 0.0_rk;
  localQ2  = 0.0_rk;

  mat_A = 0.0_rk;
  mat_B = 0.0_rk;
  mat_C = 0.0_rk;
  mat_I = 0.0_rk;

  mat_I(1,1) = 1.0_rk;
  mat_I(2,2) = 1.0_rk;
  mat_I(3,3) = 1.0_rk;
  mat_I(4,4) = 1.0_rk;
  mat_I(5,5) = 1.0_rk;

  do kc = 1-overLAP, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP,nx_local+overLAP
        ! local sound speed
        local_Cs = Cs(ic,jc,kc); 

        ! local enthalpy
        localHan = (Ucons0(5,ic,jc,kc) + P(ic,jc,kc)) / Rho(ic,jc,kc);
        
        ! local sum of square of velocity
        localQ2  = U(ic,jc,kc) * U(ic,jc,kc) &
                &+ V(ic,jc,kc) * V(ic,jc,kc) &
                &+ W(ic,jc,kc) * W(ic,jc,kc);

       ! Define the local mat_A, mat_B, mat_C
       ! Linearized Jacobian matrix of the inviscid flux along x-direction 
        mat_A(1, 2)   =  1.d0;
        
        mat_A(2, 1)   = 0.5d0 * Kc_Gamma * localQ2 - U(ic,jc,kc) * U(ic,jc,kc);
        mat_A(2, 2)   = (3.d0 - Gamma) * U(ic,jc,kc);
        mat_A(2, 3)   = (1.d0 - Gamma) * V(ic,jc,kc);
        mat_A(2, 4)   = (1.d0 - Gamma) * W(ic,jc,kc);
        mat_A(2, 5)   = Gamma - 1.d0;

        mat_A(3, 1)   = -U(ic,jc,kc) * V(ic,jc,kc);
        mat_A(3, 2)   = V(ic,jc,kc);
        mat_A(3, 3)   = U(ic,jc,kc);

        mat_A(4, 1)   = -U(ic,jc,kc) * W(ic,jc,kc);
        mat_A(4, 2)   = W(ic,jc,kc);
        mat_A(4, 4)   = U(ic,jc,kc);

        mat_A(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * U(ic,jc,kc);
        mat_A(5, 2)   = localHan - Kc_Gamma * U(ic,jc,kc) * U(ic,jc,kc);
        mat_A(5, 3)   =          - Kc_Gamma * U(ic,jc,kc) * V(ic,jc,kc);
        mat_A(5, 4)   =          - Kc_Gamma * U(ic,jc,kc) * W(ic,jc,kc);
        mat_A(5, 5)   = Kc_Gamma * U(ic,jc,kc);

        ! Linearized Jacobian matrix of the inviscid flux along y-direction
        mat_B(1, 3)   = 1.d0;

        mat_B(2, 1)   = -U(ic,jc,kc) * V(ic,jc,kc);
        mat_B(2, 2)   =  V(ic,jc,kc);
        mat_B(2, 3)   =  U(ic,jc,kc);

        mat_B(3, 1)   = 0.5d0 * Kc_Gamma * localQ2 - V(ic,jc,kc) * V(ic,jc,kc);
        mat_B(3, 2)   = (1.d0 - Gamma) * U(ic,jc,kc);
        mat_B(3, 3)   = (3.d0 - Gamma) * V(ic,jc,kc);
        mat_B(3, 4)   = (1.d0 - Gamma) * W(ic,jc,kc);
        mat_B(3, 5)   = Gamma - 1.d0;

        mat_B(4, 1)   = -V(ic,jc,kc) * W(ic,jc,kc);
        mat_B(4, 3)   = W(ic,jc,kc);
        mat_B(4, 4)   = V(ic,jc,kc);

        mat_B(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * V(ic,jc,kc);
        mat_B(5, 2)   =          - Kc_Gamma * V(ic,jc,kc) * U(ic,jc,kc);
        mat_B(5, 3)   = localHan - Kc_Gamma * V(ic,jc,kc) * V(ic,jc,kc);
        mat_B(5, 4)   =          - Kc_Gamma * V(ic,jc,kc) * W(ic,jc,kc);
        mat_B(5, 5)   = Kc_Gamma * V(ic,jc,kc);

        ! Linearized Jacobian matrix of the inviscid flux along z-direction
        mat_C(1, 4)   = 1.d0;

        mat_C(2, 1)   = -U(ic,jc,kc) * W(ic,jc,kc);
        mat_C(2, 2)   =  W(ic,jc,kc);
        mat_C(2, 4)   =  U(ic,jc,kc);

        mat_C(3, 1)   = -V(ic,jc,kc) * W(ic,jc,kc);
        mat_C(3, 3)   =  W(ic,jc,kc);
        mat_C(3, 4)   =  V(ic,jc,kc);

        mat_C(4, 1)   = 0.5d0 * Kc_Gamma * localQ2 - W(ic,jc,kc) * W(ic,jc,kc);
        mat_C(4, 2)   = (1.d0 - Gamma) * U(ic,jc,kc);
        mat_C(4, 3)   = (1.d0 - Gamma) * V(ic,jc,kc);
        mat_C(4, 4)   = (3.d0 - Gamma) * W(ic,jc,kc);
        mat_C(4, 5)   = Gamma - 1.d0;

        mat_C(5, 1)   = (0.5d0 * Kc_Gamma * localQ2 -localHan) * W(ic,jc,kc);
        mat_C(5, 2)   =          - Kc_Gamma * W(ic,jc,kc) * U(ic,jc,kc);
        mat_C(5, 3)   =          - Kc_Gamma * W(ic,jc,kc) * V(ic,jc,kc);
        mat_C(5, 4)   = localHan - Kc_Gamma * W(ic,jc,kc) * W(ic,jc,kc);
        mat_C(5, 5)   = Kc_Gamma * W(ic,jc,kc);

        DFhatDU(1:NumVar,1:NumVar,ic,jc,kc) = (mat_A *  dxidx(ic,jc,kc) &
                                            &+ mat_B *  dxidy(ic,jc,kc) &
                                            &+ mat_C *  dxidz(ic,jc,kc) &
                                            &+ mat_I *  dxidt(ic,jc,kc));

        DGhatDU(1:NumVar,1:NumVar,ic,jc,kc) = (mat_A * detadx(ic,jc,kc) &
                                            &+ mat_B * detady(ic,jc,kc) &
                                            &+ mat_C * detadz(ic,jc,kc) &
                                            &+ mat_I * detadt(ic,jc,kc));

        DHhatDU(1:NumVar,1:NumVar,ic,jc,kc) = (mat_A *dzetadx(ic,jc,kc) &
                                            &+ mat_B *dzetady(ic,jc,kc) &
                                            &+ mat_C *dzetadz(ic,jc,kc) &
                                            &+ mat_I *dzetadt(ic,jc,kc));        
      enddo
    enddo
  enddo

end subroutine JacobiMatInv_CPG

subroutine JacobiMat_TSL_CPG_approximation
  implicit none
  integer( kind = ik ):: ic,jc,kc
  ! Here the TSL approximation is used to form the Jacobian matrix for viscous flux
  ! The viscous flux is approximated as the linear combination of the primitive variables
  ! and its derivatives, for xi directions

  
  
end subroutine JacobiMat_TSL_CPG_approximation

! 注意粘性通量不能像无粘性项那样写成其次项的形式
subroutine JacobiMatVis_CPG
  implicit none
  integer( kind = ik ):: ic,jc,kc

  real( kind = rk ):: Kc_Gamma    ! Kc_Gamma = (Gamma - 1.0)
  real( kind = rk ):: gM2         ! gM2 = Gamma * Mach_Ref * Mach_Ref
  real( kind = rk ):: KgM2        ! KgM2 = Kc_Gamma * gM2
  real( kind = rk ):: Cp
  real( kind = rk ):: lclMu,lclMuT ! d mu / d T
  real( kind = rk ):: lclkap       ! Cp * Mu / Pr
  real( kind = rk ):: lclkapT      ! dKap / dT
  real( kind = rk ):: Div_vel
  real( kind = rk ):: Tauxx,Tauxy,Tauxz
  real( kind = rk ):: Tauyy,Tauyz,Tauzz
  ! Viscous Jacobian matrix is built from the viscous flux
  ! visF, visG, visH are the viscous fluxes along x-,y-,z-direction
  ! dvisF/dU0 = dvisF / dP * dP / dU0 = dvisF / dP * mat_UtoCV
  ! P is the primitive variables P = (rho, u, v, w, T)
  ! visF = dFvdP + dFvdPx * dP/dx + dFvdPy * dP/dy + dFvdPz * dP/dz
  ! visG = dGvdP + dGvdPx * dP/dx + dGvdPy * dP/dy + dGvdPz * dP/dz
  ! visH = dHvdP + dHvdPx * dP/dx + dHvdPz * dP/dy + dHvdPz * dP/dz
  ! similar processes can be used for visG and visH
  ! with respect to the primitive variables
  ! dCV/dP = mat_CVtoU, dP/dCV = mat_UtoCV
  real(kind = rk), dimension(NumVar,NumVar) :: mat_UtoCV
  real(kind = rk), dimension(NumVar,NumVar) :: dFvdP,dGvdP,dHvdP
  real(kind = rk), dimension(NumVar,NumVar) :: dFvdPx,dGvdPx,dHvdPx
  real(kind = rk), dimension(NumVar,NumVar) :: dFvdPy,dGvdPy,dHvdPy
  real(kind = rk), dimension(NumVar,NumVar) :: dFvdPz,dGvdPz,dHvdPz 
  real(kind = rk), dimension(NumVar,NumVar) :: DvisFhatDpx,DvisGhatDpx,DvisHhatDpx
  real(kind = rk), dimension(NumVar,NumVar) :: DvisFhatDpy,DvisGhatDpy,DvisHhatDpy
  real(kind = rk), dimension(NumVar,NumVar) :: DvisFhatDpz,DvisGhatDpz,DvisHhatDpz
  
  ! local mat initialization
  mat_UtoCV = 0.0_rk;
  dFvdP = 0.0_rk; dGvdP = 0.0_rk; dHvdP = 0.0_rk;
  dFvdPx = 0.0_rk; dGvdPx = 0.0_rk; dHvdPx = 0.0_rk;
  dFvdPy = 0.0_rk; dGvdPy = 0.0_rk; dHvdPy = 0.0_rk;
  dFvdPz = 0.0_rk; dGvdPz = 0.0_rk; dHvdPz = 0.0_rk;

  DvisFhatDpx = 0.0_rk; DvisGhatDpx = 0.0_rk; DvisHhatDpx = 0.0_rk;
  DvisFhatDpy = 0.0_rk; DvisGhatDpy = 0.0_rk; DvisHhatDpy = 0.0_rk;
  DvisFhatDpz = 0.0_rk; DvisGhatDpz = 0.0_rk; DvisHhatDpz = 0.0_rk;
  
  ! Some local constants
  Kc_Gamma = Gamma - 1.0_rk;

  gM2 = Gamma * Mach_Ref * Mach_Ref;

  KgM2 = Kc_Gamma * gM2;

  Cp = 1.0_rk / (Kc_Gamma * Mach_Ref * Mach_Ref);
  
  do kc = 1, nz_local
   do jc = 1, Ny
    do ic = 1, nx_local
      ! define some local variables
      lclMu = (T(ic,jc,kc)**1.5_rk) * ( 1.0_rk + Csthlnd_Ref )/ (T(ic,jc,kc) + Csthlnd_Ref)
      lclMu = lclMu / Re_Ref;

      lclMuT = sqrt(T(ic,jc,kc)) * ( 1.0_rk + Csthlnd_Ref ) * (T(ic,jc,kc) + 3.0_rk * Csthlnd_Ref)
      lclMuT = lclMuT / (T(ic,jc,kc) + Csthlnd_Ref) / Re_Ref;

      lclkap = Cp * lclMu / Pr_Ref;
      lclkapT= Cp * lclMuT / Pr_Ref
      
      Div_vel= Udx(ic,jc,kc) + Vdy(ic,jc,kc) + Wdz(ic,jc,kc);

      Tauxx = 2.0_rk * lclMu * (Udx(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauyy = 2.0_rk * lclMu * (Vdy(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauzz = 2.0_rk * lclMu * (Wdz(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauxy = lclMu * (Udy(ic,jc,kc) + Vdx(ic,jc,kc));
      Tauxz = lclMu * (Udz(ic,jc,kc) + Wdx(ic,jc,kc));
      Tauyz = lclMu * (Vdz(ic,jc,kc) + Wdy(ic,jc,kc));

      ! mat_UtoCV = dP/dCV
      mat_UtoCV(1,1) = 1.0_rk;

      mat_UtoCV(2,1) = -U(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(2,2) = 1.0_rk / Rho(ic,jc,kc);

      mat_UtoCV(3,1) = -V(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(3,3) = 1.0_rk / Rho(ic,jc,kc);

      mat_UtoCV(4,1) = -W(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(4,4) = 1.0_rk / Rho(ic,jc,kc);

      mat_UtoCV(5,1) =   KgM2 * (-Ucons0(5,ic,jc,kc)+U(ic,jc,kc)*U(ic,jc,kc)+&
                        &    V(ic,jc,kc)*V(ic,jc,kc)+W(ic,jc,kc)*W(ic,jc,kc));
      mat_UtoCV(5,2) = - KgM2 * U(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(5,3) = - KgM2 * V(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(5,4) = - KgM2 * W(ic,jc,kc) / Rho(ic,jc,kc);
      mat_UtoCV(5,5) =   KgM2               / Rho(ic,jc,kc);
      
      
      ! define the dFvdP, dGvdP, dHvdP
      ! dFvdP
      dFvdP(2,5) = lclMuT / lclMu * Tauxx
      
      dFvdP(3,5) = lclMuT / lclMu * Tauxy
      
      dFvdP(4,5) = lclMuT / lclMu * Tauxz
      
      dFvdP(5,2) = Tauxx
      dFvdP(5,3) = Tauxy
      dFvdP(5,4) = Tauxz
      dFvdP(5,5) = lclMuT / lclMu * ( U(ic,jc,kc) * Tauxx &
                                   &+ V(ic,jc,kc) * Tauxy &
                                   &+ V(ic,jc,kc) * Tauxy) &
                                   &+ lclkapT * Tdx(ic,jc,kc)
      
      ! dGvdP 
      dGvdP(2,5) = lclMuT / lclMu * Tauxy

      dGvdP(3,5) = lclMuT / lclMu * Tauyy

      dGvdP(4,5) = lclMuT / lclMu * Tauyz

      dGvdP(5,2) = Tauxy
      dGvdP(5,3) = Tauyy
      dGvdP(5,4) = Tauyz
      dGvdP(5,5) = lclMuT / lclMu * ( U(ic,jc,kc) * Tauxy &
                                   &+ V(ic,jc,kc) * Tauyy &
                                   &+ V(ic,jc,kc) * Tauyz) &
                                   &+ lclkapT * Tdy(ic,jc,kc)

      ! dHvdP
      dHvdP(2,5) = lclMuT / lclMu * Tauxz

      dHvdP(3,5) = lclMuT / lclMu * Tauyz

      dHvdP(4,5) = lclMuT / lclMu * Tauzz

      dHvdP(5,2) = Tauxz
      dHvdP(5,3) = Tauyz
      dHvdP(5,4) = Tauzz
      dHvdP(5,5) = lclMuT / lclMu * ( U(ic,jc,kc) * Tauxz &
                                   &+ V(ic,jc,kc) * Tauyz &
                                   &+ V(ic,jc,kc) * Tauzz) &
                                   &+ lclkapT * Tdz(ic,jc,kc)
      ! define the dFvdPx, dFvdPy, dFvdPz
      ! dFvdPx
      dFvdPx(2,2) = 4.0_rk * lclMu * Udx(ic,jc,kc) / 3.0_rk
      dFvdPx(3,3) = lclMu;
      dFvdPx(4,4) = lclMu;
      dFvdPx(5,2) = 4.0_rk * lclMu * U(ic,jc,kc) / 3.0_rk;
      dFvdPx(5,3) = lclMu * V(ic,jc,kc);
      dFvdPx(5,4) = lclMu * W(ic,jc,kc);
      dFvdPx(5,5) = lclkap

      ! dFvdPy
      dFvdPy(2,3) = -2.0_rk * lclMu / 3.0_rk;
      dFvdPy(3,2) = lclMu;
      dFvdPy(5,2) = lclMu * V(ic,jc,kc);
      dFvdPy(5,3) = - 2.0_rk * lclMu * U(ic,jc,kc) / 3.0_rk;
      
      ! dFvdPz
      dFvdPz(2,4) = -2.0_rk * lclMu / 3.0_rk;
      dFvdPz(4,2) = lclMu;
      dFvdPz(5,2) = lclMu * W(ic,jc,kc);
      dFvdPz(5,4) = - 2.0_rk * lclMu * U(ic,jc,kc) / 3.0_rk;

      ! define the dGvdPx, dGvdPy, dGvdPz
      ! dGvdPx
      dGvdPx(2,3) = lclMu;
      dGvdPx(3,2) = -2.0_rk * lclMu / 3.0_rk;
      dGvdPx(5,2) = 2.0_rk * lclMu * V(ic,jc,kc) / 3.0_rk;
      dGvdPx(5,3) = lclMu * U(ic,jc,kc);

      ! dGvdPy
      dGvdPy(2,2) = lclMu;
      dGvdPy(3,3) = 4.0_rk * lclMu * Vdy(ic,jc,kc) / 3.0_rk;
      dGvdPy(4,4) = lclMu;
      dGvdPy(5,2) = lclMu * U(ic,jc,kc);
      dGvdPy(5,3) = 4.0_rk * lclMu * V(ic,jc,kc) / 3.0_rk;
      dGvdPy(5,4) = lclMu * W(ic,jc,kc);
      dGvdPy(5,5) = lclkap

      ! dGvdPz
      dGvdPz(3,4) = -2.0_rk * lclMu / 3.0_rk;
      dGvdPz(4,3) = lclMu;
      dGvdPz(5,3) = lclMu * W(ic,jc,kc);
      dGvdPz(5,4) = -2.0_rk * lclMu * V(ic,jc,kc) / 3.0_rk;

      ! define the dHvdPx, dHvdPy, dHvdPz
      ! dHvdPx
      dHvdPx(2,4) = lclMu;
      dHvdPx(4,2) = -2.0_rk * lclMu / 3.0_rk;
      dHvdPx(5,2) = -2.0_rk * lclMu * W(ic,jc,kc) / 3.0_rk;
      dHvdPx(5,4) = lclMu * U(ic,jc,kc);

      ! dHvdPy
      dHvdPy(3,4) = lclMu;
      dHvdPy(4,3) = -2.0_rk * lclMu / 3.0_rk;
      dHvdPy(5,3) = -2.0_rk * lclMu * W(ic,jc,kc) / 3.0_rk;
      dHvdPy(5,4) = lclMu * V(ic,jc,kc);

      ! dHvdPz
      dHvdPz(2,2) = lclMu;
      dHvdPz(3,3) = lclMu;
      dHvdPz(4,4) = 4.0_rk * lclMu * Wdz(ic,jc,kc) / 3.0_rk;
      dHvdPz(5,2) = lclMu * U(ic,jc,kc);
      dHvdPz(5,3) = lclMu * V(ic,jc,kc);
      dHvdPz(5,4) = 4.0_rk * lclMu * W(ic,jc,kc) / 3.0_rk;
      dHvdPz(5,5) = lclkap
      
      ! p
      DvisFhatDp(1:NumVar,1:NumVar,ic,jc,kc) =(dFvdP * dxidx(ic,jc,kc) &
                                            &+ dGvdP * dxidy(ic,jc,kc) &
                                            &+ dHvdP * dxidz(ic,jc,kc)) / Jaco(ic,jc,kc);
                                            
      DvisGhatDp(1:NumVar,1:NumVar,ic,jc,kc) =(dFvdP * detadx(ic,jc,kc) &
                                            &+ dGvdP * detady(ic,jc,kc) &
                                            &+ dHvdP * detadz(ic,jc,kc)) / Jaco(ic,jc,kc);

      DvisHhatDp(1:NumVar,1:NumVar,ic,jc,kc) =(dFvdP * dzetadx(ic,jc,kc) &
                                            &+ dGvdP * dzetady(ic,jc,kc) &
                                            &+ dHvdP * dzetadz(ic,jc,kc)) / Jaco(ic,jc,kc);
      
      ! pdx
      DvisFhatDpx  =(dFvdPx * dxidx(ic,jc,kc) &
                 &+  dGvdPx * dxidy(ic,jc,kc) &
                 &+  dHvdPx * dxidz(ic,jc,kc)) / Jaco(ic,jc,kc);
      
      DvisGhatDpx  =(dFvdPx * detadx(ic,jc,kc) &
                 &+  dGvdPx * detady(ic,jc,kc) &
                 &+  dHvdPx * detadz(ic,jc,kc)) / Jaco(ic,jc,kc);
      
      DvisHhatDpx  =(dFvdPx * dzetadx(ic,jc,kc) &
                 &+  dGvdPx * dzetady(ic,jc,kc) &
                 &+  dHvdPx * dzetadz(ic,jc,kc)) / Jaco(ic,jc,kc);

      ! pdy
      DvisFhatDpy  =(dFvdPy * dxidx(ic,jc,kc) &
                 &+  dGvdPy * dxidy(ic,jc,kc) &
                 &+  dHvdPy * dxidz(ic,jc,kc)) / Jaco(ic,jc,kc);

      DvisGhatDpy  =(dFvdPy * detadx(ic,jc,kc) &
                 &+  dGvdPy * detady(ic,jc,kc) &
                 &+  dHvdPy * detadz(ic,jc,kc)) / Jaco(ic,jc,kc);

      DvisHhatDpy  =(dFvdPy * dzetadx(ic,jc,kc) &
                 &+  dGvdPy * dzetady(ic,jc,kc) &
                 &+  dHvdPy * dzetadz(ic,jc,kc)) / Jaco(ic,jc,kc);
      
      ! pdz
      DvisFhatDpz  =(dFvdPz * dxidx(ic,jc,kc) &
                 &+  dGvdPz * dxidy(ic,jc,kc) &
                 &+  dHvdPz * dxidz(ic,jc,kc)) / Jaco(ic,jc,kc);
      
      DvisGhatDpz  =(dFvdPz * detadx(ic,jc,kc) &
                 &+  dGvdPz * detady(ic,jc,kc) &
                 &+  dHvdPz * detadz(ic,jc,kc)) / Jaco(ic,jc,kc);

      DvisHhatDpz  =(dFvdPz * dzetadx(ic,jc,kc) &
                 &+  dGvdPz * dzetady(ic,jc,kc) &
                 &+  dHvdPz * dzetadz(ic,jc,kc)) / Jaco(ic,jc,kc);
      ! 计算过程中，我们是在xi,eta,zeta坐标系下进行计算的，
      ! 因此 DvisFhatDpx * d/dx + DvisFhatDpy * d/dy + DvisFhatDpz * d/dz
      ! 会转换到对应的xi,eta,zeta坐标系下进行求解
      !
      !
      ! Fvis 转换为对xi,eta,zeta的偏导数的系数
      DvisFhatDpxi(1:NumVar,1:NumVar,ic,jc,kc) = DvisFhatDpx * dxidx(ic,jc,kc) &
                                              &+ DvisFhatDpy * dxidy(ic,jc,kc) &
                                              &+ DvisFhatDpz * dxidz(ic,jc,kc);
      
      DvisFhatDpet(1:NumVar,1:NumVar,ic,jc,kc) = DvisFhatDpx * detadx(ic,jc,kc) &
                                              &+ DvisFhatDpy * detady(ic,jc,kc) &
                                              &+ DvisFhatDpz * detadz(ic,jc,kc);

      DvisFhatDpzt(1:NumVar,1:NumVar,ic,jc,kc) = DvisFhatDpx * dzetadx(ic,jc,kc) &
                                              &+ DvisFhatDpy * dzetady(ic,jc,kc) &
                                              &+ DvisFhatDpz * dzetadz(ic,jc,kc);
      ! Gvis 转化为对xi,eta,zeta的偏导数的系数
      DvisGhatDpxi(1:NumVar,1:NumVar,ic,jc,kc) = DvisGhatDpx * dxidx(ic,jc,kc) &
                                              &+ DvisGhatDpy * dxidy(ic,jc,kc) &
                                              &+ DvisGhatDpz * dxidz(ic,jc,kc);
      
      DvisGhatDpet(1:NumVar,1:NumVar,ic,jc,kc) = DvisGhatDpx * detadx(ic,jc,kc) &
                                              &+ DvisGhatDpy * detady(ic,jc,kc) &
                                              &+ DvisGhatDpz * detadz(ic,jc,kc);

      DvisGhatDpzt(1:NumVar,1:NumVar,ic,jc,kc) = DvisGhatDpx * dzetadx(ic,jc,kc) &
                                              &+ DvisGhatDpy * dzetady(ic,jc,kc) &
                                              &+ DvisGhatDpz * dzetadz(ic,jc,kc);
      ! Hvis 转化为对xi,eta,zeta的偏导数的系数
      DvisHhatDpxi(1:NumVar,1:NumVar,ic,jc,kc) = DvisHhatDpx * dxidx(ic,jc,kc) &
                                              &+ DvisHhatDpy * dxidy(ic,jc,kc) &
                                              &+ DvisHhatDpz * dxidz(ic,jc,kc);
      
      DvisHhatDpet(1:NumVar,1:NumVar,ic,jc,kc) = DvisHhatDpx * detadx(ic,jc,kc) &
                                              &+ DvisHhatDpy * detady(ic,jc,kc) &
                                              &+ DvisHhatDpz * detadz(ic,jc,kc);

      DvisHhatDpzt(1:NumVar,1:NumVar,ic,jc,kc) = DvisHhatDpx * dzetadx(ic,jc,kc) &
                                              &+ DvisHhatDpy * dzetady(ic,jc,kc) &
                                              &+ DvisHhatDpz * dzetadz(ic,jc,kc);
    enddo
   enddo
  enddo

end subroutine JacobiMatVis_CPG  

subroutine Calculate_ShkAcDeri
  ! Give the derivative of the shock acceleration
  use SF_Constant, only: ik, rk
  
  use SF_CFD_Global, only: gamma, mach_ref, nx_local, Ny, nz_local

  ! Arrays
  use SF_CFD_Global, only: ShockH, ShockHdxi, ShockHdzeta, ShockV, ShockVdxi, ShockVdzeta,&
                       &  WallSXdxi, WallSYdxi, WallSZdxi, WallSXdzeta, WallSYdzeta, WallSZdzeta,&
                       &  WallNormalX, WallNormalY, WallNormalZ, WallNormalXdxi, WallNormalYdxi, WallNormalZdxi,&
                       &  WallNormalXdzeta, WallNormalYdzeta, WallNormalZdzeta, Heta, Hetadxi, Hetadeta, Hetadzeta,&
                       &  U_INF, V_INF, W_INF, P_INF, RHO_INF, CV_INF, FINV_INF, GINV_INF, HINV_INF, dUcons0,&
                       &  ShockAcdH,ShockAcdHdxi,ShockAcdHdzeta,ShockAcdV,ShockAcdVdxi,ShockAcdVdzeta
  implicit none
  integer( kind = ik ):: ic,jc,kc

  real( kind = rk ):: lcl_shockh,lcl_shockhb,lcl_shockhdxi,lcl_shockhdxib
  real( kind = rk ):: lcl_shockhdzeta,lcl_shockhdzetab,lcl_shockv,lcl_shockvb
  real( kind = rk ):: lcl_shockvdxi,lcl_shockvdxib,lcl_shockvdzeta,lcl_shockvdzetab
  real( kind = rk ):: lcl_shockac,lcl_shockacb
  real( kind = rk ):: lcl_wallsxdxi,lcl_wallsydxi,lcl_wallszdxi,lcl_wallsxdzeta,lcl_wallsydzeta,lcl_wallszdzeta
  real( kind = rk ):: lcl_wallnormalx,lcl_wallnormaly,lcl_wallnormalz
  real( kind = rk ):: lcl_wallnormalxdxi,lcl_wallnormalydxi,lcl_wallnormalzdxi
  real( kind = rk ):: lcl_wallnormalxdzeta,lcl_wallnormalydzeta,lcl_wallnormalzdzeta
  real( kind = rk ):: lcl_heta,lcl_hetadxi,lcl_hetadeta,lcl_hetadzeta
  real( kind = rk ):: lcl_cv_inf(5),lcl_finv_inf(5),lcl_ginv_inf(5),lcl_hinv_inf(5)
  real( kind = rk ):: lcl_ducons0(5),lcl_u_inf,lcl_v_inf,lcl_w_inf,lcl_p_inf,lcl_rho_inf

  ! For steady free stream flow
  ! Those values are constant for steady flow calculation
     ! Basic variables
     lcl_rho_inf = Rho_inf
     lcl_u_inf = U_inf
     lcl_v_inf = V_inf
     lcl_w_inf = W_inf
     lcl_p_inf = P_inf
     ! COnservative variables
     lcl_cv_inf = CV_inf   
     
     ! Inviscid fluxes
     lcl_finv_inf = Finv_inf
     lcl_ginv_inf = Ginv_inf
     lcl_hinv_inf = Hinv_inf
  
  jc = Ny
  do kc = 1, nz_local
    do ic = 1, nx_local
      
      ! Define the relative geometry variables
      lcl_wallsxdxi   = WallSXdxi(ic,kc)
      lcl_wallsydxi   = WallSYdxi(ic,kc)
      lcl_wallszdxi   = WallSZdxi(ic,kc)

      lcl_wallsxdzeta = WallSXdzeta(ic,kc)
      lcl_wallsydzeta = WallSYdzeta(ic,kc)
      lcl_wallszdzeta = WallSZdzeta(ic,kc)

      lcl_wallnormalx   = WallNormalX(ic,kc)
      lcl_wallnormaly   = WallNormalY(ic,kc)
      lcl_wallnormalz   = WallNormalZ(ic,kc)

      lcl_wallnormalxdxi   = WallNormalXdxi(ic,kc)
      lcl_wallnormalydxi   = WallNormalYdxi(ic,kc)
      lcl_wallnormalzdxi   = WallNormalZdxi(ic,kc)

      lcl_wallnormalxdzeta = WallNormalXdzeta(ic,kc)
      lcl_wallnormalydzeta = WallNormalYdzeta(ic,kc)
      lcl_wallnormalzdzeta = WallNormalZdzeta(ic,kc)

      lcl_heta     = Heta(ic,jc,kc)
      lcl_hetadxi  = Hetadxi(ic,jc,kc)
      lcl_hetadeta = Hetadeta(ic,jc,kc)
      lcl_hetadzeta= Hetadzeta(ic,jc,kc)

      ! define the residual behind the shock surface
      lcl_ducons0 = dUcons0(:,ic,jc,kc)

      ! define the input the automatic differentiation tool Tapenade
      lcl_shockacb = 1.0_rk
      ! as in transverse mode AD, for Y=F(x), the output is dY = dF/dX * dX
      ! as in reverse mode AD, for Y=F(x), the output is dX = (dF/dX)^T * dY
      ! Therefore, the dF/dX , which is the jacobian matrix, should be the 
      ! output of the AD tool, with dY = 1.0_rk

      lcl_shockh      = ShockH(ic,kc)
      lcl_shockhdxi   = ShockHdxi(ic,kc)
      lcl_shockhdzeta = ShockHdzeta(ic,kc)
      lcl_shockv      = ShockV(ic,kc)
      lcl_shockvdxi   = ShockVdxi(ic,kc)
      lcl_shockvdzeta = ShockVdzeta(ic,kc)
      
      ! 这里的导数没有考虑对应的残差函数对于H和V的导数
     call CALCULATESHOCKAC_B(lcl_shockh, lcl_shockhb, lcl_shockhdxi, lcl_shockhdxib, &
          & lcl_shockhdzeta, lcl_shockhdzetab, lcl_shockv, lcl_shockvb, lcl_shockvdxi, lcl_shockvdxib, &
          & lcl_shockvdzeta, lcl_shockvdzetab, lcl_shockac, lcl_shockacb, &
          & lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, lcl_wallszdzeta, &
          & lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, &
          & lcl_wallnormalxdxi, lcl_wallnormalydxi, lcl_wallnormalzdxi, &
          & lcl_wallnormalxdzeta, lcl_wallnormalydzeta, lcl_wallnormalzdzeta, &
          & lcl_heta, lcl_hetadxi, lcl_hetadeta, lcl_hetadzeta, lcl_cv_inf, lcl_finv_inf, lcl_ginv_inf, &
          & lcl_hinv_inf, lcl_ducons0, lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_p_inf, lcl_rho_inf, gamma, &
          & mach_ref)

      ! output *b is the derivative generated by AD tool Tapenade
      !lcl_shockacb = 0.0_8
          
      ShockAcdH(ic,kc)        = lcl_shockhb
      ShockAcdHdxi(ic,kc)     = lcl_shockhdxib
      ShockAcdHdzeta(ic,kc)   = lcl_shockhdzetab

      ShockAcdV(ic,kc)        = lcl_shockvb
      ShockAcdVdxi(ic,kc)     = lcl_shockvdxib
      ShockAcdVdzeta(ic,kc)   = lcl_shockvdzetab
      
    enddo
  enddo

end subroutine Calculate_ShkAcDeri

subroutine Calculate_Us_ShockDeri
  ! Calculate the derivative of the R-H relationship with respect to the ShockH and ShockV
  ! UshkdH, UshkdHdxi, UshkdHdzeta, UshkdV
  use SF_Constant, only: ik, rk
  
  use SF_CFD_Global, only: gamma, mach_ref, nx_local, Ny, nz_local

  ! Arrays
  use SF_CFD_Global, only: ShockH, ShockHdxi, ShockHdzeta, ShockV,&
                       &  WallSXdxi, WallSYdxi, WallSZdxi, WallSXdzeta, WallSYdzeta, WallSZdzeta,&
                       &  WallNormalX, WallNormalY, WallNormalZ, WallNormalXdxi, WallNormalYdxi, WallNormalZdxi,&
                       &  WallNormalXdzeta, WallNormalYdzeta, WallNormalZdzeta, Heta, Hetadxi, Hetadeta, Hetadzeta,&
                       &  U_INF, V_INF, W_INF, P_INF, RHO_INF,&
                       &  DUsDH,DUsDHdxi,DUsDHdzeta,DUsDV

    implicit none
    integer( kind = ik ):: ic,jc,kc

    real( kind = rk ):: lcl_shockh,lcl_shockhb
    real( kind = rk ):: lcl_shockhdxi,lcl_shockhdxib
    real( kind = rk ):: lcl_shockhdzeta,lcl_shockhdzetab
    real( kind = rk ):: lcl_shockv,lcl_shockvb
    real( kind = rk ):: lcl_u_inf,lcl_v_inf,lcl_w_inf,lcl_p_inf,lcl_rho_inf
    real( kind = rk ):: lcl_heta,lcl_hetadxi,lcl_hetadeta,lcl_hetadzeta
    real( kind = rk ):: lcl_wallsxdxi,lcl_wallsydxi,lcl_wallszdxi
    real( kind = rk ):: lcl_wallsxdzeta,lcl_wallsydzeta,lcl_wallszdzeta
    real( kind = rk ):: lcl_wallnormalx,lcl_wallnormaly,lcl_wallnormalz
    real( kind = rk ):: lcl_wallnormalxdxi,lcl_wallnormalydxi,lcl_wallnormalzdxi
    real( kind = rk ):: lcl_wallnormalxdzeta,lcl_wallnormalydzeta,lcl_wallnormalzdzeta
    real( kind = rk ):: cv_1b,cv_2b,cv_3b,cv_4b,cv_5b

    lcl_u_inf = U_inf
    lcl_v_inf = V_inf
    lcl_w_inf = W_inf
    lcl_p_inf = P_inf
    lcl_rho_inf = Rho_inf
    
    cv_1b = 1.0_8
    cv_2b = 1.0_8 
    cv_3b = 1.0_8
    cv_4b = 1.0_8
    cv_5b = 1.0_8 

    jc = Ny
    do kc = 1, nz_local
      do ic = 1, nx_local
        
        ! local variables
        lcl_shockh      = ShockH(ic,kc)
        lcl_shockhdxi   = ShockHdxi(ic,kc)
        lcl_shockhdzeta = ShockHdzeta(ic,kc)

        lcl_shockv      = ShockV(ic,kc)

        lcl_heta        = Heta(ic,jc,kc)
        lcl_hetadxi     = Hetadxi(ic,jc,kc)
        lcl_hetadeta    = Hetadeta(ic,jc,kc)
        lcl_hetadzeta   = Hetadzeta(ic,jc,kc)

        lcl_wallsxdxi   = WallSXdxi(ic,kc)
        lcl_wallsydxi   = WallSYdxi(ic,kc)
        lcl_wallszdxi   = WallSZdxi(ic,kc)

        lcl_wallsxdzeta = WallSXdzeta(ic,kc)
        lcl_wallsydzeta = WallSYdzeta(ic,kc)
        lcl_wallszdzeta = WallSZdzeta(ic,kc)

        lcl_wallnormalx   = WallNormalX(ic,kc)
        lcl_wallnormaly   = WallNormalY(ic,kc)
        lcl_wallnormalz   = WallNormalZ(ic,kc)

        lcl_wallnormalxdxi   = WallNormalXdxi(ic,kc)
        lcl_wallnormalydxi   = WallNormalYdxi(ic,kc)
        lcl_wallnormalzdxi   = WallNormalZdxi(ic,kc)

        lcl_wallnormalxdzeta = WallNormalXdzeta(ic,kc)
        lcl_wallnormalydzeta = WallNormalYdzeta(ic,kc)
        lcl_wallnormalzdzeta = WallNormalZdzeta(ic,kc)

        call Calculate_UshkDeri_CV1_B(lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_p_inf, lcl_rho_inf, &
                      & lcl_shockv, lcl_shockvb, mach_ref, gamma, lcl_shockh, lcl_shockhb, lcl_shockhdxi, &
                      & lcl_shockhdxib, lcl_shockhdzeta, lcl_shockhdzetab, lcl_heta, lcl_hetadxi, lcl_hetadeta, &
                      & lcl_hetadzeta, lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, &
                      & lcl_wallszdzeta, lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, lcl_wallnormalxdxi, &
                      & lcl_wallnormalydxi, lcl_wallnormalzdxi, lcl_wallnormalxdzeta, lcl_wallnormalydzeta, &
                      & lcl_wallnormalzdzeta, cv_1b)
        ! dus1/dh, dus1/dhxi, dus1/dhzeta, dus1/dv
             DUsDH(1,ic,kc) = lcl_shockhb;
          DUsDHdxi(1,ic,kc) = lcl_shockhdxib;
        DUsDHdzeta(1,ic,kc) = lcl_shockhdzetab;
             DUsDV(1,ic,kc) = lcl_shockvb;

        call Calculate_UshkDeri_CV2_B(lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_rho_inf, &
                      & lcl_shockv, lcl_shockvb, mach_ref, gamma, lcl_shockh, lcl_shockhb, lcl_shockhdxi, &
                      & lcl_shockhdxib, lcl_shockhdzeta, lcl_shockhdzetab, lcl_heta, lcl_hetadxi, lcl_hetadeta, &
                      & lcl_hetadzeta, lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, &
                      & lcl_wallszdzeta, lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, lcl_wallnormalxdxi, &
                      & lcl_wallnormalydxi, lcl_wallnormalzdxi, lcl_wallnormalxdzeta, lcl_wallnormalydzeta, &
                      & lcl_wallnormalzdzeta, cv_2b)     
        ! dus2/dh, dus2/dhxi, dus2/dhzeta, dus2/dv
             DUsDH(2,ic,kc) = lcl_shockhb;
          DUsDHdxi(2,ic,kc) = lcl_shockhdxib;
        DUsDHdzeta(2,ic,kc) = lcl_shockhdzetab;
             DUsDV(2,ic,kc) = lcl_shockvb;

        call Calculate_UshkDeri_CV3_B(lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_rho_inf, &
                      & lcl_shockv, lcl_shockvb, mach_ref, gamma, lcl_shockh, lcl_shockhb, lcl_shockhdxi, &
                      & lcl_shockhdxib, lcl_shockhdzeta, lcl_shockhdzetab, lcl_heta, lcl_hetadxi, lcl_hetadeta, &
                      & lcl_hetadzeta, lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, &
                      & lcl_wallszdzeta, lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, lcl_wallnormalxdxi, &
                      & lcl_wallnormalydxi, lcl_wallnormalzdxi, lcl_wallnormalxdzeta, lcl_wallnormalydzeta, &
                      & lcl_wallnormalzdzeta, cv_3b)
        ! dus3/dh, dus3/dhxi, dus3/dhzeta, dus3/dv
             DUsDH(3,ic,kc) = lcl_shockhb;
          DUsDHdxi(3,ic,kc) = lcl_shockhdxib;
        DUsDHdzeta(3,ic,kc) = lcl_shockhdzetab;
             DUsDV(3,ic,kc) = lcl_shockvb;

        call Calculate_UshkDeri_CV4_B(lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_rho_inf, &
                      & lcl_shockv, lcl_shockvb, mach_ref, gamma, lcl_shockh, lcl_shockhb, lcl_shockhdxi, &
                      & lcl_shockhdxib, lcl_shockhdzeta, lcl_shockhdzetab, lcl_heta, lcl_hetadxi, lcl_hetadeta, &
                      & lcl_hetadzeta, lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, &
                      & lcl_wallszdzeta, lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, lcl_wallnormalxdxi, &
                      & lcl_wallnormalydxi, lcl_wallnormalzdxi, lcl_wallnormalxdzeta, lcl_wallnormalydzeta, &
                      & lcl_wallnormalzdzeta, cv_4b)
        ! dus4/dh, dus4/dhxi, dus4/dhzeta, dus4/dv
             DUsDH(4,ic,kc) = lcl_shockhb;
          DUsDHdxi(4,ic,kc) = lcl_shockhdxib;
        DUsDHdzeta(4,ic,kc) = lcl_shockhdzetab;
             DUsDV(4,ic,kc) = lcl_shockvb;
            
        call Calculate_UshkDeri_CV5_B(lcl_u_inf, lcl_v_inf, lcl_w_inf, lcl_p_inf, lcl_rho_inf, &
                      & lcl_shockv, lcl_shockvb, mach_ref, gamma, lcl_shockh, lcl_shockhb, lcl_shockhdxi, &
                      & lcl_shockhdxib, lcl_shockhdzeta, lcl_shockhdzetab, lcl_heta, lcl_hetadxi, lcl_hetadeta, &
                      & lcl_hetadzeta, lcl_wallsxdxi, lcl_wallsydxi, lcl_wallszdxi, lcl_wallsxdzeta, lcl_wallsydzeta, &
                      & lcl_wallszdzeta, lcl_wallnormalx, lcl_wallnormaly, lcl_wallnormalz, lcl_wallnormalxdxi, &
                      & lcl_wallnormalydxi, lcl_wallnormalzdxi, lcl_wallnormalxdzeta, lcl_wallnormalydzeta, &
                      & lcl_wallnormalzdzeta, cv_5b)
        ! dus5/dh, dus5/dhxi, dus5/dhzeta, dus5/dv
             DUsDH(5,ic,kc) = lcl_shockhb;
          DUsDHdxi(5,ic,kc) = lcl_shockhdxib;
        DUsDHdzeta(5,ic,kc) = lcl_shockhdzetab;
             DUsDV(5,ic,kc) = lcl_shockvb;
        
      enddo
    enddo

end subroutine Calculate_Us_ShockDeri

subroutine AssemblingQuisi_TridiagSystem
  ! this subroutine calculate the major matrix elements of the tridiagonal system
  ! Here, Block LU Decomposition is used to obtain the tridiagonal system
  use SF_Constant, only: ik, rk, NumVar
  use SF_CFD_Global, only: nx_local, Ny, nz_local
  
  implicit none
  ! local variables
  integer( kind = ik ):: ic,jc,kc,iVar,jVar,kVar
  real( kind = rk ),DIMENSION(NumVar,NumVar)::LMat_Local,DMat_Local,UMat_Local
  real( kind = rk ),DIMENSION(NumVar)::AcUmMat_Local,AcUsMat_Local
  real( kind = rk ):: CO_JP,CO_JM

  integer(kind = ik)::info
  
  ! Matrix Elements
  call JacobiMatInv_CPG       ! Calculate the DFhatDU, DGhatDU, DHhatDU, Here JacobiMat is not divided by Jaco
  !call Calculate_ShkAcDeri    ! Calculate the Shock Acceleration Derivative
  !call Calculate_Us_ShockDeri ! Calculate the derivative of the R-H relationship with respect to the ShockH and ShockV
  call Calculate_GeneralJacobian   ! Calculate the Jacobian Matrix

  DO kc = 1, nz_local
    DO ic = 1, nx_local
      ! jc = 1, ! Wall Surface Boundary
      jc = 1;
      UMat_Local = 0.0_rk

      ! 壁面网格上只有对应的对角元和上三角部分
      CO_JP = jaco(ic,jc,kc)/jaco(ic,jc+1,kc);
      ! 首先是和第一行相邻的上三角部分的次对角元素
      DO jVar = 1, 4 ! 注意这里我们不更新对应的能量方程部分，仅保留能量方程的对应的对角元素
        UMat_Local(1,jVar) = CO_JP * DGhatDU(1,jVar,ic,jc+1,kc);
      ENDDO
      ! 然后是对角元素(首先是谱半径部分)
      DMat_Local = 0.0_rk;
      DO iVar = 1, NumVar
        DMat_Local(iVar,iVar) = DIAG(ic,jc,kc); 
      ENDDO
      ! 考虑到壁面边界条件，仅仅存在密度和总能量的变化
      ! 密度变量的处理
        DMat_Local(1,1) = DMat_Local(1,1) - DGhatDU(1,1,ic,jc,kc);
      ! 对应的总能量的处理
        DMat_Local(5,1) = -DIAG(ic,jc,kc) * ( Ucons0(5,ic,jc,kc) / Rho(ic,jc,kc) )

      ! 这里我们调用经典的Lapack 库来求解线性系统的问题
      call dgetrf(NumVar,NumVar,DMat_Local,NumVar,ipiv(:,ic,jc,kc),info)
      call dgetrs('N',NumVar,NumVar,DMat_Local,NumVar,ipiv(:,ic,jc,kc),UMat_Local,NumVar,info)
      call dgetrs('N',NumVar,1     ,DMat_Local,NumVar,ipiv(:,ic,jc,kc),DRdH(:,ic,jc,kc),NumVar,info)
      call dgetrs('N',NumVar,1     ,DMat_Local,NumVar,ipiv(:,ic,jc,kc),DRdV(:,ic,jc,kc),NumVar,info)
      ! storage 
      DO iVar = 1, NumVar
        DO jVar = 1, NumVar
          D_LUall(iVar,jVar,ic,jc,kc) = DMat_Local(iVar,jVar); ! inv(L11)
          U_LUall(iVar,jVar,ic,jc,kc) = UMat_Local(iVar,jVar); !     U12
        ENDDO
      ENDDO
      ! DRDH, DRDV就直接存储在原始位置上

      ! jc = 2 ~ Ny - 1 ! Inner Points
      DO jc = 2, Ny - 1
        ! 下三角部分的对角块部分 
        CO_JM = - jaco(ic,jc,kc)/(2.d0 * jaco(ic,jc-1,kc));
        DO iVar = 1, NumVar
          DO jVar = 1, NumVar
          ! 下三角部分对应的Jaco矩阵部分
            LMat_Local(iVar,jVar) = CO_JM * DGhatDU(iVar,jVar,ic,jc-1,kc);
          ENDDO
          ! 谱半径对应的部分直接加在对角元素上
            LMat_Local(iVar,iVar) = LMat_Local(iVar,iVar) + CO_JM * Rds_Eta(ic,jc-1,kc);
        ENDDO
        
        ! 下三角部分的次对角元素
        ! L = D - L(m,m-1) * U(m-1,m)
        DO jVar = 1, NumVar
          ! 首先是数组运算部分
          DO iVar = 1, NumVar
            DMat_Local(iVar,jVar) = 0.0_rk;
            DO kVar = 1, NumVar
              ! LMat_Local, UMat_Local 可以直接利用之前的值      注意这里的LMat_Local都是刚更新的 
              DMat_Local(iVar,jVar) = DMat_Local(iVar,jVar) - LMat_Local(iVar,kVar) * UMat_Local(kVar,jVar);
            ENDDO
          ENDDO
          ! 其次是对应的对角部分的谱半径部分
          DMat_Local(jVar,jVar) = DMat_Local(jVar,jVar) + DIAG(ic,jc,kc);
        ENDDO

        ! 对应的上三角部分的数组
        ! JP = jaco(ic,jc,kc)/(2.0*jaco(ic,jc+1,kc))*(B(i,j+1,k) - sigma(i,j+1,k));
        CO_JP = jaco(ic,jc,kc)/(2.0_rk*jaco(ic,jc+1,kc));
        DO iVar = 1, NumVar
          DO jVar = 1, NumVar
            UMat_Local(iVar,jVar) = CO_JP * DGhatDU(iVar,jVar,ic,jc+1,kc);
          ENDDO
            UMat_Local(iVar,iVar) = UMat_Local(iVar,iVar) - CO_JP * Rds_Eta(ic,jc+1,kc);
        ENDDO
        
        ! 注意这里的dRdH(j-1,j)和dRdV(j-1,j)是上一个站位已经计算好的
        ! dRdHj = dRdH - L(j,j-1) * dRdH(j-1,j)
        ! dRdVj = dRdV - L(j,j-1) * dRdV(j-1,j)
        DO iVar = 1, NumVar
          DO jVar = 1, NumVar
            DRdH(iVar,ic,jc,kc) = DRdH(iVar,ic,jc,kc) - LMat_Local(iVar,jVar) * DRdH(jVar,ic,jc-1,kc);
            DRdV(iVar,ic,jc,kc) = DRdV(iVar,ic,jc,kc) - LMat_Local(iVar,jVar) * DRdV(jVar,ic,jc-1,kc);
          ENDDO
        ENDDO
        
        ! 这里我们调用经典的Lapack 库来求解线性系统的问题
        call dgetrf(NumVar,NumVar,DMat_Local,NumVar,ipiv(:,ic,jc,kc),info)
        call dgetrs('N',NumVar,NumVar,DMat_Local,NumVar,ipiv(:,ic,jc,kc),UMat_Local,NumVar,info)
        call dgetrs('N',NumVar,1     ,DMat_Local,NumVar,ipiv(:,ic,jc,kc),DRdH(:,ic,jc,kc),NumVar,info)
        call dgetrs('N',NumVar,1     ,DMat_Local,NumVar,ipiv(:,ic,jc,kc),DRdV(:,ic,jc,kc),NumVar,info)

        ! 存储对应的上三角元素的次对角矩阵，
        ! H和V部分相关的向量直接存储在dRdH和dRdV中
        DO iVar = 1, NumVar
          DO jVar = 1, NumVar
            D_LUall(iVar,jVar,ic,jc,kc) = DMat_Local(iVar,jVar); ! inv(L(j  ,j  )
            U_LUall(iVar,jVar,ic,jc,kc) = UMat_Local(iVar,jVar); !     U(j-1,j  )
            L_LUall(iVar,jVar,ic,jc,kc) = LMat_Local(iVar,jVar); !     L(j  ,j-1)
          ENDDO
        ENDDO

      ENDDO
      ! 这里需要注意的是对应的Ny处的一些量为0或者为单位阵，因此不需要对此进行额外的存储
      ! jc = Ny ! Shock Boundary conditions
      jc = Ny;
      
      ! 注意激波加速度部分有4块，分别是 da_s/dU(Ny-1), da_s/dU(Ny), dAs/dH, dAs/dV
      ! dAsdU(Ny-1)
      CO_JM = - jaco(ic,jc,kc)/jaco(ic,jc-1,kc);
      DO iVar = 1, NumVar
        AcUmMat_Local(iVar) = 0.0_rk;
        ! 首先考虑激波加速表达式中R前面的系数
        DO jVar = 1, NumVar
          AcUmMat_Local(iVar) = AcUmMat_Local(iVar) + CO_AC_UT(jVar,ic,kc) * DGhatDU(iVar,jVar,ic,jc-1,kc);
        ENDDO
        ! 其次考虑dR/dU部分的系数 CO_JM * dRdU(Ny-1)
        AcUmMat_Local(iVar) = CO_JM * AcUmMat_Local(iVar);
      ENDDO

      ! dAsdU(Ny)
      DO iVar = 1, NumVar
         AcUsMat_Local(iVar) = - DIAG(ic,jc,kc) * CO_AC_UT(iVar,ic,kc);  ! 这里仅仅考虑对应的R变化部分
         !AcUsMat_Local(iVar) = AcUsMat_Local(iVar) - AsUs_P1(iVar,ic,kc) - AsUs_P3(iVar,ic,kc);
        DO jVar = 1, NumVar
          AcUsMat_Local(iVar) = AcUsMat_Local(iVar) - CO_AC_UT(jVar,ic,kc) * DGhatDU(iVar,jVar,ic,jc,kc);
          ! 后续是对应的Ny-1部分变量的修正
          AcUsMat_Local(iVar) = AcUsMat_Local(iVar) - AcUmMat_Local(jVar) * UMat_Local(jVar,iVar);  
        ENDDO
      ENDDO
      
      DO iVar = 1, NumVar
        L_LUacUm(iVar,ic,kc) = AcUmMat_Local(iVar);
        L_LUacUs(iVar,ic,kc) = AcUsMat_Local(iVar);
      ENDDO
      
      ! 这里是对应的 dAs/dH位置处的下三角块矩阵对应的元素
      L_LUacH(ic,kc) = -ShockAcdH(ic,kc);
      DO iVar = 1, NumVar
        L_LUacH(ic,kc) =  L_LUacH(ic,kc) &
                        &- L_LUacUm(iVar,ic,kc) * DRdH(iVar,ic,jc-1,kc) & ! 这里的DRdH是已经添加了负号的
                        &+ L_LUacUs(iVar,ic,kc) * DUsdH(iVar,ic,kc);  
      ENDDO 

      ! 这里是对应的 dAs/dV位置处的下三角块矩阵对应的元素（就在对角阵上）
      L_LUacD(ic,kc) = 1.0_rk / DT0 - ShockAcdV(ic,kc) - DT0 * L_LUacH(ic,kc);
      DO iVar = 1, NumVar
        L_LUacD(ic,kc) = L_LUacD(ic,kc) &
                       &- L_LUacUm(iVar,ic,kc) * DRdV(iVar,ic,jc-1,kc) &
                       &+ L_LUacUs(iVar,ic,kc) * DUsdV(iVar,ic,kc);
      ENDDO

    ENDDO
  ENDDO 
  
end subroutine AssemblingQuisi_TridiagSystem

subroutine DPLUR_Line_Solver
  USE MPI_GLOBAL, only: Parallel_Exchange_NumVar, Parallel_Exchange_Surface
  USE FD5_Order,  only: Cal_Deri_SurfDxi_5th_ce,Cal_Deri_SurfDzeta_5th_ce
  implicit none
  ! 这个程序主要用于求解对应的沿着eta方向上的线性方程组
  ! 在对应的子程序 AssemblingQuisi_TridiagSystem中已经计算好了对应的系数矩阵
  ! 这里我们需要求解的就是对应的方程右端项以及对应的线性系统的解
  ! 注意这里我们仅仅在eta方向上准确求解对应的线性系统，
  ! 但是在xi和zeta方向上我们仅仅是迭代近似求解
  ! 因此这里需要额外处理对应的xi和zeta方向的边界条件，以及迭代收敛
  ! L * U * DeltaQ = f
  ! Solve processes are divided into two steps
  ! L * y = f
  ! U * DeltaQ = y
  integer( kind = ik ):: ic,jc,kc,iVar,jVar,kk
  REAL( kind = rk ):: CO_IM,CO_IP,CO_KM,CO_KP

  integer( kind = ik ):: info,dplur_iter

  call AssemblingQuisi_TridiagSystem

  DPLR_Loop: DO dplur_iter = 1, dplur_sub_iter
    
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            DU(iVar,ic,jc,kc) = dUcons0(iVar,ic,jc,kc);
          ENDDO
        ENDDO
      ENDDO
      !jc = Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
           DU(iVar,ic,Ny,kc) = 0.0_rk;
          ENDDO
           DH(ic,kc) =  ShockV(ic,kc);
           DV(ic,kc) = ShockAc(ic,kc);
        ENDDO
    ENDDO
    ! 这里我们需要区分DPLUR的第一次迭代和后续的迭代，我们假设第一次迭代的时候对应的增量为0
    IF( dplur_iter > 1 ) then
    ! 后续的更新则直接使用上一步的DPLUR_RHS的结果
      DO kc = 1, nz_local
        DO jc = 1, Ny - 1
          DO ic = 1, nx_local
            DO iVar = 1, NumVar
              DU(iVar,ic,jc,kc) = DU(iVar,ic,jc,kc) + DPLUR_RHS(iVar,ic,jc,kc);
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO ic = 1, nx_local
        DO kc = 1, nz_local
           ! Us 边界处的激波后的物理量，由RH关系式直接决定
           DO iVar = 1, NumVar
             DU(iVar,ic,Ny,kc) =   DUsDHdxi(iVar,ic,kc) *   DH_xi(ic,kc) +&
                               & DUsDHdzeta(iVar,ic,kc) * DH_zeta(ic,kc) ;
           ENDDO
           ! shock acceleration
           DV(ic,kc) = DV(ic,kc) + ShockAcdHdxi(ic,kc)   *   DH_xi(ic,kc) &
                                &+ ShockAcdVdxi(ic,kc)   *   DV_xi(ic,kc) &
                                &+ ShockAcdHdzeta(ic,kc) * DH_zeta(ic,kc) &
                                &+ ShockAcdVdzeta(ic,kc) * DV_zeta(ic,kc)
           DO kk = 1, NumVar
             DV(ic,kc) = DV(ic,kc) + CO_AC_UT(kk,ic,kc) * DPLUR_RHS(kk,ic,Ny,kc)
           ENDDO
        ENDDO
      ENDDO
      ! Surface Boundary Conditions, modifications
      DO kc = 1, nz_local
        DO ic = 1, nx_local
          DO iVar = 2, NumVar
            DU(iVar,ic,1,kc) = 0.0_rk;
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    
    ! 这里我们考虑dplur的第一次迭代,此时对应的增量部分都为0可以直接
    ! 首先求解对应的 L * y = f 的线性系统中的y
    ! 首先f部分我们分别存储在DU，DH，DV中，求解出来的y直接对对应的DU，DH，DV进行更新
    DO kc = 1, nz_local
      jc = 1;  
        ! 这时候没有对应的次对角元素，只有对应的对角元素
        DO ic = 1, nx_local
         call DGETRS('N',NumVar,1,D_LUall(:,:,ic,jc,kc),NumVar,ipiv(:,ic,jc,kc),DU(:,ic,jc,kc),NumVar,info)
        ENDDO
      DO jc = 2, Ny - 1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            DO kk = 1, NumVar
              DU(iVar,ic,jc,kc) = DU(iVar,ic,jc,kc) - L_LUall(iVar,kk,ic,jc,kc) * DU(kk,ic,jc-1,kc);
            ENDDO
          ENDDO
          call DGETRS('N',NumVar,1,D_LUall(:,:,ic,jc,kc),NumVar,ipiv(:,ic,jc,kc),DU(:,ic,jc,kc),NumVar,info)
        ENDDO
      ENDDO
    ENDDO
    ! 其次我们需要更新对应的上边界以及对应的激波相关的物理量
    DO kc = 1, nz_local
        DO ic = 1, nx_local
          ! 首先，这里的L(NY,NY) = I, 因此对应的y部分就等于对应的f部分
          ! H部分
          DH(ic,kc) = DT0 * DH(ic,kc);
          ! V部分
          DV(ic,kc) = DV(ic,kc) - L_LUacH(ic,kc) * DH(ic,kc);
          DO kk = 1, NumVar
            DV(ic,kc) = DV(ic,kc) - L_LUacUm(kk,ic,kc) * DU(kk,ic,Ny-1,kc) &
                                & - L_LUacUs(kk,ic,kc) * DU(kk,ic,Ny  ,kc) ; 
          ENDDO
          DV(ic,kc) = DV(ic,kc) / L_LUacD(ic,kc);
        ENDDO
    ENDDO

    ! 其次求解对应的 U * DQ = y 的线性系统中的DQ，此时的求解顺序是从Ny到1
    ! 首先求解对应的最后一个点的Ny处的DV，DH和DU
    DO kc = 1, nz_local
      DO ic = 1, nx_local
        ! V部分，仅有对角元素部分，不用进行更新求解
        ! DV(ic,kc) = DV(ic,kc)
        ! H部分
        DH(ic,kc) = DH(ic,kc) + DT0 * DV(ic,kc);
        ! U部分
        DO iVar = 1, NumVar
          DU(iVar,ic,Ny,kc) = DU(iVar,ic,Ny,kc) + DUsDH(iVar,ic,kc) * DH(ic,kc) &
                                               &+ DusDV(iVar,ic,kc) * DV(ic,kc);
        ENDDO
      ENDDO
    ENDDO
    ! 然后我们从Ny-1到1求解对应的线性系统
    DO kc = 1, nz_local
      DO jc = Ny - 1, 1, -1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            DO kk = 1, NumVar
              DU(iVar,ic,jc,kc) = DU(iVar,ic,jc,kc) - U_LUall(iVar,kk,ic,jc,kc) * DU(iVar,ic,jc+1,kc);
            ENDDO
            ! DRDH和DRDV部分的影响
            DU(iVar,ic,jc,kc) = DU(iVar,ic,jc,kc) - DRDH(iVar,ic,jc,kc) * DH(ic,kc) &
                                                & - DRDV(iVar,ic,jc,kc) * DV(ic,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! 这里相当于对应的线性系统部分已经求解完成了，我们需要对应的更新对应的物理量
    ! The relative linear system has been solved, we need to update 
    ! the corresponding physical quantities
    ! 在更新之前，我们可以先判断一下对应的迭代次数
    IF(dplur_iter == dplur_sub_iter) exit DPLR_Loop
    
    ! 这里需要求解对应的DH_xi,DV_xi,DH_zeta,DV_zeta的值
    ! Here, we need to solve the corresponding DH_xi,DV_xi,DH_zeta,
    ! DV_zeta values
    call Parallel_Exchange_Surface(DH)
    call Parallel_Exchange_Surface(DV)

    call Cal_Deri_SurfDxi_5th_ce(DH_xi,DH)
    call Cal_Deri_SurfDxi_5th_ce(DV_xi,DV)
    call Cal_Deri_SurfDzeta_5th_ce(DH_zeta,DH)
    call Cal_Deri_SurfDzeta_5th_ce(DV_zeta,DV)

    ! 首先，这里更新了DU部分，所以DU部分不等于0了，对应的xi方向和zeta方向的影响需要进行考虑
    ! First, the DU part has been updated, so the DU part is not equal to 0, 
    ! and the influence of the xi and zeta directions needs to be considered 
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            ADU(iVar,ic,jc,kc) = 0.d0; ! xi方向的通量Jacobian矩阵
            CDU(iVar,ic,jc,kc) = 0.d0; ! zeta方向的通量Jacobian矩阵
            DO kk = 1, NumVar
              ADU(iVar,ic,jc,kc) = ADU(iVar,ic,jc,kc) + DFhatDU(iVar,kk,ic,jc,kc) * DU(kk,ic,jc,kc);
              CDU(iVar,ic,jc,kc) = CDU(iVar,ic,jc,kc) + DHhatDU(iVar,kk,ic,jc,kc) * DU(kk,ic,jc,kc);
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    if(npx == npx0 - 1) then
      ! DU(:,nx_local+1,:,:) =  DU(:,nx_local,:,:)
      ADU(:,nx_local+1,:,:) = ADU(:,nx_local,:,:)
    endif
    if(npx == 0       ) then
      ! DU(:,0          ,:,:) = DU(:,1       ,:,:)
      ADU(:,0          ,:,:) = ADU(:,1       ,:,:)
    endif

    if(npz == npz0 - 1) then
      ! DU(:,:,:,nz_local+1) = DU(:,:,:,nz_local)
      CDU(:,:,:,nz_local+1) = CDU(:,:,:,nz_local)
    endif
    if(npz == 0       ) then
      ! DU(:,:,:,0)       = DU(:,:,:,1)
      CDU(:,:,:,0)       = CDU(:,:,:,1)
    endif

    ! 然后我们需要考虑对应的数据传输
    CALL Parallel_Exchange_NumVar(DU)
    CALL Parallel_Exchange_NumVar(ADU)
    CALL Parallel_Exchange_NumVar(CDU)

    ! 注意这些项是放在方程右端的，因此是dR/dU * DU的部分，而不是-dR/dU * DU的部分
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1,nx_local
          CO_IM = jaco(ic,jc,kc) / (2.d0 * jaco(ic-1,jc,kc));
          CO_IP =-jaco(ic,jc,kc) / (2.d0 * jaco(ic+1,jc,kc));
          CO_KM = jaco(ic,jc,kc) / (2.d0 * jaco(ic,jc,kc-1));
          CO_KP =-jaco(ic,jc,kc) / (2.d0 * jaco(ic,jc,kc+1));
          ! 这里仅仅存储对应的增量部分
          DO iVar = 1, NumVar
            DPLUR_RHS(iVar,ic,jc,kc) = CO_IM * (ADU(iVar,ic-1,jc,kc)   + Rds_Xi(ic-1,jc,  kc  )) * DU(iVar,ic-1,jc,  kc  ) &
                                   & + CO_IP * (ADU(iVar,ic+1,jc,kc)   - Rds_Xi(ic+1,jc,  kc  )) * DU(iVar,ic+1,jc,  kc  ) &
                                   & + CO_KM * (CDU(iVar,ic,  jc,kc-1) + Rds_Zeta(ic,jc  ,kc-1)) * DU(iVar,ic,  jc  ,kc-1) &
                                   & + CO_KP * (CDU(iVar,ic,  jc,kc+1) - Rds_Zeta(ic,jc  ,kc-1)) * DU(iVar,ic,  jc  ,kc+1) &
                                   & + DRDHxi(iVar,ic,jc,kc) * DH_xi(ic,kc) &
                                   & + DRDHzeta(iVar,ic,jc,kc) * DH_zeta(ic,kc);
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO kc = 1, nz_local
      jc = Ny
        DO ic = 1, nx_local
          CO_IM = jaco(ic,jc,kc) / (2.d0 * jaco(ic-1,jc,kc));
          CO_IP =-jaco(ic,jc,kc) / (2.d0 * jaco(ic+1,jc,kc));
          CO_KM = jaco(ic,jc,kc) / (2.d0 * jaco(ic,jc,kc-1));
          CO_KP =-jaco(ic,jc,kc) / (2.d0 * jaco(ic,jc,kc+1));
          DO iVar = 1, NumVar
            DPLUR_RHS(iVar,ic,jc,kc) = CO_IM * (ADU(iVar,ic-1,jc,kc)   + Rds_Xi(ic-1,jc,  kc)) * DU(iVar,ic-1,jc,  kc  ) &
                                   & + CO_IP * (ADU(iVar,ic+1,jc,kc)   - Rds_Xi(ic+1,jc,  kc)) * DU(iVar,ic+1,jc,  kc  ) &
                                   & + CO_KM * (CDU(iVar,ic,  jc,kc-1) + Rds_Zeta(ic,jc,kc-1)) * DU(iVar,ic,  jc,  kc-1) &
                                   & + CO_KP * (CDU(iVar,ic,  jc,kc+1) - Rds_Zeta(ic,jc,kc+1)) * DU(iVar,ic,  jc,  kc+1)
          ENDDO
        ENDDO
    ENDDO

  ENDDO DPLR_Loop
  
  ! 考虑到对应的壁面处的速度为0，对应增量对于壁面上三个动量方程的修正应该一直为0
  jc = 1;
  DO kc = 1, nz_local
    DO ic = 1, nx_local
      DO iVar = 2,NumVar-1
        DU(iVar,ic,jc,kc) = 0.0_rk;
      ENDDO
    ENDDO
  ENDDO

   dUcons = DU;
  dShockH = DH;
  dShockV = DV;

end subroutine DPLUR_Line_Solver

subroutine Calculate_GeneralJacobian
  use SF_Constant,     only: ik,rk,NumVar,overLAP
  use SF_CFD_Global,   only: nx_local,Ny,nz_local

  ! arrays and subs
  use SF_CFD_Global,   only: Jaco, DJDH, dUcons0, invF, visF, invG, visG, invH, visH,&
                        &  Ucons0, F_IM, G_IM, H_IM, DF_IMDxi, DF_IMDeta, DF_IMDzeta,&
                        &  DG_IMDxi, DG_IMDeta, DG_IMDzeta, DH_IMDxi, DH_IMDeta, DH_IMDzeta,&
                        &  DU_Deta,DJDHxi, DJDHzeta, DRDH, DRDHxi, DRDHzeta, DRDV, DetatJDV,&
                        &  DxixJDH, DxiyJDH, DxizJDH, DetaxJDH, DetayJDH, DetazJDH,&
                        &  DzetaxJDH, DzetayJDH, DzetazJDH, DetatJDH, DetaxJDHxi, DzetaxJDHxi,&
                        &  DetatJDHxi,DetayJDHxi,DzetayJDHxi,DetazJDHxi,DzetazJDHxi,DxixJDHzeta,&
                        &  DetaxJDHzeta,DxiyJDHzeta,DetayJDHzeta,DxizJDHzeta,DetazJDHzeta,DetatJDHzeta,&
                        &  ADUSDH,BDUSDH,CDUSDH,ADUSDHxi,BDUSDHxi,CDUSDHxi,ADUSDHzeta,BDUSDHzeta,CDUSDHzeta,&
                        &  ADUSDH_xi,CDUSDH_xi,ADUSDH_zeta,CDUSDH_zeta,ADUSDHxi_xi,CDUSDHxi_xi,&
                        &  ADUSDHxi_zeta,CDUSDHxi_zeta,ADUSDHzeta_xi,CDUSDHzeta_xi,ADUSDHzeta_zeta,CDUSDHzeta_zeta,&
                        &  ADUSDV,BDUSDV,CDUSDV,ADUSDV_xi,CDUSDV_xi,ADUSDV_zeta,CDUSDV_zeta
  use FD5_Order,       only: BC1,BC2,BC3,CO_FD_Center
  
  ! subroutines and functions
  use MPI_GLOBAL,      only: Parallel_Exchange_NumVar,Parallel_Exchange_Surface
  use FD5_Order,       only: Cal_Deri_Dxi_5th_ce_numvar, Cal_Deri_Deta_5th_ce_numvar, Cal_Deri_Dzeta_5th_ce_per_numvar,&
                            & Cal_Deri_SurfDxi_5th_ce, Cal_Deri_SurfDzeta_5th_ce
  
  implicit none


  integer( kind = ik ):: ic,jc,kc,iVar,kk
  ! Calculate dR/dU,dR/dH,dR/dV,dR/dHxi,dR/dHzeta.dR/dVxi,dR/dVzeta
  ! 实际上这里计算的就是-dR/dH,-dR/dHxi,-dR/dHzeta,-dR/dV,-dR/dVxi,-dR/dVzeta
  ! 需要额外注意的是边界条件的问题，理论上来说我们并不需要计算边界处的R和相关导数项
  ! 因为这是由边界条件单独决定的。
  ! ============================ dR/dH ============================
  ! DRDHxi,DRDHzeta
  ! R = J * (-invFlux + visFlux + U dinvJdt)
  !-R = J * (invFlux - visFlux - U dinvJdt)
  ! -DRDH_1 = - DJDH * R / J 
   do kc = 1-overLAP, nz_local+overLAP
    do jc = 1,Ny
     do ic = 1-overLAP, nx_local+overLAP
      do iVar = 1, NumVar
        ! - DRDH
           DRDH(iVar,ic,jc,kc) =     -DJDH(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
        ! - DRDHxi  
         DRDHxi(iVar,ic,jc,kc) =   -DJDHxi(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
        ! - DRDHzeta  
       DRDHzeta(iVar,ic,jc,kc) = -DJDHzeta(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
       ! As the J is not a function of V, the following terms are zero
       !     DRDV(iVar,ic,jc,kc) =     DJDV(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
       !   DRDVxi(iVar,ic,jc,kc) =   DJDVxi(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
       ! DRDVzeta(iVar,ic,jc,kc) = DJDVzeta(ic,jc,kc) * dUcons0(iVar,ic,jc,kc) / Jaco(ic,jc,kc)
      enddo
     enddo
    enddo
   enddo
   
   F_IM = invF - visF   
   G_IM = invG - visG
   H_IM = invH - visH
   
   ! Calculate the
   ! DF_IMDxi, DF_IMDeta, DF_IMDzeta
   ! DG_IMDxi, DG_IMDeta, DG_IMDzeta
   ! DH_IMDxi, DH_IMDeta, DH_IMDzeta
   ! DU_Deta 

   call Parallel_Exchange_NumVar(F_IM)
   call Parallel_Exchange_NumVar(G_IM)
   call Parallel_Exchange_NumVar(H_IM)
   call Parallel_Exchange_NumVar(Ucons0)

   call Cal_Deri_Dxi_5th_ce_numvar(      DF_IMDxi,  F_IM)
   call Cal_Deri_Deta_5th_ce_numvar(     DF_IMDeta, F_IM)
   call Cal_Deri_Dzeta_5th_ce_per_numvar(DF_IMDzeta,F_IM)

   call Cal_Deri_Dxi_5th_ce_numvar(      DG_IMDxi,  G_IM)
   call Cal_Deri_Deta_5th_ce_numvar(     DG_IMDeta, G_IM)
   call Cal_Deri_Dzeta_5th_ce_per_numvar(DG_IMDzeta,G_IM)
   
   call Cal_Deri_Dxi_5th_ce_numvar(      DH_IMDxi,  H_IM)
   call Cal_Deri_Deta_5th_ce_numvar(     DH_IMDeta, H_IM)
   call Cal_Deri_Dzeta_5th_ce_per_numvar(DH_IMDzeta,H_IM)
   
   call Cal_Deri_Deta_5th_ce_numvar(     DU_Deta, Ucons0)
   
   do kc = 1, nz_local
    do jc = 1, Ny
     do ic = 1, nx_local
      do iVar = 1, NumVar
         DRDH(iVar,ic,jc,kc) = DRDH(iVar,ic,jc,kc) + & 
        (    DxixJDH(ic,jc,kc) *   DF_IMDxi(iVar,ic,jc,kc) &
        &+  DetaxJDH(ic,jc,kc) *  DF_IMDeta(iVar,ic,jc,kc) &
        &+ DzetaxJDH(ic,jc,kc) * DF_IMDzeta(iVar,ic,jc,kc) &
        &+   DxiyJDH(ic,jc,kc) *   DG_IMDxi(iVar,ic,jc,kc) &
        &+  DetayJDH(ic,jc,kc) *  DG_IMDeta(iVar,ic,jc,kc) &
        &+ DzetayJDH(ic,jc,kc) * DG_IMDzeta(iVar,ic,jc,kc) &
        &+   DxizJDH(ic,jc,kc) *   DH_IMDxi(iVar,ic,jc,kc) &
        &+  DetazJDH(ic,jc,kc) *  DH_IMDeta(iVar,ic,jc,kc) &
        &+ DzetazJDH(ic,jc,kc) * DH_IMDzeta(iVar,ic,jc,kc) &
        &+  DetatJDH(ic,jc,kc) *    DU_Deta(iVar,ic,jc,kc)) * Jaco(ic,jc,kc)
      enddo
     enddo
    enddo
   enddo
   
   ! Surface Boundary modifications
   jc = 1;
   DO kc = 1, nz_local
     DO ic = 1, nx_local
       DO iVar = 2, NumVar ! 考虑到连续性方程的作用，壁面处的速度为0
         DRDH(iVar,ic,jc,kc) = 0.d0; 
       ENDDO
     ENDDO
   ENDDO

  ! ============================ dR/dHxi ============================
  ! DRDHxi
  ! DetaxJDHxi,  DetayJDHxi,  DetazJDHxi,  DetatJDHxi
  ! DzetaxJDHxi, DzetayJDHxi, DzetazJDHxi
  ! =================================================================
   do kc = 1, nz_local
    do jc = 2, Ny   ! Here, we do not need to calculate the surface
     do ic = 1, nx_local
      do iVar = 1, NumVar
         DRDHxi(iVar,ic,jc,kc) = DRDHxi(iVar,ic,jc,kc) + & 
        (   DetaxJDHxi(ic,jc,kc) *  DF_IMDeta(iVar,ic,jc,kc) &
        &+ DzetaxJDHxi(ic,jc,kc) * DF_IMDzeta(iVar,ic,jc,kc) &
        &+  DetayJDHxi(ic,jc,kc) *  DG_IMDeta(iVar,ic,jc,kc) &
        &+ DzetayJDHxi(ic,jc,kc) * DG_IMDzeta(iVar,ic,jc,kc) &
        &+  DetazJDHxi(ic,jc,kc) *  DH_IMDeta(iVar,ic,jc,kc) &
        &+ DzetazJDHxi(ic,jc,kc) * DH_IMDzeta(iVar,ic,jc,kc) &
        &+  DetatJDHxi(ic,jc,kc) *    DU_Deta(iVar,ic,jc,kc)) * Jaco(ic,jc,kc)
      enddo
     enddo
    enddo
   enddo
  
  ! ============================ dR/dHzeta ============================
  ! DetaxJDHzeta,DetayJDHzeta,DetazJDHzeta,DetatJDHzeta
  ! DxixJDHzeta, DxiyJDHzeta, DxizJDHzeta
  ! =================================================================
   do kc = 1, nz_local
    do jc = 2, Ny  ! Here, we do not need to calculate the surface
     do ic = 1, nx_local
      do iVar = 1, NumVar
         DRDHzeta(iVar,ic,jc,kc) = DRDHzeta(iVar,ic,jc,kc) + & 
        (    DxixJDHzeta(ic,jc,kc) *  DF_IMDxi(iVar,ic,jc,kc) &
        &+  DetaxJDHzeta(ic,jc,kc) * DF_IMDeta(iVar,ic,jc,kc) &
        &+   DxiyJDHzeta(ic,jc,kc) *  DG_IMDxi(iVar,ic,jc,kc) &
        &+  DetayJDHzeta(ic,jc,kc) * DG_IMDeta(iVar,ic,jc,kc) &
        &+   DxizJDHzeta(ic,jc,kc) *  DH_IMDxi(iVar,ic,jc,kc) &
        &+  DetazJDHzeta(ic,jc,kc) * DH_IMDeta(iVar,ic,jc,kc) &
        &+  DetatJDHzeta(ic,jc,kc) *   DU_Deta(iVar,ic,jc,kc)) * Jaco(ic,jc,kc)
      enddo
     enddo
    enddo
   enddo

  ! ============================ dR/dV ============================
   do kc = 1, nz_local
    do jc = 2, Ny
     do ic = 1, nx_local
      do iVar = 1, NumVar
        DRDV(iVar,ic,jc,kc) = - DetatJDV(ic,jc,kc) * DU_Deta(iVar,ic,jc,kc) * Jaco(ic,jc,kc)
      enddo
     enddo
    enddo
   enddo

  ! ============================ dR/dVxi ============================
  ! =========================== dR/dVzeta ===========================

  ! Rankine-Hugoniot condition can be written as
  ! Us = Us(U_inf, ShockV, ShockH, ShockHxi, ShockHzeta)
    call Calculate_Us_ShockDeri
  ! It gives us the
  ! DUsDH,DUsDHdxi,DUsDHdzeta,DUsDV
   
  ! Where the shock acceleration is dependended on more variables
  
   call Calculate_ShkAcDeri
  ! It gives us the
  ! ShockAcdH,ShockAcdHdxi,ShockAcdHdzeta,ShockAcdV,ShockAcdVdxi,ShockAcdVdzeta
    
  ! 这里有一个问题，从明面上看 R似乎仅仅是H,Hxi,Hzeta,V的函数，但实际上在激波边界处
  ! 由于R-H条件，Us可以直接由对应的条件决定。这里的Us是一个函数
  ! 最复杂的是对应的激波加速度，这个加速度不光依赖于H,Hxi,Hzeta,V，还依赖于激波的速度
  ! 对应的导数Vxi,Vzeta，因此会引入对应的导数
  ! 也就是说在激波边界有 dR/dUs*dUs/dH类似这样的项

   jc = Ny
   DO kc = 1, nz_local
     DO ic = 1, nx_local
       DO iVar = 1, NumVar
         ADUSDH(iVar,ic,kc) = 0.0_rk;
         BDUSDH(iVar,ic,kc) = 0.0_rk;
         CDUSDH(iVar,ic,kc) = 0.0_rk;
         
         ADUSDHxi(iVar,ic,kc) = 0.0_rk;
         BDUSDHxi(iVar,ic,kc) = 0.0_rk;
         CDUSDHxi(iVar,ic,kc) = 0.0_rk;

         ADUSDHzeta(iVar,ic,kc) = 0.0_rk;  
         BDUSDHzeta(iVar,ic,kc) = 0.0_rk;
         CDUSDHzeta(iVar,ic,kc) = 0.0_rk;

         ADUSDV(iVar,ic,kc) = 0.0_rk;
         BDUSDV(iVar,ic,kc) = 0.0_rk;
         CDUSDV(iVar,ic,kc) = 0.0_rk;

         DO kk = 1, NumVar
           ! dUs/dH
           ADUSDH(iVar,ic,kc) = ADUSDH(iVar,ic,kc) + DFhatDU(iVar,kk,ic,Ny,kc) * DUsDH(kk,ic,kc)
           BDUSDH(iVar,ic,kc) = BDUSDH(iVar,ic,kc) + DGhatDU(iVar,kk,ic,Ny,kc) * DUsDH(kk,ic,kc)
           CDUSDH(iVar,ic,kc) = CDUSDH(iVar,ic,kc) + DHhatDU(iVar,kk,ic,Ny,kc) * DUsDH(kk,ic,kc)
           ! dUs/dHxi
           ADUSDHxi(iVar,ic,kc) = ADUSDHxi(iVar,ic,kc) + DFhatDU(iVar,kk,ic,Ny,kc) * DUsDHdxi(kk,ic,kc)
           BDUSDHxi(iVar,ic,kc) = BDUSDHxi(iVar,ic,kc) + DGhatDU(iVar,kk,ic,Ny,kc) * DUsDHdxi(kk,ic,kc)
           CDUSDHxi(iVar,ic,kc) = CDUSDHxi(iVar,ic,kc) + DHhatDU(iVar,kk,ic,Ny,kc) * DUsDHdxi(kk,ic,kc)           
           ! dUs/dHzeta
           ADUSDHzeta(iVar,ic,kc) = ADUSDHzeta(iVar,ic,kc) + DFhatDU(iVar,kk,ic,Ny,kc) * DUsDHdzeta(kk,ic,kc)
           BDUSDHzeta(iVar,ic,kc) = BDUSDHzeta(iVar,ic,kc) + DGhatDU(iVar,kk,ic,Ny,kc) * DUsDHdzeta(kk,ic,kc)
           CDUSDHzeta(iVar,ic,kc) = CDUSDHzeta(iVar,ic,kc) + DHhatDU(iVar,kk,ic,Ny,kc) * DUsDHdzeta(kk,ic,kc)
           ! dUs/dV
           ADUSDV(iVar,ic,kc) = ADUSDV(iVar,ic,kc) + DFhatDU(iVar,kk,ic,Ny,kc) * DUsDV(kk,ic,kc)
           BDUSDV(iVar,ic,kc) = BDUSDV(iVar,ic,kc) + DGhatDU(iVar,kk,ic,Ny,kc) * DUsDV(kk,ic,kc)
           CDUSDV(iVar,ic,kc) = CDUSDV(iVar,ic,kc) + DHhatDU(iVar,kk,ic,Ny,kc) * DUsDV(kk,ic,kc)
         ENDDO 
       ENDDO
     ENDDO
   ENDDO
   DO iVar = 1, NumVar
     call Parallel_Exchange_Surface(ADUSDH(iVar,:,:))
     call Parallel_Exchange_Surface(CDUSDH(iVar,:,:))
     !
     call Cal_Deri_SurfDxi_5th_ce(  ADUSDH_xi(iVar,:,:),  ADUSDH(iVar,:,:))
     call Cal_Deri_SurfDzeta_5th_ce(CDUSDH_zeta(iVar,:,:),CDUSDH(iVar,:,:))

     call Parallel_Exchange_Surface(ADUSDHxi(iVar,:,:))
     call Parallel_Exchange_Surface(CDUSDHxi(iVar,:,:))

     call Cal_Deri_SurfDxi_5th_ce(  ADUSDHxi_xi(iVar,:,:),  ADUSDHxi(iVar,:,:))
     call Cal_Deri_SurfDzeta_5th_ce(CDUSDHxi_zeta(iVar,:,:),CDUSDHxi(iVar,:,:))

     call Parallel_Exchange_Surface(ADUSDHzeta(iVar,:,:))
     call Parallel_Exchange_Surface(CDUSDHzeta(iVar,:,:))

     call Cal_Deri_SurfDxi_5th_ce(  ADUSDHzeta_xi(iVar,:,:),  ADUSDHzeta(iVar,:,:))
     call Cal_Deri_SurfDzeta_5th_ce(CDUSDHzeta_zeta(iVar,:,:),CDUSDHzeta(iVar,:,:))
     
     call Parallel_Exchange_Surface(ADUSDV(iVar,:,:))
     call Parallel_Exchange_Surface(CDUSDV(iVar,:,:))

      call Cal_Deri_SurfDxi_5th_ce(  ADUSDV_xi(iVar,:,:),  ADUSDV(iVar,:,:))
      call Cal_Deri_SurfDzeta_5th_ce(CDUSDV_zeta(iVar,:,:),CDUSDV(iVar,:,:))
   ENDDO

   ! xi方向和zeta方向上的修正
   jc = Ny
   DO kc = 1, nz_local
     DO ic = 1, nx_local
       DO iVar = 1, NumVar
         DRDH(iVar,ic,Ny,kc) = DRDH(iVar,ic,Ny,kc) + &
                        & ADUSDH_xi(iVar,ic,kc) +     CDUSDH_zeta(iVar,ic,kc)
         DRDHxi(iVar,ic,Ny,kc) = DRDHxi(iVar,ic,Ny,kc) + &
                      & ADUSDHxi_xi(iVar,ic,kc) +   CDUSDHxi_zeta(iVar,ic,kc)
         DRDHzeta(iVar,ic,Ny,kc) = DRDHzeta(iVar,ic,Ny,kc) + &
                    & ADUSDHzeta_xi(iVar,ic,kc) + CDUSDHzeta_zeta(iVar,ic,kc)
         DRDV(iVar,ic,Ny,kc) = DRDV(iVar,ic,Ny,kc) + &
                    & ADUSDV_xi(iVar,ic,kc) +     CDUSDV_zeta(iVar,ic,kc)
       ENDDO
     ENDDO
   ENDDO

   ! Ny
   DO kc = 1, nz_local
     DO ic = 1, nx_local
       DO iVar = 1, NumVar
        ! 考虑边界条件上的差分格式的影响 
         DRDH(iVar,ic,Ny  ,kc) = DRDH(iVar,ic,Ny  ,kc) - BC1(1) * BDUSDH(iVar,ic,kc);
         DRDH(iVar,ic,Ny-1,kc) = DRDH(iVar,ic,Ny-1,kc) - BC2(1) * BDUSDH(iVar,ic,kc);
         DRDH(iVar,ic,Ny-2,kc) = DRDH(iVar,ic,Ny-2,kc) - BC3(1) * BDUSDH(iVar,ic,kc);
         DRDH(iVar,ic,Ny-3,kc) = DRDH(iVar,ic,Ny-3,kc) - CO_FD_Center(3) * BDUSDH(iVar,ic,kc);

         DRDHxi(iVar,ic,Ny  ,kc) = DRDHxi(iVar,ic,Ny  ,kc) - BC1(1) * BDUSDHxi(iVar,ic,kc);
         DRDHxi(iVar,ic,Ny-1,kc) = DRDHxi(iVar,ic,Ny-1,kc) - BC2(1) * BDUSDHxi(iVar,ic,kc);
         DRDHxi(iVar,ic,Ny-2,kc) = DRDHxi(iVar,ic,Ny-2,kc) - BC3(1) * BDUSDHxi(iVar,ic,kc);
         DRDHxi(iVar,ic,Ny-3,kc) = DRDHxi(iVar,ic,Ny-3,kc) - CO_FD_Center(3) * BDUSDHxi(iVar,ic,kc);
         
         DRDHzeta(iVar,ic,Ny  ,kc) = DRDHzeta(iVar,ic,Ny  ,kc) - BC1(1) * BDUSDHzeta(iVar,ic,kc);
         DRDHzeta(iVar,ic,Ny-1,kc) = DRDHzeta(iVar,ic,Ny-1,kc) - BC2(1) * BDUSDHzeta(iVar,ic,kc);
         DRDHzeta(iVar,ic,Ny-2,kc) = DRDHzeta(iVar,ic,Ny-2,kc) - BC3(1) * BDUSDHzeta(iVar,ic,kc);
         DRDHzeta(iVar,ic,Ny-3,kc) = DRDHzeta(iVar,ic,Ny-3,kc) - CO_FD_Center(3) * BDUSDHzeta(iVar,ic,kc);

         DRDV(iVar,ic,Ny  ,kc) = DRDV(iVar,ic,Ny  ,kc) - BC1(1) * BDUSDV(iVar,ic,kc);
         DRDV(iVar,ic,Ny-1,kc) = DRDV(iVar,ic,Ny-1,kc) - BC2(1) * BDUSDV(iVar,ic,kc);
         DRDV(iVar,ic,Ny-2,kc) = DRDV(iVar,ic,Ny-2,kc) - BC3(1) * BDUSDV(iVar,ic,kc);
         DRDV(iVar,ic,Ny-3,kc) = DRDV(iVar,ic,Ny-3,kc) - CO_FD_Center(3) * BDUSDV(iVar,ic,kc);
       ENDDO
     ENDDO
   ENDDO

   ! Shock Acceleration modifications
   DO kc = 1, nz_local
     DO ic = 1, nx_local
       DO iVar = 1, NumVar
             ShockAcdH(ic,kc) =     ShockAcdH(ic,kc) + CO_AC_UT(iVar,ic,kc) *     DRDH(iVar,ic,Ny,kc)
          ShockAcdHdxi(ic,kc) =  ShockAcdHdxi(ic,kc) + CO_AC_UT(iVar,ic,kc) *   DRDHxi(iVar,ic,Ny,kc)
        ShockAcdHdzeta(ic,kc) =ShockAcdHdzeta(ic,kc) + CO_AC_UT(iVar,ic,kc) * DRDHzeta(iVar,ic,Ny,kc)
             ShockAcdV(ic,kc) =     ShockAcdV(ic,kc) + CO_AC_UT(iVar,ic,kc) *     DRDV(iVar,ic,Ny,kc)
       ENDDO
     ENDDO
   ENDDO

end subroutine Calculate_GeneralJacobian

subroutine DPLUR_Solver
  use SFitting, only:Calculate_ShockAC3D_Implicit
  implicit none
  ! This subroutine is used to calculate the whole process of a DPLUR iteration
  ! Although we solve the normal joint equation in DPLUR, the xi and zeta directions are weakly coupled
  ! Therefore, we introduce 5 iterations to mainly update the right-hand side of the equation
  ! We call these 3-5 iterations together as a DPLUR iteration

     UconsOld = Ucons0;
    ShockH_old= ShockH;
    ShockV_old= ShockV;

  ! The first step for a DPLUR iteration is to calculate the flux 
  call Calculate_Flux                 ! Calculate the dUcons0
  call Calculate_ShockAC3D_Implicit   ! Calculate the ShockAC
  
  ! The second step is to form the matrix and solve the tridiagonal system
  call DPLUR_Line_Solver

  ! The third step is to update the variables
  ! Here we highly recommend to use the direct + in fortran without any loop
  ! The SIMD feature of the CPU can be fully utilized

  Ucons0 = UconsOld + dUcons;
  ShockH = ShockH_old + dShockH;
  ShockV = ShockV_old + dShockV;

  ! Here we recalculate the Jacobian matrix again
  if(ModelType == 1) then
    call Calculate_Jaco_Implicit  
  else    
    call Singular_Calculate_Jaco_Implicit
  endif
  
  
  call Update_Variables
   
end subroutine DPLUR_Solver

END Module DPLUR