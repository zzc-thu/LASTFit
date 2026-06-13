!SUBROUTINE Calculate_Flux_Weno5th
!  use SF_Constant,           only: ik,rk,NumVar,LaxFriSmall
!  use SF_CFD_Global,         only: nx_local,Ny,nz_local,Gamma,Pr_Ref,Mach_Ref,Re_Ref,&
!                                   Csthlnd_Ref,Cp_Ref
!
!  ! matrix and some arrays
!  use SF_CFD_Global,         only: Rho,U,V,W,P,T,Ucons0,Cs,BigU_Xi,BigU_Eta,BigU_Zeta,&
!                                   dxidx,dxidy,dxidz,dxidt,detadx,detady,detadz,detadt,&
!                                   dzetadx,dzetady,dzetadz,dzetadt,invF,invG,invH,&
!                                   Udxi,Udeta,Udzeta,Vdxi,Vdeta,Vdzeta,Wdxi,Wdeta,Wdzeta,&
!                                   Tdxi,Tdeta,Tdzeta,Mu,Udx,Udy,Udz,Vdx,Vdy,Vdz,Wdx,Wdy,Wdz,&
!                                   Tdx,Tdy,Tdz,visF,visG,visH,visFhat,visGhat,visHhat,Jaco,invJacodt,&
!                                   nablaxi,nablaeta,nablazeta,invFhat,invGhat,invHhat,Fp,Gp,Hp,&
!                                   Fm,Gm,Hm,Fpdxi,Fmdxi,Gpdeta,Gmdeta,Hpdzeta,Hmdzeta,visFcdxi,&
!                                   visGcdeta,visHcdzeta,invFlux,visFlux,dUcons0
!
!  ! subroutine and function declaration
!  use MPI_GLOBAL,            only: Parallel_Exchange,Parallel_Exchange_NumVar
!  use FD5_Order,             only: Cal_Deri_Dxi_5th_ce,Cal_Deri_Deta_5th_ce,Cal_Deri_Dzeta_5th_ce_per,& ! Central
!              Cal_Deri_Dxi_5th_ce_numvar,Cal_Deri_Deta_5th_ce_numvar,Cal_Deri_Dzeta_5th_ce_per_numvar,& ! Central
!  use FD5_Order_Weno,        only: &
!              Cal_Deri_Dxi_weno5th_up_numvar,Cal_Deri_Deta_weno5th_up_numvar,Cal_Deri_Dzeta_weno5th_up_per_numvar,& ! Upwind
!              Cal_Deri_Dxi_weno5th_do_numvar,Cal_Deri_Deta_weno5th_do_numvar,Cal_Deri_Dzeta_weno5th_do_per_numvar ! Downwind
!
!  implicit none
!  ! Local variables
!  integer( kind = ik ):: ic,jc,kc,iVar
!  real( kind = rk ):: Div_Vel,Tauxx,Tauyy,Tauzz,Tauxy,Tauxz,Tauyz
!  real( kind = rk ):: Cs_Xi,Cs_Eta,Cs_Zeta
!  real( kind = rk ):: Sigma_Xi,Sigma_Eta,Sigma_Zeta
!
!  real( kind = rk ):: local_temp
!
!  ! Get the inner variables
!  do kc = 1, nz_local
!   do jc = 1, Ny
!    do ic = 1, nx_local
!      ! conservative variables to primitive variables
!      call CVtoU_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
!      ! Calculate the temperature
!      T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc)
!      ! Calculate the sound speed
!      Cs(ic,jc,kc)= sqrt(Gamma * P(ic,jc,kc) / Rho(ic,jc,kc))
!      ! Calculate the characteristic speed
!        BigU_Xi(ic,jc,kc) = U(ic,jc,kc)*  dxidx(ic,jc,kc) &
!                         &+ V(ic,jc,kc)*  dxidy(ic,jc,kc) &
!                         &+ W(ic,jc,kc)*  dxidz(ic,jc,kc) + dxidt(ic,jc,kc)
!
!       BigU_Eta(ic,jc,kc) = U(ic,jc,kc)* detadx(ic,jc,kc) &
!                         &+ V(ic,jc,kc)* detady(ic,jc,kc) &
!                         &+ W(ic,jc,kc)* detadz(ic,jc,kc) + detadt(ic,jc,kc)
!
!      BigU_Zeta(ic,jc,kc) = U(ic,jc,kc)*dzetadx(ic,jc,kc) &
!                         &+ V(ic,jc,kc)*dzetady(ic,jc,kc) &
!                         &+ W(ic,jc,kc)*dzetadz(ic,jc,kc) + dzetadt(ic,jc,kc)
!
!      ! Calculate the flux
!      call UtoFlux_CPG3D(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),&
!                         invF(1:NumVar,ic,jc,kc),invG(1:NumVar,ic,jc,kc),invH(1:NumVar,ic,jc,kc))
!    enddo
!   enddo
!  enddo
!  
!  ! Parallel transfer
!  call Parallel_Exchange(U)
!  call Parallel_Exchange(V)
!  call Parallel_Exchange(W)
!  call Parallel_Exchange(T)
!
!  ! Calculate the viscous Flux by using the central difference
!  call Cal_Deri_Dxi_5th_ce(Udxi,U)
!  call Cal_Deri_Dxi_5th_ce(Vdxi,V)
!  call Cal_Deri_Dxi_5th_ce(Wdxi,W)
!  call Cal_Deri_Dxi_5th_ce(Tdxi,T)
!
!  call Cal_Deri_Deta_5th_ce(Udeta,U)
!  call Cal_Deri_Deta_5th_ce(Vdeta,V)
!  call Cal_Deri_Deta_5th_ce(Wdeta,W)
!  call Cal_Deri_Deta_5th_ce(Tdeta,T)
!  
!  ! Here we use the periodic boundary condition form
!  call Cal_Deri_Dzeta_5th_ce_per(Udzeta,U)
!  call Cal_Deri_Dzeta_5th_ce_per(Vdzeta,V)
!  call Cal_Deri_Dzeta_5th_ce_per(Wdzeta,W)
!  call Cal_Deri_Dzeta_5th_ce_per(Tdzeta,T)
!  
!  ! Cp_Ref = 1.0_rk / ((Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
!  local_temp = (Pr_Ref * (Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
!  ! Calculate the viscous stress
!  do kc = 1,nz_local
!   do jc = 1,Ny
!    do ic = 1,nx_local
!      ! The viscosity coefficient has been divided by Re_Ref already
!      ! T^1.5 * (1 + Csthlnd_Ref) / (T + Csthlnd_Ref) / Re_Ref 
!      Mu(ic,jc,kc) = (T(ic,jc,kc)**1.5_rk) * ( 1.0_rk + Csthlnd_Ref )/ (T(ic,jc,kc) + Csthlnd_Ref) 
!      Mu(ic,jc,kc) = Mu(ic,jc,kc) / Re_Ref;
!      ! derivatives with respect to the x direction
!      Udx(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                       Udeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                      Udzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!
!      Vdx(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                       Vdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                      Vdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);      
!
!      Wdx(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                       Wdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                      Wdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!
!      Tdx(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                       Tdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                      Tdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!
!      ! derivatives with respect to the y direction
!      Udy(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                       Udeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                      Udzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!
!      Vdy(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                       Vdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                      Vdzeta(ic,jc,kc)*dzetady(ic,jc,kc);      
!
!      Wdy(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                       Wdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                      Wdzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!
!      Tdy(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                       Tdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                      Tdzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!      
!      ! derivatives with respect to the z direction
!      Udz(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                       Udeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                      Udzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!
!      Vdz(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                       Vdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                      Vdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);      
!
!      Wdz(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                       Wdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                      Wdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!
!      Tdz(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                       Tdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                      Tdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!                      
!      Div_Vel = Udx(ic,jc,kc) + Vdy(ic,jc,kc) + Wdz(ic,jc,kc);
!
!      Tauxx = 2.0_rk * Mu(ic,jc,kc) * (Udx(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!      Tauyy = 2.0_rk * Mu(ic,jc,kc) * (Vdy(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!      Tauzz = 2.0_rk * Mu(ic,jc,kc) * (Wdz(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!      Tauxy = Mu(ic,jc,kc) * (Udy(ic,jc,kc) + Vdx(ic,jc,kc));
!      Tauxz = Mu(ic,jc,kc) * (Udz(ic,jc,kc) + Wdx(ic,jc,kc));
!      Tauyz = Mu(ic,jc,kc) * (Vdz(ic,jc,kc) + Wdy(ic,jc,kc));
!
!      ! Calculate the viscous flux in the x-,y-,z- directions
!      visF(1,ic,jc,kc) = 0.0_rk;
!      visF(2,ic,jc,kc) = Tauxx;
!      visF(3,ic,jc,kc) = Tauxy;
!      visF(4,ic,jc,kc) = Tauxz;
!      visF(5,ic,jc,kc) = U(ic,jc,kc) * Tauxx &
!                      &+ V(ic,jc,kc) * Tauxy &
!                      &+ W(ic,jc,kc) * Tauxz &
!                      &+ Mu(ic,jc,kc) * Tdx(ic,jc,kc) / local_temp;
!      
!      visG(1,ic,jc,kc) = 0.0_rk;
!      visG(2,ic,jc,kc) = Tauxy;
!      visG(3,ic,jc,kc) = Tauyy;
!      visG(4,ic,jc,kc) = Tauyz;
!      visG(5,ic,jc,kc) = U(ic,jc,kc) * Tauxy &
!                      &+ V(ic,jc,kc) * Tauyy &
!                      &+ W(ic,jc,kc) * Tauyz &
!                      &+ Mu(ic,jc,kc) * Tdy(ic,jc,kc) / local_temp;
!      
!      visH(1,ic,jc,kc) = 0.0_rk;
!      visH(2,ic,jc,kc) = Tauxz;
!      visH(3,ic,jc,kc) = Tauyz;
!      visH(4,ic,jc,kc) = Tauzz;
!      visH(5,ic,jc,kc) = U(ic,jc,kc) * Tauxz &
!                      &+ V(ic,jc,kc) * Tauyz &
!                      &+ W(ic,jc,kc) * Tauzz &
!                      &+ Mu(ic,jc,kc) * Tdz(ic,jc,kc) / local_temp;
!                      
!      do iVar = 2,NumVar
!        visFhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                                 visG(iVar,ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                                 visH(iVar,ic,jc,kc)*  dxidz(ic,jc,kc))/Jaco(ic,jc,kc);
!
!        visGhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)* detadx(ic,jc,kc) +&
!                                 visG(iVar,ic,jc,kc)* detady(ic,jc,kc) +&
!                                 visH(iVar,ic,jc,kc)* detadz(ic,jc,kc))/Jaco(ic,jc,kc);
!
!        visHhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*dzetadx(ic,jc,kc) +&
!                                 visG(iVar,ic,jc,kc)*dzetady(ic,jc,kc) +&
!                                 visH(iVar,ic,jc,kc)*dzetadz(ic,jc,kc))/Jaco(ic,jc,kc);
!      enddo
!    enddo
!   enddo
!  enddo
!
!  call Parallel_Exchange_NumVar(visFhat)
!  call Parallel_Exchange_NumVar(visGhat)
!  call Parallel_Exchange_NumVar(visHhat)
!  call Cal_Deri_Dxi_5th_ce_numvar(visFcdxi,        visFhat)
!  call Cal_Deri_Deta_5th_ce_numvar(visGcdeta,      visGhat)
!  call Cal_Deri_Dzeta_5th_ce_per_numvar(visHcdzeta,visHhat)
!
!  visFlux = visFcdxi + visGcdeta + visHcdzeta;
!
!  ! Calculate the inviscid flux
!  do kc = 1,nz_local
!   do jc = 1,Ny
!    do ic = 1,nx_local
!      Cs_Xi  = Cs(ic,jc,kc) *   nablaxi(ic,jc,kc);
!      Cs_Eta = Cs(ic,jc,kc) *  nablaeta(ic,jc,kc);
!      Cs_Zeta= Cs(ic,jc,kc) * nablazeta(ic,jc,kc);
!
!      Sigma_Xi  =(sqrt(  BigU_Xi(ic,jc,kc)*  BigU_Xi(ic,jc,kc) &
!                      &+ LaxFriSmall * LaxFriSmall * Cs_Xi  *  Cs_Xi) +  Cs_Xi)/jaco(ic,jc,kc);
!
!      Sigma_Eta =(sqrt( BigU_Eta(ic,jc,kc)* BigU_Eta(ic,jc,kc) &
!                      &+ LaxFriSmall * LaxFriSmall * Cs_Eta * Cs_Eta) + Cs_Eta)/jaco(ic,jc,kc);
!
!      Sigma_Zeta=(sqrt(BigU_Zeta(ic,jc,kc)*BigU_Zeta(ic,jc,kc) &
!                      &+ LaxFriSmall * LaxFriSmall * Cs_Zeta*Cs_Zeta) +Cs_Zeta)/jaco(ic,jc,kc);
!      do iVar = 1,NumVar
!        invFhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *   dxidx(ic,jc,kc) +&
!                                 invG(iVar,ic,jc,kc) *   dxidy(ic,jc,kc) +&
!                                 invH(iVar,ic,jc,kc) *   dxidz(ic,jc,kc) +&
!                               Ucons0(iVar,ic,jc,kc) *   dxidt(ic,jc,kc)) / Jaco(ic,jc,kc);
!        
!        invGhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *  detadx(ic,jc,kc) +&
!                                 invG(iVar,ic,jc,kc) *  detady(ic,jc,kc) +&
!                                 invH(iVar,ic,jc,kc) *  detadz(ic,jc,kc) +&
!                               Ucons0(iVar,ic,jc,kc) *  detadt(ic,jc,kc)) / Jaco(ic,jc,kc);
!        
!        invHhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) * dzetadx(ic,jc,kc) +&
!                                 invG(iVar,ic,jc,kc) * dzetady(ic,jc,kc) +&
!                                 invH(iVar,ic,jc,kc) * dzetadz(ic,jc,kc) +&
!                               Ucons0(iVar,ic,jc,kc) * dzetadt(ic,jc,kc)) / Jaco(ic,jc,kc);
!
!        Fp(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) + Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F+ in the xi direction
!        Fm(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) - Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F- in the xi direction
!
!        Gp(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) + Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G+ in the eta direction
!        Gm(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) - Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G- in the eta direction
!
!        Hp(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) + Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H+ in the zeta direction
!        Hm(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) - Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H- in the zeta direction
!      enddo
!    enddo
!   enddo
!  enddo
!  
!  call Parallel_Exchange_NumVar(Fp)
!  call Parallel_Exchange_NumVar(Fm)
!  call Parallel_Exchange_NumVar(Gp)
!  call Parallel_Exchange_NumVar(Gm)
!  call Parallel_Exchange_NumVar(Hp)
!  call Parallel_Exchange_NumVar(Hm)
!
!  call Cal_Deri_Dxi_weno5th_up_numvar(Fpdxi,Fp)
!  call Cal_Deri_Dxi_weno5th_do_numvar(Fmdxi,Fm)
!
!  call Cal_Deri_Deta_weno5th_up_numvar(Gpdeta,Gp)
!  call Cal_Deri_Deta_weno5th_do_numvar(Gmdeta,Gm)
!
!  call Cal_Deri_Dzeta_weno5th_up_per_numvar(Hpdzeta,Hp)
!  call Cal_Deri_Dzeta_weno5th_do_per_numvar(Hmdzeta,Hm)
!
!  invFlux = Fpdxi + Fmdxi + Gpdeta + Gmdeta + Hpdzeta + Hmdzeta;
!  
!  do kc = 1,nz_local
!   do jc = 1,Ny
!    do ic = 1,nx_local
!      do iVar = 1,NumVar
!        dUcons0(iVar,ic,jc,kc) = ( - invFlux(iVar,ic,jc,kc) &
!                                 & + visFlux(iVar,ic,jc,kc) &
!                                 & -  Ucons0(iVar,ic,jc,kc) * invJacodt(ic,jc,kc)) * Jaco(ic,jc,kc);
!      enddo
!    enddo
!   enddo
!  enddo
!   
!END SUBROUTINE Calculate_Flux_Weno5th