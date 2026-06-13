!SUBROUTINE Calculate_Flux_Wenochar5th
!    use SF_Constant,           only: ik,rk,NumVar,LaxFriSmall,overLap
!    use SF_CFD_Global,         only: nx_local,Ny,nz_local,Gamma,Pr_Ref,Mach_Ref,Re_Ref,&
!                                     Csthlnd_Ref,Cp_Ref
!  
!    ! matrix and some arrays
!    use SF_CFD_Global,         only: Rho,U,V,W,P,T,Ucons0,Cs,BigU_Xi,BigU_Eta,BigU_Zeta,&
!                                     dxidx,dxidy,dxidz,dxidt,detadx,detady,detadz,detadt,&
!                                     dzetadx,dzetady,dzetadz,dzetadt,invF,invG,invH,&
!                                     Udxi,Udeta,Udzeta,Vdxi,Vdeta,Vdzeta,Wdxi,Wdeta,Wdzeta,&
!                                     Tdxi,Tdeta,Tdzeta,Mu,Udx,Udy,Udz,Vdx,Vdy,Vdz,Wdx,Wdy,Wdz,&
!                                     Tdx,Tdy,Tdz,visF,visG,visH,visFhat,visGhat,visHhat,Jaco,invJacodt,&
!                                     nablaxi,nablaeta,nablazeta,invFhat,invGhat,invHhat,Fp,Gp,Hp,&
!                                     Fm,Gm,Hm,Fpdxi,Fmdxi,Gpdeta,Gmdeta,Hpdzeta,Hmdzeta,visFcdxi,&
!                                     visGcdeta,visHcdzeta,invFlux,visFlux,dUcons0
!  
!    ! subroutine and function declaration
!    use MPI_GLOBAL,            only: Parallel_Exchange,Parallel_Exchange_NumVar,npx,npx0,npz,npz0
!    use FD5_Order,             only: Cal_Deri_Dxi_5th_ce,Cal_Deri_Deta_5th_ce,Cal_Deri_Dzeta_5th_ce_per,& ! Central
!                Cal_Deri_Dxi_5th_ce_numvar,Cal_Deri_Deta_5th_ce_numvar,Cal_Deri_Dzeta_5th_ce_per_numvar,& ! Central
!    use FD5_Order_Weno,        only: &
!                Cal_Deri_Dxi_weno5th_up_numvar,Cal_Deri_Deta_weno5th_up_numvar,Cal_Deri_Dzeta_weno5th_up_per_numvar,& ! Upwind
!                Cal_Deri_Dxi_weno5th_do_numvar,Cal_Deri_Deta_weno5th_do_numvar,Cal_Deri_Dzeta_weno5th_do_per_numvar ! Downwind
!  
!    implicit none
!    ! Local variables
!    integer( kind = ik ):: ic,jc,kc,iVar,i,j,k,mm,ia,ib,ka,kb
!    real( kind = rk ):: Div_Vel,Tauxx,Tauyy,Tauzz,Tauxy,Tauxz,Tauyz
!    real( kind = rk ):: Cs_Xi,Cs_Eta,Cs_Zeta
!    real( kind = rk ):: Sigma_Xi,Sigma_Eta,Sigma_Zeta
!  
!    real( kind = rk ):: local_temp
!    ! Characteristic variables
!    real( kind = rk ):: n(3),l(3),m(3),n_amp,temp
!    real( kind = rk ):: u_ave,v_ave,w_ave,c_ave,q_ave,H_ave
!    real( kind = rk ):: un,ul,um, kdcc,kd2cc,d2c
!    real( kind = rk), dimension(5,5):: matrix_charL, matrix_charR
!    real( kind = rk ), dimension(5,5):: fpc, fmc, gpc, gmc, hpc, hmc
!    real( kind = rk ), dimension(5):: dfc, dfpc, dfmc
!    real( kind = rk ), dimension(NumVar,0:nx_local):: dfi
!    real( kind = rk ), dimension(NumVar,1:Ny):: dfj
!    real( kind = rk ), dimension(NumVar,0:nz_local):: dfk
!
!  
!    ! Get the inner variables
!    do kc = 1, nz_local
!     do jc = 1, Ny
!      do ic = 1, nx_local
!        ! conservative variables to primitive variables
!        call CVtoU_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
!        ! Calculate the temperature
!        T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc)
!        ! Calculate the sound speed
!        Cs(ic,jc,kc)= sqrt(Gamma * P(ic,jc,kc) / Rho(ic,jc,kc))
!        ! Calculate the characteristic speed
!          BigU_Xi(ic,jc,kc) = U(ic,jc,kc)*  dxidx(ic,jc,kc) &
!                           &+ V(ic,jc,kc)*  dxidy(ic,jc,kc) &
!                           &+ W(ic,jc,kc)*  dxidz(ic,jc,kc) + dxidt(ic,jc,kc)
!  
!         BigU_Eta(ic,jc,kc) = U(ic,jc,kc)* detadx(ic,jc,kc) &
!                           &+ V(ic,jc,kc)* detady(ic,jc,kc) &
!                           &+ W(ic,jc,kc)* detadz(ic,jc,kc) + detadt(ic,jc,kc)
!  
!        BigU_Zeta(ic,jc,kc) = U(ic,jc,kc)*dzetadx(ic,jc,kc) &
!                           &+ V(ic,jc,kc)*dzetady(ic,jc,kc) &
!                           &+ W(ic,jc,kc)*dzetadz(ic,jc,kc) + dzetadt(ic,jc,kc)
!  
!        ! Calculate the flux
!        call UtoFlux_CPG3D(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),&
!                           invF(1:NumVar,ic,jc,kc),invG(1:NumVar,ic,jc,kc),invH(1:NumVar,ic,jc,kc))
!      enddo
!     enddo
!    enddo
!    
!    ! Parallel transfer
!    call Parallel_Exchange(U)
!    call Parallel_Exchange(V)
!    call Parallel_Exchange(W)
!    call Parallel_Exchange(T)
!    call Parallel_Exchange(Cs)
!  
!    ! Calculate the viscous Flux by using the central difference
!    call Cal_Deri_Dxi_5th_ce(Udxi,U)
!    call Cal_Deri_Dxi_5th_ce(Vdxi,V)
!    call Cal_Deri_Dxi_5th_ce(Wdxi,W)
!    call Cal_Deri_Dxi_5th_ce(Tdxi,T)
!  
!    call Cal_Deri_Deta_5th_ce(Udeta,U)
!    call Cal_Deri_Deta_5th_ce(Vdeta,V)
!    call Cal_Deri_Deta_5th_ce(Wdeta,W)
!    call Cal_Deri_Deta_5th_ce(Tdeta,T)
!    
!    ! Here we use the periodic boundary condition form
!    call Cal_Deri_Dzeta_5th_ce_per(Udzeta,U)
!    call Cal_Deri_Dzeta_5th_ce_per(Vdzeta,V)
!    call Cal_Deri_Dzeta_5th_ce_per(Wdzeta,W)
!    call Cal_Deri_Dzeta_5th_ce_per(Tdzeta,T)
!    
!    ! Cp_Ref = 1.0_rk / ((Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
!    local_temp = (Pr_Ref * (Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
!    ! Calculate the viscous stress
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      do ic = 1,nx_local
!        ! The viscosity coefficient has been divided by Re_Ref already
!        ! T^1.5 * (1 + Csthlnd_Ref) / (T + Csthlnd_Ref) / Re_Ref 
!        Mu(ic,jc,kc) = (T(ic,jc,kc)**1.5_rk) * ( 1.0_rk + Csthlnd_Ref )/ (T(ic,jc,kc) + Csthlnd_Ref) 
!        Mu(ic,jc,kc) = Mu(ic,jc,kc) / Re_Ref;
!        ! derivatives with respect to the x direction
!        Udx(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                         Udeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                        Udzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!  
!        Vdx(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                         Vdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                        Vdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);      
!  
!        Wdx(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                         Wdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                        Wdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!  
!        Tdx(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                         Tdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
!                        Tdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);
!  
!        ! derivatives with respect to the y direction
!        Udy(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                         Udeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                        Udzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!  
!        Vdy(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                         Vdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                        Vdzeta(ic,jc,kc)*dzetady(ic,jc,kc);      
!  
!        Wdy(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                         Wdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                        Wdzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!  
!        Tdy(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                         Tdeta(ic,jc,kc)* detady(ic,jc,kc) +&
!                        Tdzeta(ic,jc,kc)*dzetady(ic,jc,kc);
!        
!        ! derivatives with respect to the z direction
!        Udz(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                         Udeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                        Udzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!  
!        Vdz(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                         Vdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                        Vdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);      
!  
!        Wdz(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                         Wdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                        Wdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!  
!        Tdz(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
!                         Tdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
!                        Tdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
!                        
!        Div_Vel = Udx(ic,jc,kc) + Vdy(ic,jc,kc) + Wdz(ic,jc,kc);
!  
!        Tauxx = 2.0_rk * Mu(ic,jc,kc) * (Udx(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!        Tauyy = 2.0_rk * Mu(ic,jc,kc) * (Vdy(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!        Tauzz = 2.0_rk * Mu(ic,jc,kc) * (Wdz(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
!        Tauxy = Mu(ic,jc,kc) * (Udy(ic,jc,kc) + Vdx(ic,jc,kc));
!        Tauxz = Mu(ic,jc,kc) * (Udz(ic,jc,kc) + Wdx(ic,jc,kc));
!        Tauyz = Mu(ic,jc,kc) * (Vdz(ic,jc,kc) + Wdy(ic,jc,kc));
!  
!        ! Calculate the viscous flux in the x-,y-,z- directions
!        visF(1,ic,jc,kc) = 0.0_rk;
!        visF(2,ic,jc,kc) = Tauxx;
!        visF(3,ic,jc,kc) = Tauxy;
!        visF(4,ic,jc,kc) = Tauxz;
!        visF(5,ic,jc,kc) = U(ic,jc,kc) * Tauxx &
!                        &+ V(ic,jc,kc) * Tauxy &
!                        &+ W(ic,jc,kc) * Tauxz &
!                        &+ Mu(ic,jc,kc) * Tdx(ic,jc,kc) / local_temp;
!        
!        visG(1,ic,jc,kc) = 0.0_rk;
!        visG(2,ic,jc,kc) = Tauxy;
!        visG(3,ic,jc,kc) = Tauyy;
!        visG(4,ic,jc,kc) = Tauyz;
!        visG(5,ic,jc,kc) = U(ic,jc,kc) * Tauxy &
!                        &+ V(ic,jc,kc) * Tauyy &
!                        &+ W(ic,jc,kc) * Tauyz &
!                        &+ Mu(ic,jc,kc) * Tdy(ic,jc,kc) / local_temp;
!        
!        visH(1,ic,jc,kc) = 0.0_rk;
!        visH(2,ic,jc,kc) = Tauxz;
!        visH(3,ic,jc,kc) = Tauyz;
!        visH(4,ic,jc,kc) = Tauzz;
!        visH(5,ic,jc,kc) = U(ic,jc,kc) * Tauxz &
!                        &+ V(ic,jc,kc) * Tauyz &
!                        &+ W(ic,jc,kc) * Tauzz &
!                        &+ Mu(ic,jc,kc) * Tdz(ic,jc,kc) / local_temp;
!                        
!        do iVar = 2,NumVar
!          visFhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*  dxidx(ic,jc,kc) +&
!                                   visG(iVar,ic,jc,kc)*  dxidy(ic,jc,kc) +&
!                                   visH(iVar,ic,jc,kc)*  dxidz(ic,jc,kc))/Jaco(ic,jc,kc);
!  
!          visGhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)* detadx(ic,jc,kc) +&
!                                   visG(iVar,ic,jc,kc)* detady(ic,jc,kc) +&
!                                   visH(iVar,ic,jc,kc)* detadz(ic,jc,kc))/Jaco(ic,jc,kc);
!  
!          visHhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*dzetadx(ic,jc,kc) +&
!                                   visG(iVar,ic,jc,kc)*dzetady(ic,jc,kc) +&
!                                   visH(iVar,ic,jc,kc)*dzetadz(ic,jc,kc))/Jaco(ic,jc,kc);
!        enddo
!      enddo
!     enddo
!    enddo
!  
!    call Parallel_Exchange_NumVar(visFhat)
!    call Parallel_Exchange_NumVar(visGhat)
!    call Parallel_Exchange_NumVar(visHhat)
!    call Cal_Deri_Dxi_5th_ce_numvar(visFcdxi,        visFhat)
!    call Cal_Deri_Deta_5th_ce_numvar(visGcdeta,      visGhat)
!    call Cal_Deri_Dzeta_5th_ce_per_numvar(visHcdzeta,visHhat)
!  
!    visFlux = visFcdxi + visGcdeta + visHcdzeta;
!  
!    ! Calculate the inviscid flux
!    invFlux  = 0.0_rk
!    ! Calculate Fp,Fm,Gp,Gm,Hp,Hm
!    do kc = 1-overLAP,nz_local+overLAP
!      do jc = 1,Ny
!        do ic = 1-overLAP,nx_local+overLAP
!          Cs_Xi  = Cs(ic,jc,kc) *   nablaxi(ic,jc,kc);
!          Cs_Eta = Cs(ic,jc,kc) *  nablaeta(ic,jc,kc);
!          Cs_Zeta= Cs(ic,jc,kc) * nablazeta(ic,jc,kc);
!
!          Sigma_Xi  =(sqrt(  BigU_Xi(ic,jc,kc)*  BigU_Xi(ic,jc,kc) &
!                          &+ LaxFriSmall * LaxFriSmall * Cs_Xi  *  Cs_Xi) +  Cs_Xi)/jaco(ic,jc,kc);
!
!          Sigma_Eta =(sqrt( BigU_Eta(ic,jc,kc)* BigU_Eta(ic,jc,kc) &
!                          &+ LaxFriSmall * LaxFriSmall * Cs_Eta * Cs_Eta) + Cs_Eta)/jaco(ic,jc,kc);
!
!          Sigma_Zeta=(sqrt(BigU_Zeta(ic,jc,kc)*BigU_Zeta(ic,jc,kc) &
!                          &+ LaxFriSmall * LaxFriSmall * Cs_Zeta*Cs_Zeta) +Cs_Zeta)/jaco(ic,jc,kc);
!          do iVar = 1,NumVar
!            invFhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *   dxidx(ic,jc,kc) +&
!                                    invG(iVar,ic,jc,kc) *   dxidy(ic,jc,kc) +&
!                                    invH(iVar,ic,jc,kc) *   dxidz(ic,jc,kc) +&
!                                  Ucons0(iVar,ic,jc,kc) *   dxidt(ic,jc,kc)) / Jaco(ic,jc,kc);
!            
!            invGhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *  detadx(ic,jc,kc) +&
!                                    invG(iVar,ic,jc,kc) *  detady(ic,jc,kc) +&
!                                    invH(iVar,ic,jc,kc) *  detadz(ic,jc,kc) +&
!                                  Ucons0(iVar,ic,jc,kc) *  detadt(ic,jc,kc)) / Jaco(ic,jc,kc);
!            
!            invHhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) * dzetadx(ic,jc,kc) +&
!                                    invG(iVar,ic,jc,kc) * dzetady(ic,jc,kc) +&
!                                    invH(iVar,ic,jc,kc) * dzetadz(ic,jc,kc) +&
!                                  Ucons0(iVar,ic,jc,kc) * dzetadt(ic,jc,kc)) / Jaco(ic,jc,kc);
!
!            Fp(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) + Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F+ in the xi direction
!            Fm(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) - Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F- in the xi direction
!
!            Gp(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) + Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G+ in the eta direction
!            Gm(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) - Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G- in the eta direction
!
!            Hp(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) + Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H+ in the zeta direction
!            Hm(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) - Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H- in the zeta direction
!          enddo
!        enddo
!      enddo
!    enddo
!    ! Xi direction
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ia = 0; ib = nx_local;
!        if (npx == 0) ia = 1
!        if (npx == npx0-1) ib = nx_local - 1
!        do ic = ia, ib
!          ! Xi direction characteristic projection
!          ! Calculation characteristic direction
!          n(1) = (dxidx(ic, jc, kc) / nablaxi(ic, jc, kc) + &
!               &  dxidx(ic+1, jc, kc) / nablaxi(ic+1, jc, kc)) / 2.0_rk
!          n(2) = (dxidy(ic, jc, kc) / nablaxi(ic, jc, kc) + &
!               &  dxidy(ic+1, jc, kc) / nablaxi(ic+1, jc, kc)) / 2.0_rk
!          n(3) = (dxidz(ic, jc, kc) / nablaxi(ic, jc, kc) + &
!               &  dxidz(ic+1, jc, kc) / nablaxi(ic+1, jc, kc)) / 2.0_rk
!          n_amp = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
!          n = n / n_amp
!          if( abs(n(2)) < abs(n(3)))then
!            temp = sqrt( n(1)**2 + n(2)**2 )
!            l(1) = -n(2) / temp;  l(2) =  n(1) / temp;  l(3) = 0.0_rk;
!          else 
!            temp = sqrt( n(1)**2 + n(3)**2 )
!            l(1) = -n(3) / temp;  l(2) = 0.0_rk;        l(3) =  n(1) / temp;
!          end if
!          m(1) = n(2)*l(3) - n(3)*l(2)
!          m(2) = n(3)*l(1) - n(1)*l(3)
!          m(3) = n(1)*l(2) - n(2)*l(1)
!          ! Calculate average speed at j+1/2
!          u_ave = ( U(ic,jc,kc) + U(ic+1,jc,kc) ) / 2.0_rk
!          v_ave = ( V(ic,jc,kc) + V(ic+1,jc,kc) ) / 2.0_rk
!          w_ave = ( W(ic,jc,kc) + W(ic+1,jc,kc) ) / 2.0_rk
!          c_ave = ( Cs(ic,jc,kc) + Cs(ic+1,jc,kc) ) / 2.0_rk
!          un = u_ave*n(1) + v_ave*n(2) + w_ave*n(3)
!          ul = u_ave*l(1) + v_ave*l(2) + w_ave*l(3)
!          um = u_ave*m(1) + v_ave*m(2) + w_ave*m(3)
!          ! Some quick variables
!          q_ave = ( u_ave**2 + v_ave**2 + w_ave**2 ) / 2.0_rk
!          kdcc = (Gamma - 1.0_rk) / ( c_ave**2 )
!          kd2cc = kcc / 2.0_rk
!          d2c = 0.5_rk / c_ave
!          H_ave = q_ave + 1.0_rk / kdcc
!          ! Construct the matrix
!          ! Left matrix
!          matrix_charL  = reshape( &
!          [ 1.0_rk-kdcc*q_ave, kdcc*u_ave, kdcc*v_ave, kdcc*w_ave, -kdcc,  &
!            -ul,               l(1),       l(2),       l(3),       0.0_rk, &
!            -um,               m(1),       m(2),       m(3),       0.0_rk, &
!            kd2cc*q_ave+d2c*un, -d2c*n(1)-kd2cc*u_ave, -d2c*n(2)-kd2cc*v_ave, -d2c*n(3)-kd2cc*w_ave, kd2cc, &
!            kd2cc*q_ave-d2c*un,  d2c*n(1)-kd2cc*u_ave,  d2c*n(2)-kd2cc*v_ave,  d2c*n(3)-kd2cc*w_ave, kd2cc], &
!           , [5,5], order = [2,1] )
!          ! Right matrix
!          matrix_charR = reshape( &
!          [ 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, &
!            u_ave,  l(1),  m(1),  u_ave-c_ave*n(1), u_ave+c_ave*n(1), &
!            v_ave,  l(2),  m(2),  v_ave-c_ave*n(2), v_ave+c_ave*n(2), &
!            w_ave,  l(3),  m(3),  w_ave-c_ave*n(3), w_ave+c_ave*n(3), &
!            q_ave,  ul,    um,    H_ave-c_ave*un,   H_ave+c_ave*un], &
!          , [5,5], order = [2,1] )
!          ! V = S * U
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                fpc(mm, i) = matrix_charL(i, j) * Fp(j, ic-3+m, jc, kc)
!              enddo
!            enddo
!          enddo
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                fmc(mm, i) = matrix_charL(i, j) * Fm(j, ic-2+mm, jc, kc)
!              enddo
!            enddo
!          enddo
!          ! Calculate WENO scheme
!          call Cal_Deri_Dxi_wenochar5thlocal_up_numvar(dfpc, fpc, ic)
!          call Cal_Deri_Dxi_wenochar5thlocal_do_numvar(dfmc, fmc, ic)
!          dfc = dfpc + dfmc
!          ! Transfer back to physical space
!          ! dfi matrix stores F(j+1/2)
!          do i = 1,5
!            do j = 1,5
!              dfi(i, ic) = dfi(i, ic) + matrix_charR(i, j) * dfc(j)
!            enddo
!          enddo
!        enddo ! Loop over Xi direction
!        ! No reflection boundary condition to F(0+1/2) and F(nx+1/2)
!        if (npx == 0) then
!          do i = 1,5
!            dfi(i, 0) = ( 2.0_rk*Fp(i,1,jc,kc) + 5.0_rk*Fp(i,2,jc,kc) - 1.0_rk*Fp(i,3,jc,kc) &
!                      & + 2.0_rk*Fm(i,3,jc,kc) - 7.0_rk*Fm(i,2,jc,kc) + 11.0_rk*Fm(i,1,jc,kc)) &
!                      & / 6.0_rk
!          enddo
!        endif 
!        if(npx == npx0-1) then
!          do i = 1,5
!            dfi(i, nx_local) = ( 2.0_rk*Fp(i,nx_local-2,jc,kc) - 7.0_rk*Fp(i,nx_local-1,jc,kc) + 11.0_rk*Fp(i,nx_local,jc,kc) &
!                             & + 2.0_rk*Fm(i,nx_local,jc,kc) + 5.0_rk*Fm(i,nx_local-1,jc,kc) - 1.0_rk*Fm(i,nx_local-2,jc,kc)) &
!                             & / 6.0_rk
!          enddo
!        enddo
!        ! Calculate inviscid flux
!        do ic = 1,nx_local
!          invFlux(:, ic, jc, kc) = invFlux(:, ic, jc, kc) + dfi(:, ic) - dfi(:, ic-1)
!        enddo
!      enddo
!    enddo   
!
!    ! Eta direction
!    do kc = 1,nz_local
!      do ic = 1,nx_local
!        ! Different from Xi and Zeta direction
!        ! Eta direction no MPI partition
!        ia = 1; ib = Ny-1;
!        do jc = ia, ib
!          ! Eta direction characteristic projection
!          ! Calculation characteristic direction
!          n(1) = (detadx(ic, jc, kc) / nablaeta(ic, jc, kc) + &
!               &  detadx(ic, jc+1, kc) / nablaeta(ic, jc+1, kc)) / 2.0_rk
!          n(2) = (detady(ic, jc, kc) / nablaeta(ic, jc, kc) + &
!               &  detady(ic, jc+1, kc) / nablaeta(ic, jc+1, kc)) / 2.0_rk
!          n(3) = (detadz(ic, jc, kc) / nablaeta(ic, jc, kc) + &
!               &  detadz(ic, jc+1, kc) / nablaeta(ic, jc+1, kc)) / 2.0_rk
!          n_amp = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
!          n = n / n_amp
!          if( abs(n(1)) < abs(n(3)))then
!            temp = sqrt( n(2)**2 + n(3)**2 )
!            l(1) = 0.0_rk;        l(2) = -n(3) / temp;  l(3) =  n(2) / temp;
!          else 
!            temp = sqrt( n(1)**2 + n(2)**2 )
!            l(1) = -n(2) / temp;  l(2) =  n(1) / temp;  l(3) = 0.0_rk;
!          end if
!          m(1) = n(2)*l(3) - n(3)*l(2)
!          m(2) = n(3)*l(1) - n(1)*l(3)
!          m(3) = n(1)*l(2) - n(2)*l(1)
!          ! Calculate average speed at j+1/2
!          u_ave = ( U(ic,jc,kc) + U(ic,jc+1,kc) ) / 2.0_rk
!          v_ave = ( V(ic,jc,kc) + V(ic,jc+1,kc) ) / 2.0_rk
!          w_ave = ( W(ic,jc,kc) + W(ic,jc+1,kc) ) / 2.0_rk
!          c_ave = ( Cs(ic,jc,kc) + Cs(ic,jc+1,kc) ) / 2.0_rk
!          un = u_ave*n(1) + v_ave*n(2) + w_ave*n(3)
!          ul = u_ave*l(1) + v_ave*l(2) + w_ave*l(3)
!          um = u_ave*m(1) + v_ave*m(2) + w_ave*m(3)
!          ! Some quick variables
!          q_ave = ( u_ave**2 + v_ave**2 + w_ave**2 ) / 2.0_rk
!          kdcc = (Gamma - 1.0_rk) / ( c_ave**2 )
!          kd2cc = kcc / 2.0_rk
!          d2c = 0.5_rk / c_ave
!          H_ave = q_ave + 1.0_rk / kdcc
!          ! Construct the matrix
!          ! Left matrix
!          matrix_charL  = reshape( &
!          [ 1.0_rk-kdcc*q_ave, kdcc*u_ave, kdcc*v_ave, kdcc*w_ave, -kdcc,  &
!            -ul,               l(1),       l(2),       l(3),       0.0_rk, &
!            -um,               m(1),       m(2),       m(3),       0.0_rk, &
!            kd2cc*q_ave+d2c*un, -d2c*n(1)-kd2cc*u_ave, -d2c*n(2)-kd2cc*v_ave, -d2c*n(3)-kd2cc*w_ave, kd2cc, &
!            kd2cc*q_ave-d2c*un,  d2c*n(1)-kd2cc*u_ave,  d2c*n(2)-kd2cc*v_ave,  d2c*n(3)-kd2cc*w_ave, kd2cc], &
!           , [5,5], order = [2,1] )
!          ! Right matrix
!          matrix_charR = reshape( &
!          [ 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, &
!            u_ave,  l(1),  m(1),  u_ave-c_ave*n(1), u_ave+c_ave*n(1), &
!            v_ave,  l(2),  m(2),  v_ave-c_ave*n(2), v_ave+c_ave*n(2), &
!            w_ave,  l(3),  m(3),  w_ave-c_ave*n(3), w_ave+c_ave*n(3), &
!            q_ave,  ul,    um,    H_ave-c_ave*un,   H_ave+c_ave*un], &
!          , [5,5], order = [2,1] )
!          ! V = S * U
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                gpc(mm, i) = matrix_charL(i, j) * Gp(j, ic, jc-3+m, kc)
!              enddo
!            enddo
!          enddo
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                gmc(mm, i) = matrix_charL(i, j) * Gm(j, ic, jc-2+mm, kc)
!              enddo
!            enddo
!          enddo
!          ! Calculate WENO scheme
!          call Cal_Deri_Deta_wenochar5thlocal_up_numvar(dgpc, gpc, jc)
!          call Cal_Deri_Deta_wenochar5thlocal_do_numvar(dgmc, gmc, jc)
!          dgc = dgpc + dgmc
!          ! Transfer back to physical space
!          ! dgi matrix stores G(j+1/2)
!          do i = 1,5
!            do j = 1,5
!              dgi(i, jc) = dgi(i, jc) + matrix_charR(i, j) * dgc(j)
!            enddo
!          enddo
!        enddo ! Loop over Eta direction
!        ! No reflection boundary condition to G(0+1/2) and G(ny+1/2)
!        do i = 1,5
!          dgi(i, 0) = ( 2.0_rk*Gp(i,ic,1,kc) + 5.0_rk*Gp(i,ic,2,kc) - 1.0_rk*Gp(i,ic,3,kc) &
!                    & + 2.0_rk*Gm(i,ic,3,kc) - 7.0_rk*Gm(i,ic,2,kc) + 11.0_rk*Gm(i,ic,1,kc)) &
!                    & / 6.0_rk
!          dgi(i, Ny) = ( 2.0_rk*Gp(i,ic,Ny-2,kc) - 7.0_rk*Gp(i,ic,Ny-1,kc) + 11.0_rk*Gp(i,ic,Ny,kc) &
!                    & + 2.0_rk*Gm(i,ic,Ny,kc) + 5.0_rk*Gm(i,ic,Ny-1,kc) - 1.0_rk*Gm(i,ic,Ny-2,kc)) &
!                    & / 6.0_rk
!        enddo
!        ! Calculate inviscid flux
!        do jc = 1,Ny
!          invFlux(:, ic, jc, kc) = invFlux(:, ic, jc, kc) + dgi(:, jc) - dgi(:, jc-1)
!        enddo
!      enddo
!    enddo
!    ! Zeta direction
!    do jc = 1,Ny
!      do ic = 1,nx_local
!        ia = 0; ib = nz_local;
!        if ( npz == 0 ) ia = 1
!        if ( npz == npz0-1 ) ib = nz_local - 1
!        do kc = ia, ib
!          ! Zeta direction characteristic projection
!          ! Calculation characteristic direction
!          n(1) = (dzetadx(ic, jc, kc) / nablazeta(ic, jc, kc) + &
!               &  dzetadx(ic, jc, kc+1) / nablazeta(ic, jc, kc+1)) / 2.0_rk
!          n(2) = (dzetady(ic, jc, kc) / nablazeta(ic, jc, kc) + &
!               &  dzetady(ic, jc, kc+1) / nablazeta(ic, jc, kc+1)) / 2.0_rk
!          n(3) = (dzetadz(ic, jc, kc) / nablazeta(ic, jc, kc) + &
!               &  dzetadz(ic, jc, kc+1) / nablazeta(ic, jc, kc+1)) / 2.0_rk
!          n_amp = sqrt(n(1)**2 + n(2)**2 + n(3)**2)
!          n = n / n_amp
!          if( abs(n(1)) < abs(n(2)))then
!            temp = sqrt( n(2)**2 + n(3)**2 )
!            l(1) = 0.0_rk;        l(2) = -n(3) / temp;  l(3) =  n(2) / temp;
!          else 
!            temp = sqrt( n(1)**2 + n(2)**2 )
!            l(1) = -n(2) / temp;  l(2) =  n(1) / temp;  l(3) = 0.0_rk;
!          end if
!          m(1) = n(2)*l(3) - n(3)*l(2)
!          m(2) = n(3)*l(1) - n(1)*l(3)
!          m(3) = n(1)*l(2) - n(2)*l(1)
!          ! Calculate average speed at j+1/2
!          u_ave = ( U(ic,jc,kc) + U(ic,jc,kc+1) ) / 2.0_rk
!          v_ave = ( V(ic,jc,kc) + V(ic,jc,kc+1) ) / 2.0_rk
!          w_ave = ( W(ic,jc,kc) + W(ic,jc,kc+1) ) / 2.0_rk
!          c_ave = ( Cs(ic,jc,kc) + Cs(ic,jc,kc+1) ) / 2.0_rk
!          un = u_ave*n(1) + v_ave*n(2) + w_ave*n(3)
!          ul = u_ave*l(1) + v_ave*l(2) + w_ave*l(3)
!          um = u_ave*m(1) + v_ave*m(2) + w_ave*m(3)
!          ! Some quick variables
!          q_ave = ( u_ave**2 + v_ave**2 + w_ave**2 ) / 2.0_rk
!          kdcc = (Gamma - 1.0_rk) / ( c_ave**2 )
!          kd2cc = kcc / 2.0_rk
!          d2c = 0.5_rk / c_ave
!          H_ave = q_ave + 1.0_rk / kdcc
!          ! Construct the matrix
!          ! Left matrix
!          matrix_charL  = reshape( &
!          [ 1.0_rk-kdcc*q_ave, kdcc*u_ave, kdcc*v_ave, kdcc*w_ave, -kdcc,  &
!            -ul,               l(1),       l(2),       l(3),       0.0_rk, &
!            -um,               m(1),       m(2),       m(3),       0.0_rk, &
!            kd2cc*q_ave+d2c*un, -d2c*n(1)-kd2cc*u_ave, -d2c*n(2)-kd2cc*v_ave, -d2c*n(3)-kd2cc*w_ave, kd2cc, &
!            kd2cc*q_ave-d2c*un,  d2c*n(1)-kd2cc*u_ave,  d2c*n(2)-kd2cc*v_ave,  d2c*n(3)-kd2cc*w_ave, kd2cc], &
!           , [5,5], order = [2,1] )
!          ! Right matrix
!          matrix_charR = reshape( &
!          [ 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, &
!            u_ave,  l(1),  m(1),  u_ave-c_ave*n(1), u_ave+c_ave*n(1), &
!            v_ave,  l(2),  m(2),  v_ave-c_ave*n(2), v_ave+c_ave*n(2), &
!            w_ave,  l(3),  m(3),  w_ave-c_ave*n(3), w_ave+c_ave*n(3), &
!            q_ave,  ul,    um,    H_ave-c_ave*un,   H_ave+c_ave*un], &
!          , [5,5], order = [2,1] )
!          ! V = S * U
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                hpc(mm, i) = matrix_charL(i, j) * Hp(j, ic, jc, kc-3+m)
!              enddo
!            enddo
!          enddo
!          do mm = 1,5
!            do i = 1,5
!              do j = 1,5
!                hmc(mm, i) = matrix_charL(i, j) * Hm(j, ic, jc, kc-2+mm)
!              enddo
!            enddo
!          enddo
!          ! Calculate WENO scheme
!          call Cal_Deri_Dzeta_wenochar5thlocal_up_numvar(dhpc, hpc, kc)
!          call Cal_Deri_Dzeta_wenochar5thlocal_do_numvar(dhmc, hmc, kc)
!          dhc = dhpc + dhmc
!          ! Transfer back to physical space
!          ! dhi matrix stores H(j+1/2)
!          do i = 1,5
!            do j = 1,5
!              dhi(i, kc) = dhi(i, kc) + matrix_charR(i, j) * dhc(j)
!            enddo
!          enddo
!        enddo ! Loop over Zeta direction
!        ! No reflection boundary condition to H(0+1/2) and H(nz+1/2)
!        if (npz == 0) then
!          do i = 1,5
!            dhi(i, 0) = ( 2.0_rk*Hp(i,ic,jc,1) + 5.0_rk*Hp(i,ic,jc,2) - 1.0_rk*Hp(i,ic,jc,3) &
!                      & + 2.0_rk*Hm(i,ic,jc,3) - 7.0_rk*Hm(i,ic,jc,2) + 11.0_rk*Hm(i,ic,jc,1)) &
!                      & / 6.0_rk
!          enddo
!        endif
!        if(npz == npz0-1) then
!          do i = 1,5
!            dhi(i, nz_local) = ( 2.0_rk*Hp(i,ic,jc,nz_local-2) - 7.0_rk*Hp(i,ic,jc,nz_local-1) + 11.0_rk*Hp(i,ic,jc,nz_local) &
!                             & + 2.0_rk*Hm(i,ic,jc,nz_local) + 5.0_rk*Hm(i,ic,jc,nz_local-1) - 1.0_rk*Hm(i,ic,jc,nz_local-2)) &
!                             & / 6.0_rk
!          enddo
!        endif 
!        ! Calculate inviscid flux
!        do kc = 1,nz_local
!          invFlux(:, ic, jc, kc) = invFlux(:, ic, jc, kc) + dhi(:, kc) - dhi(:, kc-1)
!        enddo
!      enddo
!    enddo
!
!    ! Total flux
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      do ic = 1,nx_local
!        do iVar = 1,NumVar
!          dUcons0(iVar,ic,jc,kc) = ( - invFlux(iVar,ic,jc,kc) &
!                                   & + visFlux(iVar,ic,jc,kc) &
!                                   & -  Ucons0(iVar,ic,jc,kc) * invJacodt(ic,jc,kc)) * Jaco(ic,jc,kc);
!        enddo
!      enddo
!     enddo
!    enddo
!     
!  END SUBROUTINE Calculate_Flux_Wenochar5th