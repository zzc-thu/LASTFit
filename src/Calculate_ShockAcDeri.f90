!subroutine Cal_ShockAcDeri
!  ! This subroutine calculate the shock acceleration derivatives
!  ! with respect to the geometry variables
!  !     H, Hxi, Hzeta, V, Vxi, Vzeta
!  ! dAcdH, dAcdHxi, dAcdHzeta, dAcdV, dAcdVxi, dAcdVzeta
!
!    use SF_Constant,    only:ik,rk,NumVar
!    use SF_CFD_Global,  only:nx_local,Ny,nz_local,Gamma
!    ! arrays
!    
!
!    implicit none
!
!    integer( kind = ik)::ic,jc,kc
!
!    real( kind = rk )::xdxi,xdeta,xdzeta
!    real( kind = rk )::ydxi,ydeta,ydzeta
!    real( kind = rk )::zdxi,zdeta,zdzeta
!    real( kind = rk )::xxitau,xetatau,xzetatau
!    real( kind = rk )::yxitau,yetatau,yzetatau
!    real( kind = rk )::zxitau,zetatau,zzetatau
!    real( kind = rk )::xix,etax,zetax
!    real( kind = rk )::xiy,etay,zetay
!    real( kind = rk )::xiz,etaz,zetaz
!    real( kind = rk )::jaco,temp
!    real( kind = rk )::invJacodt,nablaxi,nablaeta,nablazeta
!    real( kind = rk )::ShockNormalX,ShockNormalY,ShockNormalZ
!    real( kind = rk )::ShockXtau,ShockYtau,ShockZtau
!    real( kind = rk )::invF_IJK(5),invG_IJK(5),invH_IJK(5),U0_IJK(5)
!
!    jc = Ny ! Standing for the Shock Surface
!    do kc = 1,nz_local
!      do ic = 1,nx_local
!       ! eps0 = 1e-6
!       ! shockH = ShockH(ic,kc) + eps0
!       ! shockH = ShockH(ic,kc) - eps0
!       
!       ! In the following parts, we need to recalculate the geometry variables
!       ! 
!        xdxi = WallSXdxi(ic,kc) + (WallNormalXdxi(ic,kc)*ShockH(ic,kc) + &
!                                 & WallNormalX(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
!                                 + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
!
!        xxitau = (WallNormalXdxi(ic,kc)*ShockV(ic,kc) + &
!                  & WallNormalX(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
!                  + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
!
!        ydxi = WallSYdxi(ic,kc) + (WallNormalYdxi(ic,kc)*ShockH(ic,kc) + &
!                                 & WallNormalY(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
!                                 + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
!
!        yxitau = (WallNormalYdxi(ic,kc)*ShockV(ic,kc) + &
!                  & WallNormalY(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
!                  + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
!
!        zdxi = WallSZdxi(ic,kc) + (WallNormalZdxi(ic,kc)*ShockH(ic,kc) + &
!                                 & WallNormalZ(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
!                                 + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
!
!        zxitau = (WallNormalZdxi(ic,kc)*ShockV(ic,kc) + &
!                  & WallNormalZ(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
!                  + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
!        
!        xdeta = WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
!        ydeta = WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
!        zdeta = WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
!        
!        xetatau = WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
!        yetatau = WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
!        zetatau = WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
!
!        xdzeta = WallSXdzeta(ic,kc) + (WallNormalXdzeta(ic,kc)*ShockH(ic,kc) + &
!                                      & WallNormalX(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                                     + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
!        
!        xzetatau = (WallNormalXdzeta(ic,kc)*ShockV(ic,kc) + &
!                    & WallNormalX(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                    + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
!
!        ydzeta = WallSYdzeta(ic,kc) + (WallNormalYdzeta(ic,kc)*ShockH(ic,kc) + &
!                                      & WallNormalY(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                                     + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
!
!        yzetatau = (WallNormalYdzeta(ic,kc)*ShockV(ic,kc) + &
!                    & WallNormalY(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                    + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
!      
!        zdzeta = WallSZdzeta(ic,kc) + (WallNormalZdzeta(ic,kc)*ShockH(ic,kc) + &
!                                      & WallNormalZ(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                                     + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
!
!        zzetatau = (WallNormalZdzeta(ic,kc)*ShockV(ic,kc) + &
!                    & WallNormalZ(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
!                    + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
!        
!        temp =    xdxi * (  ydeta* zdzeta - zdeta* ydzeta) &
!             & + xdeta * ( ydzeta*   zdxi -  ydxi* zdzeta) &
!             & +xdzeta * (   ydxi * zdeta - ydeta*   zdxi)
!
!        Jaco = 1.d0/temp;
!
!        invJacodt = & 
!             &    xxitau * (  ydeta* zdzeta - zdeta* ydzeta) &
!             & + xetatau * ( ydzeta*   zdxi -  ydxi* zdzeta) &
!             & +xzetatau * (   ydxi * zdeta - ydeta*   zdxi) &
!             & +  xdxi * (  yetatau * zdzeta - zetatau* ydzeta) &
!             & + xdeta * (  yzetatau*   zdxi -  yxitau* zdzeta) &
!             & +xdzeta * (   yxitau *  zdeta - yetatau*   zdxi) &
!             & +  xdxi * ( ydeta* zzetatau - zdeta* yzetatau) &
!             & + xdeta * (ydzeta*   zxitau -  ydxi* zzetatau) &
!             & +xdzeta * (  ydxi*  zetatau - ydeta*   zxitau)            
!
!          xix=  Jaco * (ydeta*zdzeta - ydzeta* zdeta);
!          xiy= -Jaco * (xdeta*zdzeta -  zdeta*xdzeta);
!          xiz=  Jaco * (xdeta*ydzeta -  ydeta*xdzeta); 
!         etax= -Jaco * ( ydxi*zdzeta -   zdxi*ydzeta);  
!         etay=  Jaco * ( xdxi*zdzeta -   zdxi*xdzeta);
!         etaz= -Jaco * ( xdxi*ydzeta -   ydxi*xdzeta);
!        zetax=  Jaco * ( ydxi* zdeta -   zdxi* ydeta);  
!        zetay= -Jaco * ( xdxi* zdeta -   zdxi* xdeta);
!        zetaz=  Jaco * ( xdxi* ydeta -   ydxi* xdeta); 
!
!          nablaxi=sqrt(  xix *  xix +  xiy *  xiy +  xiz *  xiz)
!
!         nablaeta=sqrt( etax * etax + etay * etay + etaz * etaz)
!
!        nablazeta=sqrt(zetax *zetax +zetay *zetay +zetaz *zetaz)
!
!        ShockNormalX = etax / nablaeta
!        ShockNormalY = etay / nablaeta
!        ShockNormalZ = etaz / nablaeta
!
!        ShockXtau = WallNormalX(ic,kc) * ShockV(ic,kc)
!        ShockYtau = WallNormalY(ic,kc) * ShockV(ic,kc)
!        ShockZtau = WallNormalZ(ic,kc) * ShockV(ic,kc)
!
!        ShockXtaudxi = WallNormalXdxi(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdxi(ic,kc)
!        ShockYtaudxi = WallNormalYdxi(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdxi(ic,kc)
!        ShockZtaudxi =                                         WallNormalZ(ic,kc) * ShockVdxi(ic,kc)
!        
!        ShockXtaudzeta = WallNormalX(ic,kc) * ShockVdzeta(ic,kc)
!        ShockYtaudzeta = WallNormalY(ic,kc) * ShockVdzeta(ic,kc)
!        ShockZtaudzeta = WallNormalZ(ic,kc) * ShockVdzeta(ic,kc)
!
!        tempax = ShockYtaudzeta*zdxi + ShockZtaudxi*ydzeta - &
!               & ShockYtaudxi  *zdzeta - ShockZtaudzeta*ydxi
!        
!        tempay = ShockXtaudxi  *zdzeta + ShockZtaudzeta*xdxi - &
!               & ShockXtaudzeta*zdxi - ShockZtaudxi*xdzeta        
!        
!        tempaz = ShockXtaudzeta*ydxi + ShockYtaudxi*xdzeta - &
!               & ShockXtaudxi*ydzeta - ShockYtaudzeta*xdxi
!        
!        detadt = -( etax*ShockXtau + etay*ShockYtau + etaz*ShockZtau)
!   
!        shockVN = -detadt       /nablaeta;
!         ax_tau =    Jaco*tempax/nablaeta;
!         ay_tau =    Jaco*tempay/nablaeta;        
!         az_tau =    Jaco*tempaz/nablaeta;
!    
!        ! Obtain the variables behind the shock
!        ! Velocity ahead of the shocks with respect to the moving shock normal
!        Delta_Inf = U_inf*ShockNormalX + V_inf*ShockNormalY + W_inf*ShockNormalZ &
!                 &- ShockV(ic,kc)*(WallNormalX(ic,kc)*ShockNormalX &
!                                 &+WallNormalY(ic,kc)*ShockNormalY &
!                                 &+WallNormalZ(ic,kc)*ShockNormalZ);
!        ! Shock normal Mach number^2
!          Masn2Inf = (Mach_Ref * Delta_Inf)**2.d0;
!        ! Velocity normal to the shock after the shock
!          DeltaShkUinf = 2.d0 * Delta_Inf / (Gamma + 1.d0) * (1.d0 / Masn2Inf - 1.d0);
!          U_IJK = U_inf + DeltaShkUinf * ShockNormalX;
!          V_IJK = V_inf + DeltaShkUinf * ShockNormalY;
!          W_IJK = W_inf + DeltaShkUinf * ShockNormalZ;
!          P_IJK = P_inf * (1.d0 + 2.d0 * Gamma / (Gamma + 1.d0) * (Masn2Inf - 1.d0));
!        Rho_IJK = Rho_inf * (Gamma + 1.d0) * Masn2Inf / ((Gamma - 1.d0) * Masn2Inf + 2.d0);
!          T_IJK = P_IJK * Gamma * Mach_Ref * Mach_Ref / Rho_IJK;
!          C_IJK = sqrt(Gamma * P_IJK / Rho_IJK)
!
!        ! Update the fluxes and new conservative variables
!        ! invF_IJK(5),invG_IJK(5),invH_IJK(5),U0_IJK(5)
!        !  CALL UtoCV_CPG(Rho_IJK,U_IJK,V_IJK,W_IJK,P_IJK,U0_IJK)
!          U0_IJK(1)=Rho_IJK
!          U0_IJK(2)=Rho_IJK*U_IJK
!          U0_IJK(3)=Rho_IJK*V_IJK
!          U0_IJK(4)=Rho_IJK*W_IJK
!          U0_IJK(5)=P_IJK/(GAMMA-1.D0)+0.5D0*Rho_IJK*(U_IJK*U_IJK+V_IJK*V_IJK+W_IJK*W_IJK)
!
!        !  CALL UtoFlux_CPG3D(Rho_IJK,U_IJK,V_IJK,W_IJK,P_IJK,invF_IJK,invG_IJK,invH_IJK)
!          invF_IJK(1) = Rho_IJK * U_IJK
!          invF_IJK(2) = Rho_IJK * U_IJK * U_IJK + P_IJK
!          invF_IJK(3) = Rho_IJK * U_IJK * V_IJK
!          invF_IJK(4) = Rho_IJK * U_IJK * W_IJK
!          invF_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * U_IJK
!
!          invG_IJK(1) = Rho_IJK * V_IJK
!          invG_IJK(2) = Rho_IJK * V_IJK * U_IJK
!          invG_IJK(3) = Rho_IJK * V_IJK * V_IJK + P_IJK
!          invG_IJK(4) = Rho_IJK * V_IJK * W_IJK
!          invG_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * V_IJK
!
!          invH_IJK(1) = Rho_IJK * W_IJK
!          invH_IJK(2) = Rho_IJK * W_IJK * U_IJK
!          invH_IJK(3) = Rho_IJK * W_IJK * V_IJK
!          invH_IJK(4) = Rho_IJK * W_IJK * W_IJK + P_IJK
!          invH_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * W_IJK
!
!        ! Calculate the shock accelerations
!              ! Calculate the local normal
!        locUns = U_IJK * ShockNormalX + V_IJK * ShockNormalY + W_IJK * ShockNormalZ;
!        LAM5   = locUns - ShockVN + C_IJK;
!
!        L5(1) = 0.5d0 * (Gamma - 1.d0)*(U_IJK * U_IJK + V_IJK * V_IJK + W_IJK * W_IJK ) - locUns * C_IJK; 
!        L5(2) =        -(Gamma - 1.d0)* U_IJK + ShockNormalX * C_IJK;
!        L5(3) =        -(Gamma - 1.d0)* V_IJK + ShockNormalY * C_IJK;
!        L5(4) =        -(Gamma - 1.d0)* W_IJK + ShockNormalZ * C_IJK;
!        L5(5) =         (Gamma - 1.d0);
!    
!        SUM_L5_U  = 0.d0;
!        SUM_L5_FG = 0.d0;
!        
!        DO iVar = 1,NumVar
!          FG = LAM5 * dUcons0(iVar,ic,jc,kc) + ax_tau * (invF_IJK(iVar) - FINV_INF(iVar)) &
!                                             + ay_tau * (invG_IJK(iVar) - GINV_INF(iVar)) &
!                                             + az_tau * (invH_IJK(iVar) - HINV_INF(iVar));
!    
!          SUM_L5_U  = SUM_L5_U  + L5(iVar) * (U0_IJK(iVar) - CV_INF(iVar));
!          SUM_L5_FG = SUM_L5_FG + L5(iVar) * FG;                       
!        ENDDO
!        
!        ! Calculate the ShockAc
!        ShockAc(ic,kc) = SUM_L5_FG/SUM_L5_U - (ax_tau * ShockXtau &
!                                            &+ ay_tau * ShockYtau &
!                                            &+ az_tau * ShockZtau);
!    
!        ShockAc(ic,kc) = ShockAc(ic,kc)/( WallNormalX(ic,kc) * ShockNormalX &
!                                       &+ WallNormalY(ic,kc) * ShockNormalY &
!                                       &+ WallNormalZ(ic,kc) * ShockNormalZ);
!        
!        ! DacDH, DacDHxi, DacDHzeta, DacDV, DacDVxi, DacDVzeta
!        ! (SHockAc2 - ShockAc1) / (2 * eps0)
!      enddo
!    enddo
!
!end subroutine Cal_ShockAcDeri