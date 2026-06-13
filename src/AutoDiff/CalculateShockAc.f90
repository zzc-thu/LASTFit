subroutine CalculateShockAc(ShockH,ShockHdxi,ShockHdzeta,ShockV,ShockVdxi,ShockVdzeta,ShockAc&
                             &,WallSXdxi,WallSYdxi,WallSZdxi,WallSXdzeta,WallSYdzeta,WallSZdzeta&
                             &,WallNormalX,WallNormalY,WallNormalZ&
                             &,WallNormalXdxi,WallNormalYdxi,WallNormalZdxi&
                             &,WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta&
                             &,Heta,HetaDxi,HetaDeta,HetaDzeta&
                             &,CV_INF,FINV_INF,GINV_INF,HINV_INF,dUcons0&
                             &,U_inf,V_inf,W_inf,P_inf,Rho_inf,Gamma,Mach_Ref)
    implicit none
    integer, parameter:: NumVar = 5
    ! Input variables
    real(8),intent(in)::ShockH,ShockHdxi,ShockHdzeta,ShockV,ShockVdxi,ShockVdzeta
    ! Output variables
    real(8),intent(out)::ShockAc
    
    ! Input unchanged variables with respect to the geometry
    real(8),intent(in)::WallSXdxi,WallSYdxi,WallSZdxi,WallSXdzeta,WallSYdzeta,WallSZdzeta
    real(8),intent(in)::WallNormalX,WallNormalY,WallNormalZ
    real(8),intent(in)::WallNormalXdxi,WallNormalYdxi,WallNormalZdxi
    real(8),intent(in)::WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta
    real(8),intent(in)::Heta,HetaDxi,HetaDeta,HetaDzeta

    ! Flow related variables
    real(8),intent(in)::U_inf,V_inf,W_inf,P_inf,Rho_inf,Gamma,Mach_Ref
    real(8),intent(in)::CV_INF(5),FINV_INF(5),GINV_INF(5),HINV_INF(5),dUcons0(5)

    ! Local variables
    real(8)::xdxi,xdeta,xdzeta
    real(8)::ydxi,ydeta,ydzeta
    real(8)::zdxi,zdeta,zdzeta
    real(8)::xxitau,xetatau,xzetatau
    real(8)::yxitau,yetatau,yzetatau
    real(8)::zxitau,zetatau,zzetatau
    real(8)::temp,Jaco,invJacodt
    real(8)::xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz
    real(8)::nablaxi,nablaeta,nablazeta
    real(8)::ShockNormalX,ShockNormalY,ShockNormalZ
    real(8)::ShockXtau,ShockYtau,ShockZtau
    real(8)::ShockXtaudxi,ShockYtaudxi,ShockZtaudxi
    real(8)::ShockXtaudzeta,ShockYtaudzeta,ShockZtaudzeta
    real(8)::tempax,tempay,tempaz
    real(8)::detadt,shockVN,ax_tau,ay_tau,az_tau
    real(8)::Delta_Inf,Masn2Inf,DeltaShkUinf,U_IJK,V_IJK,W_IJK,P_IJK,Rho_IJK,T_IJK,C_IJK
    real(8)::locUns
    real(8)::LAM5,L5(5),SUM_L5_U,SUM_L5_FG,FG
    integer::iVar
    real(8)::U0_IJK(5),invF_IJK(5),invG_IJK(5),invH_IJK(5)
    
  
! In the following parts, we need to recalculate the geometry variables
       ! 
        xdxi = WallSXdxi + (WallNormalXdxi*ShockH + WallNormalX*ShockHdxi) * Heta &
                         + WallNormalX * ShockH * HetaDxi;

        xxitau = (WallNormalXdxi*ShockV + WallNormalX*ShockVdxi) * Heta &
                  + WallNormalX * ShockV * HetaDxi

        ydxi = WallSYdxi + (WallNormalYdxi*ShockH + WallNormalY*ShockHdxi) * Heta &
                         + WallNormalY *  ShockH * HetaDxi;

        yxitau = (WallNormalYdxi*  ShockV + WallNormalY*ShockVdxi) * Heta &
                  +  WallNormalY *  ShockV * HetaDxi

        zdxi = WallSZdxi + (WallNormalZdxi*ShockH + WallNormalZ*ShockHdxi) * Heta &
                         + WallNormalZ *  ShockH * HetaDxi;

        zxitau = (WallNormalZdxi*  ShockV + WallNormalZ*ShockVdxi) * Heta &
                  +  WallNormalZ *  ShockV * HetaDxi
        
        xdeta = WallNormalX * ShockH * HetaDeta;
        ydeta = WallNormalY * ShockH * HetaDeta;
        zdeta = WallNormalZ * ShockH * HetaDeta;
        
        xetatau = WallNormalX * ShockV * HetaDeta;
        yetatau = WallNormalY * ShockV * HetaDeta;
        zetatau = WallNormalZ * ShockV * HetaDeta;

        xdzeta = WallSXdzeta + (WallNormalXdzeta* ShockH + WallNormalX*ShockHdzeta) * Heta &
                             + WallNormalX * ShockH * HetaDzeta;
        
        xzetatau = (WallNormalXdzeta*  ShockV + WallNormalX*ShockVdzeta) * Heta &
                      + WallNormalX *  ShockV * HetaDzeta;

        ydzeta = WallSYdzeta + (WallNormalYdzeta*ShockH + WallNormalY*ShockHdzeta) * Heta &
                             + WallNormalY * ShockH * HetaDzeta;

        yzetatau = (WallNormalYdzeta*  ShockV + WallNormalY*ShockVdzeta) * Heta &
                      + WallNormalY *  ShockV * HetaDzeta;
      
        zdzeta = WallSZdzeta + (WallNormalZdzeta* ShockH + WallNormalZ*ShockHdzeta) * Heta &
                             + WallNormalZ * ShockH * HetaDzeta;

        zzetatau = (WallNormalZdzeta*  ShockV + WallNormalZ*ShockVdzeta) * Heta &
                      + WallNormalZ *  ShockV * HetaDzeta;
        
        temp =    xdxi * (  ydeta* zdzeta - zdeta* ydzeta) &
             & + xdeta * ( ydzeta*   zdxi -  ydxi* zdzeta) &
             & +xdzeta * (   ydxi*  zdeta - ydeta*   zdxi)

        Jaco = 1.d0/temp;

        invJacodt = & 
             &    xxitau * (  ydeta* zdzeta - zdeta* ydzeta) &
             & + xetatau * ( ydzeta*   zdxi -  ydxi* zdzeta) &
             & +xzetatau * (   ydxi * zdeta - ydeta*   zdxi) &
             & +  xdxi * (  yetatau * zdzeta - zetatau* ydzeta) &
             & + xdeta * (  yzetatau*   zdxi -  yxitau* zdzeta) &
             & +xdzeta * (   yxitau *  zdeta - yetatau*   zdxi) &
             & +  xdxi * ( ydeta* zzetatau - zdeta* yzetatau) &
             & + xdeta * (ydzeta*   zxitau -  ydxi* zzetatau) &
             & +xdzeta * (  ydxi*  zetatau - ydeta*   zxitau)            

          xix=  Jaco * (ydeta*zdzeta - ydzeta* zdeta);
          xiy= -Jaco * (xdeta*zdzeta -  zdeta*xdzeta);
          xiz=  Jaco * (xdeta*ydzeta -  ydeta*xdzeta); 
         etax= -Jaco * ( ydxi*zdzeta -   zdxi*ydzeta);  
         etay=  Jaco * ( xdxi*zdzeta -   zdxi*xdzeta);
         etaz= -Jaco * ( xdxi*ydzeta -   ydxi*xdzeta);
        zetax=  Jaco * ( ydxi* zdeta -   zdxi* ydeta);  
        zetay= -Jaco * ( xdxi* zdeta -   zdxi* xdeta);
        zetaz=  Jaco * ( xdxi* ydeta -   ydxi* xdeta); 

          nablaxi=sqrt(  xix *  xix +  xiy *  xiy +  xiz *  xiz)

         nablaeta=sqrt( etax * etax + etay * etay + etaz * etaz)

        nablazeta=sqrt(zetax *zetax +zetay *zetay +zetaz *zetaz)

        ShockNormalX = etax / nablaeta
        ShockNormalY = etay / nablaeta
        ShockNormalZ = etaz / nablaeta

        ShockXtau = WallNormalX * ShockV
        ShockYtau = WallNormalY * ShockV
        ShockZtau = WallNormalZ * ShockV

        ShockXtaudxi = WallNormalXdxi * ShockV + WallNormalX * ShockVdxi
        ShockYtaudxi = WallNormalYdxi * ShockV + WallNormalY * ShockVdxi
        ShockZtaudxi =                           WallNormalZ * ShockVdxi
        
        ShockXtaudzeta = WallNormalX * ShockVdzeta
        ShockYtaudzeta = WallNormalY * ShockVdzeta
        ShockZtaudzeta = WallNormalZ * ShockVdzeta

        tempax = ShockYtaudzeta*zdxi + ShockZtaudxi*ydzeta - &
               & ShockYtaudxi  *zdzeta - ShockZtaudzeta*ydxi
        
        tempay = ShockXtaudxi  *zdzeta + ShockZtaudzeta*xdxi - &
               & ShockXtaudzeta*zdxi - ShockZtaudxi*xdzeta        
        
        tempaz = ShockXtaudzeta*ydxi + ShockYtaudxi*xdzeta - &
               & ShockXtaudxi*ydzeta - ShockYtaudzeta*xdxi
        
        detadt = -( etax*ShockXtau + etay*ShockYtau + etaz*ShockZtau)
   
        shockVN = -detadt       /nablaeta;
         ax_tau =    Jaco*tempax/nablaeta;
         ay_tau =    Jaco*tempay/nablaeta;        
         az_tau =    Jaco*tempaz/nablaeta;
    
        ! Obtain the variables behind the shock
        ! Velocity ahead of the shocks with respect to the moving shock normal
        Delta_Inf = U_inf*ShockNormalX + V_inf*ShockNormalY + W_inf*ShockNormalZ &
                 &- ShockV*(WallNormalX*ShockNormalX &
                          &+WallNormalY*ShockNormalY &
                          &+WallNormalZ*ShockNormalZ);
        ! Shock normal Mach number^2
          Masn2Inf = (Mach_Ref * Delta_Inf)**2.d0;
        ! Velocity normal to the shock after the shock
          DeltaShkUinf = 2.d0 * Delta_Inf / (Gamma + 1.d0) * (1.d0 / Masn2Inf - 1.d0);
          U_IJK = U_inf + DeltaShkUinf * ShockNormalX;
          V_IJK = V_inf + DeltaShkUinf * ShockNormalY;
          W_IJK = W_inf + DeltaShkUinf * ShockNormalZ;
          P_IJK = P_inf * (1.d0 + 2.d0 * Gamma / (Gamma + 1.d0) * (Masn2Inf - 1.d0));
        Rho_IJK = Rho_inf * (Gamma + 1.d0) * Masn2Inf / ((Gamma - 1.d0) * Masn2Inf + 2.d0);
          T_IJK = P_IJK * Gamma * Mach_Ref * Mach_Ref / Rho_IJK;
          C_IJK = sqrt(Gamma * P_IJK / Rho_IJK)

        ! Update the fluxes and new conservative variables
        ! invF_IJK(5),invG_IJK(5),invH_IJK(5),U0_IJK(5)
        ! Update the fluxes and new conservative variables
        ! invF_IJK(5),invG_IJK(5),invH_IJK(5),U0_IJK(5)
        !  CALL UtoCV_CPG(Rho_IJK,U_IJK,V_IJK,W_IJK,P_IJK,U0_IJK)
          U0_IJK(1)=Rho_IJK
          U0_IJK(2)=Rho_IJK*U_IJK
          U0_IJK(3)=Rho_IJK*V_IJK
          U0_IJK(4)=Rho_IJK*W_IJK
          U0_IJK(5)=P_IJK/(GAMMA-1.D0)+0.5D0*Rho_IJK*(U_IJK*U_IJK+V_IJK*V_IJK+W_IJK*W_IJK)

        !  CALL UtoFlux_CPG3D(Rho_IJK,U_IJK,V_IJK,W_IJK,P_IJK,invF_IJK,invG_IJK,invH_IJK)
          invF_IJK(1) = Rho_IJK * U_IJK
          invF_IJK(2) = Rho_IJK * U_IJK * U_IJK + P_IJK
          invF_IJK(3) = Rho_IJK * U_IJK * V_IJK
          invF_IJK(4) = Rho_IJK * U_IJK * W_IJK
          invF_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * U_IJK

          invG_IJK(1) = Rho_IJK * V_IJK
          invG_IJK(2) = Rho_IJK * V_IJK * U_IJK
          invG_IJK(3) = Rho_IJK * V_IJK * V_IJK + P_IJK
          invG_IJK(4) = Rho_IJK * V_IJK * W_IJK
          invG_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * V_IJK

          invH_IJK(1) = Rho_IJK * W_IJK
          invH_IJK(2) = Rho_IJK * W_IJK * U_IJK
          invH_IJK(3) = Rho_IJK * W_IJK * V_IJK
          invH_IJK(4) = Rho_IJK * W_IJK * W_IJK + P_IJK
          invH_IJK(5) = (P_IJK/(GAMMA - 1.D0) + 0.5d0 * Rho_IJK * (U_IJK*U_IJK + V_IJK*V_IJK + W_IJK*W_IJK) + P_IJK) * W_IJK

        ! Calculate the shock accelerations
              ! Calculate the local normal
        locUns = U_IJK * ShockNormalX + V_IJK * ShockNormalY + W_IJK * ShockNormalZ;
        LAM5   = locUns - ShockVN + C_IJK;

        L5(1) = 0.5d0 * (Gamma - 1.d0)*(U_IJK * U_IJK + V_IJK * V_IJK + W_IJK * W_IJK ) - locUns * C_IJK; 
        L5(2) =        -(Gamma - 1.d0)* U_IJK + ShockNormalX * C_IJK;
        L5(3) =        -(Gamma - 1.d0)* V_IJK + ShockNormalY * C_IJK;
        L5(4) =        -(Gamma - 1.d0)* W_IJK + ShockNormalZ * C_IJK;
        L5(5) =         (Gamma - 1.d0);
    
        SUM_L5_U  = 0.d0;
        SUM_L5_FG = 0.d0;
        
        DO iVar = 1,NumVar
          FG = LAM5 * dUcons0(iVar) + ax_tau * (invF_IJK(iVar) - FINV_INF(iVar)) &
                                    + ay_tau * (invG_IJK(iVar) - GINV_INF(iVar)) &
                                    + az_tau * (invH_IJK(iVar) - HINV_INF(iVar));
    
          SUM_L5_U  = SUM_L5_U  + L5(iVar) * (U0_IJK(iVar) - CV_INF(iVar));
          SUM_L5_FG = SUM_L5_FG + L5(iVar) * FG;                       
        ENDDO
        
        ! Calculate the ShockAc
        ShockAc = SUM_L5_FG/SUM_L5_U - (ax_tau * ShockXtau &
                                     &+ ay_tau * ShockYtau &
                                     &+ az_tau * ShockZtau);
    
        ShockAc = ShockAc/( WallNormalX * ShockNormalX &
                         &+ WallNormalY * ShockNormalY &
                         &+ WallNormalZ * ShockNormalZ);

end subroutine CalculateShockAc