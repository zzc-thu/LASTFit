subroutine calculate_ShockU0(U_inf,V_inf,W_inf,P_inf,Rho_inf,ShockV,Mach_Ref,Gamma,&
                           & ShockH,ShockHdxi,ShockHdzeta,Heta,HetaDxi,HetaDeta,HetaDzeta,&
                           & WallSXdxi,WallSYdxi,WallSZdxi,WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                           & WallNormalX,WallNormalY,WallNormalZ,&
                           & WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                           & WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,CV_1,CV_2,CV_3,CV_4,CV_5)
  implicit none
  ! Input variables ShockH, ShockHdxi, ShockHdzeta
  ! Output variables: CV_1,CV_2,CV_3,CV_4,CV_5
  real( kind = 8 )::Delta_Inf,Masn2Inf,DeltaShkUinf
  real( kind = 8 )::U_inf,V_inf,W_inf,P_inf,Rho_inf
  real( kind = 8 )::ShockNormalX,ShockNormalY,ShockNormalZ
  real( kind = 8 )::WallNormalX,WallNormalY,WallNormalZ
  real( kind = 8 )::U,V,W,P,Rho,T,Mach_Ref,Gamma
  real( kind = 8 )::etax,etay,etaz,nablaeta 
  real( kind = 8 )::xdxi,ydxi,zdxi,xdzeta,ydzeta,zdzeta,Jaco
  real( kind = 8 )::xdeta,ydeta,zdeta
  real( kind = 8 )::Heta,HetaDxi,HetaDeta,HetaDzeta
  real( kind = 8 )::WallSXdxi,WallSYdxi,WallSZdxi
  real( kind = 8 )::WallSXdzeta,WallSYdzeta,WallSZdzeta
  real( kind = 8 )::WallNormalXdxi,WallNormalYdxi,WallNormalZdxi
  real( kind = 8 )::WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta
  real( kind = 8 )::temp
  real( kind = 8 ),intent(in)::ShockH,ShockHdxi,ShockHdzeta,ShockV
  real( kind = 8 ),intent(out)::CV_1,CV_2,CV_3,CV_4,CV_5
        
        ! x = WallSx + WallnormalX * ShockH * Heta
        ! y = WallSy + WallnormalY * ShockH * Heta
        ! z = WallSz + WallnormalZ * ShockH * Heta 
       
        xdxi = WallSXdxi + (WallNormalXdxi*ShockH + WallNormalX*ShockHdxi) * Heta &
                         + WallNormalX * ShockH * HetaDxi;

        ydxi = WallSYdxi + (WallNormalYdxi*ShockH + WallNormalY*ShockHdxi) * Heta &
                         + WallNormalY *  ShockH * HetaDxi;

        zdxi = WallSZdxi + (WallNormalZdxi*ShockH + WallNormalZ*ShockHdxi) * Heta &
                         + WallNormalZ *  ShockH * HetaDxi;
        
        xdeta = WallNormalX * ShockH * HetaDeta;
        ydeta = WallNormalY * ShockH * HetaDeta;
        zdeta = WallNormalZ * ShockH * HetaDeta;

        xdzeta = WallSXdzeta + (WallNormalXdzeta* ShockH + WallNormalX*ShockHdzeta) * Heta &
                             + WallNormalX * ShockH * HetaDzeta;

        ydzeta = WallSYdzeta + (WallNormalYdzeta*ShockH + WallNormalY*ShockHdzeta) * Heta &
                             + WallNormalY * ShockH * HetaDzeta;
      
        zdzeta = WallSZdzeta + (WallNormalZdzeta* ShockH + WallNormalZ*ShockHdzeta) * Heta &
                             + WallNormalZ * ShockH * HetaDzeta;
        
        temp =    xdxi * (  ydeta* zdzeta - zdeta* ydzeta) &
             & + xdeta * ( ydzeta*   zdxi -  ydxi* zdzeta) &
             & +xdzeta * (   ydxi * zdeta - ydeta*   zdxi)

        Jaco = 1.d0/temp;

    etax= -Jaco * ( ydxi*zdzeta -   zdxi*ydzeta);  
    etay=  Jaco * ( xdxi*zdzeta -   zdxi*xdzeta);
    etaz= -Jaco * ( xdxi*ydzeta -   ydxi*xdzeta);

    nablaeta=sqrt( etax * etax + etay * etay + etaz * etaz)
  ! Given ShockNormalX, ShockNormalY, ShockNormalZ 
    
    ShockNormalX = etax / nablaeta
    ShockNormalY = etay / nablaeta
    ShockNormalZ = etaz / nablaeta

  ! Velocity ahead of the shocks
    Delta_Inf = U_inf*ShockNormalX+V_inf*ShockNormalY+W_inf*ShockNormalZ &
             &- ShockV*(WallNormalX*ShockNormalX+WallNormalY*ShockNormalY+WallNormalZ*ShockNormalZ);
    ! Shock normal Mach number^2
      Masn2Inf = (Mach_Ref * Delta_Inf)**2.d0;
    ! Velocity normal to the shock after the shock
      DeltaShkUinf = 2.d0 * Delta_Inf / (Gamma + 1.d0) * (1.d0 / Masn2Inf - 1.d0);
      
      U = U_inf + DeltaShkUinf * ShockNormalX;

      V = V_inf + DeltaShkUinf * ShockNormalY;

      W = W_inf + DeltaShkUinf * ShockNormalZ;

      P = P_inf * (1.d0 + 2.d0 * Gamma / (Gamma + 1.d0) * (Masn2Inf - 1.d0));

    Rho = Rho_inf * (Gamma + 1.d0) * Masn2Inf / ((Gamma - 1.d0) * Masn2Inf + 2.d0);

      T = P * Gamma * Mach_Ref * Mach_Ref / Rho;
     
      CV_1 = Rho;
      CV_2 = Rho * U;
      CV_3 = Rho * V;
      CV_4 = Rho * W;
      CV_5 = P / (Gamma - 1.d0) + 0.5d0 * Rho * (U * U + V * V + W * W);
      
end subroutine calculate_ShockU0