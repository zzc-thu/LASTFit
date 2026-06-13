subroutine Calculate_Jaco
  use SF_Constant,  only: ik,rk,overLAP
  use SF_CFD_Global,only: nx_local,nz_local,Ny
  use MPI_GLOBAL,   only: MyID

  ! some arrays are defined in SF_CFD_Global
  use SF_CFD_Global,only: X_grid,Y_grid,Z_grid, &
                          dxdxi,dxdeta,dxdzeta,dxidx,dxidy,dxidz,nablaxi, &
                          dydxi,dydeta,dydzeta,detadx,detady,detadz,nablaeta, &
                          dzdxi,dzdeta,dzdzeta,dzetadx,dzetady,dzetadz,Jaco,nablazeta, &
                          WallSX,WallSY,WallSZ,WallSXdxi,WallSYdxi,WallSZdxi,&
                          WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                          WallNormalX,WallNormalY,WallNormalZ,&
                          WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                          WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&
                          ShockH,ShockHdxi,ShockHdzeta,ShockNormalX,ShockNormalY,ShockNormalZ,&
                          ShockV,ShockVdxi,ShockVdzeta,Heta,HetaDxi,HetaDeta,HetaDzeta,&
                          ShockXtau,ShockYtau,ShockZtau,dxidt,detadt,dzetadt,&
                          invJacodt,shockVN,ax_tau,ay_tau,az_tau,&
                          ShockNormalX_steady,ShockNormalY_steady,ShockNormalZ_steady,AnalysisType,shockXtau_steady
  ! some subroutines
  use MPI_GLOBAL,    only: Parallel_Exchange_Surface,Parallel_Exchange
  use OutputParaView,only: output_jacobian
  use FD5_Order,     only: Cal_Deri_SurfDxi_5th_ce,Cal_Deri_SurfDzeta_5th_ce,Cal_Deri_Deta_5th_ce
  implicit none
  ! local variables
  integer( kind = ik ) :: ic, jc, kc
  real( kind = rk ) :: temp,tempax,tempay,tempaz
  real( kind = rk ) :: X_tau,Y_tau,Z_tau
  real( kind = rk ) :: ShockXtaudxi,ShockYtaudxi,ShockZtaudxi
  real( kind = rk ) :: ShockXtaudzeta,ShockYtaudzeta,ShockZtaudzeta

  real( kind = rk ) :: xxitau,xetatau,xzetatau
  real( kind = rk ) :: yxitau,yetatau,yzetatau
  real( kind = rk ) :: zxitau,zetatau,zzetatau
  
  ! remeshing
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Y_grid(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Z_grid(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
  
  call Parallel_Exchange(X_grid)
  call Parallel_Exchange(Y_grid)
  call Parallel_Exchange(Z_grid)

  ! calculate the ShockHdxi and ShockHdzeta

  call Parallel_Exchange_Surface(ShockH)
  call Cal_Deri_SurfDxi_5th_ce(ShockHdxi,    ShockH)
  call Cal_Deri_SurfDzeta_5th_ce(ShockHdzeta,ShockH)
  call Parallel_Exchange_Surface(ShockHdxi)
  call Parallel_Exchange_Surface(ShockHdzeta)

  call Parallel_Exchange_Surface(ShockV)
  call Cal_Deri_SurfDxi_5th_ce(ShockVdxi,    ShockV)
  call Cal_Deri_SurfDzeta_5th_ce(ShockVdzeta,ShockV)
  call Parallel_Exchange_Surface(ShockVdxi)
  call Parallel_Exchange_Surface(ShockVdzeta)

  ! Here we do not use finite difference to calculate
  ! the dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta
  ! dzdxi, dzdeta, dzdzeta, because we have already known
  ! the analytical expression of surface grid and the 
  ! mapping function.
  ! only the derivative with respect to shock height is
  ! calculated by finite difference.
 
 
  do kc = 1-overLAP, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP, nx_local+overLAP
        dxdxi(ic,jc,kc) = WallSXdxi(ic,kc) + (WallNormalXdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalX(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        xxitau =                             (WallNormalXdxi(ic,kc)*ShockV(ic,kc) + &
                                              & WallNormalX(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                              + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        dydxi(ic,jc,kc) = WallSYdxi(ic,kc) + (WallNormalYdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        yxitau =                             (WallNormalYdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        dzdxi(ic,jc,kc) = WallSZdxi(ic,kc) + (WallNormalZdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        zxitau =                             (WallNormalZdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
        
        dxdeta(ic,jc,kc) = WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dydeta(ic,jc,kc) = WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dzdeta(ic,jc,kc) = WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        
        xetatau = WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        yetatau = WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        zetatau = WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);

        dxdzeta(ic,jc,kc) = WallSXdzeta(ic,kc) + (WallNormalXdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalX(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        xzetatau =                               (WallNormalXdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalX(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        dydzeta(ic,jc,kc) = WallSYdzeta(ic,kc) + (WallNormalYdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalY(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetatau =                              (WallNormalYdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalY(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
      
        dzdzeta(ic,jc,kc) = WallSZdzeta(ic,kc) + (WallNormalZdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalZ(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetatau =                              (WallNormalZdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalZ(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
        
        temp =    dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc))

           Jaco(ic,jc,kc)= 1.d0/temp;
        ! d(1/J)/dt = d(temp)/dt 
        invJacodt(ic,jc,kc) = & 
             &    xxitau * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + xetatau * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +xzetatau * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * (  yetatau * dzdzeta(ic,jc,kc) - zetatau* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * (  yzetatau*   dzdxi(ic,jc,kc) -  yxitau* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   yxitau *  dzdeta(ic,jc,kc) - yetatau*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * ( dydeta(ic,jc,kc)* zzetatau - dzdeta(ic,jc,kc)* yzetatau) &
             & + dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)*   zxitau -  dydxi(ic,jc,kc)* zzetatau) &
             & +dxdzeta(ic,jc,kc) * (  dydxi(ic,jc,kc)*  zetatau - dydeta(ic,jc,kc)*   zxitau)            

          dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
          dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
          dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc)); 
         detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));  
         detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
         detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));  
        dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
        dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc)); 

          nablaxi(ic,jc,kc)=sqrt(  dxidx(ic,jc,kc) *  dxidx(ic,jc,kc) &
                               &+  dxidy(ic,jc,kc) *  dxidy(ic,jc,kc) &
                               &+  dxidz(ic,jc,kc) *  dxidz(ic,jc,kc))

         nablaeta(ic,jc,kc)=sqrt( detadx(ic,jc,kc) * detadx(ic,jc,kc) &
                               &+ detady(ic,jc,kc) * detady(ic,jc,kc) &
                               &+ detadz(ic,jc,kc) * detadz(ic,jc,kc))

        nablazeta(ic,jc,kc)=sqrt(dzetadx(ic,jc,kc) *dzetadx(ic,jc,kc) &
                               &+dzetady(ic,jc,kc) *dzetady(ic,jc,kc) &
                               &+dzetadz(ic,jc,kc) *dzetadz(ic,jc,kc))

      end do
    end do
  end do   
 
  ! calculate the Shock Normal directions
  ! only along the jc = Ny surface
  
  
  jc = Ny
  do kc = 1-overLAP, nz_local+overLAP
    do ic = 1-overLAP, nx_local+overLAP
      ShockNormalX(ic,kc) = detadx(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalY(ic,kc) = detady(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalZ(ic,kc) = detadz(ic,jc,kc) / nablaeta(ic,jc,kc)
    enddo
  enddo

  if (AnalysisType == 3) then    
      ShockNormalX_steady = ShockNormalX
      ShockNormalY_steady = ShockNormalY
      ShockNormalZ_steady = ShockNormalZ
  endif
  

  do kc = 1-overLAP,nz_local+overLAP
    do ic = 1-overLAP,nx_local+overLAP
      do jc = 1,Ny
       X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       !Z_tau = 0.d0;

        dxidt(ic,jc,kc) =  -(  dxidx(ic,jc,kc)*X_tau +  dxidy(ic,jc,kc)*Y_tau +  dxidz(ic,jc,kc)*Z_tau);
         !dxidt(ic,jc,kc) =   0.d0;
        detadt(ic,jc,kc) =  -( detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

        dzetadt(ic,jc,kc) =  -(dzetadx(ic,jc,kc)*X_tau +dzetady(ic,jc,kc)*Y_tau +dzetadz(ic,jc,kc)*Z_tau);
       !dzetadt(ic,jc,kc) =   0.d0;
       
       !detadt_invJ(ic,jc,kc) = detadt(ic,jc,kc)/Jaco(ic,jc,kc);
      enddo
      ShockXtau(ic,kc) = WallNormalX(ic,kc) * ShockV(ic,kc)
      ShockYtau(ic,kc) = WallNormalY(ic,kc) * ShockV(ic,kc)
      ShockZtau(ic,kc) = WallNormalZ(ic,kc) * ShockV(ic,kc)
      !ShockZtau(ic,kc) = 0.d0
    enddo
  enddo

  if (AnalysisType == 3) then    
      ShockXtau_steady = ShockXtau
  endif
  ! version 1
  !CALL Cal_Deri_Deta_5th_ce(detadt_invJdeta,detadt_invJ)
  !invJacodt = - detadt_invJdeta;

  jc = Ny
   do ic = 1-overLAP,nx_local+overLAP
    do kc = 1-overLAP,nz_local+overLAP
      ShockXtaudxi = WallNormalXdxi(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdxi(ic,kc)
      ShockYtaudxi = WallNormalYdxi(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdxi(ic,kc)
      ShockZtaudxi =                                         WallNormalZ(ic,kc) * ShockVdxi(ic,kc)

      ShockXtaudzeta = WallNormalX(ic,kc) * ShockVdzeta(ic,kc)
      ShockYtaudzeta = WallNormalY(ic,kc) * ShockVdzeta(ic,kc)
      ShockZtaudzeta = WallNormalZ(ic,kc) * ShockVdzeta(ic,kc)

      tempax = ShockYtaudzeta*dzdxi(ic,jc,kc) + ShockZtaudxi*dydzeta(ic,jc,kc) - &
             & ShockYtaudxi  *dzdzeta(ic,jc,kc) - ShockZtaudzeta*dydxi(ic,jc,kc)
      
      tempay = ShockXtaudxi  *dzdzeta(ic,jc,kc) + ShockZtaudzeta*dxdxi(ic,jc,kc) - &
             & ShockXtaudzeta*dzdxi(ic,jc,kc) - ShockZtaudxi*dxdzeta(ic,jc,kc)        
      
      tempaz = ShockXtaudzeta*dydxi(ic,jc,kc) + ShockYtaudxi*dxdzeta(ic,jc,kc) - &
             & ShockXtaudxi*dydzeta(ic,jc,kc) - ShockYtaudzeta*dxdxi(ic,jc,kc)         

     shockVN(ic,kc) = -detadt(ic,jc,kc)       /nablaeta(ic,jc,kc);
      ax_tau(ic,kc) =    Jaco(ic,jc,kc)*tempax/nablaeta(ic,jc,kc);
      ay_tau(ic,kc) =    Jaco(ic,jc,kc)*tempay/nablaeta(ic,jc,kc);        
      az_tau(ic,kc) =    Jaco(ic,jc,kc)*tempaz/nablaeta(ic,jc,kc);
    enddo
   enddo
   

   
  call Parallel_Exchange(Jaco)
  call Parallel_Exchange(detadx)
  call Parallel_Exchange(detady)   
  call Parallel_Exchange(detadz)
  call Parallel_Exchange(dxidx)
  call Parallel_Exchange(dxidy)
  call Parallel_Exchange(dxidz)
  call Parallel_Exchange(dzetadx)
  call Parallel_Exchange(dzetady)
  call Parallel_Exchange(dzetadz)
  call Parallel_Exchange(dxdxi)
  call Parallel_Exchange(dydxi)
  call Parallel_Exchange(dzdxi)
  call Parallel_Exchange(dxdzeta)
  call Parallel_Exchange(dydzeta)
  call Parallel_Exchange(dzdzeta)
  call Parallel_Exchange(dxdeta)
  call Parallel_Exchange(dydeta)
  call Parallel_Exchange(dzdeta)

!write(*,*)ShockNormalX(:,1)
! We need to check the Jacobian to see if every thing is correct
!#ifdef DEBUG
!  if(MyID == 0) then
!    write(*,*)"You are debugging the Jacobian"
!    write(*,*)"Output the Jacobian info for testing"
!  endif
!  ! We need to output the Jacobian to see if every thing is correct
!  call output_jacobian
!#endif

end subroutine Calculate_Jaco

subroutine Calculate_Jaco_Implicit
  ! This is the subroutine that calculate the derivatives with respect to the shock Height
  ! This is used in the implicit method
  use SF_Constant,  only: ik,rk,overLAP
  use SF_CFD_Global,only: nx_local,nz_local,Ny
  use MPI_GLOBAL,   only: MyID

  ! some arrays are defined in SF_CFD_Global
  use SF_CFD_Global,only: X_grid,Y_grid,Z_grid, &
                          dxdxi,dxdeta,dxdzeta,dxidx,dxidy,dxidz,nablaxi, &
                          dydxi,dydeta,dydzeta,detadx,detady,detadz,nablaeta, &
                          dzdxi,dzdeta,dzdzeta,dzetadx,dzetady,dzetadz,Jaco,nablazeta, &
                          WallSX,WallSY,WallSZ,WallSXdxi,WallSYdxi,WallSZdxi,&
                          WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                          WallNormalX,WallNormalY,WallNormalZ,&
                          WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                          WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&
                          ShockH,ShockHdxi,ShockHdzeta,ShockNormalX,ShockNormalY,ShockNormalZ,&
                          ShockV,ShockVdxi,ShockVdzeta,Heta,HetaDxi,HetaDeta,HetaDzeta,&
                          ShockXtau,ShockYtau,ShockZtau,dxidt,detadt,dzetadt,&
                          invJacodt,shockVN,ax_tau,ay_tau,az_tau,DxixJDH,DxiyJDH,DxizJDH,&
                          DetaxJDH,DetayJDH,DetazJDH,DzetaxJDH,DzetayJDH,DzetazJDH,DJDH,DetatDH,DetatJDH,&
                          DetaxJDHxi,  DetayJDHxi,  DetazJDHxi,  DetatJDHxi,DetatJDV,&
                          DetaxJDHzeta,DetayJDHzeta,DetazJDHzeta,DetatJDHzeta,&
                          DzetaxJDHxi,   DzetayJDHxi,   DzetazJDHxi,&
                          DxixJDHzeta, DxiyJDHzeta, DxizJDHzeta,DJDHxi,DJDHzeta
  ! some subroutines
  use MPI_GLOBAL,    only: Parallel_Exchange_Surface,Parallel_Exchange
  use OutputParaView,only: output_jacobian
  use FD5_Order,     only: Cal_Deri_SurfDxi_5th_ce,Cal_Deri_SurfDzeta_5th_ce,Cal_Deri_Deta_5th_ce
  implicit none
  ! local variables
  integer( kind = ik ) :: ic, jc, kc
  real( kind = rk ) :: temp,tempax,tempay,tempaz
  real( kind = rk ) :: X_tau,Y_tau,Z_tau
  real( kind = rk ) :: ShockXtaudxi,ShockYtaudxi,ShockZtaudxi
  real( kind = rk ) :: ShockXtaudzeta,ShockYtaudzeta,ShockZtaudzeta

  real( kind = rk ) :: xxiH,xetaH,xzetaH
  real( kind = rk ) :: yxiH,yetaH,yzetaH
  real( kind = rk ) :: zxiH,zetaH,zzetaH
  
  real( kind = rk ) :: xxiHxi,yxiHxi,zxiHxi
  real( kind = rk ) :: xzetaHzeta,yzetaHzeta,zzetaHzeta
  real( kind = rk ) :: xxitau,xetatau,xzetatau
  real( kind = rk ) :: yxitau,yetatau,yzetatau
  real( kind = rk ) :: zxitau,zetatau,zzetatau
  
    ! remeshing
    do kc = 1-overLAP,nz_local+overLAP
      do jc = 1,Ny
        do ic = 1-overLAP,nx_local+overLAP
          X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
          Y_grid(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
          Z_grid(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        end do
      end do
    end do
    ! 
    call Parallel_Exchange(X_grid)
    call Parallel_Exchange(Y_grid)
    call Parallel_Exchange(Z_grid)
  
    ! calculate the ShockHdxi and ShockHdzeta
    call Parallel_Exchange_Surface(ShockH)
    call Cal_Deri_SurfDxi_5th_ce(ShockHdxi,    ShockH)
    call Cal_Deri_SurfDzeta_5th_ce(ShockHdzeta,ShockH)
    call Parallel_Exchange_Surface(ShockHdxi)
    call Parallel_Exchange_Surface(ShockHdzeta)
  
    call Parallel_Exchange_Surface(ShockV)
    call Cal_Deri_SurfDxi_5th_ce(ShockVdxi,    ShockV)
    call Cal_Deri_SurfDzeta_5th_ce(ShockVdzeta,ShockV)
    call Parallel_Exchange_Surface(ShockVdxi)
    call Parallel_Exchange_Surface(ShockVdzeta)

  ! Here we do not use finite difference to calculate
  ! the dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta
  ! dzdxi, dzdeta, dzdzeta, because we have already known
  ! the analytical expression of surface grid and the 
  ! mapping function.
  ! only the derivative with respect to shock height is
  ! calculated by finite difference.
  
    
  do kc = 1-overLap, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP, nx_local+overLAP
        dxdxi(ic,jc,kc) = WallSXdxi(ic,kc) + (WallNormalXdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalX(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        xxitau =                             (WallNormalXdxi(ic,kc)*ShockV(ic,kc) + &
                                              & WallNormalX(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                              + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        xxiH = WallNormalXdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalX(ic,kc) * HetaDxi(ic,jc,kc); 

        xxiHxi = WallNormalX(ic,kc) * Heta(ic,jc,kc);

        dydxi(ic,jc,kc) = WallSYdxi(ic,kc) + (WallNormalYdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
        
        yxitau =                             (WallNormalYdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
                                           
        yxiH = WallNormalYdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalY(ic,kc) * HetaDxi(ic,jc,kc);

        yxiHxi = WallNormalY(ic,kc) * Heta(ic,jc,kc);
        
        dzdxi(ic,jc,kc) = WallSZdxi(ic,kc) + (WallNormalZdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
        
        zxitau =                             (WallNormalZdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        zxiH = WallNormalZdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalZ(ic,kc) * HetaDxi(ic,jc,kc);

        zxiHxi = WallNormalZ(ic,kc) * Heta(ic,jc,kc);
        
        dxdeta(ic,jc,kc) = WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dydeta(ic,jc,kc) = WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dzdeta(ic,jc,kc) = WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);

        xetaH = WallNormalX(ic,kc) * HetaDeta(ic,jc,kc);
        yetaH = WallNormalY(ic,kc) * HetaDeta(ic,jc,kc);
        zetaH = WallNormalZ(ic,kc) * HetaDeta(ic,jc,kc);

        xetatau = WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        yetatau = WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        zetatau = WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);

        dxdzeta(ic,jc,kc) = WallSXdzeta(ic,kc) + (WallNormalXdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalX(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        xzetatau =                               (WallNormalXdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalX(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        xzetaH = WallNormalXdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalX(ic,kc) * HetaDzeta(ic,jc,kc);

        xzetaHzeta = WallNormalX(ic,kc) * Heta(ic,jc,kc);

        dydzeta(ic,jc,kc) = WallSYdzeta(ic,kc) + (WallNormalYdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalY(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        yzetatau =                              (WallNormalYdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalY(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                                + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetaH = WallNormalYdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalY(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetaHzeta = WallNormalY(ic,kc) * Heta(ic,jc,kc);

        dzdzeta(ic,jc,kc) = WallSZdzeta(ic,kc) + (WallNormalZdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalZ(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        zzetatau =                              (WallNormalZdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalZ(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                                + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetaH = WallNormalZdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalZ(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetaHzeta = WallNormalZ(ic,kc) * Heta(ic,jc,kc);

        temp =    dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc))

        Jaco(ic,jc,kc)= 1.d0/temp;          

          dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
          dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
          dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc)); 
         detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));  
         detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
         detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));  
        dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
        dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));

          nablaxi(ic,jc,kc)=sqrt(  dxidx(ic,jc,kc) *  dxidx(ic,jc,kc) &
                               &+  dxidy(ic,jc,kc) *  dxidy(ic,jc,kc) &
                               &+  dxidz(ic,jc,kc) *  dxidz(ic,jc,kc))

         nablaeta(ic,jc,kc)=sqrt( detadx(ic,jc,kc) * detadx(ic,jc,kc) &
                               &+ detady(ic,jc,kc) * detady(ic,jc,kc) &
                               &+ detadz(ic,jc,kc) * detadz(ic,jc,kc))

        nablazeta(ic,jc,kc)=sqrt(dzetadx(ic,jc,kc) *dzetadx(ic,jc,kc) &
                               &+dzetady(ic,jc,kc) *dzetady(ic,jc,kc) &
                               &+dzetadz(ic,jc,kc) *dzetadz(ic,jc,kc))

        ! d(1/J)/dt = d(temp)/dt 
        invJacodt(ic,jc,kc) = & 
             &    xxitau * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + xetatau * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +xzetatau * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * (  yetatau * dzdzeta(ic,jc,kc) - zetatau* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * (  yzetatau*   dzdxi(ic,jc,kc) -  yxitau* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   yxitau *  dzdeta(ic,jc,kc) - yetatau*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * ( dydeta(ic,jc,kc)* zzetatau - dzdeta(ic,jc,kc)* yzetatau) &
             & + dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)*   zxitau -  dydxi(ic,jc,kc)* zzetatau) &
             & +dxdzeta(ic,jc,kc) * (  dydxi(ic,jc,kc)*  zetatau - dydeta(ic,jc,kc)*   zxitau)     
        
        ! This part calculate the   d(xix/J)/dH,   d(xiy/J)/dH,   d(xiz/J)/dH
        !                          d(etax/J)/dH,  d(etay/J)/dH,  d(etaz/J)/dH
        !                         d(zetax/J)/dH, d(zetay/J)/dH, d(zetaz/J)/dH 
        !dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
         DxixJDH(ic,jc,kc) = (  yetaH*dzdzeta(ic,jc,kc) +  dydeta(ic,jc,kc)*  zzetaH &
                           -   yzetaH* dzdeta(ic,jc,kc) - dydzeta(ic,jc,kc)*   zetaH)
        
        !dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
         DxiyJDH(ic,jc,kc) =-(  xetaH*dzdzeta(ic,jc,kc) +  dxdeta(ic,jc,kc)*  zzetaH &
                           -    zetaH*dxdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*  xzetaH)
        
        !dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
         DxizJDH(ic,jc,kc) = (  xetaH*dydzeta(ic,jc,kc) +  dxdeta(ic,jc,kc)*  yzetaH &
                           -    yetaH*dxdzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*  xzetaH)
        
        !detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));
        DetaxJDH(ic,jc,kc) =-(   yxiH*dzdzeta(ic,jc,kc) +   dydxi(ic,jc,kc)*  zzetaH &
                           -     zxiH*dydzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*  yzetaH)
        
        !detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        DetayJDH(ic,jc,kc) = (   xxiH*dzdzeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*  zzetaH &
                           -     zxiH*dxdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*  xzetaH)
        
        !detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        DetazJDH(ic,jc,kc) =-(   xxiH*dydzeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*  yzetaH &
                           -     yxiH*dxdzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*  xzetaH)
       
       !dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));
       DzetaxJDH(ic,jc,kc) = (   yxiH* dzdeta(ic,jc,kc) +   dydxi(ic,jc,kc)*   zetaH &
                           -     zxiH* dydeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*   yetaH)
       
       !dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
       DzetayJDH(ic,jc,kc) =-(   xxiH* dzdeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*   zetaH &
                           -     zxiH* dxdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*   xetaH) 
       
       !dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));
       DzetazJDH(ic,jc,kc) = (   xxiH* dydeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*   yetaH &
                           -     yxiH* dxdeta(ic,jc,kc) -   dydxi(ic,jc,kc)*   xetaH)
            
            ! J = 1/temp, DJ/DH = -Dtemp/DH/(temp*temp)
            !             DJ/DHxi = -Dtemp/DHxi/(temp*temp)
            !             DJ/DHzeta = -Dtemp/DHzeta/(temp*temp)
            DJDH(ic,jc,kc) = xxiH * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
                        & + xetaH * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
                        & +xzetaH * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
                        & +  dxdxi(ic,jc,kc) * (  yetaH* dzdzeta(ic,jc,kc) - zetaH* dydzeta(ic,jc,kc)) &
                        & + dxdeta(ic,jc,kc) * ( yzetaH*   dzdxi(ic,jc,kc) -  yxiH* dzdzeta(ic,jc,kc)) &
                        & +dxdzeta(ic,jc,kc) * (   yxiH * dzdeta(ic,jc,kc) - yetaH*   dzdxi(ic,jc,kc)) &
                        & +  dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* zzetaH - dzdeta(ic,jc,kc)* yzetaH) &
                        & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   zxiH -  dydxi(ic,jc,kc)* zzetaH) &
                        & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * zetaH - dydeta(ic,jc,kc)*   zxiH)
            
            DJDH(ic,jc,kc) = - DJDH(ic,jc,kc) / (temp * temp);

            DJDHxi(ic,jc,kc) = xxiHxi * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
                          & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*  zxiHxi -  yxiHxi* dzdzeta(ic,jc,kc)) &
                          & +dxdzeta(ic,jc,kc) * (   yxiHxi * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   zxiHxi)
            
            DJDHxi(ic,jc,kc) = - DJDHxi(ic,jc,kc) / (temp * temp);

            DJDHzeta(ic,jc,kc) = dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* zzetaHzeta - dzdeta(ic,jc,kc)* yzetaHzeta) &
                          & +  dxdeta(ic,jc,kc) * ( yzetaHzeta *   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* zzetaHzeta) &
                          & +xzetaHzeta * ( dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*dzdxi(ic,jc,kc)) 
            
            DJDHzeta(ic,jc,kc) = - DJDHzeta(ic,jc,kc) / (temp * temp);

            X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
            Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
            Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);

            !
            detadt(ic,jc,kc) = -(detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

            DetatJDH(ic,jc,kc) = -(DetaxJDH(ic,jc,kc)*X_tau + DetayJDH(ic,jc,kc)*Y_tau + DetazJDH(ic,jc,kc)*Z_tau);

            DetatJDV(ic,jc,kc) = -(detadx(ic,jc,kc) * WallNormalX(ic,kc) * Heta(ic,jc,kc) +&
                                 & detady(ic,jc,kc) * WallNormalY(ic,kc) * Heta(ic,jc,kc) +&
                                 & detadz(ic,jc,kc) * WallNormalZ(ic,kc) * Heta(ic,jc,kc))/jaco(ic,jc,kc);

           !dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
           !dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
           !dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
            DxixJDHzeta(ic,jc,kc) = ( dydeta(ic,jc,kc)*zzetaHzeta - yzetaHzeta* dzdeta(ic,jc,kc) )

            DxiyJDHzeta(ic,jc,kc) =-( dxdeta(ic,jc,kc)*zzetaHzeta -  dzdeta(ic,jc,kc)*xzetaHzeta )

            DxizJDHzeta(ic,jc,kc) = ( dxdeta(ic,jc,kc)*yzetaHzeta -  dydeta(ic,jc,kc)*xzetaHzeta )

           !detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));
           !detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
           !detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
              DetaxJDHxi(ic,jc,kc) =-(   yxiHxi * dzdzeta(ic,jc,kc) - zxiHxi * dydzeta(ic,jc,kc))
            DetaxJDHzeta(ic,jc,kc) =-(   dydxi(ic,jc,kc)*  zzetaHzeta - dzdxi(ic,jc,kc) * yzetaHzeta)
            
              DetayJDHxi(ic,jc,kc) = (   xxiHxi * dzdzeta(ic,jc,kc) - zxiHxi * dxdzeta(ic,jc,kc))
            DetayJDHzeta(ic,jc,kc) = (   dxdxi(ic,jc,kc)*  zzetaHzeta - dzdxi(ic,jc,kc) * xzetaHzeta)
            
              DetazJDHxi(ic,jc,kc) =-(   xxiHxi * dydzeta(ic,jc,kc) - yxiHxi * dxdzeta(ic,jc,kc))
            DetazJDHzeta(ic,jc,kc) =-(   dxdxi(ic,jc,kc)*  yzetaHzeta - dydxi(ic,jc,kc) * xzetaHzeta)
           
           !dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));
           !dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
           !dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));
              DzetaxJDHxi(ic,jc,kc) = (   yxiHxi * dzdeta(ic,jc,kc) - zxiHxi * dydeta(ic,jc,kc))

              DzetayJDHxi(ic,jc,kc) =-(   xxiHxi * dzdeta(ic,jc,kc) + zxiHxi * dxdeta(ic,jc,kc))
            
              DzetazJDHxi(ic,jc,kc) = (   xxiHxi * dydeta(ic,jc,kc) + yxiHxi * dxdeta(ic,jc,kc))
           
           !etat related variables 
              DetatJDHxi(ic,jc,kc) = -(  DetaxJDHxi(ic,jc,kc)*X_tau +   DetayJDHxi(ic,jc,kc)*Y_tau +   DetazJDHxi(ic,jc,kc)*Z_tau);

            DetatJDHzeta(ic,jc,kc) = -(DetaxJDHzeta(ic,jc,kc)*X_tau + DetayJDHzeta(ic,jc,kc)*Y_tau + DetazJDHzeta(ic,jc,kc)*Z_tau);

           !Then the etat related variables with respect to V and Vxi, Vzeta
           ! Vxi and Vzeta term are only appear at the shock boundary conditions
          
      end do
    end do
  end do
! calculate the Shock Normal directions
  ! only along the jc = Ny surface
  jc = Ny
  do kc = 1-overLAP, nz_local+overLAP
    do ic = 1-overLAP, nx_local+overLAP
      ShockNormalX(ic,kc) = detadx(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalY(ic,kc) = detady(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalZ(ic,kc) = detadz(ic,jc,kc) / nablaeta(ic,jc,kc)
    enddo
  enddo

  do kc = 1-overLAP,nz_local+overLAP
    do ic = 1-overLAP,nx_local+overLAP
      do jc = 1,Ny
       X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);

        dxidt(ic,jc,kc) =  -(  dxidx(ic,jc,kc)*X_tau +  dxidy(ic,jc,kc)*Y_tau +  dxidz(ic,jc,kc)*Z_tau);
         !dxidt(ic,jc,kc) =   0.d0;

        detadt(ic,jc,kc) =  -( detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

        dzetadt(ic,jc,kc) =  -(dzetadx(ic,jc,kc)*X_tau +dzetady(ic,jc,kc)*Y_tau +dzetadz(ic,jc,kc)*Z_tau);
       !dzetadt(ic,jc,kc) =   0.d0;
       
       !detadt_invJ(ic,jc,kc) = detadt(ic,jc,kc)/Jaco(ic,jc,kc);
      enddo
      ShockXtau(ic,kc) = WallNormalX(ic,kc) * ShockV(ic,kc)
      ShockYtau(ic,kc) = WallNormalY(ic,kc) * ShockV(ic,kc)
      ShockZtau(ic,kc) = WallNormalZ(ic,kc) * ShockV(ic,kc)
    enddo
  enddo

  jc = Ny
   do ic = 1-overLAP,nx_local+overLAP
    do kc = 1-overLAP,nz_local+overLAP
      ShockXtaudxi = WallNormalXdxi(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdxi(ic,kc)
      ShockYtaudxi = WallNormalYdxi(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdxi(ic,kc)
      ShockZtaudxi = WallNormalZdxi(ic,kc) * ShockV(ic,kc) + WallNormalZ(ic,kc) * ShockVdxi(ic,kc)

      ShockXtaudzeta = WallNormalX(ic,kc) * ShockVdzeta(ic,kc)
      ShockYtaudzeta = WallNormalY(ic,kc) * ShockVdzeta(ic,kc)
      ShockZtaudzeta = WallNormalZ(ic,kc) * ShockVdzeta(ic,kc)

      tempax = ShockYtaudzeta*dzdxi(ic,jc,kc) + ShockZtaudxi*dydzeta(ic,jc,kc) - &
             & ShockYtaudxi  *dzdzeta(ic,jc,kc) - ShockZtaudzeta*dydxi(ic,jc,kc)
      
      tempay = ShockXtaudxi  *dzdzeta(ic,jc,kc) + ShockZtaudzeta*dxdxi(ic,jc,kc) - &
             & ShockXtaudzeta*dzdxi(ic,jc,kc) - ShockZtaudxi*dxdzeta(ic,jc,kc)        
      
      tempaz = ShockXtaudzeta*dydxi(ic,jc,kc) + ShockYtaudxi*dxdzeta(ic,jc,kc) - &
             & ShockXtaudxi*dydzeta(ic,jc,kc) - ShockYtaudzeta*dxdxi(ic,jc,kc)         

     shockVN(ic,kc) = -detadt(ic,jc,kc)       /nablaeta(ic,jc,kc);
      ax_tau(ic,kc) =    Jaco(ic,jc,kc)*tempax/nablaeta(ic,jc,kc);
      ay_tau(ic,kc) =    Jaco(ic,jc,kc)*tempay/nablaeta(ic,jc,kc);        
      az_tau(ic,kc) =    Jaco(ic,jc,kc)*tempaz/nablaeta(ic,jc,kc);
    enddo
   enddo

  call Parallel_Exchange(Jaco)
  call Parallel_Exchange(detadx)
  call Parallel_Exchange(detady)   
  call Parallel_Exchange(detadz)
  call Parallel_Exchange(dxidx)
  call Parallel_Exchange(dxidy)
  call Parallel_Exchange(dxidz)
  call Parallel_Exchange(dzetadx)
  call Parallel_Exchange(dzetady)
  call Parallel_Exchange(dzetadz)
  call Parallel_Exchange(dxdxi)
  call Parallel_Exchange(dydxi)
  call Parallel_Exchange(dzdxi)
  call Parallel_Exchange(dxdzeta)
  call Parallel_Exchange(dydzeta)
  call Parallel_Exchange(dzdzeta)
  call Parallel_Exchange(dxdeta)
  call Parallel_Exchange(dydeta)
  call Parallel_Exchange(dzdeta)

! We need to check the Jacobian to see if every thing is correct
!#ifdef DEBUG
!  if(MyID == 0) then
!    write(*,*)"You are debugging the Jacobian"
!    write(*,*)"Output the Jacobian info for testing"
!  endif
!  ! We need to output the Jacobian to see if every thing is correct
!  call output_jacobian
!#endif

end subroutine Calculate_Jaco_Implicit

subroutine Singular_Calculate_Jaco
  use SF_Constant,  only: ik,rk,overLAP
  use SF_CFD_Global,only: nx_local,nz_local,Ny,Nz
  use MPI_GLOBAL,   only: MyID,npx

  ! some arrays are defined in SF_CFD_Global
  use SF_CFD_Global,only: X_grid,Y_grid,Z_grid, &
                          dxdxi,dxdeta,dxdzeta,dxidx,dxidy,dxidz,nablaxi, &
                          dydxi,dydeta,dydzeta,detadx,detady,detadz,nablaeta, &
                          dzdxi,dzdeta,dzdzeta,dzetadx,dzetady,dzetadz,Jaco,nablazeta, &
                          WallSX,WallSY,WallSZ,WallSXdxi,WallSYdxi,WallSZdxi,&
                          WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                          WallNormalX,WallNormalY,WallNormalZ,&
                          WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                          WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&
                          ShockH,ShockHdxi,ShockHdzeta,ShockNormalX,ShockNormalY,ShockNormalZ,&
                          ShockV,ShockVdxi,ShockVdzeta,Heta,HetaDxi,HetaDeta,HetaDzeta,&
                          ShockXtau,ShockYtau,ShockZtau,dxidt,detadt,dzetadt,&
                          invJacodt,shockVN,ax_tau,ay_tau,az_tau,&
                          X_grid_total,Y_grid_total,Z_grid_total,&
                          ShockH_total,ShockHdxi_total,ShockHdzeta_total,&
                          ShockV_total,ShockVdxi_total,ShockVdzeta_total
  ! some subroutines
  use MPI_GLOBAL,    only: Parallel_Exchange_Surface,Parallel_Exchange
  use OutputParaView,only: output_jacobian
  use FD5_Order,     only: Cal_Deri_SurfDxi_5th_ce,Cal_Deri_SurfDzeta_5th_ce,Cal_Deri_Deta_5th_ce
  implicit none
  ! local variables
  integer( kind = ik ) :: ic, jc, kc
  real( kind = rk ) :: temp,tempax,tempay,tempaz
  real( kind = rk ) :: X_tau,Y_tau,Z_tau
  real( kind = rk ) :: ShockXtaudxi,ShockYtaudxi,ShockZtaudxi
  real( kind = rk ) :: ShockXtaudzeta,ShockYtaudzeta,ShockZtaudzeta

  real( kind = rk ) :: xxitau,xetatau,xzetatau
  real( kind = rk ) :: yxitau,yetatau,yzetatau
  real( kind = rk ) :: zxitau,zetatau,zzetatau
  
  ! remeshing
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSX(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Y_grid(ic,jc,kc) = WallSY(ic,kc) + WallNormalY(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Z_grid(ic,jc,kc) = WallSZ(ic,kc) + WallNormalZ(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
  
  call Parallel_Exchange_singularity(X_grid,X_grid_total)
  call Parallel_Exchange_singularity(Y_grid,Y_grid_total)
  call Parallel_Exchange_singularity(Z_grid,Z_grid_total)
  
  call Parallel_Exchange(X_grid)
  call Parallel_Exchange(Y_grid)
  call Parallel_Exchange(Z_grid)  
  
  ! calculate the ShockHdxi and ShockHdzeta
  call Parallel_Exchange_singularity_surface(shockH,ShockH_total) 
  call Parallel_Exchange_Surface(ShockH)
  call Cal_Deri_SurfDxi_5th_ce(ShockHdxi,    ShockH)
  call Cal_Deri_SurfDzeta_5th_ce(ShockHdzeta,ShockH)
  
  call Parallel_Exchange_singularity_surface_negative(shockHdxi,ShockHdxi_total)
  call Parallel_Exchange_singularity_surface(shockHdzeta,ShockHdzeta_total)
  call Parallel_Exchange_Surface(ShockHdxi)
  call Parallel_Exchange_Surface(ShockHdzeta)
  
  call Parallel_Exchange_singularity_surface(shockV,ShockV_total)
  call Parallel_Exchange_Surface(ShockV)
  call Cal_Deri_SurfDxi_5th_ce(ShockVdxi,    ShockV)
  call Cal_Deri_SurfDzeta_5th_ce(ShockVdzeta,ShockV)
  
  call Parallel_Exchange_singularity_surface_negative(shockVdxi,ShockVdxi_total)
  call Parallel_Exchange_singularity_surface(shockVdzeta,ShockVdzeta_total)
  call Parallel_Exchange_Surface(ShockVdxi)
  call Parallel_Exchange_Surface(ShockVdzeta)


  ! calculated by finite difference.
do kc = 1-overLAP, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP, nx_local+overLAP
  !do kc = 1, nz_local
  !  do jc = 1, Ny
  !    do ic = 1, nx_local  
  
        dxdxi(ic,jc,kc) = WallSXdxi(ic,kc) + (WallNormalXdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalX(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        xxitau =                             (WallNormalXdxi(ic,kc)*ShockV(ic,kc) + &
                                              & WallNormalX(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                              + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        dydxi(ic,jc,kc) = WallSYdxi(ic,kc) + (WallNormalYdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        yxitau =                             (WallNormalYdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        dzdxi(ic,jc,kc) = WallSZdxi(ic,kc) + (WallNormalZdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        zxitau =                             (WallNormalZdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
        
        dxdeta(ic,jc,kc) = WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dydeta(ic,jc,kc) = WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dzdeta(ic,jc,kc) = WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        
        xetatau = WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        yetatau = WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        zetatau = WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);

        dxdzeta(ic,jc,kc) = WallSXdzeta(ic,kc) + (WallNormalXdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalX(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        xzetatau =                               (WallNormalXdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalX(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        dydzeta(ic,jc,kc) = WallSYdzeta(ic,kc) + (WallNormalYdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalY(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetatau =                              (WallNormalYdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalY(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
      
        dzdzeta(ic,jc,kc) = WallSZdzeta(ic,kc) + (WallNormalZdzeta(ic,kc)*ShockH(ic,kc) + &
                                               & WallNormalZ(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetatau =                              (WallNormalZdzeta(ic,kc)*ShockV(ic,kc) + &
                                               & WallNormalZ(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);
        
        temp =    dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc))

           Jaco(ic,jc,kc)= 1.d0/temp;
        ! d(1/J)/dt = d(temp)/dt 
        invJacodt(ic,jc,kc) = & 
             &    xxitau * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + xetatau * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +xzetatau * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * (  yetatau * dzdzeta(ic,jc,kc) - zetatau* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * (  yzetatau*   dzdxi(ic,jc,kc) -  yxitau* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   yxitau *  dzdeta(ic,jc,kc) - yetatau*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * ( dydeta(ic,jc,kc)* zzetatau - dzdeta(ic,jc,kc)* yzetatau) &
             & + dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)*   zxitau -  dydxi(ic,jc,kc)* zzetatau) &
             & +dxdzeta(ic,jc,kc) * (  dydxi(ic,jc,kc)*  zetatau - dydeta(ic,jc,kc)*   zxitau)            

          dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
          dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
          dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc)); 
         detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));  
         detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
         detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));  
        dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
        dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc)); 

          nablaxi(ic,jc,kc)=sqrt(  dxidx(ic,jc,kc) *  dxidx(ic,jc,kc) &
                               &+  dxidy(ic,jc,kc) *  dxidy(ic,jc,kc) &
                               &+  dxidz(ic,jc,kc) *  dxidz(ic,jc,kc))

         nablaeta(ic,jc,kc)=sqrt( detadx(ic,jc,kc) * detadx(ic,jc,kc) &
                               &+ detady(ic,jc,kc) * detady(ic,jc,kc) &
                               &+ detadz(ic,jc,kc) * detadz(ic,jc,kc))

        nablazeta(ic,jc,kc)=sqrt(dzetadx(ic,jc,kc) *dzetadx(ic,jc,kc) &
                               &+dzetady(ic,jc,kc) *dzetady(ic,jc,kc) &
                               &+dzetadz(ic,jc,kc) *dzetadz(ic,jc,kc))

      end do
    end do
  end do   
  
  ! calculate the Shock Normal directions
  ! only along the jc = Ny surface
  jc = Ny
  do kc = 1-overLAP, nz_local+overLAP
    do ic = 1-overLAP, nx_local+overLAP
  !jc = Ny
  !do kc = 1, nz_local
  !  do ic = 1, nx_local
      ShockNormalX(ic,kc) = detadx(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalY(ic,kc) = detady(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalZ(ic,kc) = detadz(ic,jc,kc) / nablaeta(ic,jc,kc)
    enddo
  enddo

  !write(*,*),'shockNormalX(:,1)',shockNormalX(:,1)
  !write(*,*),'shockNormalY(:,1)',shockNormalY(:,1)
  !write(*,*),'shockNormalZ(:,1)',shockNormalZ(:,1)
  !
  do kc = 1-overLAP,nz_local+overLAP
    do ic = 1-overLAP,nx_local+overLAP
      do jc = 1,Ny

  !do kc = 1,nz_local
    !do ic = 1,nx_local     
      ! do jc = 1,Ny     
       X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);

       dxidt(ic,jc,kc) =  -(  dxidx(ic,jc,kc)*X_tau +  dxidy(ic,jc,kc)*Y_tau +  dxidz(ic,jc,kc)*Z_tau);

       detadt(ic,jc,kc) =  -( detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

       dzetadt(ic,jc,kc) =  -(dzetadx(ic,jc,kc)*X_tau +dzetady(ic,jc,kc)*Y_tau +dzetadz(ic,jc,kc)*Z_tau);

      enddo
      ShockXtau(ic,kc) = WallNormalX(ic,kc) * ShockV(ic,kc)
      ShockYtau(ic,kc) = WallNormalY(ic,kc) * ShockV(ic,kc)
      ShockZtau(ic,kc) = WallNormalZ(ic,kc) * ShockV(ic,kc)

    enddo
  enddo

  jc = Ny
   do ic = 1-overLAP,nx_local+overLAP
    do kc = 1-overLAP,nz_local+overLAP
  !jc = Ny
   !do ic = 1,nx_local
    !do kc = 1,nz_local
      ShockXtaudxi = WallNormalXdxi(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdxi(ic,kc)
      ShockYtaudxi = WallNormalYdxi(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdxi(ic,kc)
      ShockZtaudxi = WallNormalZdxi(ic,kc) * ShockV(ic,kc) + WallNormalZ(ic,kc) * ShockVdxi(ic,kc)

      ShockXtaudzeta = WallNormalXdzeta(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdzeta(ic,kc)
      ShockYtaudzeta = WallNormalYdzeta(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdzeta(ic,kc)
      ShockZtaudzeta = WallNormalZdzeta(ic,kc) * ShockV(ic,kc) + WallNormalZ(ic,kc) * ShockVdzeta(ic,kc)

      tempax = ShockYtaudzeta*dzdxi(ic,jc,kc) + ShockZtaudxi*dydzeta(ic,jc,kc) - &
             & ShockYtaudxi  *dzdzeta(ic,jc,kc) - ShockZtaudzeta*dydxi(ic,jc,kc)
      
      tempay = ShockXtaudxi  *dzdzeta(ic,jc,kc) + ShockZtaudzeta*dxdxi(ic,jc,kc) - &
             & ShockXtaudzeta*dzdxi(ic,jc,kc) - ShockZtaudxi*dxdzeta(ic,jc,kc)        
      
      tempaz = ShockXtaudzeta*dydxi(ic,jc,kc) + ShockYtaudxi*dxdzeta(ic,jc,kc) - &
             & ShockXtaudxi*dydzeta(ic,jc,kc) - ShockYtaudzeta*dxdxi(ic,jc,kc)         

     shockVN(ic,kc) = -detadt(ic,jc,kc)       /nablaeta(ic,jc,kc);
      ax_tau(ic,kc) =    Jaco(ic,jc,kc)*tempax/nablaeta(ic,jc,kc);
      ay_tau(ic,kc) =    Jaco(ic,jc,kc)*tempay/nablaeta(ic,jc,kc);        
      az_tau(ic,kc) =    Jaco(ic,jc,kc)*tempaz/nablaeta(ic,jc,kc);
    enddo
   enddo

   

  call Parallel_Exchange(Jaco)
  call Parallel_Exchange(detadx)
  call Parallel_Exchange(detady)   
  call Parallel_Exchange(detadz)
  call Parallel_Exchange(dxidx)
  call Parallel_Exchange(dxidy)
  call Parallel_Exchange(dxidz)
  call Parallel_Exchange(dzetadx)
  call Parallel_Exchange(dzetady)
  call Parallel_Exchange(dzetadz)
  call Parallel_Exchange(dxdxi)
  call Parallel_Exchange(dydxi)
  call Parallel_Exchange(dzdxi)
  call Parallel_Exchange(dxdzeta)
  call Parallel_Exchange(dydzeta)
  call Parallel_Exchange(dzdzeta)
  call Parallel_Exchange(dxdeta)
  call Parallel_Exchange(dydeta)
  call Parallel_Exchange(dzdeta)
  call Parallel_Exchange(invJacodt)
  call Parallel_Exchange(nablaxi)
  call Parallel_Exchange(nablazeta)
  call Parallel_Exchange(nablaeta)
  call Parallel_Exchange_surface(ShockNormalX)
  call Parallel_Exchange_surface(ShockNormalY)
  call Parallel_Exchange_surface(ShockNormalZ)
  call Parallel_Exchange(dxidt)
  call Parallel_Exchange(detadt)
  call Parallel_Exchange(dzetadt)
  call Parallel_Exchange_surface(ShockXtau)
  call Parallel_Exchange_surface(ShockYtau)
  call Parallel_Exchange_surface(ShockZtau)
  call Parallel_Exchange_surface(ax_tau)
  call Parallel_Exchange_surface(ay_tau)
  call Parallel_Exchange_surface(az_tau)
  call Parallel_Exchange_surface(shockVN)

! We need to check the Jacobian to see if every thing is correct
!#ifdef DEBUG
!  if(MyID == 0) then
!    write(*,*)"You are debugging the Jacobian"
!    write(*,*)"Output the Jacobian info for testing"
!  endif
   !We need to output the Jacobian to see if every thing is correct
!  call output_jacobian
!#endif

    end subroutine Singular_Calculate_Jaco
    
SUBROUTINE Parallel_Exchange_singularity_surface(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,ic,local_z,local_z_opposite
  real( kind = rk )::f(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1-overLap:nx_local+overLap,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1-overLap:nx_local+overLap,1-overLap:Nz+overLap)
  character(len=8)  :: date 
  character(len=10) :: time
  !integer :: color,MPI_comm_npx0
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, MPI_comm_npx0, ierr) 
 
  f_local_extended = 0.0d0
  
  if(npx == 0) then
    
    do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1     
       do ic = 1,overLap    
           f_local_extended(ic,local_z) = f(ic,kc)
       enddo
    enddo
        !call DATE_AND_TIME(date, time)
        !WRITE(*,*) 'before', MyId+1, NumProcess, time
        call MPI_ALLREDUCE(f_local_extended(1:overlap,1:Nz),f_total(1:overlap,1:Nz),overLap*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)
        ! call DATE_AND_TIME(date, time)
        !WRITE(*,*) 'after', MyId+1, NumProcess, time  
        !pause

     
     do kc = 1,Nz_local
         local_z = K_OFFSET(npz)+ kc - 1
         local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
         do ic = 1,overLap
             f(1-ic,kc) = f_total(ic,local_z_opposite)
         enddo
     enddo
     
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL      
  endif
  
    END SUBROUTINE Parallel_Exchange_singularity_surface 
        
SUBROUTINE Parallel_Exchange_singularity_surface_negative(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,ic,local_z,local_z_opposite
  real( kind = rk )::f(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1-overLap:nx_local+overLap,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1-overLap:nx_local+overLap,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)   
 
  f_local_extended = 0.0d0
  
  if(npx == 0) then

    do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1     
       do ic = 1,overLap           
           f_local_extended(ic,local_z) = f(ic,kc)
       enddo
    enddo   
 
 call MPI_ALLREDUCE(f_local_extended(1:overlap,1:Nz),f_total(1:overlap,1:Nz),overLap*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)


     do kc = 1,Nz_local
             local_z = K_OFFSET(npz)+ kc - 1
             local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
         do ic = 1,overLap
             f(1-ic,kc) = -f_total(ic,local_z_opposite)
         enddo
     enddo
     
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL    
      
  endif
  
END SUBROUTINE Parallel_Exchange_singularity_surface_negative
 
 
    
SUBROUTINE Parallel_Exchange_singularity(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,jc,ic,local_z,local_z_opposite
  real( kind = rk )::f(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)   
  
   f_local_extended = 0.0d0
  
  if(npx == 0) then

    do kc = 1,Nz_local
       local_z = K_OFFSET(npz)+ kc - 1  
       do jc = 1,Ny
         do ic = 1,overLap           
            f_local_extended(ic,jc,local_z) = f(ic,jc,kc)
         enddo
       enddo
    enddo 
 
call MPI_ALLREDUCE(f_local_extended(1:overlap,:,1:Nz),f_total(1:overlap,:,1:Nz),overLap*Ny*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)
  

     do kc = 1,Nz_local
          local_z = K_OFFSET(npz) + kc - 1
          local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
       do jc = 1,Ny     
         do ic = 1,overLap         
             f(1-ic,jc,kc) = f_total(ic,jc,local_z_opposite)
         enddo
       enddo
     enddo
     
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL 
  endif

  END SUBROUTINE Parallel_Exchange_singularity
    
    
SUBROUTINE Parallel_Exchange_singularity_negative(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,jc,ic,local_z,local_z_opposite
  real( kind = rk )::f(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)   
  
   f_local_extended = 0.0d0
  
  if(npx == 0) then

       do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1   
          do jc = 1,Ny
              do ic = 1,overLap
                 f_local_extended(ic,jc,local_z) = f(ic,jc,kc)
              enddo
          enddo
       enddo

    call MPI_ALLREDUCE(f_local_extended(1:overlap,:,1:Nz),f_total(1:overlap,:,1:Nz),overLap*Ny*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)
 
    do kc = 1,Nz_local
          local_z = K_OFFSET(npz) + kc - 1
          local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
          do jc = 1,Ny
              do ic = 1,overLap             
                 f(1-ic,jc,kc) = -f_total(ic,jc,local_z_opposite)
              enddo
          enddo
    enddo
     
     
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL 
  endif
END SUBROUTINE Parallel_Exchange_singularity_negative    
    

SUBROUTINE Parallel_Exchange_singularity_numvar(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,jc,ic,iVar,local_z,local_z_opposite
  real( kind = rk )::f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)   
 
   f_local_extended = 0.0d0
  
  if(npx == 0) then
    do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1 
       do jc = 1,Ny
         do ic = 1,overLap
             do iVar = 1,NumVar
                 f_local_extended(iVar,ic,jc,local_z) = f(iVar,ic,jc,kc)
             enddo
         enddo
       enddo
    enddo  
 
    call MPI_ALLREDUCE(f_local_extended(:,1:overlap,:,1:Nz),f_total(:,1:overlap,:,1:Nz),NumVar*overLap*Ny*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)

     do kc = 1,Nz_local
             local_z = K_OFFSET(npz) + kc - 1
             local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
       do jc = 1,Ny
         do ic = 1,overLap
             do iVar = 1,NumVar             
                 f(iVar,1-ic,jc,kc) = f_total(iVar,ic,jc,local_z_opposite)
             enddo   
         enddo
       enddo  
     enddo          
      
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL      
  endif
  
END SUBROUTINE Parallel_Exchange_singularity_numvar    
    
SUBROUTINE Parallel_Exchange_singularity_numvar_negative(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,jc,ic,iVar,local_z,local_z_opposite
  real( kind = rk )::f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)  
  
   f_local_extended = 0.0d0
  
  if(npx == 0) then
       do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1     
        do jc = 1,Ny
          do ic = 1,overLap
            do iVar = 1,NumVar
                 f_local_extended(iVar,ic,jc,local_z) = f(iVar,ic,jc,kc)
            enddo
          enddo 
        enddo
       enddo  
 
  call MPI_ALLREDUCE(f_local_extended(:,1:overlap,:,1:Nz),f_total(:,1:overlap,:,1:Nz),NumVar*overLap*Ny*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)
  

           do kc = 1,Nz_local
               local_z = K_OFFSET(npz) + kc - 1
               local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
              do jc = 1,Ny
                  do ic = 1,overLap
                      do iVar = 1,2                     
                         f(iVar,1-ic,jc,kc) = -f_total(iVar,ic,jc,local_z_opposite)
                      enddo
                  enddo 
              enddo
           enddo     
           
           iVar = 3  
           do kc = 1,Nz_local
               local_z = K_OFFSET(npz) + kc - 1
               local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
              do jc = 1,Ny
                  do ic = 1,overLap
                         f(iVar,1-ic,jc,kc) = f_total(iVar,ic,jc,local_z_opposite)
                  enddo 
              enddo
           enddo              
     
           do kc = 1,Nz_local
               local_z = K_OFFSET(npz) + kc - 1
               local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
              do jc = 1,Ny
                  do ic = 1,overLap
                      do iVar = 4,5                     
                         f(iVar,1-ic,jc,kc) = -f_total(iVar,ic,jc,local_z_opposite)
                      enddo
                  enddo 
              enddo
           enddo    


      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL  
      
  endif
  
  END SUBROUTINE Parallel_Exchange_singularity_numvar_negative
    
  SUBROUTINE Parallel_Exchange_singularity_numvar_negative_1(f,f_total)
  use SF_Constant,  only:ik,rk,overLAP,NumVar
  use SF_Constant, only: ik, rk, overLap
  use SF_CFD_Global, only: nx_local, nz_local, Ny, Nz,NumVar
  use MPI_Global,   only:npx,npz,i_offset,k_offset,MPI_COMM_WORLD,ierr,npz0,MPI_DOUBLE_PRECISION,MPI_SUM,NumProcess,MyId,MPI_UNDEFINED,MPI_COMM_NULL,MPI_comm_npx0
  !use mpi
  implicit none
  integer:: kc,jc,ic,iVar,local_z,local_z_opposite
  real( kind = rk )::f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  real( kind = rk )::f_total(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  real( kind = rk )::f_local_extended(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:Nz+overLap)
  !integer :: color, comm_npx0 
  !color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0�Ľ���color=0������ΪMPI_UNDEFINED 
  !call MPI_Comm_split(MPI_COMM_WORLD, color, MyId, comm_npx0, ierr)   
 
   f_local_extended = 0.0d0
  
  if(npx == 0) then

    do kc = 1,Nz_local
           local_z = K_OFFSET(npz)+ kc - 1     
        do jc = 1,Ny
          do ic = 1,overLap
            do iVar = 1,NumVar
                f_local_extended(iVar,ic,jc,local_z) = f(iVar,ic,jc,kc)
            enddo
          enddo 
        enddo
    enddo        

 
  call MPI_ALLREDUCE(f_local_extended(:,1:overlap,:,1:Nz),f_total(:,1:overlap,:,1:Nz),NumVar*overLap*Ny*Nz,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_comm_npx0,ierr)
  
           do kc = 1,Nz_local
               local_z = K_OFFSET(npz) + kc - 1
               local_z_opposite = MOD(local_z+Nz/2-1,Nz)+1
              do jc = 1,Ny
                  do ic = 1,overLap
                      do iVar = 1,NumVar                     
                         f(iVar,1-ic,jc,kc) = -f_total(iVar,ic,jc,local_z_opposite)
                      enddo
                  enddo 
              enddo
           enddo        
  
      !call MPI_Comm_free(MPI_comm_npx0, ierr)
      !MPI_comm_npx0 = MPI_COMM_NULL      
  endif
  
    END SUBROUTINE Parallel_Exchange_singularity_numvar_negative_1
    
    
subroutine Singular_Calculate_Jaco_Implicit
  ! This is the subroutine that calculate the derivatives with respect to the shock Height
  ! This is used in the implicit method
  use SF_Constant,  only: ik,rk,overLAP
  use SF_CFD_Global,only: nx_local,nz_local,Ny
  use MPI_GLOBAL,   only: MyID

  ! some arrays are defined in SF_CFD_Global
  use SF_CFD_Global,only: X_grid,Y_grid,Z_grid, &
                          dxdxi,dxdeta,dxdzeta,dxidx,dxidy,dxidz,nablaxi, &
                          dydxi,dydeta,dydzeta,detadx,detady,detadz,nablaeta, &
                          dzdxi,dzdeta,dzdzeta,dzetadx,dzetady,dzetadz,Jaco,nablazeta, &
                          WallSX,WallSY,WallSZ,WallSXdxi,WallSYdxi,WallSZdxi,&
                          WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                          WallNormalX,WallNormalY,WallNormalZ,&
                          WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                          WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&
                          ShockH,ShockHdxi,ShockHdzeta,ShockNormalX,ShockNormalY,ShockNormalZ,&
                          ShockV,ShockVdxi,ShockVdzeta,Heta,HetaDxi,HetaDeta,HetaDzeta,&
                          ShockXtau,ShockYtau,ShockZtau,dxidt,detadt,dzetadt,&
                          invJacodt,shockVN,ax_tau,ay_tau,az_tau,DxixJDH,DxiyJDH,DxizJDH,&
                          DetaxJDH,DetayJDH,DetazJDH,DzetaxJDH,DzetayJDH,DzetazJDH,DJDH,DetatDH,DetatJDH,&
                          DetaxJDHxi,  DetayJDHxi,  DetazJDHxi,  DetatJDHxi,DetatJDV,&
                          DetaxJDHzeta,DetayJDHzeta,DetazJDHzeta,DetatJDHzeta,&
                          DzetaxJDHxi,   DzetayJDHxi,   DzetazJDHxi,&
                          DxixJDHzeta, DxiyJDHzeta, DxizJDHzeta,DJDHxi,DJDHzeta,&
                          ShockNormalX_steady,ShockNormalY_steady,ShockNormalZ_steady,AnalysisType,shockXtau_steady,&
                          X_grid_total,Y_grid_total,Z_grid_total,&
                          ShockH_total,ShockHdxi_total,ShockHdzeta_total,&
                          ShockV_total,ShockVdxi_total,ShockVdzeta_total                          
 
  ! some subroutines
  use MPI_GLOBAL,    only: Parallel_Exchange_Surface,Parallel_Exchange
  use OutputParaView,only: output_jacobian
  use FD5_Order,     only: Cal_Deri_SurfDxi_5th_ce,Cal_Deri_SurfDzeta_5th_ce,Cal_Deri_Deta_5th_ce
  implicit none
  ! local variables
  integer( kind = ik ) :: ic, jc, kc
  real( kind = rk ) :: temp,tempax,tempay,tempaz
  real( kind = rk ) :: X_tau,Y_tau,Z_tau
  real( kind = rk ) :: ShockXtaudxi,ShockYtaudxi,ShockZtaudxi
  real( kind = rk ) :: ShockXtaudzeta,ShockYtaudzeta,ShockZtaudzeta

  real( kind = rk ) :: xxiH,xetaH,xzetaH
  real( kind = rk ) :: yxiH,yetaH,yzetaH
  real( kind = rk ) :: zxiH,zetaH,zzetaH
  
  real( kind = rk ) :: xxiHxi,yxiHxi,zxiHxi
  real( kind = rk ) :: xzetaHzeta,yzetaHzeta,zzetaHzeta
  real( kind = rk ) :: xxitau,xetatau,xzetatau
  real( kind = rk ) :: yxitau,yetatau,yzetatau
  real( kind = rk ) :: zxitau,zetatau,zzetatau
  
    ! remeshing
    do kc = 1-overLAP,nz_local+overLAP
      do jc = 1,Ny
        do ic = 1-overLAP,nx_local+overLAP
          X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
          Y_grid(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
          Z_grid(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        end do
      end do
    end do
    ! 
    call Parallel_Exchange_singularity(X_grid,X_grid_total)
    call Parallel_Exchange_singularity(Y_grid,Y_grid_total)
    call Parallel_Exchange_singularity(Z_grid,Z_grid_total)
    
    call Parallel_Exchange(X_grid)
    call Parallel_Exchange(Y_grid)
    call Parallel_Exchange(Z_grid)
  
    ! calculate the ShockHdxi and ShockHdzeta
  call Parallel_Exchange_singularity_surface(shockH,ShockH_total) 
  call Parallel_Exchange_Surface(ShockH)
  call Cal_Deri_SurfDxi_5th_ce(ShockHdxi,    ShockH)
  call Cal_Deri_SurfDzeta_5th_ce(ShockHdzeta,ShockH)
  
  call Parallel_Exchange_singularity_surface_negative(shockHdxi,ShockHdxi_total)
  call Parallel_Exchange_singularity_surface(shockHdzeta,ShockHdzeta_total)
  call Parallel_Exchange_Surface(ShockHdxi)
  call Parallel_Exchange_Surface(ShockHdzeta)
  
  call Parallel_Exchange_singularity_surface(shockV,ShockV_total)
  call Parallel_Exchange_Surface(ShockV)
  call Cal_Deri_SurfDxi_5th_ce(ShockVdxi,    ShockV)
  call Cal_Deri_SurfDzeta_5th_ce(ShockVdzeta,ShockV)
  
  call Parallel_Exchange_singularity_surface_negative(shockVdxi,ShockVdxi_total)
  call Parallel_Exchange_singularity_surface(shockVdzeta,ShockVdzeta_total)
  call Parallel_Exchange_Surface(ShockVdxi)
  call Parallel_Exchange_Surface(ShockVdzeta)

  ! Here we do not use finite difference to calculate
  ! the dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta
  ! dzdxi, dzdeta, dzdzeta, because we have already known
  ! the analytical expression of surface grid and the 
  ! mapping function.
  ! only the derivative with respect to shock height is
  ! calculated by finite difference.
  
    
  do kc = 1-overLap, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP, nx_local+overLAP
        dxdxi(ic,jc,kc) = WallSXdxi(ic,kc) + (WallNormalXdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalX(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);

        xxitau =                             (WallNormalXdxi(ic,kc)*ShockV(ic,kc) + &
                                              & WallNormalX(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                              + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        xxiH = WallNormalXdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalX(ic,kc) * HetaDxi(ic,jc,kc); 

        xxiHxi = WallNormalX(ic,kc) * Heta(ic,jc,kc);

        dydxi(ic,jc,kc) = WallSYdxi(ic,kc) + (WallNormalYdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
        
        yxitau =                             (WallNormalYdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalY(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)
                                           
        yxiH = WallNormalYdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalY(ic,kc) * HetaDxi(ic,jc,kc);

        yxiHxi = WallNormalY(ic,kc) * Heta(ic,jc,kc);
        
        dzdxi(ic,jc,kc) = WallSZdxi(ic,kc) + (WallNormalZdxi(ic,kc)*ShockH(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockHdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDxi(ic,jc,kc);
        
        zxitau =                             (WallNormalZdxi(ic,kc)*ShockV(ic,kc) + &
                                            & WallNormalZ(ic,kc)*ShockVdxi(ic,kc)) * Heta(ic,jc,kc) &
                                           + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDxi(ic,jc,kc)

        zxiH = WallNormalZdxi(ic,kc) * Heta(ic,jc,kc) + WallNormalZ(ic,kc) * HetaDxi(ic,jc,kc);

        zxiHxi = WallNormalZ(ic,kc) * Heta(ic,jc,kc);
        
        dxdeta(ic,jc,kc) = WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dydeta(ic,jc,kc) = WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);
        dzdeta(ic,jc,kc) = WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDeta(ic,jc,kc);

        xetaH = WallNormalX(ic,kc) * HetaDeta(ic,jc,kc);
        yetaH = WallNormalY(ic,kc) * HetaDeta(ic,jc,kc);
        zetaH = WallNormalZ(ic,kc) * HetaDeta(ic,jc,kc);

        xetatau = WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        yetatau = WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);
        zetatau = WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDeta(ic,jc,kc);

        dxdzeta(ic,jc,kc) = WallSXdzeta(ic,kc) + (WallNormalXdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalX(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        xzetatau =                               (WallNormalXdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalX(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalX(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        xzetaH = WallNormalXdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalX(ic,kc) * HetaDzeta(ic,jc,kc);

        xzetaHzeta = WallNormalX(ic,kc) * Heta(ic,jc,kc);

        dydzeta(ic,jc,kc) = WallSYdzeta(ic,kc) + (WallNormalYdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalY(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalY(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        yzetatau =                              (WallNormalYdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalY(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                                + WallNormalY(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetaH = WallNormalYdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalY(ic,kc) * HetaDzeta(ic,jc,kc);

        yzetaHzeta = WallNormalY(ic,kc) * Heta(ic,jc,kc);

        dzdzeta(ic,jc,kc) = WallSZdzeta(ic,kc) + (WallNormalZdzeta(ic,kc)*ShockH(ic,kc) + &
                                                & WallNormalZ(ic,kc)*ShockHdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                               + WallNormalZ(ic,kc) * ShockH(ic,kc) * HetaDzeta(ic,jc,kc);
        
        zzetatau =                              (WallNormalZdzeta(ic,kc)*ShockV(ic,kc) + &
                                                & WallNormalZ(ic,kc)*ShockVdzeta(ic,kc)) * Heta(ic,jc,kc) &
                                                + WallNormalZ(ic,kc) * ShockV(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetaH = WallNormalZdzeta(ic,kc) * Heta(ic,jc,kc) + WallNormalZ(ic,kc) * HetaDzeta(ic,jc,kc);

        zzetaHzeta = WallNormalZ(ic,kc) * Heta(ic,jc,kc);

        temp =    dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc))

        Jaco(ic,jc,kc)= 1.d0/temp;          

          dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
          dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
          dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc)); 
         detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));  
         detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
         detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));  
        dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
        dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));

          nablaxi(ic,jc,kc)=sqrt(  dxidx(ic,jc,kc) *  dxidx(ic,jc,kc) &
                               &+  dxidy(ic,jc,kc) *  dxidy(ic,jc,kc) &
                               &+  dxidz(ic,jc,kc) *  dxidz(ic,jc,kc))

         nablaeta(ic,jc,kc)=sqrt( detadx(ic,jc,kc) * detadx(ic,jc,kc) &
                               &+ detady(ic,jc,kc) * detady(ic,jc,kc) &
                               &+ detadz(ic,jc,kc) * detadz(ic,jc,kc))

        nablazeta(ic,jc,kc)=sqrt(dzetadx(ic,jc,kc) *dzetadx(ic,jc,kc) &
                               &+dzetady(ic,jc,kc) *dzetady(ic,jc,kc) &
                               &+dzetadz(ic,jc,kc) *dzetadz(ic,jc,kc))

        ! d(1/J)/dt = d(temp)/dt 
        invJacodt(ic,jc,kc) = & 
             &    xxitau * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
             & + xetatau * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
             & +xzetatau * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * (  yetatau * dzdzeta(ic,jc,kc) - zetatau* dydzeta(ic,jc,kc)) &
             & + dxdeta(ic,jc,kc) * (  yzetatau*   dzdxi(ic,jc,kc) -  yxitau* dzdzeta(ic,jc,kc)) &
             & +dxdzeta(ic,jc,kc) * (   yxitau *  dzdeta(ic,jc,kc) - yetatau*   dzdxi(ic,jc,kc)) &
             & +  dxdxi(ic,jc,kc) * ( dydeta(ic,jc,kc)* zzetatau - dzdeta(ic,jc,kc)* yzetatau) &
             & + dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)*   zxitau -  dydxi(ic,jc,kc)* zzetatau) &
             & +dxdzeta(ic,jc,kc) * (  dydxi(ic,jc,kc)*  zetatau - dydeta(ic,jc,kc)*   zxitau)     
        
        ! This part calculate the   d(xix/J)/dH,   d(xiy/J)/dH,   d(xiz/J)/dH
        !                          d(etax/J)/dH,  d(etay/J)/dH,  d(etaz/J)/dH
        !                         d(zetax/J)/dH, d(zetay/J)/dH, d(zetaz/J)/dH 
        !dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
         DxixJDH(ic,jc,kc) = (  yetaH*dzdzeta(ic,jc,kc) +  dydeta(ic,jc,kc)*  zzetaH &
                           -   yzetaH* dzdeta(ic,jc,kc) - dydzeta(ic,jc,kc)*   zetaH)
        
        !dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
         DxiyJDH(ic,jc,kc) =-(  xetaH*dzdzeta(ic,jc,kc) +  dxdeta(ic,jc,kc)*  zzetaH &
                           -    zetaH*dxdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*  xzetaH)
        
        !dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
         DxizJDH(ic,jc,kc) = (  xetaH*dydzeta(ic,jc,kc) +  dxdeta(ic,jc,kc)*  yzetaH &
                           -    yetaH*dxdzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*  xzetaH)
        
        !detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));
        DetaxJDH(ic,jc,kc) =-(   yxiH*dzdzeta(ic,jc,kc) +   dydxi(ic,jc,kc)*  zzetaH &
                           -     zxiH*dydzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*  yzetaH)
        
        !detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        DetayJDH(ic,jc,kc) = (   xxiH*dzdzeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*  zzetaH &
                           -     zxiH*dxdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*  xzetaH)
        
        !detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
        DetazJDH(ic,jc,kc) =-(   xxiH*dydzeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*  yzetaH &
                           -     yxiH*dxdzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*  xzetaH)
       
       !dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));
       DzetaxJDH(ic,jc,kc) = (   yxiH* dzdeta(ic,jc,kc) +   dydxi(ic,jc,kc)*   zetaH &
                           -     zxiH* dydeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*   yetaH)
       
       !dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
       DzetayJDH(ic,jc,kc) =-(   xxiH* dzdeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*   zetaH &
                           -     zxiH* dxdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*   xetaH) 
       
       !dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));
       DzetazJDH(ic,jc,kc) = (   xxiH* dydeta(ic,jc,kc) +   dxdxi(ic,jc,kc)*   yetaH &
                           -     yxiH* dxdeta(ic,jc,kc) -   dydxi(ic,jc,kc)*   xetaH)
            
            ! J = 1/temp, DJ/DH = -Dtemp/DH/(temp*temp)
            !             DJ/DHxi = -Dtemp/DHxi/(temp*temp)
            !             DJ/DHzeta = -Dtemp/DHzeta/(temp*temp)
            DJDH(ic,jc,kc) = xxiH * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
                        & + xetaH * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
                        & +xzetaH * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
                        & +  dxdxi(ic,jc,kc) * (  yetaH* dzdzeta(ic,jc,kc) - zetaH* dydzeta(ic,jc,kc)) &
                        & + dxdeta(ic,jc,kc) * ( yzetaH*   dzdxi(ic,jc,kc) -  yxiH* dzdzeta(ic,jc,kc)) &
                        & +dxdzeta(ic,jc,kc) * (   yxiH * dzdeta(ic,jc,kc) - yetaH*   dzdxi(ic,jc,kc)) &
                        & +  dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* zzetaH - dzdeta(ic,jc,kc)* yzetaH) &
                        & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*   zxiH -  dydxi(ic,jc,kc)* zzetaH) &
                        & +dxdzeta(ic,jc,kc) * (   dydxi(ic,jc,kc) * zetaH - dydeta(ic,jc,kc)*   zxiH)
            
            DJDH(ic,jc,kc) = - DJDH(ic,jc,kc) / (temp * temp);

            DJDHxi(ic,jc,kc) = xxiHxi * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
                          & + dxdeta(ic,jc,kc) * ( dydzeta(ic,jc,kc)*  zxiHxi -  yxiHxi* dzdzeta(ic,jc,kc)) &
                          & +dxdzeta(ic,jc,kc) * (   yxiHxi * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   zxiHxi)
            
            DJDHxi(ic,jc,kc) = - DJDHxi(ic,jc,kc) / (temp * temp);

            DJDHzeta(ic,jc,kc) = dxdxi(ic,jc,kc) * (  dydeta(ic,jc,kc)* zzetaHzeta - dzdeta(ic,jc,kc)* yzetaHzeta) &
                          & +  dxdeta(ic,jc,kc) * ( yzetaHzeta *   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* zzetaHzeta) &
                          & +xzetaHzeta * ( dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*dzdxi(ic,jc,kc)) 
            
            DJDHzeta(ic,jc,kc) = - DJDHzeta(ic,jc,kc) / (temp * temp);

            X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
            Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
            Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);

            !
            detadt(ic,jc,kc) = -(detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

            DetatJDH(ic,jc,kc) = -(DetaxJDH(ic,jc,kc)*X_tau + DetayJDH(ic,jc,kc)*Y_tau + DetazJDH(ic,jc,kc)*Z_tau);

            DetatJDV(ic,jc,kc) = -(detadx(ic,jc,kc) * WallNormalX(ic,kc) * Heta(ic,jc,kc) +&
                                 & detady(ic,jc,kc) * WallNormalY(ic,kc) * Heta(ic,jc,kc) +&
                                 & detadz(ic,jc,kc) * WallNormalZ(ic,kc) * Heta(ic,jc,kc))/jaco(ic,jc,kc);

           !dxidx(ic,jc,kc)=  Jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)* dzdeta(ic,jc,kc));
           !dxidy(ic,jc,kc)= -Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) -  dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
           !dxidz(ic,jc,kc)=  Jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) -  dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
            DxixJDHzeta(ic,jc,kc) = ( dydeta(ic,jc,kc)*zzetaHzeta - yzetaHzeta* dzdeta(ic,jc,kc) )

            DxiyJDHzeta(ic,jc,kc) =-( dxdeta(ic,jc,kc)*zzetaHzeta -  dzdeta(ic,jc,kc)*xzetaHzeta )

            DxizJDHzeta(ic,jc,kc) = ( dxdeta(ic,jc,kc)*yzetaHzeta -  dydeta(ic,jc,kc)*xzetaHzeta )

           !detadx(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));
           !detady(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) -   dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
           !detadz(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) -   dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
              DetaxJDHxi(ic,jc,kc) =-(   yxiHxi * dzdzeta(ic,jc,kc) - zxiHxi * dydzeta(ic,jc,kc))
            DetaxJDHzeta(ic,jc,kc) =-(   dydxi(ic,jc,kc)*  zzetaHzeta - dzdxi(ic,jc,kc) * yzetaHzeta)
            
              DetayJDHxi(ic,jc,kc) = (   xxiHxi * dzdzeta(ic,jc,kc) - zxiHxi * dxdzeta(ic,jc,kc))
            DetayJDHzeta(ic,jc,kc) = (   dxdxi(ic,jc,kc)*  zzetaHzeta - dzdxi(ic,jc,kc) * xzetaHzeta)
            
              DetazJDHxi(ic,jc,kc) =-(   xxiHxi * dydzeta(ic,jc,kc) - yxiHxi * dxdzeta(ic,jc,kc))
            DetazJDHzeta(ic,jc,kc) =-(   dxdxi(ic,jc,kc)*  yzetaHzeta - dydxi(ic,jc,kc) * xzetaHzeta)
           
           !dzetadx(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dydeta(ic,jc,kc));
           !dzetady(ic,jc,kc)= -Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dzdeta(ic,jc,kc) -   dzdxi(ic,jc,kc)* dxdeta(ic,jc,kc));
           !dzetadz(ic,jc,kc)=  Jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)* dydeta(ic,jc,kc) -   dydxi(ic,jc,kc)* dxdeta(ic,jc,kc));
              DzetaxJDHxi(ic,jc,kc) = (   yxiHxi * dzdeta(ic,jc,kc) - zxiHxi * dydeta(ic,jc,kc))

              DzetayJDHxi(ic,jc,kc) =-(   xxiHxi * dzdeta(ic,jc,kc) + zxiHxi * dxdeta(ic,jc,kc))
            
              DzetazJDHxi(ic,jc,kc) = (   xxiHxi * dydeta(ic,jc,kc) + yxiHxi * dxdeta(ic,jc,kc))
           
           !etat related variables 
              DetatJDHxi(ic,jc,kc) = -(  DetaxJDHxi(ic,jc,kc)*X_tau +   DetayJDHxi(ic,jc,kc)*Y_tau +   DetazJDHxi(ic,jc,kc)*Z_tau);

            DetatJDHzeta(ic,jc,kc) = -(DetaxJDHzeta(ic,jc,kc)*X_tau + DetayJDHzeta(ic,jc,kc)*Y_tau + DetazJDHzeta(ic,jc,kc)*Z_tau);

           !Then the etat related variables with respect to V and Vxi, Vzeta
           ! Vxi and Vzeta term are only appear at the shock boundary conditions
          
      end do
    end do
  end do
! calculate the Shock Normal directions
  ! only along the jc = Ny surface
  jc = Ny
  do kc = 1-overLAP, nz_local+overLAP
    do ic = 1-overLAP, nx_local+overLAP
      ShockNormalX(ic,kc) = detadx(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalY(ic,kc) = detady(ic,jc,kc) / nablaeta(ic,jc,kc)
      ShockNormalZ(ic,kc) = detadz(ic,jc,kc) / nablaeta(ic,jc,kc)
    enddo
  enddo

    if (AnalysisType == 3) then    
      ShockNormalX_steady = ShockNormalX
      ShockNormalY_steady = ShockNormalY
      ShockNormalZ_steady = ShockNormalZ
  endif
  
  do kc = 1-overLAP,nz_local+overLAP
    do ic = 1-overLAP,nx_local+overLAP
      do jc = 1,Ny
       X_tau = WallNormalX(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Y_tau = WallNormalY(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);
       Z_tau = WallNormalZ(ic,kc) * ShockV(ic,kc) * Heta(ic,jc,kc);

        dxidt(ic,jc,kc) =  -(  dxidx(ic,jc,kc)*X_tau +  dxidy(ic,jc,kc)*Y_tau +  dxidz(ic,jc,kc)*Z_tau);

        detadt(ic,jc,kc) =  -( detadx(ic,jc,kc)*X_tau + detady(ic,jc,kc)*Y_tau + detadz(ic,jc,kc)*Z_tau);

        dzetadt(ic,jc,kc) =  -(dzetadx(ic,jc,kc)*X_tau +dzetady(ic,jc,kc)*Y_tau +dzetadz(ic,jc,kc)*Z_tau);
       
       !detadt_invJ(ic,jc,kc) = detadt(ic,jc,kc)/Jaco(ic,jc,kc);
      enddo
      ShockXtau(ic,kc) = WallNormalX(ic,kc) * ShockV(ic,kc)
      ShockYtau(ic,kc) = WallNormalY(ic,kc) * ShockV(ic,kc)
      ShockZtau(ic,kc) = WallNormalZ(ic,kc) * ShockV(ic,kc)
    enddo
  enddo

  if (AnalysisType == 3) then    
      ShockXtau_steady = ShockXtau
  endif
  
  jc = Ny
   do ic = 1-overLAP,nx_local+overLAP
    do kc = 1-overLAP,nz_local+overLAP
      ShockXtaudxi = WallNormalXdxi(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdxi(ic,kc)
      ShockYtaudxi = WallNormalYdxi(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdxi(ic,kc)
      ShockZtaudxi = WallNormalZdxi(ic,kc) * ShockV(ic,kc) + WallNormalZ(ic,kc) * ShockVdxi(ic,kc)

      ShockXtaudzeta = WallNormalXdzeta(ic,kc) * ShockV(ic,kc) + WallNormalX(ic,kc) * ShockVdzeta(ic,kc)
      ShockYtaudzeta = WallNormalYdzeta(ic,kc) * ShockV(ic,kc) + WallNormalY(ic,kc) * ShockVdzeta(ic,kc)
      ShockZtaudzeta = WallNormalZdzeta(ic,kc) * ShockV(ic,kc) + WallNormalZ(ic,kc) * ShockVdzeta(ic,kc)

      tempax = ShockYtaudzeta*dzdxi(ic,jc,kc) + ShockZtaudxi*dydzeta(ic,jc,kc) - &
             & ShockYtaudxi  *dzdzeta(ic,jc,kc) - ShockZtaudzeta*dydxi(ic,jc,kc)
      
      tempay = ShockXtaudxi  *dzdzeta(ic,jc,kc) + ShockZtaudzeta*dxdxi(ic,jc,kc) - &
             & ShockXtaudzeta*dzdxi(ic,jc,kc) - ShockZtaudxi*dxdzeta(ic,jc,kc)        
      
      tempaz = ShockXtaudzeta*dydxi(ic,jc,kc) + ShockYtaudxi*dxdzeta(ic,jc,kc) - &
             & ShockXtaudxi*dydzeta(ic,jc,kc) - ShockYtaudzeta*dxdxi(ic,jc,kc)         

     shockVN(ic,kc) = -detadt(ic,jc,kc)       /nablaeta(ic,jc,kc);
      ax_tau(ic,kc) =    Jaco(ic,jc,kc)*tempax/nablaeta(ic,jc,kc);
      ay_tau(ic,kc) =    Jaco(ic,jc,kc)*tempay/nablaeta(ic,jc,kc);        
      az_tau(ic,kc) =    Jaco(ic,jc,kc)*tempaz/nablaeta(ic,jc,kc);
    enddo
   enddo

  call Parallel_Exchange(Jaco)
  call Parallel_Exchange(detadx)
  call Parallel_Exchange(detady)   
  call Parallel_Exchange(detadz)
  call Parallel_Exchange(dxidx)
  call Parallel_Exchange(dxidy)
  call Parallel_Exchange(dxidz)
  call Parallel_Exchange(dzetadx)
  call Parallel_Exchange(dzetady)
  call Parallel_Exchange(dzetadz)
  call Parallel_Exchange(dxdxi)
  call Parallel_Exchange(dydxi)
  call Parallel_Exchange(dzdxi)
  call Parallel_Exchange(dxdzeta)
  call Parallel_Exchange(dydzeta)
  call Parallel_Exchange(dzdzeta)
  call Parallel_Exchange(dxdeta)
  call Parallel_Exchange(dydeta)
  call Parallel_Exchange(dzdeta)
  call Parallel_Exchange(invJacodt)
  call Parallel_Exchange(nablaxi)
  call Parallel_Exchange(nablazeta)
  call Parallel_Exchange(nablaeta)
  call Parallel_Exchange_surface(ShockNormalX)
  call Parallel_Exchange_surface(ShockNormalY)
  call Parallel_Exchange_surface(ShockNormalZ)
  call Parallel_Exchange(dxidt)
  call Parallel_Exchange(detadt)
  call Parallel_Exchange(dzetadt)
  call Parallel_Exchange_surface(ShockXtau)
  call Parallel_Exchange_surface(ShockYtau)
  call Parallel_Exchange_surface(ShockZtau)
  call Parallel_Exchange_surface(ax_tau)
  call Parallel_Exchange_surface(ay_tau)
  call Parallel_Exchange_surface(az_tau)
  call Parallel_Exchange_surface(shockVN)

! We need to check the Jacobian to see if every thing is correct
!#ifdef DEBUG
!  if(MyID == 0) then
!    write(*,*)"You are debugging the Jacobian"
!    write(*,*)"Output the Jacobian info for testing"
!  endif
!  ! We need to output the Jacobian to see if every thing is correct
!  call output_jacobian
!#endif

end subroutine Singular_Calculate_Jaco_Implicit    