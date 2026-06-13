SUBROUTINE Linearized_shockBC(shockH,shockV,&
                      &Rho_free,U_free,V_free,W_free,P_free,&
                      &nx_local,ny_local,nz_local,overlap,Gamma,Mach_Ref,Re_Ref,&
                      &Wallsx,Wallsy,Wallsz,k_infty,epsilon,pert_type,&
                      &Wallsxdxi,Wallsydxi,Wallszdxi,&
                      &Wallsxdzeta,Wallsydzeta,Wallszdzeta,&
                      &WallNormalX,WallNormalY,WallNormalZ,&
                      &WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                      &WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&
                      &Heta,HetaDxi,HetaDzeta,HetaDeta,&
                      &sRho,sU,sV,sW,sP,sT)
 !use MPI_GLOBAL,   only: Parallel_Exchange_Surface, Parallel_Exchange,&
 !                     & Parallel_Exchange_NumVar
    
 implicit none
    integer(4),intent(in) :: nx_local,ny_local,nz_local,overlap,pert_type
    real(8),intent(in) :: gamma,mach_ref,Re_ref,epsilon,k_infty
    real(8),intent(in) :: Rho_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: U_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: V_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: W_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: P_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: shockh(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: shockv(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap) 
    real(8),intent(in) :: wallsx(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsy(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsz(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsxdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsydxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallszdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsxdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallsydzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallszdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalx(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormaly(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalz(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalxdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalydxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalzdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalxdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalydzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: wallnormalzdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: heta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: hetadxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: hetadeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: hetadzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sRho(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sU(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sV(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sW(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sP(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out) :: sT(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    
    integer::ic,jc,kc
    
    real(8),dimension(:,:,:),  allocatable::x_grid,y_grid,z_grid
    real(8),dimension(:,:),    allocatable::shockhdzeta,shockhdxi,shockvdzeta,shockvdxi
    real(8),dimension(:,:,:),  allocatable::dxdxi,dydxi,dzdxi,dxdeta,dydeta,dzdeta,dxdzeta,dydzeta,dzdzeta
    real(8),dimension(:,:,:),  allocatable::jaco,invjacodt
    real(8),dimension(:,:,:),  allocatable::dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz
    real(8),dimension(:,:,:),  allocatable::nablaxi,nablaeta,nablazeta
    real(8),dimension(:,:),    allocatable::shocknormalx,shocknormaly,shocknormalz
      
    real(8)::xxitau,yxitau,zxitau,xetatau,yetatau,zetatau,xzetatau,yzetatau,zzetatau,temp
    real(8)::Delta_Inf,Delta2
    
    allocate(x_grid(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(y_grid(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(z_grid(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(shockhdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockhdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockvdzeta(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockvdxi(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(dxdxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dydxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzdxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dxdeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dydeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzdeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dxdzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dydzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzdzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(jaco(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invjacodt(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dxidx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dxidy(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dxidz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(detadx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(detady(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(detadz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzetadx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzetady(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzetadz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(nablaxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(nablaeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(nablazeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(shocknormalx(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shocknormaly(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shocknormalz(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny_local
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Y_grid(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
        Z_grid(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
  
  !CALL Parallel_Exchange(X_grid)
  !CALL Parallel_Exchange(Y_grid)
  !CALL Parallel_Exchange(Z_grid)  
  !
  !CALL Parallel_Exchange_surface(shockH)
  !CALL Parallel_Exchange_surface(shockV)

  CALL calculate_fdxi_surf_ce_BC(shockhdxi,shockh,nx_local,nz_local,overlap)
  CALL calculate_fdxi_surf_ce_BC(shockvdxi,shockv,nx_local,nz_local,overlap)
  
  CALL calculate_fdzeta_surf_ce_BC(shockhdzeta,shockh,nx_local,nz_local,overlap)
  CALL calculate_fdzeta_surf_ce_BC(shockvdzeta,shockv,nx_local,nz_local,overlap)  

  !CALL Parallel_Exchange_surface(shockhdxi)
  !CALL Parallel_Exchange_surface(shockvdxi)
  !CALL Parallel_Exchange_surface(shockhdzeta)
  !CALL Parallel_Exchange_surface(shockvdzeta)  
  
  do kc=1-overlap, nz_local+overlap
     do jc=1, ny_local
         do ic=1-overlap, nx_local+overlap
         dxdxi(ic,jc,kc) = wallsxdxi(ic,kc) + (wallnormalxdxi(ic,kc)*shockh(ic,kc) + &
                         & wallnormalx(ic,kc)*shockhdxi(ic,kc)) * heta(ic,jc,kc) &
                         + wallnormalx(ic,kc) * shockh(ic,kc) * hetadxi(ic,jc,kc);
         
         xxitau =          (wallnormalxdxi(ic,kc)*shockv(ic,kc) + &
                          & wallnormalx(ic,kc)*shockvdxi(ic,kc)) * heta(ic,jc,kc) &
                          + wallnormalx(ic,kc) * shockv(ic,kc) * hetadxi(ic,jc,kc)
         
         dydxi(ic,jc,kc) = wallsydxi(ic,kc) + (wallnormalydxi(ic,kc)*shockh(ic,kc) + &
                         & wallnormaly(ic,kc)*shockhdxi(ic,kc)) * heta(ic,jc,kc) &
                         + wallnormaly(ic,kc) * shockh(ic,kc) * hetadxi(ic,jc,kc);
         
         yxitau =          (wallnormalydxi(ic,kc)*shockv(ic,kc) + &
                          & wallnormaly(ic,kc)*shockvdxi(ic,kc)) * heta(ic,jc,kc) &
                          + wallnormaly(ic,kc) * shockv(ic,kc) * hetadxi(ic,jc,kc)
         
         dzdxi(ic,jc,kc) = wallszdxi(ic,kc) + (wallnormalzdxi(ic,kc)*shockh(ic,kc) + &
                         & wallnormalz(ic,kc)*shockhdxi(ic,kc)) * heta(ic,jc,kc) &
                         + wallnormalz(ic,kc) * shockh(ic,kc) * hetadxi(ic,jc,kc); 
         
         zxitau =         (wallnormalzdxi(ic,kc)*shockv(ic,kc) + &
                         & wallnormalz(ic,kc)*shockvdxi(ic,kc)) * heta(ic,jc,kc) &
                         + wallnormalz(ic,kc) * shockv(ic,kc) * hetadxi(ic,jc,kc)
         
         dxdeta(ic,jc,kc) = wallnormalx(ic,kc) * shockh(ic,kc) * hetadeta(ic,jc,kc);
         dydeta(ic,jc,kc) = wallnormaly(ic,kc) * shockh(ic,kc) * hetadeta(ic,jc,kc);
         dzdeta(ic,jc,kc) = wallnormalz(ic,kc) * shockh(ic,kc) * hetadeta(ic,jc,kc);
         
         xetatau = wallnormalx(ic,kc) * shockv(ic,kc) * hetadeta(ic,jc,kc);
         yetatau = wallnormaly(ic,kc) * shockv(ic,kc) * hetadeta(ic,jc,kc);
         zetatau = wallnormalz(ic,kc) * shockv(ic,kc) * hetadeta(ic,jc,kc);
         
         dxdzeta(ic,jc,kc) = wallsxdzeta(ic,kc) + (wallnormalxdzeta(ic,kc)*shockh(ic,kc) + &
                           & wallnormalx(ic,kc)*shockhdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormalx(ic,kc) * shockh(ic,kc) * hetadzeta(ic,jc,kc);
         
         xzetatau =         (wallnormalxdzeta(ic,kc)*shockv(ic,kc) + &
                           & wallnormalx(ic,kc)*shockvdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormalx(ic,kc) * shockv(ic,kc) * hetadzeta(ic,jc,kc);
         
         dydzeta(ic,jc,kc) = wallsydzeta(ic,kc) + (wallnormalydzeta(ic,kc)*shockh(ic,kc) + &
                           & wallnormaly(ic,kc)*shockhdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormaly(ic,kc) * shockh(ic,kc) * hetadzeta(ic,jc,kc);
         
         yzetatau =         (wallnormalydzeta(ic,kc)*shockv(ic,kc) + &
                           & wallnormaly(ic,kc)*shockvdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormaly(ic,kc) * shockv(ic,kc) * hetadzeta(ic,jc,kc);
         
         dzdzeta(ic,jc,kc) = wallszdzeta(ic,kc) + (wallnormalzdzeta(ic,kc)*shockh(ic,kc) + &
                           & wallnormalz(ic,kc)*shockhdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormalz(ic,kc) * shockh(ic,kc) * hetadzeta(ic,jc,kc);
         
         zzetatau =         (wallnormalzdzeta(ic,kc)*shockv(ic,kc) + &
                           & wallnormalz(ic,kc)*shockvdzeta(ic,kc)) * heta(ic,jc,kc) &
                           + wallnormalz(ic,kc) * shockv(ic,kc) * hetadzeta(ic,jc,kc);
         
         temp =      dxdxi(ic,jc,kc) * (dydeta(ic,jc,kc) * dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc) * dydzeta(ic,jc,kc)) &
                & + dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)* dzdxi(ic,jc,kc)   - dydxi(ic,jc,kc)  * dzdzeta(ic,jc,kc)) &
               & + dxdzeta(ic,jc,kc) * (dydxi(ic,jc,kc)  * dzdeta(ic,jc,kc)  - dydeta(ic,jc,kc) * dzdxi(ic,jc,kc))
         
         jaco(ic,jc,kc) = 1.0d0/temp;
         
         invjacodt(ic,jc,kc) = & 
                            &     xxitau * (  dydeta(ic,jc,kc)* dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)* dydzeta(ic,jc,kc)) &
                            & +  xetatau * ( dydzeta(ic,jc,kc)*   dzdxi(ic,jc,kc) -  dydxi(ic,jc,kc)* dzdzeta(ic,jc,kc)) &
                            & + xzetatau * (   dydxi(ic,jc,kc) * dzdeta(ic,jc,kc) - dydeta(ic,jc,kc)*   dzdxi(ic,jc,kc)) &
                            & +   dxdxi(ic,jc,kc) * (  yetatau * dzdzeta(ic,jc,kc) - zetatau* dydzeta(ic,jc,kc)) &
                            & +  dxdeta(ic,jc,kc) * (  yzetatau*   dzdxi(ic,jc,kc) -  yxitau* dzdzeta(ic,jc,kc)) &
                            & + dxdzeta(ic,jc,kc) * (   yxitau *  dzdeta(ic,jc,kc) - yetatau*   dzdxi(ic,jc,kc)) &
                            & +   dxdxi(ic,jc,kc) * ( dydeta(ic,jc,kc)* zzetatau - dzdeta(ic,jc,kc)* yzetatau) &
                            & +  dxdeta(ic,jc,kc) * (dydzeta(ic,jc,kc)*   zxitau -  dydxi(ic,jc,kc)* zzetatau) &
                            & + dxdzeta(ic,jc,kc) * (  dydxi(ic,jc,kc)*  zetatau - dydeta(ic,jc,kc)*   zxitau)        
         
           dxidx(ic,jc,kc) =  jaco(ic,jc,kc) * (dydeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dydzeta(ic,jc,kc)*dzdeta(ic,jc,kc));
           dxidy(ic,jc,kc) = -jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dzdzeta(ic,jc,kc) - dzdeta(ic,jc,kc)*dxdzeta(ic,jc,kc));
           dxidz(ic,jc,kc) =  jaco(ic,jc,kc) * (dxdeta(ic,jc,kc)*dydzeta(ic,jc,kc) - dydeta(ic,jc,kc)*dxdzeta(ic,jc,kc)); 
          detadx(ic,jc,kc) = -jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdzeta(ic,jc,kc) - dzdxi(ic,jc,kc)*dydzeta(ic,jc,kc));  
          detady(ic,jc,kc) =  jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdzeta(ic,jc,kc) - dzdxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
          detadz(ic,jc,kc) = -jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydzeta(ic,jc,kc) - dydxi(ic,jc,kc)*dxdzeta(ic,jc,kc));
         dzetadx(ic,jc,kc) =  jaco(ic,jc,kc) * ( dydxi(ic,jc,kc)*dzdeta(ic,jc,kc) - dzdxi(ic,jc,kc)*dydeta(ic,jc,kc));  
         dzetady(ic,jc,kc) = -jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dzdeta(ic,jc,kc) - dzdxi(ic,jc,kc)*dxdeta(ic,jc,kc));
         dzetadz(ic,jc,kc) =  jaco(ic,jc,kc) * ( dxdxi(ic,jc,kc)*dydeta(ic,jc,kc) - dydxi(ic,jc,kc)*dxdeta(ic,jc,kc));  
         
         nablaxi(ic,jc,kc)=sqrt(  dxidx(ic,jc,kc)  *  dxidx(ic,jc,kc) &
                           &+  dxidy(ic,jc,kc)  *  dxidy(ic,jc,kc) &
                           &+  dxidz(ic,jc,kc)  *  dxidz(ic,jc,kc))
 
        nablaeta(ic,jc,kc)=sqrt( detadx(ic,jc,kc) *  detadx(ic,jc,kc) &
                           &+ detady(ic,jc,kc) *  detady(ic,jc,kc) &
                           &+ detadz(ic,jc,kc) *  detadz(ic,jc,kc))
 
       nablazeta(ic,jc,kc)=sqrt( dzetadx(ic,jc,kc) * dzetadx(ic,jc,kc) &
                           &+ dzetady(ic,jc,kc) * dzetady(ic,jc,kc) &
                           &+ dzetadz(ic,jc,kc) * dzetadz(ic,jc,kc))
         
      enddo
    enddo
enddo
      
  jc = ny_local
  do kc = 1-overlap, nz_local+overlap
     do ic = 1-overlap, nx_local+overlap
      shocknormalx(ic,kc) = detadx(ic,jc,kc) / nablaeta(ic,jc,kc)
      shocknormaly(ic,kc) = detady(ic,jc,kc) / nablaeta(ic,jc,kc)
      shocknormalz(ic,kc) = detadz(ic,jc,kc) / nablaeta(ic,jc,kc)
     enddo
  enddo      

  jc = ny_local ! Standing for the Shock Surface
  do kc = 1-overLap,nz_local+overLap
   do ic = 1-overLap,nx_local+overLap
    !Velocity ahead of the shocks
    Delta_Inf = U_free(ic,jc,kc) * shockNormalx(ic,kc) + &
              & V_free(ic,jc,kc) * shockNormaly(ic,kc) + &
              & W_free(ic,jc,kc) * shockNormalz(ic,kc) - &
              & ShockV(ic,kc) * (WallNormalX(ic,kc)*ShockNormalX(ic,kc)+&
                              &  WallNormalY(ic,kc)*ShockNormalY(ic,kc)+&
                              &  WallNormalZ(ic,kc)*ShockNormalZ(ic,kc));
    Delta2 = (Gamma - 1.d0) * Delta_Inf/(Gamma + 1.d0) + 2.d0 * Gamma * P_free(ic,jc,kc)/ &
          & ((Gamma + 1.d0) * Rho_free(ic,jc,kc) * Delta_Inf);
    !Velocity normal to the shock after the shock
      sU(ic,jc,kc) = U_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalX(ic,kc);
      sV(ic,jc,kc) = V_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalY(ic,kc);
      sW(ic,jc,kc) = W_free(ic,jc,kc) + (Delta2 - Delta_Inf) * ShockNormalZ(ic,kc);
      sP(ic,jc,kc) = 2.0 * Rho_free(ic,jc,kc) * delta_inf * delta_inf/(gamma + 1.0) - &
                   &(gamma - 1.0) * P_free(ic,jc,kc) / (gamma + 1.0);
    sRho(ic,jc,kc) = Delta_inf * Rho_free(ic,jc,kc) / delta2;
      sT(ic,jc,kc) = sP(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / sRho(ic,jc,kc);
   enddo
  enddo
   
  
  deallocate(x_grid,y_grid,z_grid)
  deallocate(shockhdxi,shockhdzeta,shockvdxi,shockvdzeta)
  deallocate(dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta)
  deallocate(jaco,invjacodt,dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz)
  deallocate(nablaxi,nablaeta,nablazeta)
  deallocate(shocknormalx,shocknormaly,shocknormalz)    
  
END SUBROUTINE Linearized_shockBC    
    
subroutine calculate_fdxi_surf_ce_BC(df,f,nx,nz,overlap)
implicit none
integer,intent(in)::nx,nz,overlap
real(8),intent(in)::  f(1-overlap:nx+overlap,1-overlap:nz+overlap)
real(8),intent(out)::df(1-overlap:nx+overlap,1-overlap:nz+overlap)

integer(4)::ic,kc

df = 0.d0;
do kc=1,nz
   do ic=1,nx
       df(ic,kc) = -1.d0*f(ic-3,kc)&
                &+9.d0*f(ic-2,kc)&
               &-45.d0*f(ic-1,kc)&
               &+45.d0*f(ic+1,kc)&
                &-9.d0*f(ic+2,kc)&
                     &+f(ic+3,kc)
       df(ic,kc) = df(ic,kc)/(60.d0)
   enddo
enddo

do kc = 1,nz
       ic = 1
       df(ic,kc) = (-125.d0*f(ic,kc)+240*f(ic+1,kc)-180.d0*f(ic+2,kc)+80.d0*f(ic+3,kc)-15.d0*f(ic+4,kc))/(60.d0)
       ic = 2 
       df(ic,kc) = (-15.d0*f(ic-1,kc)-50.d0*f(ic,kc)+90.d0*f(ic+1,kc)-30.d0*f(ic+2,kc)+5.d0*f(ic+3,kc))/(60.d0)
       ic = 3
       df(ic,kc) = (5.d0*f(ic-2,kc)-40.d0*f(ic-1,kc)+0.d0*f(ic,kc)+40.d0*f(ic+1,kc)-5.d0*f(ic+2,kc))/(60.d0)
enddo


do kc = 1,nz
       ic = nx-2
       df(ic,kc) = (-5.d0*f(ic+2,kc)+40.d0*f(ic+1,kc)+0.d0*f(ic,kc)-40.d0*f(ic-1,kc)+5.d0*f(ic-2,kc))/(60.d0)
       ic = nx-1
       df(ic,kc) = (15.d0*f(ic+1,kc)+50.d0*f(ic,kc)-90.d0*f(ic-1,kc)+30.d0*f(ic-2,kc)-5.d0*f(ic-3,kc))/(60.d0)
       ic = nx           
       df(ic,kc) = (125.d0*f(ic,kc)-240.d0*f(ic-1,kc)+180.d0*f(ic-2,kc)-80.d0*f(ic-3,kc)+15.d0*f(ic-4,kc))/(60.d0)
enddo
 
end subroutine calculate_fdxi_surf_ce_BC
    
subroutine calculate_fdzeta_surf_ce_BC(df,f,nx,nz,overlap)
implicit none
integer,intent(in)::nx,nz,overlap
real(8),intent(in)::f(1-overlap:nx+overlap,1-overlap:nz+overlap)
real(8),intent(out)::df(1-overlap:nx+overlap,1-overlap:nz+overlap)

integer::ic,kc
df = 0.d0;

    do kc=1,nz
        do ic=1,nx
            df(ic,kc) = -1.d0*f(ic,kc-3)+9.d0*f(ic,kc-2)-45.d0*f(ic,kc-1)+45.d0*f(ic,kc+1)-9.d0*f(ic,kc+2)+f(ic,kc+3)
            df(ic,kc) = df(ic,kc)/(60.d0)
        enddo
    enddo

    do ic = 1,nx
        kc = 1
        df(ic,kc) = (-125.d0*f(ic,kc)+240*f(ic,kc+1)-180.d0*f(ic,kc+2)+80.d0*f(ic,kc+3)-15.d0*f(ic,kc+4))/(60.d0)
        kc = 2 
        df(ic,kc) = (-15.d0*f(ic,kc-1)-50.d0*f(ic,kc)+90.d0*f(ic,kc+1)-30.d0*f(ic,kc+2)+5.d0*f(ic,kc+3))/(60.d0)
        kc = 3
        df(ic,kc) = (5.d0*f(ic,kc-2)-40.d0*f(ic,kc-1)+0.d0*f(ic,kc)+40.d0*f(ic,kc+1)-5.d0*f(ic,kc+2))/(60.d0)
    enddo


    do ic = 1,nx
        kc = nz-2
        df(ic,kc) = (-5.d0*f(ic,kc+2)+40.d0*f(ic,kc+1)+0.d0*f(ic,kc)-40.d0*f(ic,kc-1)+5.d0*f(ic,kc-2))/(60.d0)
        kc = nz-1
        df(ic,kc) = (15.d0*f(ic,kc+1)+50.d0*f(ic,kc)-90.d0*f(ic,kc-1)+30.d0*f(ic,kc-2)-5.d0*f(ic,kc-3))/(60.d0)
        kc = nz           
        df(ic,kc) = (125.d0*f(ic,kc)-240.d0*f(ic,kc-1)+180.d0*f(ic,kc-2)-80.d0*f(ic,kc-3)+15.d0*f(ic,kc-4))/(60.d0)
    enddo
end subroutine calculate_fdzeta_surf_ce_BC    