subroutine completeflux_un(ucons0,shockh,shockv,&
                        &rho_free,u_free,v_free,w_free,p_free,&
                        &Rho_pert_tau,U_pert_tau,V_pert_tau,W_pert_tau,P_pert_tau,&
                        &nx_local,ny_local,nz_local,overlap,time,&
                        &gamma,mach_ref,re_ref,csthlnd_ref,pr_ref,laxfrismall,&
                        &finv_inf,ginv_inf,hinv_inf,cv_inf,c_inf,k_infty,epsilon,pert_type,&
                        &wallsx,wallsy,wallsz,&
                        &wallsxdxi,wallsydxi,wallszdxi,wallsxdzeta,wallsydzeta,wallszdzeta,&
                        &wallnormalx,wallnormaly,wallnormalz,&
                        &wallnormalxdxi,wallnormalydxi,wallnormalzdxi,wallnormalxdzeta,wallnormalydzeta,wallnormalzdzeta,&    
                        &heta,hetadxi,hetadeta,hetadzeta,&
                        &r,shockac)
    implicit none
    integer,intent(in) :: nx_local,ny_local,nz_local,overlap,pert_type
    real(8),intent(in) :: rho_inf,u_inf,v_inf,w_inf,p_inf,t_inf,epsilon,k_infty,c_inf
    real(8),intent(in) :: ucons0(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: shockh(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    real(8),intent(in) :: shockv(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)    
    real(8),intent(in) :: gamma,mach_ref,re_ref,csthlnd_ref,pr_ref,laxfrismall    
    real(8),intent(in) :: finv_inf(1:5),ginv_inf(1:5),hinv_inf(1:5),cv_inf(1:5) 
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
    real(8),intent(in) :: Rho_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: U_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: V_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: W_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: P_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: Rho_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: U_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: V_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: W_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(in) :: P_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out):: r(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap)
    real(8),intent(out):: shockac(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap)
    
    integer::i,j,k,l
    real(8),dimension(:,:,:),  allocatable::x_grid,y_grid,z_grid
    real(8),dimension(:,:),    allocatable::shockhdzeta,shockhdxi,shockvdzeta,shockvdxi
    real(8),dimension(:,:,:),  allocatable::dxdxi,dydxi,dzdxi,dxdeta,dydeta,dzdeta,dxdzeta,dydzeta,dzdzeta
    real(8),dimension(:,:,:),  allocatable::jaco,invjacodt
    real(8),dimension(:,:,:),  allocatable::dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz
    real(8),dimension(:,:,:),  allocatable::nablaxi,nablaeta,nablazeta
    real(8),dimension(:,:),    allocatable::shocknormalx,shocknormaly,shocknormalz
    real(8),dimension(:,:,:),  allocatable::dxidt,detadt,dzetadt
    real(8),dimension(:,:),    allocatable::shockxtau,shockytau,shockztau
    real(8),dimension(:,:),    allocatable::shockvn
    real(8),dimension(:,:),    allocatable::ax_tau,ay_tau,az_tau
    real(8),dimension(:,:,:),  allocatable::rho,u,v,w,p,t,cs
    real(8),dimension(:,:,:),  allocatable::bigu_xi,bigu_eta,bigu_zeta
    real(8),dimension(:,:,:,:),allocatable::invf,invg,invh,invfhat,invghat,invhhat
    real(8),dimension(:,:,:,:),allocatable::fp,gp,hp,fm,gm,hm 
    real(8),dimension(:,:,:,:),allocatable::fpdxi,fmdxi,gpdeta,gmdeta,hpdzeta,hmdzeta
    real(8),dimension(:,:,:,:),allocatable::invflux
    real(8),dimension(:,:,:),  allocatable::udxi,vdxi,wdxi,tdxi,udeta,vdeta,wdeta,tdeta,udzeta,vdzeta,wdzeta,tdzeta,mu
    real(8),dimension(:,:,:),  allocatable::udx,vdx,wdx,tdx,udy,vdy,wdy,tdy,udz,vdz,wdz,tdz
    real(8),dimension(:,:,:,:),allocatable::visf,visg,vish,visfhat,visghat,vishhat
    real(8),dimension(:,:,:,:),allocatable::visfhatdxi,visghatdeta,vishhatdzeta
    real(8),dimension(:,:,:,:),allocatable::visflux  
    real(8),dimension(:),      allocatable::l5
    real(8),dimension(:,:,:,:),allocatable::cv_free,fcv_free,gcv_free,hcv_free
    real(8),dimension(:,:,:,:),allocatable::cv_free_tau,fcv_free_tau,gcv_free_tau,hcv_free_tau
    real(8),dimension(:,:,:,:),allocatable::fgcv_free_tau
    !real(8),dimension(:,:,:),  allocatable::rho_free,u_free,v_free,w_free,p_free,t_free
    !real(8),dimension(:,:,:),  allocatable::rho_pert_tau,u_pert_tau,v_pert_tau,w_pert_tau,p_pert_tau
    !real(8),dimension(:,:,:),  allocatable::rho_per,u_per,v_per,w_per,p_per
    !real(8),dimension(:,:,:),  allocatable::i_t,i_t_tau
    
    real(8)::xxitau,yxitau,zxitau,xetatau,yetatau,zetatau,xzetatau,yzetatau,zzetatau
    real(8)::temp,x_tau,y_tau,z_tau
    real(8)::shockxtaudxi,shockytaudxi,shockztaudxi,shockxtaudzeta,shockytaudzeta,shockztaudzeta
    real(8)::tempax,tempay,tempaz
    real(8)::cs_xi,cs_eta,cs_zeta,sigma_xi,sigma_eta,sigma_zeta
    real(8)::local_temp,div_vel,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    real(8)::locrho,locu,locv,locw,locp,locc,locuns,lam5,sum_l5_u,sum_l5_fg,fg,sum_l5_fg_f,fg_f

    !allocate(rho_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(u_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(v_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(w_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(p_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(t_free(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !
    !allocate(rho_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(u_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(v_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(w_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(p_pert_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(rho_per(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(u_per(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(v_per(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(w_per(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(p_per(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !
    !allocate(i_t(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    !allocate(i_t_tau(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    
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
    allocate(dxidt(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(detadt(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(dzetadt(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(shockxtau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockytau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockztau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(shockvn(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(ax_tau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(ay_tau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(az_tau(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap))
    allocate(rho(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(u(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(v(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(w(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(p(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(t(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(cs(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(bigu_xi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(bigu_eta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(bigu_zeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invf(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invg(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invh(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invfhat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invghat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invhhat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fp(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gp(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hp(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fm(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gm(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hm(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fpdxi(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fmdxi(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gpdeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gmdeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hpdzeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hmdzeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(invflux(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdxi(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdzeta(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(mu(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udy(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(udz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdy(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vdz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdy(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(wdz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdx(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdy(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(tdz(1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visf(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visg(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vish(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visfhat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visghat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vishhat(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visfhatdxi(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visghatdeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(vishhatdzeta(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(visflux(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(l5(1:5))
    allocate(cv_free(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fcv_free(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gcv_free(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hcv_free(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(cv_free_tau(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fcv_free_tau(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(gcv_free_tau(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(hcv_free_tau(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    allocate(fgcv_free_tau(1:5,1-overlap:nx_local+overlap,1:ny_local,1-overlap:nz_local+overlap))
    
    do k = 1-overlap,nz_local+overlap
      do j = 1,ny_local
        do i = 1-overlap,nx_local+overlap
        x_grid(i,j,k) = wallsx(i,k) + wallnormalx(i,k) * shockh(i,k) * heta(i,j,k);
        y_grid(i,j,k) = wallsy(i,k) + wallnormaly(i,k) * shockh(i,k) * heta(i,j,k);
        z_grid(i,j,k) = wallsz(i,k) + wallnormalz(i,k) * shockh(i,k) * heta(i,j,k);
      end do
    end do
    end do
    
    !CALL Parallel_Exchange(X_grid)
    !CALL Parallel_Exchange(Y_grid)
    !CALL Parallel_Exchange(Z_grid)
    !
    !CALL Parallel_Exchange_surface(shockH)
    !CALL Parallel_Exchange_surface(shockV)

    call calculate_fdxi_surf_ce(shockhdxi,shockh,nx_local,nz_local,overlap)
    call calculate_fdxi_surf_ce(shockvdxi,shockv,nx_local,nz_local,overlap)
    
    call calculate_fdzeta_surf_ce(shockhdzeta,shockh,nx_local,nz_local,overlap)
    call calculate_fdzeta_surf_ce(shockvdzeta,shockv,nx_local,nz_local,overlap)
   
    !CALL Parallel_Exchange_surface(shockhdxi)
    !CALL Parallel_Exchange_surface(shockvdxi)
    !CALL Parallel_Exchange_surface(shockhdzeta)
    !CALL Parallel_Exchange_surface(shockvdzeta)
 
    do k=1-overlap, nz_local+overlap
        do j=1, ny_local
             do i=1-overlap, nx_local+overlap
            dxdxi(i,j,k) = wallsxdxi(i,k) + (wallnormalxdxi(i,k)*shockh(i,k) + &
                                          & wallnormalx(i,k)*shockhdxi(i,k)) * heta(i,j,k) &
                                          + wallnormalx(i,k) * shockh(i,k) * hetadxi(i,j,k);
            
            xxitau =                       (wallnormalxdxi(i,k)*shockv(i,k) + &
                                          & wallnormalx(i,k)*shockvdxi(i,k)) * heta(i,j,k) &
                                          + wallnormalx(i,k) * shockv(i,k) * hetadxi(i,j,k)
            
            dydxi(i,j,k) = wallsydxi(i,k) + (wallnormalydxi(i,k)*shockh(i,k) + &
                                          & wallnormaly(i,k)*shockhdxi(i,k)) * heta(i,j,k) &
                                          + wallnormaly(i,k) * shockh(i,k) * hetadxi(i,j,k);
            
            yxitau =                       (wallnormalydxi(i,k)*shockv(i,k) + &
                                          & wallnormaly(i,k)*shockvdxi(i,k)) * heta(i,j,k) &
                                          + wallnormaly(i,k) * shockv(i,k) * hetadxi(i,j,k)
            
            dzdxi(i,j,k) = wallszdxi(i,k) + (wallnormalzdxi(i,k)*shockh(i,k) + &
                                          & wallnormalz(i,k)*shockhdxi(i,k)) * heta(i,j,k) &
                                          + wallnormalz(i,k) * shockh(i,k) * hetadxi(i,j,k); 
            
            zxitau =                       (wallnormalzdxi(i,k)*shockv(i,k) + &
                                          & wallnormalz(i,k)*shockvdxi(i,k)) * heta(i,j,k) &
                                          + wallnormalz(i,k) * shockv(i,k) * hetadxi(i,j,k)
            
            dxdeta(i,j,k) = wallnormalx(i,k) * shockh(i,k) * hetadeta(i,j,k);
            dydeta(i,j,k) = wallnormaly(i,k) * shockh(i,k) * hetadeta(i,j,k);
            dzdeta(i,j,k) = wallnormalz(i,k) * shockh(i,k) * hetadeta(i,j,k);
            
            xetatau = wallnormalx(i,k) * shockv(i,k) * hetadeta(i,j,k);
            yetatau = wallnormaly(i,k) * shockv(i,k) * hetadeta(i,j,k);
            zetatau = wallnormalz(i,k) * shockv(i,k) * hetadeta(i,j,k);
            
            dxdzeta(i,j,k) = wallsxdzeta(i,k) + (wallnormalxdzeta(i,k)*shockh(i,k) + &
                                          & wallnormalx(i,k)*shockhdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormalx(i,k) * shockh(i,k) * hetadzeta(i,j,k);
            
            xzetatau =                     (wallnormalxdzeta(i,k)*shockv(i,k) + &
                                          & wallnormalx(i,k)*shockvdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormalx(i,k) * shockv(i,k) * hetadzeta(i,j,k);
            
            dydzeta(i,j,k) = wallsydzeta(i,k) + (wallnormalydzeta(i,k)*shockh(i,k) + &
                                          & wallnormaly(i,k)*shockhdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormaly(i,k) * shockh(i,k) * hetadzeta(i,j,k);
            
            yzetatau =                     (wallnormalydzeta(i,k)*shockv(i,k) + &
                                          & wallnormaly(i,k)*shockvdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormaly(i,k) * shockv(i,k) * hetadzeta(i,j,k);
            
            dzdzeta(i,j,k) = wallszdzeta(i,k) + (wallnormalzdzeta(i,k)*shockh(i,k) + &
                                          & wallnormalz(i,k)*shockhdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormalz(i,k) * shockh(i,k) * hetadzeta(i,j,k);
            
            zzetatau =                     (wallnormalzdzeta(i,k)*shockv(i,k) + &
                                          & wallnormalz(i,k)*shockvdzeta(i,k)) * heta(i,j,k) &
                                          + wallnormalz(i,k) * shockv(i,k) * hetadzeta(i,j,k);
            
            temp =      dxdxi(i,j,k) * (dydeta(i,j,k) * dzdzeta(i,j,k) - dzdeta(i,j,k) * dydzeta(i,j,k)) &
                   & + dxdeta(i,j,k) * (dydzeta(i,j,k)* dzdxi(i,j,k)   - dydxi(i,j,k)  * dzdzeta(i,j,k)) &
                  & + dxdzeta(i,j,k) * (dydxi(i,j,k)  * dzdeta(i,j,k)  - dydeta(i,j,k) * dzdxi(i,j,k))
            
            jaco(i,j,k) = 1.0d0/temp;
            
            invjacodt(i,j,k) = & 
                               &     xxitau * (  dydeta(i,j,k)* dzdzeta(i,j,k) - dzdeta(i,j,k)* dydzeta(i,j,k)) &
                               & +  xetatau * ( dydzeta(i,j,k)*   dzdxi(i,j,k) -  dydxi(i,j,k)* dzdzeta(i,j,k)) &
                               & + xzetatau * (   dydxi(i,j,k) * dzdeta(i,j,k) - dydeta(i,j,k)*   dzdxi(i,j,k)) &
                               & +   dxdxi(i,j,k) * (  yetatau * dzdzeta(i,j,k) - zetatau* dydzeta(i,j,k)) &
                               & +  dxdeta(i,j,k) * (  yzetatau*   dzdxi(i,j,k) -  yxitau* dzdzeta(i,j,k)) &
                               & + dxdzeta(i,j,k) * (   yxitau *  dzdeta(i,j,k) - yetatau*   dzdxi(i,j,k)) &
                               & +   dxdxi(i,j,k) * ( dydeta(i,j,k)* zzetatau - dzdeta(i,j,k)* yzetatau) &
                               & +  dxdeta(i,j,k) * (dydzeta(i,j,k)*   zxitau -  dydxi(i,j,k)* zzetatau) &
                               & + dxdzeta(i,j,k) * (  dydxi(i,j,k)*  zetatau - dydeta(i,j,k)*   zxitau)        
            
              dxidx(i,j,k) =  jaco(i,j,k) * (dydeta(i,j,k)*dzdzeta(i,j,k) - dydzeta(i,j,k)*dzdeta(i,j,k));
              dxidy(i,j,k) = -jaco(i,j,k) * (dxdeta(i,j,k)*dzdzeta(i,j,k) - dzdeta(i,j,k)*dxdzeta(i,j,k));
              dxidz(i,j,k) =  jaco(i,j,k) * (dxdeta(i,j,k)*dydzeta(i,j,k) - dydeta(i,j,k)*dxdzeta(i,j,k)); 
             detadx(i,j,k) = -jaco(i,j,k) * ( dydxi(i,j,k)*dzdzeta(i,j,k) - dzdxi(i,j,k)*dydzeta(i,j,k));  
             detady(i,j,k) =  jaco(i,j,k) * ( dxdxi(i,j,k)*dzdzeta(i,j,k) - dzdxi(i,j,k)*dxdzeta(i,j,k));
             detadz(i,j,k) = -jaco(i,j,k) * ( dxdxi(i,j,k)*dydzeta(i,j,k) - dydxi(i,j,k)*dxdzeta(i,j,k));
            dzetadx(i,j,k) =  jaco(i,j,k) * ( dydxi(i,j,k)*dzdeta(i,j,k) - dzdxi(i,j,k)*dydeta(i,j,k));  
            dzetady(i,j,k) = -jaco(i,j,k) * ( dxdxi(i,j,k)*dzdeta(i,j,k) - dzdxi(i,j,k)*dxdeta(i,j,k));
            dzetadz(i,j,k) =  jaco(i,j,k) * ( dxdxi(i,j,k)*dydeta(i,j,k) - dydxi(i,j,k)*dxdeta(i,j,k));  
            
           nablaxi(i,j,k)=sqrt(  dxidx(i,j,k)  *  dxidx(i,j,k) &
                             &+  dxidy(i,j,k)  *  dxidy(i,j,k) &
                             &+  dxidz(i,j,k)  *  dxidz(i,j,k))
    
           nablaeta(i,j,k)=sqrt( detadx(i,j,k) *  detadx(i,j,k) &
                              &+ detady(i,j,k) *  detady(i,j,k) &
                              &+ detadz(i,j,k) *  detadz(i,j,k))
    
          nablazeta(i,j,k)=sqrt( dzetadx(i,j,k) * dzetadx(i,j,k) &
                              &+ dzetady(i,j,k) * dzetady(i,j,k) &
                              &+ dzetadz(i,j,k) * dzetadz(i,j,k))
            
         enddo
       enddo
   enddo
      
  j = ny_local
  do k = 1-overlap, nz_local+overlap
     do i = 1-overlap, nx_local+overlap
      shocknormalx(i,k) = detadx(i,j,k) / nablaeta(i,j,k)
      shocknormaly(i,k) = detady(i,j,k) / nablaeta(i,j,k)
      shocknormalz(i,k) = detadz(i,j,k) / nablaeta(i,j,k)
     enddo
  enddo         
   
   do k=1-overlap, nz_local+overlap
       do i=1-overlap, nx_local+overlap
           do j=1, ny_local
               x_tau = wallnormalx(i,k) * shockv(i,k) * heta(i,j,k);
               y_tau = wallnormaly(i,k) * shockv(i,k) * heta(i,j,k);
               z_tau = wallnormalz(i,k) * shockv(i,k) * heta(i,j,k);
               
               dxidt(i,j,k) = 0.d0
              detadt(i,j,k) = -( detadx(i,j,k)*x_tau + detady(i,j,k)*y_tau + detadz(i,j,k)*z_tau);
             dzetadt(i,j,k) = 0.d0
           enddo
           shockxtau(i,k) = wallnormalx(i,k) * shockv(i,k);
           shockytau(i,k) = wallnormaly(i,k) * shockv(i,k);
           shockztau(i,k) = wallnormalz(i,k) * shockv(i,k);
        enddo
   enddo
   
   j = ny_local
   do k = 1-overlap, nz_local+overlap
      do i = 1-overlap, nx_local+overlap
   
               shockxtaudxi = wallnormalxdxi(i,k) * shockv(i,k) + wallnormalx(i,k) * shockvdxi(i,k);
               shockytaudxi = wallnormalydxi(i,k) * shockv(i,k) + wallnormaly(i,k) * shockvdxi(i,k);
               shockztaudxi =                                     wallnormalz(i,k) * shockvdxi(i,k);
               
               shockxtaudzeta = wallnormalx(i,k) * shockvdzeta(i,k);
               shockytaudzeta = wallnormaly(i,k) * shockvdzeta(i,k);
               shockztaudzeta = wallnormalz(i,k) * shockvdzeta(i,k);
               
               tempax = shockytaudzeta * dzdxi(i,j,k) + shockztaudxi * dydzeta(i,j,k) - &
                      & shockytaudxi   * dzdzeta(i,j,k) - shockztaudzeta * dydxi(i,j,k);
               tempay = shockxtaudxi   * dzdzeta(i,j,k) + shockztaudzeta * dxdxi(i,j,k) - &
                      & shockxtaudzeta * dzdxi(i,j,k) - shockztaudxi * dxdzeta(i,j,k);
               tempaz = shockxtaudzeta * dydxi(i,j,k) + shockytaudxi * dxdzeta(i,j,k) - &
                      & shockxtaudxi   * dydzeta(i,j,k) - shockytaudzeta * dxdxi(i,j,k);
               
              shockvn(i,k) = -detadt(i,j,k)      /nablaeta(i,j,k);
               ax_tau(i,k) = jaco(i,j,k) * tempax/nablaeta(i,j,k);
               ay_tau(i,k) = jaco(i,j,k) * tempay/nablaeta(i,j,k);
               az_tau(i,k) = jaco(i,j,k) * tempaz/nablaeta(i,j,k);
      enddo
   enddo
   
    ! inviscous term
    do k=1-overlap, nz_local+overlap
        do j=1,ny_local
            do i=1-overlap, nx_local+overlap
               rho(i,j,k) = ucons0(1,i,j,k)
                 u(i,j,k) = ucons0(2,i,j,k)/ucons0(1,i,j,k)
                 v(i,j,k) = ucons0(3,i,j,k)/ucons0(1,i,j,k)
                 w(i,j,k) = ucons0(4,i,j,k)/ucons0(1,i,j,k)
                 p(i,j,k) = (gamma-1.d0)*(ucons0(5,i,j,k)-0.5d0*rho(i,j,k)*(u(i,j,k)**2.d0+v(i,j,k)**2.d0+w(i,j,k)**2.d0))
                 t(i,j,k) = gamma * mach_ref * mach_ref * p(i,j,k) / rho(i,j,k)
                cs(i,j,k) = sqrt(gamma * p(i,j,k) / rho(i,j,k))
                bigu_xi(i,j,k) = u(i,j,k)*  dxidx(i,j,k) &
                              &+ v(i,j,k)*  dxidy(i,j,k) &
                              &+ w(i,j,k)*  dxidz(i,j,k) + dxidt(i,j,k)
               bigu_eta(i,j,k) = u(i,j,k)* detadx(i,j,k) &
                              &+ v(i,j,k)* detady(i,j,k) &
                              &+ w(i,j,k)* detadz(i,j,k) + detadt(i,j,k)
              bigu_zeta(i,j,k) = u(i,j,k)* dzetadx(i,j,k) &
                              &+ v(i,j,k)* dzetady(i,j,k) &
                              &+ w(i,j,k)* dzetadz(i,j,k) + dzetadt(i,j,k)
                
               call utoflux_cpg3d1(rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),&
                                  &invf(1:5,i,j,k),invg(1:5,i,j,k),invh(1:5,i,j,k),gamma)

                cs_xi   = cs(i,j,k) * nablaxi(i,j,k)
                cs_eta  = cs(i,j,k) * nablaeta(i,j,k)
                cs_zeta = cs(i,j,k) * nablazeta(i,j,k)
                
                sigma_xi  =(sqrt( bigu_xi(i,j,k)*  bigu_xi(i,j,k) &
                      &+ laxfrismall * laxfrismall * cs_xi  *  cs_xi) + cs_xi)/jaco(i,j,k);
                
                sigma_eta =(sqrt( bigu_eta(i,j,k)* bigu_eta(i,j,k) &
                      &+ laxfrismall * laxfrismall * cs_eta * cs_eta) + cs_eta)/jaco(i,j,k);

                sigma_zeta=(sqrt(bigu_zeta(i,j,k)*bigu_zeta(i,j,k) &
                      &+ laxfrismall * laxfrismall * cs_zeta*cs_zeta) + cs_zeta)/jaco(i,j,k);
                              
               do l=1,5
                   invfhat(l,i,j,k) = (invf(l,i,j,k) * dxidx(i,j,k) + &
                                       invg(l,i,j,k) * dxidy(i,j,k) + &
                                       invh(l,i,j,k) * dxidz(i,j,k) + &
                                     ucons0(l,i,j,k) * dxidt(i,j,k))/ jaco(i,j,k)
                   
                   invghat(l,i,j,k) = (invf(l,i,j,k) * detadx(i,j,k) + &
                                       invg(l,i,j,k) * detady(i,j,k) + &
                                       invh(l,i,j,k) * detadz(i,j,k) + &
                                     ucons0(l,i,j,k) * detadt(i,j,k))/jaco(i,j,k)
                   
                   invhhat(l,i,j,k) = (invf(l,i,j,k) * dzetadx(i,j,k) + &
                                       invg(l,i,j,k) * dzetady(i,j,k) + &
                                       invh(l,i,j,k) * dzetadz(i,j,k) + &
                                     ucons0(l,i,j,k) * dzetadt(i,j,k))/jaco(i,j,k)
                   
                   fp(l,i,j,k) = 0.5d0 * (invfhat(l,i,j,k) + sigma_xi   * ucons0(l,i,j,k));
                   gp(l,i,j,k) = 0.5d0 * (invghat(l,i,j,k) + sigma_eta  * ucons0(l,i,j,k));
                   hp(l,i,j,k) = 0.5d0 * (invhhat(l,i,j,k) + sigma_zeta * ucons0(l,i,j,k));
                   
                   fm(l,i,j,k) = 0.5d0 * (invfhat(l,i,j,k) - sigma_xi   * ucons0(l,i,j,k));
                   gm(l,i,j,k) = 0.5d0 * (invghat(l,i,j,k) - sigma_eta  * ucons0(l,i,j,k));  
                   hm(l,i,j,k) = 0.5d0 * (invhhat(l,i,j,k) - sigma_zeta * ucons0(l,i,j,k));
               
               enddo
            enddo
        enddo
    enddo

    !call Parallel_Exchange(U)
    !call Parallel_Exchange(V)
    !call Parallel_Exchange(W)
    !call Parallel_Exchange(T)  
    !call Parallel_Exchange_NumVar(Fp)
    !call Parallel_Exchange_NumVar(Fm)
    !call Parallel_Exchange_NumVar(Gp)
    !call Parallel_Exchange_NumVar(Gm)
    !call Parallel_Exchange_NumVar(Hp)
    !call Parallel_Exchange_NumVar(Hm) 
    
    
    call calculate_fdxi_numvar_up(fpdxi,fp,nx_local,ny_local,nz_local,overlap)
    call calculate_fdxi_numvar_do(fmdxi,fm,nx_local,ny_local,nz_local,overlap)
    
    call calculate_fdeta_numvar_up(gpdeta,gp,nx_local,ny_local,nz_local,overlap)
    call calculate_fdeta_numvar_do(gmdeta,gm,nx_local,ny_local,nz_local,overlap)
    
    call calculate_fdzeta_numvar_up_per(hpdzeta,hp,nx_local,ny_local,nz_local,overlap)
    call calculate_fdzeta_numvar_do_per(hmdzeta,hm,nx_local,ny_local,nz_local,overlap)
    
    
    invflux = fpdxi + gpdeta + hpdzeta + fmdxi + gmdeta + hmdzeta;
    !call parallel_exchange_numvar(invflux)

  ! viscous term
    
    call calculate_fdxi_ce(udxi,u,nx_local,ny_local,nz_local,overlap)
    call calculate_fdxi_ce(vdxi,v,nx_local,ny_local,nz_local,overlap)
    call calculate_fdxi_ce(wdxi,w,nx_local,ny_local,nz_local,overlap)
    call calculate_fdxi_ce(tdxi,t,nx_local,ny_local,nz_local,overlap)
    
    call calculate_fdeta_ce(udeta,u,nx_local,ny_local,nz_local,overlap)
    call calculate_fdeta_ce(vdeta,v,nx_local,ny_local,nz_local,overlap)
    call calculate_fdeta_ce(wdeta,w,nx_local,ny_local,nz_local,overlap)
    call calculate_fdeta_ce(tdeta,t,nx_local,ny_local,nz_local,overlap)
    
    call calculate_fdzeta_ce_per(udzeta,u,nx_local,ny_local,nz_local,overlap)
    call calculate_fdzeta_ce_per(vdzeta,v,nx_local,ny_local,nz_local,overlap)
    call calculate_fdzeta_ce_per(wdzeta,w,nx_local,ny_local,nz_local,overlap)
    call calculate_fdzeta_ce_per(tdzeta,t,nx_local,ny_local,nz_local,overlap)
    
    
    !CALL Parallel_Exchange(udxi);   
    !CALL Parallel_Exchange(vdxi);   
    !CALL Parallel_Exchange(wdxi);   
    !CALL Parallel_Exchange(tdxi);   
    !CALL Parallel_Exchange(udeta); 
    !CALL Parallel_Exchange(vdeta);  
    !CALL Parallel_Exchange(wdeta); 
    !CALL Parallel_Exchange(tdeta); 
    !CALL Parallel_Exchange(udzeta); 
    !CALL Parallel_Exchange(vdzeta);
    !CALL Parallel_Exchange(wdzeta); 
    !CALL Parallel_Exchange(tdzeta); 
    
    local_temp = (pr_ref * (gamma - 1.d0) * mach_ref * mach_ref);
    
     do k=1-overlap,nz_local+overlap
        do j=1,ny_local
            do i=1-overlap,nx_local+overlap
               mu(i,j,k) = (t(i,j,k)**1.5d0)*(1.d0+csthlnd_ref)/(t(i,j,k) + csthlnd_ref);
               mu(i,j,k) =  mu(i,j,k)/re_ref;
               
               udx(i,j,k) = udxi(i,j,k)   * dxidx(i,j,k)  + &
                            udeta(i,j,k)  * detadx(i,j,k) + &
                            udzeta(i,j,k) * dzetadx(i,j,k);
               
               vdx(i,j,k) = vdxi(i,j,k)   * dxidx(i,j,k)  + &
                            vdeta(i,j,k)  * detadx(i,j,k) + &
                            vdzeta(i,j,k) * dzetadx(i,j,k);     
               
               wdx(i,j,k) = wdxi(i,j,k)   * dxidx(i,j,k)  + &
                            wdeta(i,j,k)  * detadx(i,j,k) + &
                            wdzeta(i,j,k) * dzetadx(i,j,k);
               
               tdx(i,j,k) = tdxi(i,j,k)   * dxidx(i,j,k)  + &
                            tdeta(i,j,k)  * detadx(i,j,k) + &
                            tdzeta(i,j,k) * dzetadx(i,j,k);
               
               udy(i,j,k) = udxi(i,j,k)   * dxidy(i,j,k)  + &
                            udeta(i,j,k)  * detady(i,j,k) + &
                            udzeta(i,j,k) * dzetady(i,j,k);
               
               vdy(i,j,k) = vdxi(i,j,k)   * dxidy(i,j,k)  + &
                            vdeta(i,j,k)  * detady(i,j,k) + &
                            vdzeta(i,j,k) * dzetady(i,j,k);      
               
               wdy(i,j,k) = wdxi(i,j,k)   * dxidy(i,j,k)  + &
                            wdeta(i,j,k)  * detady(i,j,k) + &
                            wdzeta(i,j,k) * dzetady(i,j,k);
               
               tdy(i,j,k) = tdxi(i,j,k)   * dxidy(i,j,k)  + &
                            tdeta(i,j,k)  * detady(i,j,k) + &
                            tdzeta(i,j,k) * dzetady(i,j,k);
               
               udz(i,j,k) = udxi(i,j,k)   * dxidz(i,j,k)  + &
                            udeta(i,j,k)  * detadz(i,j,k) + &
                            udzeta(i,j,k) * dzetadz(i,j,k);
               
               vdz(i,j,k) = vdxi(i,j,k)   * dxidz(i,j,k)  + &
                            vdeta(i,j,k)  * detadz(i,j,k) + &
                            vdzeta(i,j,k) * dzetadz(i,j,k);     
               
               wdz(i,j,k) = wdxi(i,j,k)   * dxidz(i,j,k)  + &
                            wdeta(i,j,k)  * detadz(i,j,k) + &
                            wdzeta(i,j,k) * dzetadz(i,j,k);
               
               tdz(i,j,k) = tdxi(i,j,k)   * dxidz(i,j,k)  + &
                            tdeta(i,j,k)  * detadz(i,j,k) + &
                            tdzeta(i,j,k) * dzetadz(i,j,k);
               
               div_vel = udx(i,j,k) + vdy(i,j,k) + wdz(i,j,k);
               
               tauxx = 2.d0 * mu(i,j,k) * (udx(i,j,k) - 1.d0 / 3.d0 * div_vel);
               tauyy = 2.d0 * mu(i,j,k) * (vdy(i,j,k) - 1.d0 / 3.d0 * div_vel);
               tauzz = 2.d0 * mu(i,j,k) * (wdz(i,j,k) - 1.d0 / 3.d0 * div_vel);
               tauxy = mu(i,j,k) * (udy(i,j,k) + vdx(i,j,k));
               tauxz = mu(i,j,k) * (udz(i,j,k) + wdx(i,j,k));
               tauyz = mu(i,j,k) * (vdz(i,j,k) + wdy(i,j,k));
               
               visf(1,i,j,k) = 0.d0;
               visf(2,i,j,k) = tauxx;
               visf(3,i,j,k) = tauxy;
               visf(4,i,j,k) = tauxz;
               visf(5,i,j,k) = u(i,j,k) * tauxx &
                            &+ v(i,j,k) * tauxy &
                            &+ w(i,j,k) * tauxz &
                            &+ mu(i,j,k) * tdx(i,j,k) / local_temp;
               
               visg(1,i,j,k) = 0.d0;
               visg(2,i,j,k) = tauxy;
               visg(3,i,j,k) = tauyy;
               visg(4,i,j,k) = tauyz;
               visg(5,i,j,k) = u(i,j,k) * tauxy &
                            &+ v(i,j,k) * tauyy &
                            &+ w(i,j,k) * tauyz &
                            &+ mu(i,j,k) * tdy(i,j,k) / local_temp;
               
               vish(1,i,j,k) = 0.d0;
               vish(2,i,j,k) = tauxz;
               vish(3,i,j,k) = tauyz;
               vish(4,i,j,k) = tauzz;
               vish(5,i,j,k) = u(i,j,k) * tauxz &
                            &+ v(i,j,k) * tauyz &
                            &+ w(i,j,k) * tauzz &
                            &+ mu(i,j,k) * tdz(i,j,k) / local_temp;
              
              do l = 2,5
               visfhat(l,i,j,k) =(visf(l,i,j,k)*  dxidx(i,j,k) +&
                                  visg(l,i,j,k)*  dxidy(i,j,k) +&
                                  vish(l,i,j,k)*  dxidz(i,j,k))/jaco(i,j,k);

               visghat(l,i,j,k) =(visf(l,i,j,k)* detadx(i,j,k) +&
                                  visg(l,i,j,k)* detady(i,j,k) +&
                                  vish(l,i,j,k)* detadz(i,j,k))/jaco(i,j,k);

               vishhat(l,i,j,k) =(visf(l,i,j,k)*dzetadx(i,j,k) +&
                                  visg(l,i,j,k)*dzetady(i,j,k) +&
                                  vish(l,i,j,k)*dzetadz(i,j,k))/jaco(i,j,k);
              enddo
            enddo
        enddo
     enddo

      !call parallel_exchange_numvar(visfhat)   
      !call parallel_exchange_numvar(visghat)
      !call parallel_exchange_numvar(vishhat)
      
     call calculate_fdxi_numvar_ce(visfhatdxi,visfhat,nx_local,ny_local,nz_local,overlap)
     call calculate_fdeta_numvar_ce(visghatdeta,visghat,nx_local,ny_local,nz_local,overlap)
     call calculate_fdzeta_numvar_ce_per(vishhatdzeta,vishhat,nx_local,ny_local,nz_local,overlap)
     
     visflux = visfhatdxi + visghatdeta + vishhatdzeta;
     !call parallel_exchange_numvar(visflux)
     
   do k = 1-overlap,nz_local+overlap
       do j = 1,ny_local
         do i = 1-overlap,nx_local+overlap
             do l = 1,5
                 r(l,i,j,k) = (-invflux(l,i,j,k) + visflux(l,i,j,k) - ucons0(l,i,j,k) * invjacodt(i,j,k)) * jaco(i,j,k);
             enddo
         enddo
       enddo
   enddo

   
 !j = ny_local;
 !     do k = 1-overlap,nz_local+overlap
 !       do i = 1-overlap,nx_local+overlapLN
 !
 !   if (pert_type .eq. 1)then !fast 
 !         i_t(i,j,k) = cos(k_infty*((x_grid(i,j,k)-1.d0) - c_inf * time));
 !     i_t_tau(i,j,k) = -k_infty*(shockxtau(i,k)-c_inf)*sin(k_infty*((x_grid(i,j,k)-1.d0)-c_inf * time));
 !      
 !         rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);  
 !           u_per(i,j,k)    = epsilon * i_t(i,j,k);
 !           v_per(i,j,k)    = 0.d0;
 !           w_per(i,j,k)    = 0.d0;
 !           p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
 !          
 !       rho_pert_tau(i,j,k) = epsilon * mach_ref * i_t_tau(i,j,k);
 !         u_pert_tau(i,j,k) = epsilon * i_t_tau(i,j,k);
 !         v_pert_tau(i,j,k) = 0.d0;
 !         w_pert_tau(i,j,k) = 0.d0;
 !         p_pert_tau(i,j,k) = epsilon * mach_ref * gamma * i_t_tau(i,j,k)/ ( gamma * mach_ref * mach_ref );
 !           
 !          rho_free(i,j,k)  = rho_inf + rho_per(i,j,k);
 !            u_free(i,j,k)  = u_inf   + u_per(i,j,k);
 !            v_free(i,j,k)  = v_inf   + v_per(i,j,k);
 !            w_free(i,j,k)  = w_inf   + w_per(i,j,k);
 !            p_free(i,j,k)  = p_inf   + p_per(i,j,k);
 !            t_free(i,j,k)  = gamma * mach_ref * mach_ref * p_free(i,j,k) / rho_free(i,j,k);
 !         elseif (pert_type .eq. 2)then !slow
 !         rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
 !           u_per(i,j,k)    = - epsilon * i_t(i,j,k);
 !           v_per(i,j,k)    = 0.d0;
 !           w_per(i,j,k)    = 0.d0;            
 !           p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
 !         elseif (pert_type .eq. 3)then
 !         rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
 !           u_per(i,j,k)    = 0.d0;
 !           v_per(i,j,k)    = 0.d0;
 !           w_per(i,j,k)    = 0.d0;   
 !           p_per(i,j,k)    = 0.d0;
 !         !if (pert_type .eq. 4)then
 !         else
 !         rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
 !           u_per(i,j,k)    = 0.d0;
 !           v_per(i,j,k)    = epsilon * i_t(i,j,k);
 !           w_per(i,j,k)    = 0.d0;            
 !           p_per(i,j,k)    = 0.d0;                
 !         endif    
 !       enddo
 !     enddo  

  j = ny_local
  do k = 1-overlap,nz_local+overlap
   do i = 1-overlap,nx_local+overlap
    locrho = rho(i,j,k);
    locu   =   u(i,j,k);
    locv   =   v(i,j,k);
    locw   =   w(i,j,k);
    locp   =   p(i,j,k);
    locc   = sqrt(gamma * locp / locrho);
    ! calculate the local normal
    locuns = locu * shocknormalx(i,k) + locv * shocknormaly(i,k) + locw * shocknormalz(i,k);
    lam5   = locuns - shockvn(i,k) + locc;

    l5(1) = 0.5d0 * (gamma - 1.d0)*( locu * locu + locv * locv + locw * locw ) - locuns * locc; 
    l5(2) =        -(gamma - 1.d0)* locu + shocknormalx(i,k) * locc;
    l5(3) =        -(gamma - 1.d0)* locv + shocknormaly(i,k) * locc;
    l5(4) =        -(gamma - 1.d0)* locw + shocknormalz(i,k) * locc;
    l5(5) =         (gamma - 1.d0);

    sum_l5_u  = 0.d0;
    sum_l5_fg = 0.d0;
    sum_l5_fg_f = 0.d0;

    call utocv_cpg1(rho_free(i,j,k),u_free(i,j,k),v_free(i,j,k),w_free(i,j,k),p_free(i,j,k),cv_free(:,i,j,k),gamma)
    call utoflux_cpg3d1(rho_free(i,j,k),u_free(i,j,k),v_free(i,j,k),w_free(i,j,k),p_free(i,j,k),&
         &fcv_free(:,i,j,k),gcv_free(:,i,j,k),hcv_free(:,i,j,k),gamma)
    call utocv_tau1(rho_free(i,j,k),u_free(i,j,k),v_free(i,j,k),w_free(i,j,k),&
           &rho_pert_tau(i,j,k),u_pert_tau(i,j,k),v_pert_tau(i,j,k),w_pert_tau(i,j,k),p_pert_tau(i,j,k),cv_free_tau(:,i,j,k),gamma)
    call utoflux_tau1(rho_free(i,j,k),u_free(i,j,k),v_free(i,j,k),w_free(i,j,k),p_free(i,j,k),&
           &rho_pert_tau(i,j,k),u_pert_tau(i,j,k),v_pert_tau(i,j,k),w_pert_tau(i,j,k),p_pert_tau(i,j,k),&
           &fcv_free_tau(:,i,j,k),gcv_free_tau(:,i,j,k),hcv_free_tau(:,i,j,k),gamma)
    
    do l = 1,5
      fg = lam5 * r(l,i,j,k) + ax_tau(i,k) * (invf(l,i,j,k) - fcv_free(l,i,j,k)) &
                             + ay_tau(i,k) * (invg(l,i,j,k) - gcv_free(l,i,j,k)) &
                             + az_tau(i,k) * (invh(l,i,j,k) - hcv_free(l,i,j,k));

      fgcv_free_tau(l,i,j,k) = (shocknormalx(i,k)* fcv_free_tau(l,i,j,k)) + &
                             & (shocknormaly(i,k)* gcv_free_tau(l,i,j,k)) + &
                             & (shocknormalz(i,k)* hcv_free_tau(l,i,j,k));
      fg_f = fgcv_free_tau(l,i,j,k) - shockvn(i,k) * cv_free_tau(l,i,j,k);
      
      sum_l5_u    = sum_l5_u    + l5(l) * (ucons0(l,i,j,k) - cv_free(l,i,j,k));
      sum_l5_fg   = sum_l5_fg   + l5(l) * fg;     
      sum_l5_fg_f = sum_l5_fg_f + l5(l) * fg_f;
    enddo

    shockac(i,k) = sum_l5_fg/sum_l5_u - (ax_tau(i,k) * shockxtau(i,k) &
                                     &+  ay_tau(i,k) * shockytau(i,k) &
                                     &+  az_tau(i,k) * shockztau(i,k))&
                                     &-  (sum_l5_fg_f/sum_l5_u);

    shockac(i,k) = shockac(i,k)/( wallnormalx(i,k) * shocknormalx(i,k) &
                               &+ wallnormaly(i,k) * shocknormaly(i,k) &
                               &+ wallnormalz(i,k) * shocknormalz(i,k));
   enddo
  enddo

  
  deallocate(x_grid,y_grid,z_grid)
  deallocate(shockhdxi,shockhdzeta,shockvdxi,shockvdzeta)
  deallocate(dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta)
  deallocate(jaco,invjacodt,dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz)
  deallocate(nablaxi,nablaeta,nablazeta)
  deallocate(shocknormalx,shocknormaly,shocknormalz)
  deallocate(dxidt,detadt,dzetadt)
  deallocate(shockxtau,shockytau,shockztau)
  deallocate(shockvn,ax_tau,ay_tau,az_tau)
  deallocate(rho,u,v,w,p,t,cs)
  deallocate(bigu_xi,bigu_eta,bigu_zeta)
  deallocate(invf,invg,invh,invfhat,invghat,invhhat)
  deallocate(fp,fm,gp,gm,hp,hm)
  deallocate(fpdxi,fmdxi,gpdeta,gmdeta,hpdzeta,hmdzeta)
  deallocate(invflux)
  deallocate(udxi,udeta,udzeta,vdxi,vdeta,vdzeta,wdxi,wdeta,wdzeta,tdxi,tdeta,tdzeta,mu)
  deallocate(udx,udy,udz,vdx,vdy,vdz,wdx,wdy,wdz,tdx,tdy,tdz,visf,visg,vish,visfhat,visghat,vishhat)
  deallocate(visfhatdxi,visghatdeta,vishhatdzeta,visflux)
  deallocate(l5)
  deallocate(cv_free,fcv_free,gcv_free,hcv_free)
  deallocate(cv_free_tau,fcv_free_tau,gcv_free_tau,hcv_free_tau)
  deallocate(fgcv_free_tau)
  
end subroutine completeflux_un
 
    
subroutine utoflux_cpg3d1(rho,u,v,w,p,f,g,h,gamma)
    implicit none
    real(8),intent(in)::rho,u,v,w,p,gamma
    real(8),intent(out)::f(1:5),g(1:5),h(1:5)

    f(1) = rho * u
    f(2) = rho * u * u + p
    f(3) = rho * u * v
    f(4) = rho * u * w
    f(5) = (p/(gamma - 1.d0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * u

    g(1) = rho * v
    g(2) = rho * v * u
    g(3) = rho * v * v + p
    g(4) = rho * v * w
    g(5) = (p/(gamma - 1.d0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * v

    h(1) = rho * w
    h(2) = rho * w * u
    h(3) = rho * w * v
    h(4) = rho * w * w + p
    h(5) = (p/(gamma - 1.d0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * w

end subroutine utoflux_cpg3d1        
    
subroutine calculate_fdxi_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k
    real(8),parameter::alpha=0.d0
    
    df = 0.d0;
    do k=1,nz
        do j=1,ny
            do i=1,nx
                df(i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(i-3,j,k)+&
                             &(9.d0-1.d0/2.d0*alpha)*f(i-2,j,k)+&
                           &(-45.d0+5.d0/4.d0*alpha)*f(i-1,j,k)+&
                                 &(-5.d0/3.d0*alpha)*f(i,j,k)  +&
                            &(45.d0+5.d0/4.d0*alpha)*f(i+1,j,k)+&
                            &(-9.d0-1.d0/2.d0*alpha)*f(i+2,j,k)+&
                            &(1.d0+1.d0/12.d0*alpha)*f(i+3,j,k)
                df(i,j,k) = df(i,j,k)/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = 1
            df(i,j,k) = (-125.d0*f(i,j,k)+240*f(i+1,j,k)-180.d0*f(i+2,j,k)+80.d0*f(i+3,j,k)-15.d0*f(i+4,j,k))/(60.d0)
            i = 2 
            df(i,j,k) = (-15.d0*f(i-1,j,k)-50.d0*f(i,j,k)+90.d0*f(i+1,j,k)-30.d0*f(i+2,j,k)+5.d0*f(i+3,j,k))/(60.d0)
            i = 3
            df(i,j,k) = (5.d0*f(i-2,j,k)-40.d0*f(i-1,j,k)+0.d0*f(i,j,k)+40.d0*f(i+1,j,k)-5.d0*f(i+2,j,k))/(60.d0)
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = nx-2
            df(i,j,k) = (-5.d0*f(i+2,j,k)+40.d0*f(i+1,j,k)+0.d0*f(i,j,k)-40.d0*f(i-1,j,k)+5.d0*f(i-2,j,k))/(60.d0)
            i = nx-1
            df(i,j,k) = (15.d0*f(i+1,j,k)+50.d0*f(i,j,k)-90.d0*f(i-1,j,k)+30.d0*f(i-2,j,k)-5.d0*f(i-3,j,k))/(60.d0)
            i = nx           
            df(i,j,k) = (125.d0*f(i,j,k)-240.d0*f(i-1,j,k)+180.d0*f(i-2,j,k)-80.d0*f(i-3,j,k)+15.d0*f(i-4,j,k))/(60.d0)
        enddo
    enddo

end subroutine calculate_fdxi_ce
    
subroutine calculate_fdeta_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k
    real(8),parameter::alpha=0.d0
    df = 0.d0;    
    do k=1,nz
        do j=4,ny-3
            do i=1,nx
            
                df(i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(i,j-3,k)&
                            &+(9.d0-1.d0/2.d0*alpha)*f(i,j-2,k)+&
                           &(-45.d0+5.d0/4.d0*alpha)*f(i,j-1,k)+&
                                 &(-5.d0/3.d0*alpha)*f(i,j,k)+&
                            &(45.d0+5.d0/4.d0*alpha)*f(i,j+1,k)+&
                            &(-9.d0-1.d0/2.d0*alpha)*f(i,j+2,k)+&
                            &(1.d0+1.d0/12.d0*alpha)*f(i,j+3,k)
                df(i,j,k) = df(i,j,k)/(60.d0)
                
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = 1
            df(i,j,k) = (-125.d0*f(i,j,k)+240*f(i,j+1,k)-180.d0*f(i,j+2,k)+80.d0*f(i,j+3,k)-15.d0*f(i,j+4,k))/(60.d0)
            j = 2 
            df(i,j,k) = (-15.d0*f(i,j-1,k)-50.d0*f(i,j,k)+90.d0*f(i,j+1,k)-30.d0*f(i,j+2,k)+5.d0*f(i,j+3,k))/(60.d0)
            j = 3
            df(i,j,k) = (5.d0*f(i,j-2,k)-40.d0*f(i,j-1,k)+0.d0*f(i,j,k)+40.d0*f(i,j+1,k)-5.d0*f(i,j+2,k))/(60.d0)
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = ny-2
            df(i,j,k) = (-5.d0*f(i,j+2,k)+40.d0*f(i,j+1,k)+0.d0*f(i,j,k)-40.d0*f(i,j-1,k)+5.d0*f(i,j-2,k))/(60.d0)
            j = ny-1
            df(i,j,k) = (15.d0*f(i,j+1,k)+50.d0*f(i,j,k)-90.d0*f(i,j-1,k)+30.d0*f(i,j-2,k)-5.d0*f(i,j-3,k))/(60.d0)
            j = ny           
            df(i,j,k) = (125.d0*f(i,j,k)-240.d0*f(i,j-1,k)+180.d0*f(i,j-2,k)-80.d0*f(i,j-3,k)+15.d0*f(i,j-4,k))/(60.d0)
        enddo
    enddo

end subroutine calculate_fdeta_ce
    
    
subroutine calculate_fdzeta_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::  f(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k
    real(8),parameter::alpha=0.d0
    df = 0.d0;
    do k=1,nz
       do j=1,ny
            do i=1,nx
            
                df(i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(i,j,k-3)+&
                             &(9.d0-1.d0/2.d0*alpha)*f(i,j,k-2)+&
                           &(-45.d0+5.d0/4.d0*alpha)*f(i,j,k-1)+&
                                 &(-5.d0/3.d0*alpha)*f(i,j,k)+&
                            &(45.d0+5.d0/4.d0*alpha)*f(i,j,k+1)+&
                            &(-9.d0-1.d0/2.d0*alpha)*f(i,j,k+2)+&
                            &(1.d0+1.d0/12.d0*alpha)*f(i,j,k+3)
                df(i,j,k) = df(i,j,k)/(60.d0)
                
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = 1
            df(i,j,k) = (-125.d0*f(i,j,k)+240*f(i,j,k+1)-180.d0*f(i,j,k+2)+80.d0*f(i,j,k+3)-15.d0*f(i,j,k+4))/(60.d0)
            k = 2 
            df(i,j,k) = (-15.d0*f(i,j,k-1)-50.d0*f(i,j,k)+90.d0*f(i,j,k+1)-30.d0*f(i,j,k+2)+5.d0*f(i,j,k+3))/(60.d0)
            k = 3
            df(i,j,k) = (5.d0*f(i,j,k-2)-40.d0*f(i,j,k-1)+0.d0*f(i,j,k)+40.d0*f(i,j,k+1)-5.d0*f(i,j,k+2))/(60.d0)
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = nz-2
            df(i,j,k) = (-5.d0*f(i,j,k+2)+40.d0*f(i,j,k+1)+0.d0*f(i,j,k)-40.d0*f(i,j,k-1)+5.d0*f(i,j,k-2))/(60.d0)
            k = nz-1
            df(i,j,k) = (15.d0*f(i,j,k+1)+50.d0*f(i,j,k)-90.d0*f(i,j,k-1)+30.d0*f(i,j,k-2)-5.d0*f(i,j,k-3))/(60.d0)
            k = nz           
            df(i,j,k) = (125.d0*f(i,j,k)-240.d0*f(i,j,k-1)+180.d0*f(i,j,k-2)-80.d0*f(i,j,k-3)+15.d0*f(i,j,k-4))/(60.d0)
        enddo
    enddo

end subroutine calculate_fdzeta_ce
    
subroutine calculate_fdxi_surf_ce(df,f,nx,nz,overlap)
    implicit none
    integer,intent(in)::nx,nz,overlap
    real(8),intent(in)::  f(1-overlap:nx+overlap,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1-overlap:nz+overlap)
 
    integer(4)::i,k
    
    df = 0.d0;
    do k=1,nz
       do i=1,nx
           df(i,k) = -1.d0*f(i-3,k)&
                    &+9.d0*f(i-2,k)&
                   &-45.d0*f(i-1,k)&
                   &+45.d0*f(i+1,k)&
                    &-9.d0*f(i+2,k)&
                         &+f(i+3,k)
           df(i,k) = df(i,k)/(60.d0)
       enddo
    enddo
    
    do k = 1,nz
           i = 1
           df(i,k) = (-125.d0*f(i,k)+240*f(i+1,k)-180.d0*f(i+2,k)+80.d0*f(i+3,k)-15.d0*f(i+4,k))/(60.d0)
           i = 2 
           df(i,k) = (-15.d0*f(i-1,k)-50.d0*f(i,k)+90.d0*f(i+1,k)-30.d0*f(i+2,k)+5.d0*f(i+3,k))/(60.d0)
           i = 3
           df(i,k) = (5.d0*f(i-2,k)-40.d0*f(i-1,k)+0.d0*f(i,k)+40.d0*f(i+1,k)-5.d0*f(i+2,k))/(60.d0)
    enddo

    
    do k = 1,nz
           i = nx-2
           df(i,k) = (-5.d0*f(i+2,k)+40.d0*f(i+1,k)+0.d0*f(i,k)-40.d0*f(i-1,k)+5.d0*f(i-2,k))/(60.d0)
           i = nx-1
           df(i,k) = (15.d0*f(i+1,k)+50.d0*f(i,k)-90.d0*f(i-1,k)+30.d0*f(i-2,k)-5.d0*f(i-3,k))/(60.d0)
           i = nx           
           df(i,k) = (125.d0*f(i,k)-240.d0*f(i-1,k)+180.d0*f(i-2,k)-80.d0*f(i-3,k)+15.d0*f(i-4,k))/(60.d0)
    enddo
 
end subroutine calculate_fdxi_surf_ce
    
subroutine calculate_fdzeta_surf_ce(df,f,nx,nz,overlap)
    implicit none
    integer,intent(in)::nx,nz,overlap
    real(8),intent(in)::f(1-overlap:nx+overlap,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1-overlap:nz+overlap)
 
    integer::i,k
    df = 0.d0;
    
        do k=1,nz
            do i=1,nx
                df(i,k) = -1.d0*f(i,k-3)+9.d0*f(i,k-2)-45.d0*f(i,k-1)+45.d0*f(i,k+1)-9.d0*f(i,k+2)+f(i,k+3)
                df(i,k) = df(i,k)/(60.d0)
            enddo
        enddo
    
        do i = 1,nx
            k = 1
            df(i,k) = (-125.d0*f(i,k)+240*f(i,k+1)-180.d0*f(i,k+2)+80.d0*f(i,k+3)-15.d0*f(i,k+4))/(60.d0)
            k = 2 
            df(i,k) = (-15.d0*f(i,k-1)-50.d0*f(i,k)+90.d0*f(i,k+1)-30.d0*f(i,k+2)+5.d0*f(i,k+3))/(60.d0)
            k = 3
            df(i,k) = (5.d0*f(i,k-2)-40.d0*f(i,k-1)+0.d0*f(i,k)+40.d0*f(i,k+1)-5.d0*f(i,k+2))/(60.d0)
        enddo

    
        do i = 1,nx
            k = nz-2
            df(i,k) = (-5.d0*f(i,k+2)+40.d0*f(i,k+1)+0.d0*f(i,k)-40.d0*f(i,k-1)+5.d0*f(i,k-2))/(60.d0)
            k = nz-1
            df(i,k) = (15.d0*f(i,k+1)+50.d0*f(i,k)-90.d0*f(i,k-1)+30.d0*f(i,k-2)-5.d0*f(i,k-3))/(60.d0)
            k = nz           
            df(i,k) = (125.d0*f(i,k)-240.d0*f(i,k-1)+180.d0*f(i,k-2)-80.d0*f(i,k-3)+15.d0*f(i,k-4))/(60.d0)
        enddo


end subroutine calculate_fdzeta_surf_ce
    
subroutine calculate_fdxi_numvar_up(df,f,nx,ny,nz,overlap)
     implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=-6.d0
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i-3,j,k)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i-2,j,k)+&
                                &(-45.d0+1.25d0*alpha)*f(l,i-1,j,k)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                                 &(45.d0+1.25d0*alpha)*f(l,i+1,j,k)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i+2,j,k)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i+3,j,k)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i+1,j,k)-180.d0*f(l,i+2,j,k)+80.d0*f(l,i+3,j,k)-15.d0*f(l,i+4,j,k))/(60.d0)
            enddo
            i = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i-1,j,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i+1,j,k)-30.d0*f(l,i+2,j,k)+5.d0*f(l,i+3,j,k))/(60.d0)
            enddo
            i = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i-2,j,k)-40.d0*f(l,i-1,j,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i+1,j,k)-5.d0*f(l,i+2,j,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = nx-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i+2,j,k)+40.d0*f(l,i+1,j,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i-1,j,k)+5.d0*f(l,i-2,j,k))/(60.d0)
            enddo
            i = nx-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i+1,j,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i-1,j,k)+30.d0*f(l,i-2,j,k)-5.d0*f(l,i-3,j,k))/(60.d0)
            enddo
            i = nx  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i-1,j,k)+180.d0*f(l,i-2,j,k)-80.d0*f(l,i-3,j,k)+15.d0*f(l,i-4,j,k))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdxi_numvar_up
    
subroutine calculate_fdxi_numvar_do(df,f,nx,ny,nz,overlap)
     implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=6.d0
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                  df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i-3,j,k)+&
                                     &(9.d0-0.5d0*alpha)*f(l,i-2,j,k)+&
                                  &(-45.d0+1.25d0*alpha)*f(l,i-1,j,k)+&
                                     &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                                   &(45.d0+1.25d0*alpha)*f(l,i+1,j,k)+&
                                    &(-9.d0-0.5d0*alpha)*f(l,i+2,j,k)+&
                                &(1.d0+1.d0/12.d0*alpha)*f(l,i+3,j,k)
                  df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i+1,j,k)-180.d0*f(l,i+2,j,k)+80.d0*f(l,i+3,j,k)-15.d0*f(l,i+4,j,k))/(60.d0)
            enddo
            i = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i-1,j,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i+1,j,k)-30.d0*f(l,i+2,j,k)+5.d0*f(l,i+3,j,k))/(60.d0)
            enddo
            i = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i-2,j,k)-40.d0*f(l,i-1,j,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i+1,j,k)-5.d0*f(l,i+2,j,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = nx-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i+2,j,k)+40.d0*f(l,i+1,j,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i-1,j,k)+5.d0*f(l,i-2,j,k))/(60.d0)
            enddo
            i = nx-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i+1,j,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i-1,j,k)+30.d0*f(l,i-2,j,k)-5.d0*f(l,i-3,j,k))/(60.d0)
            enddo
            i = nx  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i-1,j,k)+180.d0*f(l,i-2,j,k)-80.d0*f(l,i-3,j,k)+15.d0*f(l,i-4,j,k))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdxi_numvar_do    
    
subroutine calculate_fdeta_numvar_up(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::  f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=-6.d0
    df = 0.d0;
    
    do k=1,nz
        do j=4,ny-3
           do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j-3,k)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j-2,k)+&
                                &(-45.d0+1.25d0*alpha)*f(l,i,j-1,k)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                                 &(45.d0+1.25d0*alpha)*f(l,i,j+1,k)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i,j+2,k)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j+3,k)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j+1,k)-180.d0*f(l,i,j+2,k)+80.d0*f(l,i,j+3,k)-15.d0*f(l,i,j+4,k))/(60.d0)
            enddo
            j = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j-1,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j+1,k)-30.d0*f(l,i,j+2,k)+5.d0*f(l,i,j+3,k))/(60.d0)
            enddo
            j = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j-2,k)-40.d0*f(l,i,j-1,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j+1,k)-5.d0*f(l,i,j+2,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = ny-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j+2,k)+40.d0*f(l,i,j+1,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j-1,k)+5.d0*f(l,i,j-2,k))/(60.d0)
            enddo
            j = ny-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j+1,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j-1,k)+30.d0*f(l,i,j-2,k)-5.d0*f(l,i,j-3,k))/(60.d0)
            enddo
            j = ny  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j-1,k)+180.d0*f(l,i,j-2,k)-80.d0*f(l,i,j-3,k)+15.d0*f(l,i,j-4,k))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdeta_numvar_up 
    
subroutine calculate_fdeta_numvar_do(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=6.d0
    df = 0.d0;
    do k=1,nz
        do j=4,ny-3
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j-3,k)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j-2,k)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j-1,k)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j+1,k)+&
                              &(-9.d0-1.d0/2.d0*alpha)*f(l,i,j+2,k)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j+3,k)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j+1,k)-180.d0*f(l,i,j+2,k)+80.d0*f(l,i,j+3,k)-15.d0*f(l,i,j+4,k))/(60.d0)
            enddo
            j = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j-1,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j+1,k)-30.d0*f(l,i,j+2,k)+5.d0*f(l,i,j+3,k))/(60.d0)
            enddo
            j = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j-2,k)-40.d0*f(l,i,j-1,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j+1,k)-5.d0*f(l,i,j+2,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = ny-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j+2,k)+40.d0*f(l,i,j+1,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j-1,k)+5.d0*f(l,i,j-2,k))/(60.d0)
            enddo
            j = ny-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j+1,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j-1,k)+30.d0*f(l,i,j-2,k)-5.d0*f(l,i,j-3,k))/(60.d0)
            enddo
            j = ny  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j-1,k)+180.d0*f(l,i,j-2,k)-80.d0*f(l,i,j-3,k)+15.d0*f(l,i,j-4,k))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdeta_numvar_do 
    
subroutine calculate_fdzeta_numvar_up(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=-6.d0
    df = 0.d0;
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                                &(-45.d0+1.25d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                                 &(45.d0+1.25d0*alpha)*f(l,i,j,k+1)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j,k+1)-180.d0*f(l,i,j,k+2)+80.d0*f(l,i,j,k+3)-15.d0*f(l,i,j,k+4))/(60.d0)
            enddo
            k = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j,k-1)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j,k+1)-30.d0*f(l,i,j,k+2)+5.d0*f(l,i,j,k+3))/(60.d0)
            enddo
            k = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j,k-2)-40.d0*f(l,i,j,k-1)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j,k+1)-5.d0*f(l,i,j,k+2))/(60.d0)
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = nz-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j,k+2)+40.d0*f(l,i,j,k+1)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j,k-1)+5.d0*f(l,i,j,k-2))/(60.d0)
            enddo
            k = nz-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j,k+1)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j,k-1)+30.d0*f(l,i,j,k-2)-5.d0*f(l,i,j,k-3))/(60.d0)
            enddo
            k = nz  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j,k-1)+180.d0*f(l,i,j,k-2)-80.d0*f(l,i,j,k-3)+15.d0*f(l,i,j,k-4))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdzeta_numvar_up 
    
subroutine calculate_fdzeta_numvar_do(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=6.d0
    
    df = 0.d0;
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j,k+1)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j,k+1)-180.d0*f(l,i,j,k+2)+80.d0*f(l,i,j,k+3)-15.d0*f(l,i,j,k+4))/(60.d0)
            enddo
            k = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j,k-1)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j,k+1)-30.d0*f(l,i,j,k+2)+5.d0*f(l,i,j,k+3))/(60.d0)
            enddo
            k = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j,k-2)-40.d0*f(l,i,j,k-1)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j,k+1)-5.d0*f(l,i,j,k+2))/(60.d0)
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = nz-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j,k+2)+40.d0*f(l,i,j,k+1)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j,k-1)+5.d0*f(l,i,j,k-2))/(60.d0)
            enddo
            k = nz-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j,k+1)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j,k-1)+30.d0*f(l,i,j,k-2)-5.d0*f(l,i,j,k-3))/(60.d0)
            enddo
            k = nz  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j,k-1)+180.d0*f(l,i,j,k-2)-80.d0*f(l,i,j,k-3)+15.d0*f(l,i,j,k-4))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdzeta_numvar_do
    
    
subroutine calculate_fdxi_numvar_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=0.d0
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i-3,j,k)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i-2,j,k)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i-1,j,k)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i+1,j,k)+&
                              &(-9.d0-1.d0/2.d0*alpha)*f(l,i+2,j,k)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i+3,j,k)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i+1,j,k)-180.d0*f(l,i+2,j,k)+80.d0*f(l,i+3,j,k)-15.d0*f(l,i+4,j,k))/(60.d0)
            enddo
            i = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i-1,j,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i+1,j,k)-30.d0*f(l,i+2,j,k)+5.d0*f(l,i+3,j,k))/(60.d0)
            enddo
            i = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i-2,j,k)-40.d0*f(l,i-1,j,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i+1,j,k)-5.d0*f(l,i+2,j,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do j = 1,ny
            i = nx-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i+2,j,k)+40.d0*f(l,i+1,j,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i-1,j,k)+5.d0*f(l,i-2,j,k))/(60.d0)
            enddo
            i = nx-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i+1,j,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i-1,j,k)+30.d0*f(l,i-2,j,k)-5.d0*f(l,i-3,j,k))/(60.d0)
            enddo
            i = nx  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i-1,j,k)+180.d0*f(l,i-2,j,k)-80.d0*f(l,i-3,j,k)+15.d0*f(l,i-4,j,k))/(60.d0)
            enddo
        enddo
    enddo
end subroutine calculate_fdxi_numvar_ce
    
subroutine calculate_fdeta_numvar_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=0.d0
    df = 0.d0;
    do k=1,nz
        do j=4,ny-3
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j-3,k)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j-2,k)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j-1,k)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j+1,k)+&
                              &(-9.d0-1.d0/2.d0*alpha)*f(l,i,j+2,k)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j+3,k)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j+1,k)-180.d0*f(l,i,j+2,k)+80.d0*f(l,i,j+3,k)-15.d0*f(l,i,j+4,k))/(60.d0)
            enddo
            j = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j-1,k)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j+1,k)-30.d0*f(l,i,j+2,k)+5.d0*f(l,i,j+3,k))/(60.d0)
            enddo
            j = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j-2,k)-40.d0*f(l,i,j-1,k)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j+1,k)-5.d0*f(l,i,j+2,k))/(60.d0)
            enddo
        enddo
    enddo
    
    do k = 1,nz
        do i = 1,nx
            j = ny-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j+2,k)+40.d0*f(l,i,j+1,k)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j-1,k)+5.d0*f(l,i,j-2,k))/(60.d0)
            enddo
            j = ny-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j+1,k)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j-1,k)+30.d0*f(l,i,j-2,k)-5.d0*f(l,i,j-3,k))/(60.d0)
            enddo
            j = ny  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j-1,k)+180.d0*f(l,i,j-2,k)-80.d0*f(l,i,j-3,k)+15.d0*f(l,i,j-4,k))/(60.d0)
            enddo
        enddo
    enddo

end subroutine calculate_fdeta_numvar_ce
    

    subroutine calculate_fdzeta_ce_per(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k
    real(8),parameter::alpha=0.d0
    df = 0.d0;
    do k=1,nz
       do j=1,ny
            do i=1,nx
            
                df(i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(i,j,k-3)+&
                             &(9.d0-1.d0/2.d0*alpha)*f(i,j,k-2)+&
                           &(-45.d0+5.d0/4.d0*alpha)*f(i,j,k-1)+&
                                 &(-5.d0/3.d0*alpha)*f(i,j,k  )+&
                            &(45.d0+5.d0/4.d0*alpha)*f(i,j,k+1)+&
                            &(-9.d0-1.d0/2.d0*alpha)*f(i,j,k+2)+&
                            &(1.d0+1.d0/12.d0*alpha)*f(i,j,k+3)
                df(i,j,k) = df(i,j,k)/(60.d0)
                
            enddo
        enddo
    enddo

    end subroutine calculate_fdzeta_ce_per
    
    subroutine calculate_fdzeta_up_per(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k
    real(8),parameter::alpha=-6.d0
    df = 0.d0;
    do k=1,nz
       do j=1,ny
            do i=1,nx
            
                df(i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(i,j,k-3)+&
                             &(9.d0-1.d0/2.d0*alpha)*f(i,j,k-2)+&
                           &(-45.d0+5.d0/4.d0*alpha)*f(i,j,k-1)+&
                                 &(-5.d0/3.d0*alpha)*f(i,j,k  )+&
                            &(45.d0+5.d0/4.d0*alpha)*f(i,j,k+1)+&
                            &(-9.d0-1.d0/2.d0*alpha)*f(i,j,k+2)+&
                            &(1.d0+1.d0/12.d0*alpha)*f(i,j,k+3)
                df(i,j,k) = df(i,j,k)/(60.d0)
                
            enddo
        enddo
    enddo
    end subroutine calculate_fdzeta_up_per 
    
    subroutine calculate_fdzeta_numvar_ce_per(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=0.d0
    df = 0.d0;
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k  )+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j,k+1)+&
                              &(-9.d0-1.d0/2.d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    end subroutine calculate_fdzeta_numvar_ce_per
    
    subroutine calculate_fdzeta_numvar_up_per(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=-6.d0
    df = 0.d0;
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k  )+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j,k+1)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    end subroutine calculate_fdzeta_numvar_up_per
    
    
    subroutine calculate_fdzeta_numvar_do_per(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=6.d0
    df = 0.d0;
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k  )+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j,k+1)+&
                              &(-9.d0-1.d0/2.d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    end subroutine calculate_fdzeta_numvar_do_per
    
    
    subroutine calculate_fdzeta_numvar_ce(df,f,nx,ny,nz,overlap)
    implicit none
    integer,intent(in)::nx,ny,nz,overlap
    real(8),intent(in)::f(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
    real(8),intent(out)::df(1:5,1-overlap:nx+overlap,1:ny,1-overlap:nz+overlap)
 
    integer::i,j,k,l
    real(8),parameter::alpha=0.d0
    
    df = 0.d0;
    do k=1,nz
        do j=1,ny
            do i=1,nx
                do l=1,5
                df(l,i,j,k) = (-1.d0+1.d0/12.d0*alpha)*f(l,i,j,k-3)+&
                                   &(9.d0-0.5d0*alpha)*f(l,i,j,k-2)+&
                             &(-45.d0+5.d0/4.d0*alpha)*f(l,i,j,k-1)+&
                                   &(-5.d0/3.d0*alpha)*f(l,i,j,k)+&
                              &(45.d0+5.d0/4.d0*alpha)*f(l,i,j,k+1)+&
                                  &(-9.d0-0.5d0*alpha)*f(l,i,j,k+2)+&
                              &(1.d0+1.d0/12.d0*alpha)*f(l,i,j,k+3)
                df(l,i,j,k) = df(l,i,j,k)/(60.d0)
                enddo
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = 1
            do l=1,5
            df(l,i,j,k) = (-125.d0*f(l,i,j,k)+240*f(l,i,j,k+1)-180.d0*f(l,i,j,k+2)+80.d0*f(l,i,j,k+3)-15.d0*f(l,i,j,k+4))/(60.d0)
            enddo
            k = 2 
            do l=1,5
            df(l,i,j,k) = (-15.d0*f(l,i,j,k-1)-50.d0*f(l,i,j,k)+90.d0*f(l,i,j,k+1)-30.d0*f(l,i,j,k+2)+5.d0*f(l,i,j,k+3))/(60.d0)
            enddo
            k = 3
            do l=1,5
            df(l,i,j,k) = (5.d0*f(l,i,j,k-2)-40.d0*f(l,i,j,k-1)+0.d0*f(l,i,j,k)+40.d0*f(l,i,j,k+1)-5.d0*f(l,i,j,k+2))/(60.d0)
            enddo
        enddo
    enddo
    
    do j = 1,ny
        do i = 1,nx
            k = nz-2
            do l=1,5
            df(l,i,j,k) = (-5.d0*f(l,i,j,k+2)+40.d0*f(l,i,j,k+1)+0.d0*f(l,i,j,k)-40.d0*f(l,i,j,k-1)+5.d0*f(l,i,j,k-2))/(60.d0)
            enddo
            k = nz-1
            do l=1,5
            df(l,i,j,k) = (15.d0*f(l,i,j,k+1)+50.d0*f(l,i,j,k)-90.d0*f(l,i,j,k-1)+30.d0*f(l,i,j,k-2)-5.d0*f(l,i,j,k-3))/(60.d0)
            enddo
            k = nz  
            do l=1,5
            df(l,i,j,k) = (125.d0*f(l,i,j,k)-240.d0*f(l,i,j,k-1)+180.d0*f(l,i,j,k-2)-80.d0*f(l,i,j,k-3)+15.d0*f(l,i,j,k-4))/(60.d0)
            enddo
        enddo
    enddo
    end subroutine calculate_fdzeta_numvar_ce
    
    subroutine utocv_cpg1(rho,u,v,w,p,cv,gamma)
    ! transform primitive variables to consevative variables
    implicit none
    real( kind = 8 ),intent(in)::rho,u,v,w,p,gamma
    real( kind = 8 ),intent(out)::cv(1:5)

    cv(1)=rho
    cv(2)=rho*u
    cv(3)=rho*v
    cv(4)=rho*w
    cv(5)=p/(gamma-1.d0)+0.5d0*rho*(u*u+v*v+w*w)

end subroutine utocv_cpg1

! this subroutine is used to calculate the flux in 3d
    subroutine utocv_tau1(rho,u,v,w,rho_tau,u_tau,v_tau,w_tau,p_tau,cv,gamma)
    implicit none
    real(8),intent(in)::rho,u,v,w,rho_tau,u_tau,v_tau,w_tau,p_tau,gamma
    real(8),intent(out)::cv(5)
    
     cv(1)=rho_tau
     cv(2)=rho * u_tau + u * rho_tau
     cv(3)=rho * v_tau + v * rho_tau
     cv(4)=rho * w_tau + w * rho_tau
     cv(5)=p_tau/(gamma-1.d0)+0.5d0*rho_tau*(u*u+v*v+w*w)+rho*(u_tau*u+v_tau*v+w_tau*w)
     
    end subroutine utocv_tau1
    
    subroutine utoflux_tau1(rho,u,v,w,p,rho_tau,u_tau,v_tau,w_tau,p_tau,f,g,h,gamma)
    implicit none
    real(8),intent(in)::rho,u,v,w,p,rho_tau,u_tau,v_tau,w_tau,p_tau,gamma
    real(8),intent(out)::f(5),g(5),h(5)
  
     f(1) = rho * u_tau + u * rho_tau
     f(2) = rho_tau * u * u + rho * 2 * u * u_tau + p_tau
     f(3) = rho_tau * u * v + rho * u_tau * v + rho * u * v_tau
     f(4) = rho_tau * u * w + rho * u_tau * w + rho * u * w_tau
     f(5) = (p_tau/(gamma - 1.0) + 0.5d0 * rho_tau * (u*u + v*v + w*w) + rho *(u_tau*u+v_tau*v+w_tau*w) + p_tau) * u &
         &+ (p/(gamma - 1.0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * u_tau
     
     g(1) = rho * v_tau + v * rho_tau
     g(2) = rho_tau * v * u + rho * v_tau * u + rho * v * u_tau
     g(3) = rho_tau * v * v + rho * 2 * v * v_tau + p_tau
     g(4) = rho_tau * v * w + rho * v_tau * w + rho * v * w_tau
     g(5) = (p_tau/(gamma - 1.0) + 0.5d0 * rho_tau * (u*u + v*v + w*w) + rho *(u_tau*u+v_tau*v+w_tau*w) + p_tau) * v &
         &+ (p/(gamma - 1.0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * v_tau
  
     h(1) = rho * w_tau + w * rho_tau
     h(2) = rho_tau * w * u + rho * w_tau * u + rho * w * u_tau
     h(3) = rho_tau * w * v + rho * w_tau * v + rho * w * v_tau
     h(4) = rho_tau * w * w + rho * 2 * w * w_tau + p_tau
     h(5) = (p_tau/(gamma - 1.0) + 0.5d0 * rho_tau * (u*u + v*v + w*w) + rho *(u_tau*u+v_tau*v+w_tau*w) + p_tau) * w &
         &+ (p/(gamma - 1.0) + 0.5d0 * rho * (u*u + v*v + w*w) + p) * w_tau
     
     end subroutine utoflux_tau1