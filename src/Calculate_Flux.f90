SUBROUTINE Calculate_Flux
  use SF_Constant,           only: ik,rk,NumVar,LaxFriSmall,overLAP,fifth_order,first_order
  use SF_CFD_Global,         only: nx_local,Ny,nz_local,Gamma,Pr_Ref,Mach_Ref,Re_Ref,&
                                   Csthlnd_Ref,Cp_Ref,Select_Scheme,ModelType
                                   

  ! matrix and some arrays
  use SF_CFD_Global,         only: Rho,U,V,W,P,T,Ucons0,Cs,BigU_Xi,BigU_Eta,BigU_Zeta,&
                                   dxidx,dxidy,dxidz,dxidt,detadx,detady,detadz,detadt,&
                                   dzetadx,dzetady,dzetadz,dzetadt,invF,invG,invH,&
                                   Udxi,Udeta,Udzeta,Vdxi,Vdeta,Vdzeta,Wdxi,Wdeta,Wdzeta,&
                                   Tdxi,Tdeta,Tdzeta,Mu,Udx,Udy,Udz,Vdx,Vdy,Vdz,Wdx,Wdy,Wdz,&
                                   Tdx,Tdy,Tdz,visF,visG,visH,visFhat,visGhat,visHhat,Jaco,invJacodt,&
                                   nablaxi,nablaeta,nablazeta,invFhat,invGhat,invHhat,Fp,Gp,Hp,&
                                   Fm,Gm,Hm,Fpdxi,Fmdxi,Gpdeta,Gmdeta,Hpdzeta,Hmdzeta,visFcdxi,&
                                   visGcdeta,visHcdzeta,invFlux,visFlux,dUcons0,Cs_Xi,Cs_Eta,Cs_Zeta,&
                                   U_total,V_total,W_total,T_total,&
                                   Udxi_total,Udeta_total,Udzeta_total,Vdxi_total,Vdeta_total,Vdzeta_total,&
                                   Wdxi_total,Wdeta_total,Wdzeta_total,Tdxi_total,Tdeta_total,Tdzeta_total,&
                                   visFhat_total,visGhat_total,visHhat_total,visFcdxi_total,visGcdeta_total,visHcdzeta_total,&
                                   visFlux_total,dUcons0_total,invFhatdxi_new

  ! subroutine and function declaration
  use MPI_GLOBAL,            only: Parallel_Exchange,Parallel_Exchange_NumVar,npx,npx0
  use FD5_Order,             only: Cal_Deri_Dxi_5th_ce,Cal_Deri_Deta_5th_ce,Cal_Deri_Dzeta_5th_ce_per,& ! Central
              Cal_Deri_Dxi_5th_ce_numvar,Cal_Deri_Deta_5th_ce_numvar,Cal_Deri_Dzeta_5th_ce_per_numvar,& ! Central
              Cal_Deri_Dxi_5th_up_numvar,Cal_Deri_Deta_5th_up_numvar,Cal_Deri_Dzeta_5th_up_per_numvar,& ! Upwind
              Cal_Deri_Dxi_5th_do_numvar,Cal_Deri_Deta_5th_do_numvar,Cal_Deri_Dzeta_5th_do_per_numvar ! Downwind

  use FD1st_order,           only: Cal_Deri_Dxi_2nd_ce,Cal_Deri_Deta_2nd_ce,Cal_Deri_Dzeta_2nd_ce_per,&
              Cal_Deri_Dxi_1st_up_numvar,Cal_Deri_Deta_1st_up_numvar,Cal_Deri_Dzeta_1st_up_per_numvar,&
              Cal_Deri_Dxi_2nd_ce_numvar,Cal_Deri_Deta_2nd_ce_numvar,Cal_Deri_Dzeta_2nd_ce_per_numvar,&
              Cal_Deri_Dxi_1st_do_numvar,Cal_Deri_Deta_1st_do_numvar,Cal_Deri_Dzeta_1st_do_per_numvar
  implicit none
  ! Local variables
  integer( kind = ik ):: ic,jc,kc,iVar
  real( kind = rk ):: Div_Vel,Tauxx,Tauyy,Tauzz,Tauxy,Tauxz,Tauyz
  real( kind = rk ):: Sigma_Xi,Sigma_Eta,Sigma_Zeta

  real( kind = rk ):: local_temp

  ! Get the inner variables
  do kc = 1-overLAP, nz_local+overLAP
   do jc = 1, Ny
    do ic = 1-overLAP, nx_local+overLAP
  !do kc = 1, nz_local
  ! do jc = 1, ny
  !  do ic = 1, nx_local
      ! conservative variables to primitive variables
      call CVtoU_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
      ! Calculate the temperature
      T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc)
      ! Calculate the sound speed
      Cs(ic,jc,kc)= sqrt(Gamma * P(ic,jc,kc) / Rho(ic,jc,kc))
      !! Calculate the characteristic speed
      !  BigU_Xi(ic,jc,kc) = U(ic,jc,kc)*  dxidx(ic,jc,kc) &
      !                   &+ V(ic,jc,kc)*  dxidy(ic,jc,kc) &
      !                   &+ W(ic,jc,kc)*  dxidz(ic,jc,kc) + dxidt(ic,jc,kc)

      ! BigU_Eta(ic,jc,kc) = U(ic,jc,kc)* detadx(ic,jc,kc) &
      !                   &+ V(ic,jc,kc)* detady(ic,jc,kc) &
      !                   &+ W(ic,jc,kc)* detadz(ic,jc,kc) + detadt(ic,jc,kc)

      !BigU_Zeta(ic,jc,kc) = U(ic,jc,kc)*dzetadx(ic,jc,kc) &
      !                   &+ V(ic,jc,kc)*dzetady(ic,jc,kc) &
      !                   &+ W(ic,jc,kc)*dzetadz(ic,jc,kc) + dzetadt(ic,jc,kc)

      ! Calculate the flux
      call UtoFlux_CPG3D(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),&
                         invF(1:NumVar,ic,jc,kc),invG(1:NumVar,ic,jc,kc),invH(1:NumVar,ic,jc,kc))
      
    enddo
   enddo
  enddo
  
  ! Parallel transfer
if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then
 call Parallel_Exchange_singularity(U,U_total)
 call Parallel_Exchange_singularity(V,V_total)
 call Parallel_Exchange_singularity(W,W_total)
 call Parallel_Exchange_singularity(T,T_total)
endif    
  

  call Parallel_Exchange(U)
  call Parallel_Exchange(V)
  call Parallel_Exchange(W)
  call Parallel_Exchange(T)
  
  select case(Select_Scheme)
  case(fifth_order)  
    ! Calculate the viscous Flux by using the central difference
      
    call Cal_Deri_Dxi_5th_ce(Udxi,U)
    call Cal_Deri_Dxi_5th_ce(Vdxi,V)
    call Cal_Deri_Dxi_5th_ce(Wdxi,W)
    call Cal_Deri_Dxi_5th_ce(Tdxi,T)
    
    call Cal_Deri_Deta_5th_ce(Udeta,U)
    call Cal_Deri_Deta_5th_ce(Vdeta,V)
    call Cal_Deri_Deta_5th_ce(Wdeta,W)
    call Cal_Deri_Deta_5th_ce(Tdeta,T)
    
    call Cal_Deri_Dzeta_5th_ce_per(Udzeta,U)
    call Cal_Deri_Dzeta_5th_ce_per(Vdzeta,V)
    call Cal_Deri_Dzeta_5th_ce_per(Wdzeta,W)
    call Cal_Deri_Dzeta_5th_ce_per(Tdzeta,T)

   case(first_order)
        ! Calculate the viscous Flux by using the central difference
    call Cal_Deri_Dxi_2nd_ce(Udxi,U)
    call Cal_Deri_Dxi_2nd_ce(Vdxi,V)
    call Cal_Deri_Dxi_2nd_ce(Wdxi,W)
    call Cal_Deri_Dxi_2nd_ce(Tdxi,T)
  
    call Cal_Deri_Deta_2nd_ce(Udeta,U)
    call Cal_Deri_Deta_2nd_ce(Vdeta,V)
    call Cal_Deri_Deta_2nd_ce(Wdeta,W)
    call Cal_Deri_Deta_2nd_ce(Tdeta,T)
    ! Here we use the periodic boundary condition form
    call Cal_Deri_Dzeta_2nd_ce_per(Udzeta,U)
    call Cal_Deri_Dzeta_2nd_ce_per(Vdzeta,V)
    call Cal_Deri_Dzeta_2nd_ce_per(Wdzeta,W)
    call Cal_Deri_Dzeta_2nd_ce_per(Tdzeta,T)
  
    
    
   end select
if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then
  call Parallel_Exchange_singularity_negative(Udxi,Udxi_total); call Parallel_Exchange_singularity(Udeta,Udeta_total); call Parallel_Exchange_singularity(Udzeta,Udzeta_total)
  call Parallel_Exchange_singularity_negative(Vdxi,Vdxi_total); call Parallel_Exchange_singularity(Vdeta,Vdeta_total); call Parallel_Exchange_singularity(Vdzeta,Vdzeta_total)
  call Parallel_Exchange_singularity_negative(Wdxi,Wdxi_total); call Parallel_Exchange_singularity(Wdeta,Wdeta_total); call Parallel_Exchange_singularity(Wdzeta,Wdzeta_total)
  call Parallel_Exchange_singularity_negative(Tdxi,Tdxi_total); call Parallel_Exchange_singularity(Tdeta,Tdeta_total); call Parallel_Exchange_singularity(Tdzeta,Tdzeta_total)
endif  
  call Parallel_Exchange(Udxi); call Parallel_Exchange(Udeta); call Parallel_Exchange(Udzeta)
  call Parallel_Exchange(Vdxi); call Parallel_Exchange(Vdeta); call Parallel_Exchange(Vdzeta)
  call Parallel_Exchange(Wdxi); call Parallel_Exchange(Wdeta); call Parallel_Exchange(Wdzeta)
  call Parallel_Exchange(Tdxi); call Parallel_Exchange(Tdeta); call Parallel_Exchange(Tdzeta)      
  
  ! Cp_Ref = 1.0_rk / ((Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
  local_temp = (Pr_Ref * (Gamma - 1.0_rk) * Mach_Ref * Mach_Ref);
  ! Calculate the viscous stress
  do kc = 1-overLAP,nz_local+overLAP
   do jc = 1,Ny
    do ic = 1-overLAP,nx_local+overLAP
  !do kc = 1,nz_local
  ! do jc = 1,Ny
   ! do ic = 1,nx_local       
     ! The viscosity coefficient has been divided by Re_Ref already
     ! T^1.5 * (1 + Csthlnd_Ref) / (T + Csthlnd_Ref) / Re_Ref 
      Mu(ic,jc,kc) = (T(ic,jc,kc)**1.5_rk) * ( 1.0_rk + Csthlnd_Ref )/ (T(ic,jc,kc) + Csthlnd_Ref) 
      Mu(ic,jc,kc) = Mu(ic,jc,kc) / Re_Ref;
     ! derivatives with respect to the x direction
      Udx(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
                       Udeta(ic,jc,kc)* detadx(ic,jc,kc) +&
                      Udzeta(ic,jc,kc)*dzetadx(ic,jc,kc);

      Vdx(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
                       Vdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
                      Vdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);      

      Wdx(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
                       Wdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
                      Wdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);

      Tdx(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidx(ic,jc,kc) +&
                       Tdeta(ic,jc,kc)* detadx(ic,jc,kc) +&
                      Tdzeta(ic,jc,kc)*dzetadx(ic,jc,kc);

      ! derivatives with respect to the y direction
      Udy(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
                       Udeta(ic,jc,kc)* detady(ic,jc,kc) +&
                      Udzeta(ic,jc,kc)*dzetady(ic,jc,kc);

      Vdy(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
                       Vdeta(ic,jc,kc)* detady(ic,jc,kc) +&
                      Vdzeta(ic,jc,kc)*dzetady(ic,jc,kc);      

      Wdy(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
                       Wdeta(ic,jc,kc)* detady(ic,jc,kc) +&
                      Wdzeta(ic,jc,kc)*dzetady(ic,jc,kc);

      Tdy(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidy(ic,jc,kc) +&
                       Tdeta(ic,jc,kc)* detady(ic,jc,kc) +&
                      Tdzeta(ic,jc,kc)*dzetady(ic,jc,kc);
      
      ! derivatives with respect to the z direction
      Udz(ic,jc,kc) =   Udxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
                       Udeta(ic,jc,kc)* detadz(ic,jc,kc) +&
                      Udzeta(ic,jc,kc)*dzetadz(ic,jc,kc);

      Vdz(ic,jc,kc) =   Vdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
                       Vdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
                      Vdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);      

      Wdz(ic,jc,kc) =   Wdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
                       Wdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
                      Wdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);

      Tdz(ic,jc,kc) =   Tdxi(ic,jc,kc)*  dxidz(ic,jc,kc) +&
                       Tdeta(ic,jc,kc)* detadz(ic,jc,kc) +&
                      Tdzeta(ic,jc,kc)*dzetadz(ic,jc,kc);
                      
      Div_Vel = Udx(ic,jc,kc) + Vdy(ic,jc,kc) + Wdz(ic,jc,kc);

      Tauxx = 2.0_rk * Mu(ic,jc,kc) * (Udx(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauyy = 2.0_rk * Mu(ic,jc,kc) * (Vdy(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauzz = 2.0_rk * Mu(ic,jc,kc) * (Wdz(ic,jc,kc) - 1.0_rk / 3.0_rk * Div_Vel);
      Tauxy = Mu(ic,jc,kc) * (Udy(ic,jc,kc) + Vdx(ic,jc,kc));
      Tauxz = Mu(ic,jc,kc) * (Udz(ic,jc,kc) + Wdx(ic,jc,kc));
      Tauyz = Mu(ic,jc,kc) * (Vdz(ic,jc,kc) + Wdy(ic,jc,kc));

      ! Calculate the viscous flux in the x-,y-,z- directions
      visF(1,ic,jc,kc) = 0.0_rk;
      visF(2,ic,jc,kc) = Tauxx;
      visF(3,ic,jc,kc) = Tauxy;
      visF(4,ic,jc,kc) = Tauxz;
      visF(5,ic,jc,kc) = U(ic,jc,kc) * Tauxx &
                      &+ V(ic,jc,kc) * Tauxy &
                      &+ W(ic,jc,kc) * Tauxz &
                      &+ Mu(ic,jc,kc) * Tdx(ic,jc,kc) / local_temp;
      
      visG(1,ic,jc,kc) = 0.0_rk;
      visG(2,ic,jc,kc) = Tauxy;
      visG(3,ic,jc,kc) = Tauyy;
      visG(4,ic,jc,kc) = Tauyz;
      visG(5,ic,jc,kc) = U(ic,jc,kc) * Tauxy &
                      &+ V(ic,jc,kc) * Tauyy &
                      &+ W(ic,jc,kc) * Tauyz &
                      &+ Mu(ic,jc,kc) * Tdy(ic,jc,kc) / local_temp;
      
      visH(1,ic,jc,kc) = 0.0_rk;
      visH(2,ic,jc,kc) = Tauxz;
      visH(3,ic,jc,kc) = Tauyz;
      visH(4,ic,jc,kc) = Tauzz;
      visH(5,ic,jc,kc) = U(ic,jc,kc) * Tauxz &
                      &+ V(ic,jc,kc) * Tauyz &
                      &+ W(ic,jc,kc) * Tauzz &
                      &+ Mu(ic,jc,kc) * Tdz(ic,jc,kc) / local_temp;
      enddo
   enddo
  enddo    

      
 do kc = 1,nz_local
   do jc = 1,Ny
    do ic = 1,nx_local 
 !  do kc = 1-overLAP,nz_local+overLAP
 !  do jc = 1,Ny
 !   do ic = 1-overLAP,nx_local+overLAP 
      do iVar = 2,NumVar
        visFhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*  dxidx(ic,jc,kc) +&
                                 visG(iVar,ic,jc,kc)*  dxidy(ic,jc,kc) +&
                                 visH(iVar,ic,jc,kc)*  dxidz(ic,jc,kc))/Jaco(ic,jc,kc);

        visGhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)* detadx(ic,jc,kc) +&
                                 visG(iVar,ic,jc,kc)* detady(ic,jc,kc) +&
                                 visH(iVar,ic,jc,kc)* detadz(ic,jc,kc))/Jaco(ic,jc,kc);

        visHhat(iVar,ic,jc,kc) =(visF(iVar,ic,jc,kc)*dzetadx(ic,jc,kc) +&
                                 visG(iVar,ic,jc,kc)*dzetady(ic,jc,kc) +&
                                 visH(iVar,ic,jc,kc)*dzetadz(ic,jc,kc))/Jaco(ic,jc,kc);
      enddo
    enddo
   enddo
  enddo

if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then
   call Parallel_Exchange_singularity_numvar(visFhat,visFhat_total)
   call Parallel_Exchange_singularity_numvar_negative_1(visGhat,visGhat_total)
   call Parallel_Exchange_singularity_numvar_negative_1(visHhat,visHhat_total)
endif

   call Parallel_Exchange_NumVar(visFhat)
   call Parallel_Exchange_NumVar(visGhat)
   call Parallel_Exchange_NumVar(visHhat) 
   
   select case(Select_Scheme)
   case(fifth_order)
 
     call Cal_Deri_Dxi_5th_ce_numvar(visFcdxi,        visFhat)   
     call Cal_Deri_Deta_5th_ce_numvar(visGcdeta,      visGhat)
     call Cal_Deri_Dzeta_5th_ce_per_numvar(visHcdzeta,visHhat)
   case(first_order)
   
     call Cal_Deri_Dxi_2nd_ce_numvar(visFcdxi,        visFhat)
     call Cal_Deri_Deta_2nd_ce_numvar(visGcdeta,      visGhat)
     call Cal_Deri_Dzeta_2nd_ce_per_numvar(visHcdzeta,visHhat)
    end select
 
    
    call Parallel_Exchange_NumVar(visFcdxi)
    call Parallel_Exchange_NumVar(visGcdeta)
    call Parallel_Exchange_NumVar(visHcdzeta)

   visFlux = visFcdxi + visGcdeta + visHcdzeta;

   
 !if (ModelType == 2) then 
 !  call Parallel_Exchange_singularity_numvar_negative_1(visFlux,visFlux_total)
 !endif     
 
   call Parallel_Exchange_NumVar(visFlux)
   
  !==============================================================================================================================
  ! Calculate the inviscid flux
   
  do kc = 1-overLAP,nz_local+overLAP
   do jc = 1,Ny
    do ic = 1-overLAP,nx_local+overLAP
  !do kc = 1,nz_local
  ! do jc = 1,Ny
  !  do ic = 1,nx_local 
       do iVar = 1,NumVar
        invFhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *   dxidx(ic,jc,kc) +&
                                 invG(iVar,ic,jc,kc) *   dxidy(ic,jc,kc) +&
                                 invH(iVar,ic,jc,kc) *   dxidz(ic,jc,kc) +&
                               Ucons0(iVar,ic,jc,kc) *   dxidt(ic,jc,kc)) / Jaco(ic,jc,kc);
        
        invGhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *  detadx(ic,jc,kc) +&
                                 invG(iVar,ic,jc,kc) *  detady(ic,jc,kc) +&
                                 invH(iVar,ic,jc,kc) *  detadz(ic,jc,kc) +&
                               Ucons0(iVar,ic,jc,kc) *  detadt(ic,jc,kc)) / Jaco(ic,jc,kc);
        
        invHhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) * dzetadx(ic,jc,kc) +&
                                 invG(iVar,ic,jc,kc) * dzetady(ic,jc,kc) +&
                                 invH(iVar,ic,jc,kc) * dzetadz(ic,jc,kc) +&
                               Ucons0(iVar,ic,jc,kc) * dzetadt(ic,jc,kc)) / Jaco(ic,jc,kc);
       enddo
    enddo
   enddo
  enddo

  call Parallel_Exchange_NumVar(invFhat)
  call Parallel_Exchange_NumVar(invGhat)
  call Parallel_Exchange_NumVar(invHhat)
  
  do kc = 1-overLAP,nz_local+overLAP
   do jc = 1,Ny
    do ic = 1-overLAP,nx_local+overLAP   
  !do kc = 1,nz_local
  ! do jc = 1,Ny
  !  do ic = 1,nx_local   
        
      Cs_Xi(ic,jc,kc)   = Cs(ic,jc,kc) *   nablaxi(ic,jc,kc);
      Cs_Eta(ic,jc,kc)  = Cs(ic,jc,kc) *  nablaeta(ic,jc,kc);
      Cs_Zeta(ic,jc,kc) = Cs(ic,jc,kc) * nablazeta(ic,jc,kc);  
        
        
      Sigma_Xi  =(sqrt(  BigU_Xi(ic,jc,kc)*  BigU_Xi(ic,jc,kc) &
                      &+ LaxFriSmall * LaxFriSmall * Cs_Xi(ic,jc,kc)   *  Cs_Xi(ic,jc,kc) ) +  Cs_Xi(ic,jc,kc) )/jaco(ic,jc,kc);

      Sigma_Eta =(sqrt( BigU_Eta(ic,jc,kc)* BigU_Eta(ic,jc,kc) &
                      &+ LaxFriSmall * LaxFriSmall * Cs_Eta(ic,jc,kc)  * Cs_Eta(ic,jc,kc) ) + Cs_Eta(ic,jc,kc) )/jaco(ic,jc,kc);

      Sigma_Zeta=(sqrt(BigU_Zeta(ic,jc,kc)*BigU_Zeta(ic,jc,kc) &
                      &+ LaxFriSmall * LaxFriSmall * Cs_Zeta(ic,jc,kc) *Cs_Zeta(ic,jc,kc) ) +Cs_Zeta(ic,jc,kc) )/jaco(ic,jc,kc);
      do iVar = 1,NumVar

        !invFhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *   dxidx(ic,jc,kc) +&
        !                         invG(iVar,ic,jc,kc) *   dxidy(ic,jc,kc) +&
        !                         invH(iVar,ic,jc,kc) *   dxidz(ic,jc,kc) +&
        !                       Ucons0(iVar,ic,jc,kc) *   dxidt(ic,jc,kc)) / Jaco(ic,jc,kc);
        !
        !invGhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) *  detadx(ic,jc,kc) +&
        !                         invG(iVar,ic,jc,kc) *  detady(ic,jc,kc) +&
        !                         invH(iVar,ic,jc,kc) *  detadz(ic,jc,kc) +&
        !                       Ucons0(iVar,ic,jc,kc) *  detadt(ic,jc,kc)) / Jaco(ic,jc,kc);
        !
        !invHhat(iVar,ic,jc,kc) =(invF(iVar,ic,jc,kc) * dzetadx(ic,jc,kc) +&
        !                         invG(iVar,ic,jc,kc) * dzetady(ic,jc,kc) +&
        !                         invH(iVar,ic,jc,kc) * dzetadz(ic,jc,kc) +&
        !                       Ucons0(iVar,ic,jc,kc) * dzetadt(ic,jc,kc)) / Jaco(ic,jc,kc);          
          
        Fp(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) + Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F+ in the xi direction
        Fm(iVar,ic,jc,kc) = 0.5_rk * (invFhat(iVar,ic,jc,kc) - Sigma_Xi   * Ucons0(iVar,ic,jc,kc)); ! F- in the xi direction

        Gp(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) + Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G+ in the eta direction
        Gm(iVar,ic,jc,kc) = 0.5_rk * (invGhat(iVar,ic,jc,kc) - Sigma_Eta  * Ucons0(iVar,ic,jc,kc)); ! G- in the eta direction

        Hp(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) + Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H+ in the zeta direction
        Hm(iVar,ic,jc,kc) = 0.5_rk * (invHhat(iVar,ic,jc,kc) - Sigma_Zeta * Ucons0(iVar,ic,jc,kc)); ! H- in the zeta direction
   
      enddo
      
    enddo
   enddo
  enddo
  
    call Parallel_Exchange_NumVar(Fp)
    call Parallel_Exchange_NumVar(Fm)
    call Parallel_Exchange_NumVar(Gp)
    call Parallel_Exchange_NumVar(Gm)
    call Parallel_Exchange_NumVar(Hp)
    call Parallel_Exchange_NumVar(Hm)   

    
    Select case(Select_Scheme)
    case(fifth_order)
  
      call Cal_Deri_Dxi_5th_up_numvar(Fpdxi,Fp)
      call Cal_Deri_Dxi_5th_do_numvar(Fmdxi,Fm)  
      
      call Cal_Deri_Deta_5th_up_numvar(Gpdeta,Gp)
      call Cal_Deri_Deta_5th_do_numvar(Gmdeta,Gm)
  
      call Cal_Deri_Dzeta_5th_up_per_numvar(Hpdzeta,Hp)
      call Cal_Deri_Dzeta_5th_do_per_numvar(Hmdzeta,Hm)

     case(first_order)
      call Cal_Deri_Dxi_1st_up_numvar(Fpdxi,Fp)
      call Cal_Deri_Dxi_1st_do_numvar(Fmdxi,Fm) 
      
      call Cal_Deri_Deta_1st_up_numvar(Gpdeta,Gp)
      call Cal_Deri_Deta_1st_do_numvar(Gmdeta,Gm)

      call Cal_Deri_Dzeta_1st_up_per_numvar(Hpdzeta,Hp)
      call Cal_Deri_Dzeta_1st_do_per_numvar(Hmdzeta,Hm)
     end Select
    
   
     call Parallel_Exchange_NumVar(Fpdxi)
     call Parallel_Exchange_NumVar(Fmdxi)
     call Parallel_Exchange_NumVar(Gpdeta)
     call Parallel_Exchange_NumVar(Gmdeta)
     call Parallel_Exchange_NumVar(Hpdzeta)
     call Parallel_Exchange_NumVar(Hmdzeta)
     
  invFlux = Fpdxi + Fmdxi + Gpdeta + Gmdeta + Hpdzeta + Hmdzeta;    
   
  !write(*,*),'Fpdxi(:,nx_local,1,1)+Fmdxi(:,nx_local,1,1)',invFlux(:,nx_local,Ny-1,1)
  
  call SOLVE_LD_outflow

  
  if(npx == npx0 - 1)then
    ic = nx_local 
     do kc = 1, nz_local
        do jc = 2, Ny-1
            do iVar = 1, NumVar
             invFlux(iVar,ic,jc,kc) = invFhatdxi_new(iVar,ic,jc,kc) + Gpdeta(iVar,ic,jc,kc) + Gmdeta(iVar,ic,jc,kc) + Hpdzeta(iVar,ic,jc,kc) + Hmdzeta(iVar,ic,jc,kc);
            enddo
        enddo
     enddo
  endif

  !if (npx == npx0 - 1) then
  !    write(*,*),'Fpdxi(:,nx_local,1,1)+Fmdxi(:,nx_local,1,1)',Fpdxi(1,nx_local,:,1)+Fmdxi(1,nx_local,:,1)
  !    write(*,*),'Fpdxi(:,nx_local,1,1)+Fmdxi(:,nx_local,1,1)',invFlux(:,nx_local,Ny-1,1)
  !endif
  
  !if (ModelType == 2) then
  !  call Parallel_Exchange_singularity_numvar_negative(invFlux)
  !endif   
  
    call Parallel_Exchange_NumVar(invFlux)
 
do kc = 1-overLAP,nz_local+overLAP
   do jc = 1,Ny
    do ic = 1-overLAP,nx_local+overLAP
!do kc = 1,nz_local
!   do jc = 1,Ny
!    do ic = 1,nx_local       
      do iVar = 1,NumVar
        dUcons0(iVar,ic,jc,kc) = ( - invFlux(iVar,ic,jc,kc) &
                                 & + visFlux(iVar,ic,jc,kc) &
                                 & -  Ucons0(iVar,ic,jc,kc) * invJacodt(ic,jc,kc)) * Jaco(ic,jc,kc);
      enddo
      
    enddo
   enddo
enddo

    
if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then
    call Parallel_Exchange_singularity_numvar(dUcons0,dUcons0_total)
endif

    call Parallel_Exchange_NumVar(dUcons0)

END SUBROUTINE Calculate_Flux
    
subroutine SOLVE_LD_outflow
  use SF_Constant,           only: ik,rk,NumVar,LaxFriSmall,overLAP,fifth_order,first_order
  use SF_CFD_Global,         only: nx_local,Ny,nz_local,Gamma,Pr_Ref,Mach_Ref,Re_Ref,&
                                   Csthlnd_Ref,Cp_Ref,Select_Scheme,ModelType,P_inf
  use SF_CFD_Global,         only: jaco,invF,invG,invH,dxidx,dxidy,dxidz,dxidt,Ucons0,Fpdxi,Fmdxi,Rho,U,V,W,T,P,Cs,&
                                   invFhatdxi_new
  use MPI_GLOBAL,            only: Parallel_Exchange,Parallel_Exchange_NumVar,npx,npx0

 implicit none
 integer( kind = ik ) :: ic,jc,kc,iVar,m,l
 real( kind = rk ):: P_Jinv(NumVar,NumVar),P_J(NumVar,NumVar)
 real( kind = rk ):: l_xi1(nx_local,Ny,nz_local),l_xi2(nx_local,Ny,nz_local),l_xi3(nx_local,Ny,nz_local),&
                     v_cdot_I(nx_local,Ny,nz_local),v_cross_Ix(nx_local,Ny,nz_local),v_cross_Iy(nx_local,Ny,nz_local),&
                     v_cross_Iz(nx_local,Ny,nz_local),l_xi_abs(nx_local,Ny,nz_local),&
                     dxidx_J(nx_local,Ny,nz_local),dxidy_J(nx_local,Ny,nz_local),dxidz_J(nx_local,Ny,nz_local),&
                     dxidt_J(nx_local,Ny,nz_local),dxidx_Jdxi(nx_local,Ny,nz_local),dxidy_Jdxi(nx_local,Ny,nz_local),&
                     dxidz_Jdxi(nx_local,Ny,nz_local),dxidt_Jdxi(nx_local,Ny,nz_local),LODI_L(NumVar,nx_local,Ny,nz_local),&
                     invfluxdxi(NumVar,nx_local,Ny,nz_local),invFhatdxi_old(NumVar,nx_local,Ny,nz_local),&
                     Mach2_local(nx_local,Ny,nz_local),H_LODI(nx_local,Ny,nz_local),V_temp(NumVar),Mach_local(nx_local,Ny,nz_local)
 real( kind = rk ):: temp,sigma_out, maxMach2
 
  P_Jinv = 0.0_rk
  P_J = 0.0_rk
  invFluxdxi = 0.0_rk
  sigma_out = 0.25_rk
  H_LODI = 0.0_rk
  LODI_L = 0.0_rk   
  Mach2_local = 0.0_rk
  invFhatdxi_old = 0.0_rk
  V_temp = 0.0_rk
  Mach_local = 0.0_rk

  do kc = 1,nz_local
   do jc = 1,Ny
    do ic = 1,nx_local
      do iVar = 1,NumVar
        invFhatdxi_old(iVar,ic,jc,kc) = Fpdxi(iVar,ic,jc,kc) + Fmdxi(iVar,ic,jc,kc)
      enddo
    enddo
   enddo
  enddo
  
  if(npx == npx0 - 1) then
   do kc = 1, nz_local
    do jc = 1, Ny
     do ic = 1, nx_local
        dxidx_J(ic,jc,kc) = dxidx(ic,jc,kc) / Jaco(ic,jc,kc)
        dxidy_J(ic,jc,kc) = dxidy(ic,jc,kc) / Jaco(ic,jc,kc)
        dxidz_J(ic,jc,kc) = dxidz(ic,jc,kc) / Jaco(ic,jc,kc)
        dxidt_J(ic,jc,kc) = dxidt(ic,jc,kc) / Jaco(ic,jc,kc)
        
     enddo
    enddo
   enddo
  endif 

  
  if(npx == npx0 - 1) then   
  ic = nx_local 
  do kc = 1, nz_local
   do jc = 1, Ny

      l_xi_abs(ic,jc,kc) = 1.0_rk/sqrt(dxidx(ic,jc,kc)*dxidx(ic,jc,kc) + &
                                      &dxidy(ic,jc,kc)*dxidy(ic,jc,kc) + &
                                      &dxidz(ic,jc,kc)*dxidz(ic,jc,kc))
      
      l_xi1(ic,jc,kc) = l_xi_abs(ic,jc,kc)*dxidx(ic,jc,kc)
      l_xi2(ic,jc,kc) = l_xi_abs(ic,jc,kc)*dxidy(ic,jc,kc)
      l_xi3(ic,jc,kc) = l_xi_abs(ic,jc,kc)*dxidz(ic,jc,kc)
      
      v_cdot_I(ic,jc,kc)   = U(ic,jc,kc)*l_xi1(ic,jc,kc) + &
                            &V(ic,jc,kc)*l_xi2(ic,jc,kc) + &
                            &W(ic,jc,kc)*l_xi3(ic,jc,kc) + dxidt(ic,jc,kc)
       
      Mach_local(ic,jc,kc) = v_cdot_I(ic,jc,kc) / Cs(ic,jc,kc)
      
      Mach2_local(ic,jc,kc) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2) / (Cs(ic,jc,kc)**2)

      enddo
  enddo
  ! 
  !  !maxMach2 = maxval(Mach2_local)
  
  ic = nx_local 
  do kc = 1, nz_local
   do jc = 1, Ny
      
     if (Mach_local(ic,jc,kc) < 1.0_rk )then
      
      v_cross_Ix(ic,jc,kc) = V(ic,jc,kc)*l_xi3(ic,jc,kc)-W(ic,jc,kc)*l_xi2(ic,jc,kc)
      v_cross_Iy(ic,jc,kc) = W(ic,jc,kc)*l_xi1(ic,jc,kc)-U(ic,jc,kc)*l_xi3(ic,jc,kc)
      v_cross_Iz(ic,jc,kc) = U(ic,jc,kc)*l_xi2(ic,jc,kc)-V(ic,jc,kc)*l_xi1(ic,jc,kc)
      
      
      
      dxidx_Jdxi(ic,jc,kc) = 25.d0/12.d0*dxidx_J(ic,jc,kc)-&
                            &       4.d0*dxidx_J(ic-1,jc,kc)+&
                            &       3.d0*dxidx_J(ic-2,jc,kc)-&
                            &  4.d0/3.d0*dxidx_J(ic-3,jc,kc)+&
                            &     0.25d0*dxidx_J(ic-4,jc,kc)
      dxidy_Jdxi(ic,jc,kc) = 25.d0/12.d0*dxidy_J(ic,jc,kc)-&
                            &       4.d0*dxidy_J(ic-1,jc,kc)+&
                            &       3.d0*dxidy_J(ic-2,jc,kc)-&
                            &  4.d0/3.d0*dxidy_J(ic-3,jc,kc)+&
                            &     0.25d0*dxidy_J(ic-4,jc,kc)
      dxidz_Jdxi(ic,jc,kc) = 25.d0/12.d0*dxidz_J(ic,jc,kc)-&
                            &       4.d0*dxidz_J(ic-1,jc,kc)+&
                            &       3.d0*dxidz_J(ic-2,jc,kc)-&
                            &  4.d0/3.d0*dxidz_J(ic-3,jc,kc)+&
                            &     0.25d0*dxidz_J(ic-4,jc,kc)
      dxidt_Jdxi(ic,jc,kc) = 25.d0/12.d0*dxidt_J(ic,jc,kc)-&
                            &       4.d0*dxidt_J(ic-1,jc,kc)+&
                            &       3.d0*dxidt_J(ic-2,jc,kc)-&
                            &  4.d0/3.d0*dxidt_J(ic-3,jc,kc)+&
                            &0.25d0     *dxidt_J(ic-4,jc,kc)
       
       
       P_Jinv = 0.0_rk
       P_J = 0.0_rk
      
       P_Jinv(1,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_xi1(ic,jc,kc) - 1/Rho(ic,jc,kc) * v_cross_Ix(ic,jc,kc)
       P_Jinv(2,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_xi2(ic,jc,kc) - 1/Rho(ic,jc,kc) * v_cross_Iy(ic,jc,kc)
       P_Jinv(3,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_xi3(ic,jc,kc) - 1/Rho(ic,jc,kc) * v_cross_Iz(ic,jc,kc)
       P_Jinv(4,1) = Cs(ic,jc,kc)/Rho(ic,jc,kc)*((gamma-1.d0)/2.d0*Mach2_local(ic,jc,kc) - v_cdot_I(ic,jc,kc)/Cs(ic,jc,kc))
       P_Jinv(5,1) = Cs(ic,jc,kc)/Rho(ic,jc,kc)*((gamma-1.d0)/2.d0*Mach2_local(ic,jc,kc) + v_cdot_I(ic,jc,kc)/Cs(ic,jc,kc))
       
       P_Jinv(1,2) = (gamma -1.d0)*U(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi1(ic,jc,kc)
       P_Jinv(2,2) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi2(ic,jc,kc) - l_xi3(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(3,2) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi3(ic,jc,kc) + l_xi2(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(4,2) = l_xi1(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*U(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,2) = -l_xi1(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*U(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi1(ic,jc,kc) + l_xi3(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(2,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi2(ic,jc,kc)
       P_Jinv(3,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi3(ic,jc,kc) - l_xi1(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(4,3) = l_xi2(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*V(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,3) = -l_xi2(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*V(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi1(ic,jc,kc) - l_xi2(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(2,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi2(ic,jc,kc) + l_xi1(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(3,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi3(ic,jc,kc)
       P_Jinv(4,4) = l_xi3(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*W(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,4) = -l_xi3(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*W(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi1(ic,jc,kc)
       P_Jinv(2,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi2(ic,jc,kc)
       P_Jinv(3,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_xi3(ic,jc,kc)
       P_Jinv(4,5) = (gamma - 1.d0)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,5) = (gamma - 1.d0)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       do l = 1, NumVar
           invfluxdxi(l,ic,jc,kc) = invF(l,ic,jc,kc)*dxidx_Jdxi(ic,jc,kc) + &
                                    invG(l,ic,jc,kc)*dxidy_Jdxi(ic,jc,kc) + &
                                    invH(l,ic,jc,kc)*dxidz_Jdxi(ic,jc,kc) + &
                                    Ucons0(l,ic,jc,kc)*dxidt_Jdxi(ic,jc,kc); 
           
           V_temp(l) = invFhatdxi_old(l,ic,jc,kc) - invfluxdxi(l,ic,jc,kc)
       enddo
       
           
          do m = 1, 4
               LODI_L(m,ic,jc,kc) = 0.0_rk
            do l = 1, NumVar 
               LODI_L(m,ic,jc,kc) = LODI_L(m,ic,jc,kc) + P_Jinv(m,l) * Jaco(ic,jc,kc) * V_temp(l)
            enddo
          enddo
          
               LODI_L(5,ic,jc,kc) = sigma_out*(1.d0 - Mach_local(ic,jc,kc)**2)*(P(ic,jc,kc)-P_inf*(Gamma * Mach_Ref * Mach_Ref))/Rho(ic,jc,kc)

         
               
               
               
       H_LODI(ic,jc,kc) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0 + Cs(ic,jc,kc)**2 / (gamma - 1.d0)   
               
       P_J(1,1) = l_xi1(ic,jc,kc)
       P_J(2,1) = U(ic,jc,kc)*l_xi1(ic,jc,kc)
       P_J(3,1) = V(ic,jc,kc)*l_xi1(ic,jc,kc) + Rho(ic,jc,kc)*l_xi3(ic,jc,kc)
       P_J(4,1) = W(ic,jc,kc)*l_xi1(ic,jc,kc) - Rho(ic,jc,kc)*l_xi2(ic,jc,kc)
       P_J(5,1) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_xi1(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Ix(ic,jc,kc)
       
       P_J(1,2) = l_xi2(ic,jc,kc)
       P_J(2,2) = U(ic,jc,kc)*l_xi2(ic,jc,kc) - Rho(ic,jc,kc)*l_xi3(ic,jc,kc)
       P_J(3,2) = V(ic,jc,kc)*l_xi2(ic,jc,kc)
       P_J(4,2) = W(ic,jc,kc)*l_xi2(ic,jc,kc) + Rho(ic,jc,kc)*l_xi1(ic,jc,kc)
       P_J(5,2) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_xi2(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Iy(ic,jc,kc)
       
       P_J(1,3) = l_xi3(ic,jc,kc) 
       P_J(2,3) = U(ic,jc,kc)*l_xi3(ic,jc,kc) + Rho(ic,jc,kc)*l_xi2(ic,jc,kc)
       P_J(3,3) = V(ic,jc,kc)*l_xi3(ic,jc,kc) - Rho(ic,jc,kc)*l_xi1(ic,jc,kc)
       P_J(4,3) = W(ic,jc,kc)*l_xi3(ic,jc,kc)
       P_J(5,3) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_xi3(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Iz(ic,jc,kc)
              
       P_J(1,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)
       P_J(2,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(U(ic,jc,kc) + l_xi1(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(3,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(V(ic,jc,kc) + l_xi2(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(4,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(W(ic,jc,kc) + l_xi3(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(5,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(H_LODI(ic,jc,kc) + Cs(ic,jc,kc) * v_cdot_I(ic,jc,kc))
       
       P_J(1,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)
       P_J(2,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(U(ic,jc,kc) - l_xi1(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(3,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(V(ic,jc,kc) - l_xi2(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(4,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(W(ic,jc,kc) - l_xi3(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(5,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(H_LODI(ic,jc,kc) - Cs(ic,jc,kc) * v_cdot_I(ic,jc,kc))
       
      
      do l = 1, NumVar
             temp = 0.0_rk
             do m = 1, NumVar
                 temp = temp + P_J(l,m)*LODI_L(m,ic,jc,kc)
             enddo
             
             invFhatdxi_new(l,ic,jc,kc) =  temp/Jaco(ic,jc,kc) + invFluxdxi(l,ic,jc,kc)
      enddo
    
      else
           do l = 1, NumVar
             invFhatdxi_new(l,ic,jc,kc) = invFhatdxi_old(l,ic,jc,kc)
           enddo
      endif  
      
   enddo
  enddo  
  
  endif
  
end subroutine SOLVE_LD_outflow  
    
    
subroutine SOLVE_LD_wall
  use SF_Constant,           only: ik,rk,NumVar,LaxFriSmall,overLAP,fifth_order,first_order
  use SF_CFD_Global,         only: nx_local,Ny,nz_local,Gamma,Pr_Ref,Mach_Ref,Re_Ref,&
                                   Csthlnd_Ref,Cp_Ref,Select_Scheme,ModelType,P_inf
  use SF_CFD_Global,         only: jaco,invF,invG,invH,detadx,detady,detadz,detadt,Ucons0,Gpdeta,Gmdeta,Rho,U,V,W,T,P,Cs,&
                                   invGhatdeta_new
  use MPI_GLOBAL,            only: Parallel_Exchange,Parallel_Exchange_NumVar,npx,npx0

 implicit none
 integer( kind = ik ) :: ic,jc,kc,iVar,m,l
 real( kind = rk ):: P_Jinv(NumVar,NumVar),P_J(NumVar,NumVar)
 real( kind = rk ):: l_eta1(nx_local,Ny,nz_local),l_eta2(nx_local,Ny,nz_local),l_eta3(nx_local,Ny,nz_local),&
                     v_cdot_I(nx_local,Ny,nz_local),v_cross_Ix(nx_local,Ny,nz_local),v_cross_Iy(nx_local,Ny,nz_local),&
                     v_cross_Iz(nx_local,Ny,nz_local),l_eta_abs(nx_local,Ny,nz_local),&
                     detadx_J(nx_local,Ny,nz_local),detady_J(nx_local,Ny,nz_local),detadz_J(nx_local,Ny,nz_local),&
                     detadt_J(nx_local,Ny,nz_local),detadx_Jdeta(nx_local,Ny,nz_local),detady_Jdeta(nx_local,Ny,nz_local),&
                     detadz_Jdeta(nx_local,Ny,nz_local),detadt_Jdeta(nx_local,Ny,nz_local),LODI_L(NumVar,nx_local,Ny,nz_local),&
                     invfluxdeta(NumVar,nx_local,Ny,nz_local),invFhatdeta_old(NumVar,nx_local,Ny,nz_local),&
                     Mach2_local(nx_local,Ny,nz_local),H_LODI(nx_local,Ny,nz_local),V_temp(NumVar),Mach_local(nx_local,Ny,nz_local)
 real( kind = rk ):: temp,sigma_out, maxMach2
 
  P_Jinv = 0.0_rk
  P_J = 0.0_rk
  invFluxdeta = 0.0_rk
  sigma_out = 0.25_rk
  H_LODI = 0.0_rk
  LODI_L = 0.0_rk   
  Mach2_local = 0.0_rk
  invFhatdeta_old = 0.0_rk
  V_temp = 0.0_rk
  Mach_local = 0.0_rk
  
  if(npx == npx0 - 1) then
   do kc = 1, nz_local
    do jc = 1, Ny
     do ic = 1, nx_local
        detadx_J(ic,jc,kc) = detadx(ic,jc,kc) / Jaco(ic,jc,kc)
        detady_J(ic,jc,kc) = detady(ic,jc,kc) / Jaco(ic,jc,kc)
        detadz_J(ic,jc,kc) = detadz(ic,jc,kc) / Jaco(ic,jc,kc)
        detadt_J(ic,jc,kc) = detadt(ic,jc,kc) / Jaco(ic,jc,kc)
        
     enddo
    enddo
   enddo
  endif 

  
  if(npx == npx0 - 1) then   
  jc = 1 
  do kc = 1, Nz_local
   do ic = 1, Nx_local

      l_eta_abs(ic,jc,kc) = 1.0_rk/sqrt(detadx(ic,jc,kc)*detadx(ic,jc,kc) + &
                                       &detady(ic,jc,kc)*detady(ic,jc,kc) + &
                                       &detadz(ic,jc,kc)*detadz(ic,jc,kc))
      
      l_eta1(ic,jc,kc) = l_eta_abs(ic,jc,kc)*detadx(ic,jc,kc)
      l_eta2(ic,jc,kc) = l_eta_abs(ic,jc,kc)*detady(ic,jc,kc)
      l_eta3(ic,jc,kc) = l_eta_abs(ic,jc,kc)*detadz(ic,jc,kc)
      
      v_cdot_I(ic,jc,kc)   = U(ic,jc,kc)*l_eta1(ic,jc,kc) + &
                            &V(ic,jc,kc)*l_eta2(ic,jc,kc) + &
                            &W(ic,jc,kc)*l_eta3(ic,jc,kc) + detadt(ic,jc,kc)
       
      Mach_local(ic,jc,kc) = v_cdot_I(ic,jc,kc) / Cs(ic,jc,kc)
      
      Mach2_local(ic,jc,kc) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2) / (Cs(ic,jc,kc)**2)

      enddo
  enddo
  ! 
  !  !maxMach2 = maxval(Mach2_local)
  
  jc = 1 
  do kc = 1, Nz_local
   do ic = 1, Nx_local
      
      v_cross_Ix(ic,jc,kc) = V(ic,jc,kc)*l_eta3(ic,jc,kc)-W(ic,jc,kc)*l_eta2(ic,jc,kc)
      v_cross_Iy(ic,jc,kc) = W(ic,jc,kc)*l_eta1(ic,jc,kc)-U(ic,jc,kc)*l_eta3(ic,jc,kc)
      v_cross_Iz(ic,jc,kc) = U(ic,jc,kc)*l_eta2(ic,jc,kc)-V(ic,jc,kc)*l_eta1(ic,jc,kc)
      
      
      
      detadx_Jdeta(ic,jc,kc) = -25.d0/12.d0*detadx_J(ic,jc,  kc)+&
                             &         4.d0*detadx_J(ic,jc+1,kc)-&
                             &         3.d0*detadx_J(ic,jc+2,kc)+&
                             &    4.d0/3.d0*detadx_J(ic,jc+3,kc)-&
                             &       0.25d0*detadx_J(ic,jc+4,kc)
      
      detady_Jdeta(ic,jc,kc) = -25.d0/12.d0*detady_J(ic,jc,  kc)+&
                             &         4.d0*detady_J(ic,jc+1,kc)-&
                             &         3.d0*detady_J(ic,jc+2,kc)+&
                             &    4.d0/3.d0*detady_J(ic,jc+3,kc)-&
                             &       0.25d0*detady_J(ic,jc+4,kc)
      
      detadz_Jdeta(ic,jc,kc) = -25.d0/12.d0*detadz_J(ic,jc,  kc)+&
                             &         4.d0*detadz_J(ic,jc+1,kc)-&
                             &         3.d0*detadz_J(ic,jc+2,kc)+&
                             &    4.d0/3.d0*detadz_J(ic,jc+3,kc)-&
                             &       0.25d0*detadz_J(ic,jc+4,kc)
      
      detadt_Jdeta(ic,jc,kc) = -25.d0/12.d0*detadt_J(ic,jc,  kc)+&
                             &         4.d0*detadt_J(ic,jc+1,kc)-&
                             &         3.d0*detadt_J(ic,jc+2,kc)+&
                             &    4.d0/3.d0*detadt_J(ic,jc+3,kc)-&
                             &       0.25d0*detadt_J(ic,jc+4,kc)
       
       
       P_Jinv = 0.0_rk
       P_J = 0.0_rk
      
       P_Jinv(1,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_eta1(ic,jc,kc) - v_cross_Ix(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(2,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_eta2(ic,jc,kc) - v_cross_Iy(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(3,1) = (1-(gamma - 1)/2.d0 * Mach2_local(ic,jc,kc))*l_eta3(ic,jc,kc) - v_cross_Iz(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(4,1) = Cs(ic,jc,kc)/Rho(ic,jc,kc)*((gamma-1.d0)/2.d0*Mach2_local(ic,jc,kc) - v_cdot_I(ic,jc,kc)/Cs(ic,jc,kc))
       P_Jinv(5,1) = Cs(ic,jc,kc)/Rho(ic,jc,kc)*((gamma-1.d0)/2.d0*Mach2_local(ic,jc,kc) + v_cdot_I(ic,jc,kc)/Cs(ic,jc,kc))
       
       P_Jinv(1,2) = (gamma -1.d0)*U(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta1(ic,jc,kc)
       P_Jinv(2,2) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta2(ic,jc,kc) - l_eta3(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(3,2) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta3(ic,jc,kc) + l_eta2(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(4,2) = l_eta1(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*U(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,2) = -l_eta1(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*U(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta1(ic,jc,kc) + l_eta3(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(2,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta2(ic,jc,kc)
       P_Jinv(3,3) = (gamma -1.d0)*V(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta3(ic,jc,kc) - l_eta1(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(4,3) = l_eta2(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*V(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,3) = -l_eta2(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*V(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta1(ic,jc,kc) - l_eta2(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(2,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta2(ic,jc,kc) + l_eta1(ic,jc,kc)/Rho(ic,jc,kc)
       P_Jinv(3,4) = (gamma -1.d0)*W(ic,jc,kc)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta3(ic,jc,kc)
       P_Jinv(4,4) = l_eta3(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*W(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,4) = -l_eta3(ic,jc,kc)/Rho(ic,jc,kc) - (gamma - 1.d0)*W(ic,jc,kc)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       P_Jinv(1,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta1(ic,jc,kc)
       P_Jinv(2,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta2(ic,jc,kc)
       P_Jinv(3,5) = -(gamma -1.d0)/Cs(ic,jc,kc)/Cs(ic,jc,kc)*l_eta3(ic,jc,kc)
       P_Jinv(4,5) = (gamma - 1.d0)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       P_Jinv(5,5) = (gamma - 1.d0)/Rho(ic,jc,kc)/Cs(ic,jc,kc)
       
       do l = 1, NumVar
           invfluxdeta(l,ic,jc,kc) = invF(l,ic,jc,kc)*detadx_Jdeta(ic,jc,kc) + &
                                     invG(l,ic,jc,kc)*detady_Jdeta(ic,jc,kc) + &
                                     invH(l,ic,jc,kc)*detadz_Jdeta(ic,jc,kc) + &
                                   Ucons0(l,ic,jc,kc)*detadt_Jdeta(ic,jc,kc); 
           
           V_temp(l) = invFhatdeta_old(l,ic,jc,kc) - invfluxdeta(l,ic,jc,kc)
       enddo
       
           
          do m = 1, 3
            do l = 1, NumVar 
               LODI_L(m,ic,jc,kc) = 0.0_rk
            enddo
          enddo
          
            m = 4
               LODI_L(m,ic,jc,kc) = 0.0_rk
            do l = 1, NumVar 
               LODI_L(m,ic,jc,kc) = LODI_L(m,ic,jc,kc) + P_Jinv(m,l) * Jaco(ic,jc,kc) * V_temp(l)
            enddo
                   
          
               LODI_L(5,ic,jc,kc) = LODI_L(4,ic,jc,kc)
           
               
       H_LODI(ic,jc,kc) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0 + Cs(ic,jc,kc)**2 / (gamma - 1.d0)   
               
       P_J(1,1) = l_eta1(ic,jc,kc)
       P_J(2,1) = U(ic,jc,kc)*l_eta1(ic,jc,kc)
       P_J(3,1) = V(ic,jc,kc)*l_eta1(ic,jc,kc) + Rho(ic,jc,kc)*l_eta3(ic,jc,kc)
       P_J(4,1) = W(ic,jc,kc)*l_eta1(ic,jc,kc) - Rho(ic,jc,kc)*l_eta2(ic,jc,kc)
       P_J(5,1) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_eta1(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Ix(ic,jc,kc)
       
       P_J(1,2) = l_eta2(ic,jc,kc)
       P_J(2,2) = U(ic,jc,kc)*l_eta2(ic,jc,kc) - Rho(ic,jc,kc)*l_eta3(ic,jc,kc)
       P_J(3,2) = V(ic,jc,kc)*l_eta2(ic,jc,kc)
       P_J(4,2) = W(ic,jc,kc)*l_eta2(ic,jc,kc) + Rho(ic,jc,kc)*l_eta1(ic,jc,kc)
       P_J(5,2) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_eta2(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Iy(ic,jc,kc)
       
       P_J(1,3) = l_eta3(ic,jc,kc) 
       P_J(2,3) = U(ic,jc,kc)*l_eta3(ic,jc,kc) + Rho(ic,jc,kc)*l_eta2(ic,jc,kc)
       P_J(3,3) = V(ic,jc,kc)*l_eta3(ic,jc,kc) - Rho(ic,jc,kc)*l_eta1(ic,jc,kc)
       P_J(4,3) = W(ic,jc,kc)*l_eta3(ic,jc,kc)
       P_J(5,3) = (U(ic,jc,kc)**2 + V(ic,jc,kc)**2 + W(ic,jc,kc)**2)/2.d0*l_eta3(ic,jc,kc) + Rho(ic,jc,kc)*v_cross_Iz(ic,jc,kc)
              
       P_J(1,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)
       P_J(2,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(U(ic,jc,kc) + l_eta1(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(3,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(V(ic,jc,kc) + l_eta2(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(4,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(W(ic,jc,kc) + l_eta3(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(5,4) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(H_LODI(ic,jc,kc) + Cs(ic,jc,kc) * v_cdot_I(ic,jc,kc))
       
       P_J(1,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)
       P_J(2,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(U(ic,jc,kc) - l_eta1(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(3,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(V(ic,jc,kc) - l_eta2(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(4,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(W(ic,jc,kc) - l_eta3(ic,jc,kc) * Cs(ic,jc,kc))
       P_J(5,5) = Rho(ic,jc,kc)/2.d0/Cs(ic,jc,kc)*(H_LODI(ic,jc,kc) - Cs(ic,jc,kc) * v_cdot_I(ic,jc,kc))
       
      
      do l = 1, NumVar
             temp = 0.0_rk
             do m = 1, NumVar
                 temp = temp + P_J(l,m)*LODI_L(m,ic,jc,kc)
             enddo
             
             invGhatdeta_new(l,ic,jc,kc) =  temp/Jaco(ic,jc,kc) + invFluxdeta(l,ic,jc,kc)
      enddo

      
   enddo
  enddo  
  
  endif
  
end subroutine SOLVE_LD_wall   