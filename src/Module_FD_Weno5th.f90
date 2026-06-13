!MODULE FD5_Order_Weno
!  ! This is the module that contains all the subs and dates
!  ! used for 5th order upwind/downwind and 6th order finite
!  ! difference methods
!  use SF_Constant,  only: ik, rk, NumVar, overLap
!  use SF_CFD_Global,only: nx_local,Ny,nz_local
!  use MPI_Global,   only: if_parallel, MyId, NumProcess, MPI_COMM_WORLD, ierr,&
!                        & npx,npx0,npz,npz0
!  use mpi
!  implicit none
!  ! in order to improve the calculation speed, all the coefficients for finite
!  ! difference calculation are used as constants
!  
!  ! Boundary: the same as Module FD5_Order
!  real( kind = rk ), parameter:: BC1(1:5) =(/ -25.d0/12.d0,      4.d0,  -3.d0, 4.d0/3.d0,    -0.25d0 /) 
!  real( kind = rk ), parameter:: BC2(1:5) =(/      -0.25d0,-5.d0/6.d0,  1.5d0,    -0.5d0, 1.d0/12.d0 /) 
!  real( kind = rk ), parameter:: BC3(1:5) =(/   1.d0/12.d0,-2.d0/3.d0,   0.d0, 2.d0/3.d0,-1.d0/12.d0 /)
!  real( kind = rk ), parameter:: BCchar1(1:3) = (/ 2.d0/6.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!  real( kind = rk ), parameter:: eps = 1.d-6
!
!  real( kind = rk ), parameter:: CO_FD_WENO_H1(-1:1) = (/ 1.d0/3.d0,  -7.d0/6.d0, 11.d0/6.d0  /)
!  real( kind = rk ), parameter:: CO_FD_WENO_H2(-1:1) = (/ -1.d0/6.d0, 5.d0/6.d0,  1.d0/3.d0   /)
!  real( kind = rk ), parameter:: CO_FD_WENO_H3(-1:1) = (/ 1.d0/3.d0,  5.d0/6.d0,  -1.d0/6.d0  /)
!
!    contains
!
!! WENO f(i-1/2),f(i+1/2) calculation template
!function weno_halfvalue(varloc) result(f_half)
!  implicit none
!  real( kind = rk ), intent(in)::varloc(-2:2)
!  real( kind = rk ):: f_half
!  real( kind = rk ):: h(3),beta(3),w(3)
!
!  h(1) = CO_FD_WENO_H1(-1) * varloc(-2) &
!      &+ CO_FD_WENO_H1( 0) * varloc(-1) &
!      &+ CO_FD_WENO_H1( 1) * varloc( 0)
!  h(2) = CO_FD_WENO_H2(-1) * varloc(-1) &
!      &+ CO_FD_WENO_H2( 0) * varloc( 0) &
!      &+ CO_FD_WENO_H2( 1) * varloc( 1)
!  h(3) = CO_FD_WENO_H3(-1) * varloc( 0) &
!      &+ CO_FD_WENO_H3( 0) * varloc( 1) &
!      &+ CO_FD_WENO_H3( 1) * varloc( 2)
!
!  beta(1) = 1.d0 / 4.d0 * (       varloc(-2) - 4.d0 * varloc(-1) + 3.d0 * varloc( 0))**2 &
!      &+  13.d0 / 12.d0 * (       varloc(-2) - 2.d0 * varloc(-1) +        varloc( 0))**2
!  beta(2) = 1.d0 / 4.d0 * (       varloc(-1)                     -        varloc( 1))**2 &
!      &+  13.d0/12.d0   * (       varloc(-1) - 2.d0 * varloc( 0) +        varloc( 1))**2
!  beta(3) = 1.d0 / 4.d0 * (3.d0 * varloc( 0) - 4.d0 * varloc( 1) +        varloc( 2))**2 &
!      &+  13.d0/12.d0   * (       varloc( 0) - 2.d0 * varloc( 1) +        varloc( 2))**2
!
!  w(1) = 1.d0 / 10.d0 / ( (eps + beta(1))**2 )
!  w(2) = 6.d0 / 10.d0 / ( (eps + beta(2))**2 )
!  w(3) = 3.d0 / 10.d0 / ( (eps + beta(3))**2 )
!
!  f_half = dot_product(w, h) / sum(w)
!
!end function weno_halfvalue
!
!function weno_halfvalue_del_left(varloc) result(f_half)
!  implicit none
!  real( kind = rk ), intent(in)::varloc(-2:2)
!  real( kind = rk ):: f_half
!  real( kind = rk ):: h(3),beta(3),w(3)
!
!  h(1) = 0.d0
!  h(2) = CO_FD_WENO_H2(-1) * varloc(-1) &
!      &+ CO_FD_WENO_H2( 0) * varloc( 0) &
!      &+ CO_FD_WENO_H2( 1) * varloc( 1)
!  h(3) = CO_FD_WENO_H3(-1) * varloc( 0) &
!      &+ CO_FD_WENO_H3( 0) * varloc( 1) &
!      &+ CO_FD_WENO_H3( 1) * varloc( 2)
!
!  beta(2) = 1.d0 / 4.d0 * (       varloc(-1)                     -        varloc( 1))**2 &
!      &+  13.d0/12.d0   * (       varloc(-1) - 2.d0 * varloc( 0) +        varloc( 1))**2
!  beta(3) = 1.d0 / 4.d0 * (3.d0 * varloc( 0) - 4.d0 * varloc( 1) +        varloc( 2))**2 &
!      &+  13.d0/12.d0   * (       varloc( 0) - 2.d0 * varloc( 1) +        varloc( 2))**2
!
!  w(1) = 0.d0
!  w(2) = 6.d0 / 10.d0 / ( (eps + beta(2))**2 )
!  w(3) = 3.d0 / 10.d0 / ( (eps + beta(3))**2 )
!
!  f_half = dot_product(w, h) / sum(w)
!
!end function weno_halfvalue_del_left
!
!function weno_halfvalue_del_right(varloc) result(f_half)
!  implicit none
!  real( kind = rk ), intent(in)::varloc(-2:2)
!  real( kind = rk ):: f_half
!  real( kind = rk ):: h(3),beta(3),w(3)
!
!  h(1) = CO_FD_WENO_H1(-1) * varloc(-2) &
!      &+ CO_FD_WENO_H1( 0) * varloc(-1) &
!      &+ CO_FD_WENO_H1( 1) * varloc( 0)
!  h(2) = CO_FD_WENO_H2(-1) * varloc(-1) &
!      &+ CO_FD_WENO_H2( 0) * varloc( 0) &
!      &+ CO_FD_WENO_H2( 1) * varloc( 1)
!  h(3) = 0.d0
!
!  beta(1) = 1.d0 / 4.d0 * (       varloc(-2) - 4.d0 * varloc(-1) + 3.d0 * varloc( 0))**2 &
!      &+  13.d0 / 12.d0 * (       varloc(-2) - 2.d0 * varloc(-1) +        varloc( 0))**2
!  beta(2) = 1.d0 / 4.d0 * (       varloc(-1)                     -        varloc( 1))**2 &
!      &+  13.d0/12.d0   * (       varloc(-1) - 2.d0 * varloc( 0) +        varloc( 1))**2
!
!  w(1) = 1.d0 / 10.d0 / ( (eps + beta(1))**2 )
!  w(2) = 6.d0 / 10.d0 / ( (eps + beta(2))**2 )
!  w(3) = 0.d0
!
!  f_half = dot_product(w, h) / sum(w)
!
!end function weno_halfvalue_del_right
!
!! 3D Block Derivatives
!! Xi direction 5th upwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Dxi_weno5th_up_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        do iVar = 1, NumVar
!          varloc(-2:2) = f_loc(iVar, ic-2:ic+2, jc, kc)
!          fp = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic-3:ic+1, jc, kc)
!          fm = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npx == 0) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = 1
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC1(2) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC1(3) * f_loc(iVar,ic+2,jc,kc) &
!                           &+ BC1(4) * f_loc(iVar,ic+3,jc,kc) &
!                           &+ BC1(5) * f_loc(iVar,ic+4,jc,kc)
!        enddo
!                  
!        ic = 2
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic-1,jc,kc) &
!                           &+ BC2(2) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC2(3) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC2(4) * f_loc(iVar,ic+2,jc,kc) &
!                           &+ BC2(5) * f_loc(iVar,ic+3,jc,kc)
!        enddo
!    
!        ic = 3
!        do iVar = 1, NumVar                                
!          df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic-2,jc,kc) &
!                           &+ BC3(2) * f_loc(iVar,ic-1,jc,kc) &
!                           &+ BC3(3) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC3(4) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC3(5) * f_loc(iVar,ic+2,jc,kc)
!        enddo
!      enddo 
!    enddo
!  endif
!
!  if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = nx_local
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC1(2) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC1(3) * f_loc(iVar,ic-2,jc,kc) &
!                            &- BC1(4) * f_loc(iVar,ic-3,jc,kc) &
!                            &- BC1(5) * f_loc(iVar,ic-4,jc,kc)
!        enddo
!
!        ic = nx_local-1
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic+1,jc,kc) &
!                            &- BC2(2) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC2(3) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC2(4) * f_loc(iVar,ic-2,jc,kc) &
!                            &- BC2(5) * f_loc(iVar,ic-3,jc,kc)
!        enddo
!
!        ic = nx_local-2
!        do iVar = 1, NumVar                                
!          df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic+2,jc,kc) &
!                            &- BC3(2) * f_loc(iVar,ic+1,jc,kc) &
!                            &- BC3(3) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC3(4) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC3(5) * f_loc(iVar,ic-2,jc,kc)
!        enddo
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dxi_weno5th_up_numvar
!
!! Xi direction 5th downwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Dxi_weno5th_do_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        do iVar = 1, NumVar
!          varloc(-2:2) = f_loc(iVar, ic+2:ic-2:-1, jc, kc)
!          fm = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic+3:ic-1:-1, jc, kc)
!          fp = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npx == 0) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = 1
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC1(2) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC1(3) * f_loc(iVar,ic+2,jc,kc) &
!                           &+ BC1(4) * f_loc(iVar,ic+3,jc,kc) &
!                           &+ BC1(5) * f_loc(iVar,ic+4,jc,kc)
!        enddo
!                  
!        ic = 2
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic-1,jc,kc) &
!                           &+ BC2(2) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC2(3) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC2(4) * f_loc(iVar,ic+2,jc,kc) &
!                           &+ BC2(5) * f_loc(iVar,ic+3,jc,kc)
!        enddo
!    
!        ic = 3
!        do iVar = 1, NumVar                                
!          df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic-2,jc,kc) &
!                           &+ BC3(2) * f_loc(iVar,ic-1,jc,kc) &
!                           &+ BC3(3) * f_loc(iVar,ic  ,jc,kc) &
!                           &+ BC3(4) * f_loc(iVar,ic+1,jc,kc) &
!                           &+ BC3(5) * f_loc(iVar,ic+2,jc,kc)
!        enddo
!      enddo 
!    enddo
!  endif
!
!  if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = nx_local
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC1(2) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC1(3) * f_loc(iVar,ic-2,jc,kc) &
!                            &- BC1(4) * f_loc(iVar,ic-3,jc,kc) &
!                            &- BC1(5) * f_loc(iVar,ic-4,jc,kc)
!        enddo
!
!        ic = nx_local-1
!        do iVar = 1, NumVar
!          df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic+1,jc,kc) &
!                            &- BC2(2) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC2(3) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC2(4) * f_loc(iVar,ic-2,jc,kc) &
!                            &- BC2(5) * f_loc(iVar,ic-3,jc,kc)
!        enddo
!
!        ic = nx_local-2
!        do iVar = 1, NumVar                                
!          df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic+2,jc,kc) &
!                            &- BC3(2) * f_loc(iVar,ic+1,jc,kc) &
!                            &- BC3(3) * f_loc(iVar,ic  ,jc,kc) &
!                            &- BC3(4) * f_loc(iVar,ic-1,jc,kc) &
!                            &- BC3(5) * f_loc(iVar,ic-2,jc,kc)
!        enddo
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dxi_weno5th_do_numvar
!
!
!! Xi direction 5th upwind WENO normal derivatives f(:,:,:) 
!subroutine Cal_Deri_Dxi_weno5th_up(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic-2:ic+2, jc, kc)
!        fp = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic-3:ic+1, jc, kc)
!        fm = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npx == 0) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = 1
!        df(ic,jc,kc) = BC1(1) * f_loc(ic  ,jc,kc) &
!                    &+ BC1(2) * f_loc(ic+1,jc,kc) &
!                    &+ BC1(3) * f_loc(ic+2,jc,kc) &
!                    &+ BC1(4) * f_loc(ic+3,jc,kc) &
!                    &+ BC1(5) * f_loc(ic+4,jc,kc)
!                    
!        ic = 2
!        df(ic,jc,kc) = BC2(1) * f_loc(ic-1,jc,kc) &
!                    &+ BC2(2) * f_loc(ic  ,jc,kc) &
!                    &+ BC2(3) * f_loc(ic+1,jc,kc) &
!                    &+ BC2(4) * f_loc(ic+2,jc,kc) &
!                    &+ BC2(5) * f_loc(ic+3,jc,kc)
!      
!        ic = 3                                
!        df(ic,jc,kc) = BC3(1) * f_loc(ic-2,jc,kc) &
!                    &+ BC3(2) * f_loc(ic-1,jc,kc) &
!                    &+ BC3(3) * f_loc(ic  ,jc,kc) &
!                    &+ BC3(4) * f_loc(ic+1,jc,kc) &
!                    &+ BC3(5) * f_loc(ic+2,jc,kc)
!      enddo 
!    enddo
!  endif
!
!  if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = nx_local
!        df(ic,jc,kc) =- BC1(1) * f_loc(ic  ,jc,kc) &
!                    &- BC1(2) * f_loc(ic-1,jc,kc) &
!                    &- BC1(3) * f_loc(ic-2,jc,kc) &
!                    &- BC1(4) * f_loc(ic-3,jc,kc) &
!                    &- BC1(5) * f_loc(ic-4,jc,kc)
!
!        ic = nx_local-1
!        df(ic,jc,kc) =- BC2(1) * f_loc(ic+1,jc,kc) &
!                    &- BC2(2) * f_loc(ic  ,jc,kc) &
!                    &- BC2(3) * f_loc(ic-1,jc,kc) &
!                    &- BC2(4) * f_loc(ic-2,jc,kc) &
!                    &- BC2(5) * f_loc(ic-3,jc,kc)
!
!        ic = nx_local-2                                
!        df(ic,jc,kc) = -BC3(1) * f_loc(ic+2,jc,kc) &
!                    &- BC3(2) * f_loc(ic+1,jc,kc) &
!                    &- BC3(3) * f_loc(ic  ,jc,kc) &
!                    &- BC3(4) * f_loc(ic-1,jc,kc) &
!                    &- BC3(5) * f_loc(ic-2,jc,kc)
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dxi_weno5th_up
!
!
!! Xi direction 5th downwind WENO normal derivatives f(:,:,:) 
!subroutine Cal_Deri_Dxi_weno5th_do(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic+2:ic-2:-1, jc, kc)
!        fm = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic+3:ic-1:-1, jc, kc)
!        fp = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npx == 0) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = 1
!        df(ic,jc,kc) = BC1(1) * f_loc(ic  ,jc,kc) &
!                    &+ BC1(2) * f_loc(ic+1,jc,kc) &
!                    &+ BC1(3) * f_loc(ic+2,jc,kc) &
!                    &+ BC1(4) * f_loc(ic+3,jc,kc) &
!                    &+ BC1(5) * f_loc(ic+4,jc,kc)
!                    
!        ic = 2
!        df(ic,jc,kc) = BC2(1) * f_loc(ic-1,jc,kc) &
!                    &+ BC2(2) * f_loc(ic  ,jc,kc) &
!                    &+ BC2(3) * f_loc(ic+1,jc,kc) &
!                    &+ BC2(4) * f_loc(ic+2,jc,kc) &
!                    &+ BC2(5) * f_loc(ic+3,jc,kc)
!      
!        ic = 3                                
!        df(ic,jc,kc) = BC3(1) * f_loc(ic-2,jc,kc) &
!                    &+ BC3(2) * f_loc(ic-1,jc,kc) &
!                    &+ BC3(3) * f_loc(ic  ,jc,kc) &
!                    &+ BC3(4) * f_loc(ic+1,jc,kc) &
!                    &+ BC3(5) * f_loc(ic+2,jc,kc)
!      enddo 
!    enddo
!  endif
!
!  if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!      do jc = 1,Ny
!        ic = nx_local
!        df(ic,jc,kc) =- BC1(1) * f_loc(ic  ,jc,kc) &
!                    &- BC1(2) * f_loc(ic-1,jc,kc) &
!                    &- BC1(3) * f_loc(ic-2,jc,kc) &
!                    &- BC1(4) * f_loc(ic-3,jc,kc) &
!                    &- BC1(5) * f_loc(ic-4,jc,kc)
!
!        ic = nx_local-1
!        df(ic,jc,kc) =- BC2(1) * f_loc(ic+1,jc,kc) &
!                    &- BC2(2) * f_loc(ic  ,jc,kc) &
!                    &- BC2(3) * f_loc(ic-1,jc,kc) &
!                    &- BC2(4) * f_loc(ic-2,jc,kc) &
!                    &- BC2(5) * f_loc(ic-3,jc,kc)
!
!        ic = nx_local-2                                
!        df(ic,jc,kc) = -BC3(1) * f_loc(ic+2,jc,kc) &
!                    &- BC3(2) * f_loc(ic+1,jc,kc) &
!                    &- BC3(3) * f_loc(ic  ,jc,kc) &
!                    &- BC3(4) * f_loc(ic-1,jc,kc) &
!                    &- BC3(5) * f_loc(ic-2,jc,kc)
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dxi_weno5th_do
!
!! Eta direction 5th upwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Deta_weno5th_up_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 4, Ny - 3
!      do ic = 1, nx_local
!        do iVar = 1, NumVar
!          varloc(-2:2) = f_loc(iVar, ic, jc-2:jc+2, kc)
!          fp = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic, jc-3:jc+1, kc)
!          fm = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = 1
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC1(2) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC1(3) * f_loc(iVar,ic,jc+2,kc) &
!                         &+ BC1(4) * f_loc(iVar,ic,jc+3,kc) &
!                         &+ BC1(5) * f_loc(iVar,ic,jc+4,kc)
!      enddo 
!      jc = 2
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic,jc-1,kc) &
!                         &+ BC2(2) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC2(3) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC2(4) * f_loc(iVar,ic,jc+2,kc) &
!                         &+ BC2(5) * f_loc(iVar,ic,jc+3,kc)
!      enddo
!      jc = 3
!      do iVar = 1, NumVar                                
!        df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic,jc-2,kc) &
!                         &+ BC3(2) * f_loc(iVar,ic,jc-1,kc) &
!                         &+ BC3(3) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC3(4) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC3(5) * f_loc(iVar,ic,jc+2,kc)
!      enddo
!    enddo 
!  enddo
!
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = Ny
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC1(2) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC1(3) * f_loc(iVar,ic,jc-2,kc) &
!                          &- BC1(4) * f_loc(iVar,ic,jc-3,kc) &
!                          &- BC1(5) * f_loc(iVar,ic,jc-4,kc)
!      enddo
!      jc = Ny - 1
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic,jc+1,kc) &
!                          &- BC2(2) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC2(3) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC2(4) * f_loc(iVar,ic,jc-2,kc) &
!                          &- BC2(5) * f_loc(iVar,ic,jc-3,kc)
!      enddo
!      jc = Ny - 2
!      do iVar = 1, NumVar                                
!        df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic,jc+2,kc) &
!                          &- BC3(2) * f_loc(iVar,ic,jc+1,kc) &
!                          &- BC3(3) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC3(4) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC3(5) * f_loc(iVar,ic,jc-2,kc)
!      enddo
!    enddo 
!  enddo
!end subroutine Cal_Deri_Deta_weno5th_up_numvar
!
!
!! Eta direction 5th downwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Deta_weno5th_do_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 4, Ny-3
!      do ic = 1, nx_local
!        do iVar = 1, NumVar
!          varloc(-2:2) = f_loc(iVar, ic, jc+2:jc-2:-1, kc)
!          fm = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic, jc+3:jc-1:-1, kc)
!          fp = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = 1
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC1(2) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC1(3) * f_loc(iVar,ic,jc+2,kc) &
!                         &+ BC1(4) * f_loc(iVar,ic,jc+3,kc) &
!                         &+ BC1(5) * f_loc(iVar,ic,jc+4,kc)
!      enddo 
!      jc = 2
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic,jc-1,kc) &
!                         &+ BC2(2) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC2(3) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC2(4) * f_loc(iVar,ic,jc+2,kc) &
!                         &+ BC2(5) * f_loc(iVar,ic,jc+3,kc)
!      enddo
!      jc = 3
!      do iVar = 1, NumVar                                
!        df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic,jc-2,kc) &
!                         &+ BC3(2) * f_loc(iVar,ic,jc-1,kc) &
!                         &+ BC3(3) * f_loc(iVar,ic,jc  ,kc) &
!                         &+ BC3(4) * f_loc(iVar,ic,jc+1,kc) &
!                         &+ BC3(5) * f_loc(iVar,ic,jc+2,kc)
!      enddo
!    enddo 
!  enddo
!
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = Ny
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC1(2) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC1(3) * f_loc(iVar,ic,jc-2,kc) &
!                          &- BC1(4) * f_loc(iVar,ic,jc-3,kc) &
!                          &- BC1(5) * f_loc(iVar,ic,jc-4,kc)
!      enddo
!      jc = Ny - 1
!      do iVar = 1, NumVar
!        df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic,jc+1,kc) &
!                          &- BC2(2) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC2(3) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC2(4) * f_loc(iVar,ic,jc-2,kc) &
!                          &- BC2(5) * f_loc(iVar,ic,jc-3,kc)
!      enddo
!      jc = Ny - 2
!      do iVar = 1, NumVar                                
!        df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic,jc+2,kc) &
!                          &- BC3(2) * f_loc(iVar,ic,jc+1,kc) &
!                          &- BC3(3) * f_loc(iVar,ic,jc  ,kc) &
!                          &- BC3(4) * f_loc(iVar,ic,jc-1,kc) &
!                          &- BC3(5) * f_loc(iVar,ic,jc-2,kc)
!      enddo
!    enddo 
!  enddo
!end subroutine Cal_Deri_Deta_weno5th_do_numvar
!
!
!! Eta direction 5th upwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Deta_weno5th_up(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 4, Ny-3
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc-2:jc+2, kc)
!        fp = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc-3:jc+1, kc)
!        fm = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = 1
!      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc  ,kc) &
!                  &+ BC1(2) * f_loc(ic,jc+1,kc) &
!                  &+ BC1(3) * f_loc(ic,jc+2,kc) &
!                  &+ BC1(4) * f_loc(ic,jc+3,kc) &
!                  &+ BC1(5) * f_loc(ic,jc+4,kc)
!                  
!      jc = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc-1,kc) &
!                  &+ BC2(2) * f_loc(ic,jc  ,kc) &
!                  &+ BC2(3) * f_loc(ic,jc+1,kc) &
!                  &+ BC2(4) * f_loc(ic,jc+2,kc) &
!                  &+ BC2(5) * f_loc(ic,jc+3,kc)
!
!      jc = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc-2,kc) &
!                  &+ BC3(2) * f_loc(ic,jc-1,kc) &
!                  &+ BC3(3) * f_loc(ic,jc  ,kc) &
!                  &+ BC3(4) * f_loc(ic,jc+1,kc) &
!                  &+ BC3(5) * f_loc(ic,jc+2,kc)
!    enddo 
!  enddo
!
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = Ny
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc  ,kc) &
!                   &- BC1(2) * f_loc(ic,jc-1,kc) &
!                   &- BC1(3) * f_loc(ic,jc-2,kc) &
!                   &- BC1(4) * f_loc(ic,jc-3,kc) &
!                   &- BC1(5) * f_loc(ic,jc-4,kc)
!
!      jc = Ny - 1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc+1,kc) &
!                   &- BC2(2) * f_loc(ic,jc  ,kc) &
!                   &- BC2(3) * f_loc(ic,jc-1,kc) &
!                   &- BC2(4) * f_loc(ic,jc-2,kc) &
!                   &- BC2(5) * f_loc(ic,jc-3,kc)
!
!      jc = Ny - 2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc+2,kc) &
!                   &- BC3(2) * f_loc(ic,jc+1,kc) &
!                   &- BC3(3) * f_loc(ic,jc  ,kc) &
!                   &- BC3(4) * f_loc(ic,jc-1,kc) &
!                   &- BC3(5) * f_loc(ic,jc-2,kc)
!    enddo 
!  enddo
!end subroutine Cal_Deri_Deta_weno5th_up
!
!
!! Eta direction 5th downwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Deta_weno5th_do(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 4, Ny-3
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc+2:jc-2:-1, kc)
!        fm = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc+3:jc-1:-1, kc)
!        fp = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = 1
!      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc  ,kc) &
!                  &+ BC1(2) * f_loc(ic,jc+1,kc) &
!                  &+ BC1(3) * f_loc(ic,jc+2,kc) &
!                  &+ BC1(4) * f_loc(ic,jc+3,kc) &
!                  &+ BC1(5) * f_loc(ic,jc+4,kc)
!                  
!      jc = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc-1,kc) &
!                  &+ BC2(2) * f_loc(ic,jc  ,kc) &
!                  &+ BC2(3) * f_loc(ic,jc+1,kc) &
!                  &+ BC2(4) * f_loc(ic,jc+2,kc) &
!                  &+ BC2(5) * f_loc(ic,jc+3,kc)
!
!      jc = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc-2,kc) &
!                  &+ BC3(2) * f_loc(ic,jc-1,kc) &
!                  &+ BC3(3) * f_loc(ic,jc  ,kc) &
!                  &+ BC3(4) * f_loc(ic,jc+1,kc) &
!                  &+ BC3(5) * f_loc(ic,jc+2,kc)
!    enddo 
!  enddo
!
!  do kc = 1,nz_local
!    do ic = 1,nx_local
!      jc = Ny
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc  ,kc) &
!                   &- BC1(2) * f_loc(ic,jc-1,kc) &
!                   &- BC1(3) * f_loc(ic,jc-2,kc) &
!                   &- BC1(4) * f_loc(ic,jc-3,kc) &
!                   &- BC1(5) * f_loc(ic,jc-4,kc)
!
!      jc = Ny - 1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc+1,kc) &
!                   &- BC2(2) * f_loc(ic,jc  ,kc) &
!                   &- BC2(3) * f_loc(ic,jc-1,kc) &
!                   &- BC2(4) * f_loc(ic,jc-2,kc) &
!                   &- BC2(5) * f_loc(ic,jc-3,kc)
!
!      jc = Ny - 2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc+2,kc) &
!                   &- BC3(2) * f_loc(ic,jc+1,kc) &
!                   &- BC3(3) * f_loc(ic,jc  ,kc) &
!                   &- BC3(4) * f_loc(ic,jc-1,kc) &
!                   &- BC3(5) * f_loc(ic,jc-2,kc)
!    enddo 
!  enddo
!end subroutine Cal_Deri_Deta_weno5th_do
!
!
!! Zeta direction 5th upwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_up(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc, kc-2:kc+2)
!        fp = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc, kc-3:kc+1)
!        fm = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npz == 0) then
!    do jc = 1,Ny
!      do ic = 1,nx_local
!        kc = 1
!        df(ic,jc,kc) = BC1(1) * f_loc(ic,jc,kc  ) &
!                    &+ BC1(2) * f_loc(ic,jc,kc+1) &
!                    &+ BC1(3) * f_loc(ic,jc,kc+2) &
!                    &+ BC1(4) * f_loc(ic,jc,kc+3) &
!                    &+ BC1(5) * f_loc(ic,jc,kc+4)
!                    
!        kc = 2
!        df(ic,jc,kc) = BC2(1) * f_loc(ic,jc,kc-1) &
!                    &+ BC2(2) * f_loc(ic,jc,kc  ) &
!                    &+ BC2(3) * f_loc(ic,jc,kc+1) &
!                    &+ BC2(4) * f_loc(ic,jc,kc+2) &
!                    &+ BC2(5) * f_loc(ic,jc,kc+3)
!
!        kc = 3                                
!        df(ic,jc,kc) = BC3(1) * f_loc(ic,jc,kc-2) &
!                    &+ BC3(2) * f_loc(ic,jc,kc-1) &
!                    &+ BC3(3) * f_loc(ic,jc,kc  ) &
!                    &+ BC3(4) * f_loc(ic,jc,kc+1) &
!                    &+ BC3(5) * f_loc(ic,jc,kc+2)
!      enddo 
!    enddo
!  endif
!
!  if(npz == npz0 - 1) then
!    do jc = 1,Ny
!      do ic = 1,nx_local
!        kc = nz_local
!        df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc,kc  ) &
!                     &- BC1(2) * f_loc(ic,jc,kc-1) &
!                     &- BC1(3) * f_loc(ic,jc,kc-2) &
!                     &- BC1(4) * f_loc(ic,jc,kc-3) &
!                     &- BC1(5) * f_loc(ic,jc,kc-4)
!
!        kc = nz_local-1
!        df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc,kc+1) &
!                     &- BC2(2) * f_loc(ic,jc,kc  ) &
!                     &- BC2(3) * f_loc(ic,jc,kc-1) &
!                     &- BC2(4) * f_loc(ic,jc,kc-2) &
!                     &- BC2(5) * f_loc(ic,jc,kc-3)
!
!        kc = nz_local-2                                
!        df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc,kc+2) &
!                     &- BC3(2) * f_loc(ic,jc,kc+1) &
!                     &- BC3(3) * f_loc(ic,jc,kc  ) &
!                     &- BC3(4) * f_loc(ic,jc,kc-1) &
!                     &- BC3(5) * f_loc(ic,jc,kc-2)
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dzeta_weno5th_up
!
!
!! Zeta direction 5th downwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_do(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc, kc+2:kc-2:-1)
!        fm = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc, kc+3:kc-1:-1)
!        fp = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!
!  ! Boundarys
!  if(npz == 0) then
!    do jc = 1,Ny
!      do ic = 1,nx_local
!        kc = 1
!        df(ic,jc,kc) = BC1(1) * f_loc(ic,jc,kc  ) &
!                    &+ BC1(2) * f_loc(ic,jc,kc+1) &
!                    &+ BC1(3) * f_loc(ic,jc,kc+2) &
!                    &+ BC1(4) * f_loc(ic,jc,kc+3) &
!                    &+ BC1(5) * f_loc(ic,jc,kc+4)
!                    
!        kc = 2
!        df(ic,jc,kc) = BC2(1) * f_loc(ic,jc,kc-1) &
!                    &+ BC2(2) * f_loc(ic,jc,kc  ) &
!                    &+ BC2(3) * f_loc(ic,jc,kc+1) &
!                    &+ BC2(4) * f_loc(ic,jc,kc+2) &
!                    &+ BC2(5) * f_loc(ic,jc,kc+3)
!
!        kc = 3                                
!        df(ic,jc,kc) = BC3(1) * f_loc(ic,jc,kc-2) &
!                    &+ BC3(2) * f_loc(ic,jc,kc-1) &
!                    &+ BC3(3) * f_loc(ic,jc,kc  ) &
!                    &+ BC3(4) * f_loc(ic,jc,kc+1) &
!                    &+ BC3(5) * f_loc(ic,jc,kc+2)
!      enddo 
!    enddo
!  endif
!
!  if(npz == npz0 - 1) then
!    do jc = 1,Ny
!      do ic = 1,nx_local
!        kc = nz_local
!        df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc,kc  ) &
!                     &- BC1(2) * f_loc(ic,jc,kc-1) &
!                     &- BC1(3) * f_loc(ic,jc,kc-2) &
!                     &- BC1(4) * f_loc(ic,jc,kc-3) &
!                     &- BC1(5) * f_loc(ic,jc,kc-4)
!
!        kc = nz_local-1
!        df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc,kc+1) &
!                     &- BC2(2) * f_loc(ic,jc,kc  ) &
!                     &- BC2(3) * f_loc(ic,jc,kc-1) &
!                     &- BC2(4) * f_loc(ic,jc,kc-2) &
!                     &- BC2(5) * f_loc(ic,jc,kc-3)
!
!        kc = nz_local-2                                
!        df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc,kc+2) &
!                     &- BC3(2) * f_loc(ic,jc,kc+1) &
!                     &- BC3(3) * f_loc(ic,jc,kc  ) &
!                     &- BC3(4) * f_loc(ic,jc,kc-1) &
!                     &- BC3(5) * f_loc(ic,jc,kc-2)
!      enddo 
!    enddo
!  endif
!end subroutine Cal_Deri_Dzeta_weno5th_do
!
!
!! Zeta direction(periodic) 5th upwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_up_per_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        do iVar = 1,NumVar
!          varloc(-2:2) = f_loc(iVar, ic, jc, kc-2:kc+2)
!          fp = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic, jc, kc-3:kc+1)
!          fm = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!end subroutine Cal_Deri_Dzeta_weno5th_up_per_numvar
!
!
!! Zeta direction(periodic) 5th downwind WENO normal derivatives f(:,:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_do_per_numvar(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc,iVar 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        do iVar = 1,NumVar
!          varloc(-2:2) = f_loc(iVar, ic, jc, kc+2:kc-2:-1)
!          fm = weno_halfvalue(varloc)
!          varloc(-2:2) = f_loc(iVar, ic, jc, kc+3:kc-1:-1)
!          fp = weno_halfvalue(varloc)
!          df(iVar, ic, jc, kc) = fp - fm
!        enddo
!      enddo
!    enddo
!  enddo
!end subroutine Cal_Deri_Dzeta_weno5th_do_per_numvar
!
!
!! Zeta direction(periodic) 5th upwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_up_per(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc, kc-2:kc+2)
!        fp = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc, kc-3:kc+1)
!        fm = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!end subroutine Cal_Deri_Dzeta_weno5th_up_per
!
!
!! Zeta direction(periodic) 5th downwind WENO normal derivatives f(:,:,:)
!subroutine Cal_Deri_Dzeta_weno5th_do_per(df,f_loc)
!  implicit none
!  real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!  
!  real( kind = rk ):: varloc(-2:2),fp,fm
!  integer( kind = ik )::ic,jc,kc 
!
!  df = 0.d0;
!  ! To ensure the memory consist 
!  do kc = 1, nz_local
!    do jc = 1, Ny
!      do ic = 1, nx_local
!        varloc(-2:2) = f_loc(ic, jc, kc+2:kc-2:-1)
!        fm = weno_halfvalue(varloc)
!        varloc(-2:2) = f_loc(ic, jc, kc+3:kc-1:-1)
!        fp = weno_halfvalue(varloc)
!        df(ic, jc, kc) = fp - fm
!      enddo
!    enddo
!  enddo
!end subroutine Cal_Deri_Dzeta_weno5th_do_per
!
!! 2D Surface Derivatives: use Module_FD_5th.f90
!! 1D Line Derivatives: use Module_FD_5th.f90
!
!!===============================================================================
!! Characteristic Reconstruction, local, including boundarys
!!===============================================================================
!! Xi direction 5th upwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Dxi_wenochar5thlocal_up_numvar(df,f_loc,i)
!  implicit none
!  integer( kind = ik ), intent(in):: i
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(npx == 0 .and. i == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(3, iVar) &
!              & + BCchar1(2) * f_loc(4, iVar) &
!              & + BCchar1(3) * f_loc(5, iVar)
!    enddo
!  elseif(npx == 0 .and. i == 2)then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(:,iVar))
!    enddo
!  elseif(npx == npx0-1 .and. i == nx_local - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(:,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(:,iVar))
!    enddo
!  endif    
!end subroutine Cal_Deri_Dxi_wenochar5th_up_numvar
!
!! Xi direction 5th downwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Dxi_wenochar5th_do_numvar(df,f_loc,i)
!  implicit none
!  integer( kind = ik ), intent(in):: i
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(npx == 0 .and. i == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(5:1:-1,iVar))
!    enddo
!  elseif(npx == npx0-1 .and. i == nx_local - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(4, iVar) &
!              & + BCchar1(2) * f_loc(3, iVar) &
!              & + BCchar1(3) * f_loc(2, iVar)
!    enddo
!  elseif(npx == npx0-1 .and. i == nx_local - 2) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(5:1:-1,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(5:1:-1,iVar))
!    enddo
!  endif 
!end subroutine Cal_Deri_Dxi_wenochar5th_do_numvar
!
!! Eta direction 5th upwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Deta_wenochar5th_up_numvar(df,f_loc,j)
!  implicit none
!  integer( kind = ik ), intent(in):: j
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(j == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(3, iVar) &
!              & + BCchar1(2) * f_loc(4, iVar) &
!              & + BCchar1(3) * f_loc(5, iVar)
!    enddo
!  elseif(j == 2)then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(:,iVar))
!    enddo
!  elseif(j == Ny - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(:,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(:,iVar))
!    enddo
!  endif    
!end subroutine Cal_Deri_Deta_wenochar5th_up_numvar
!
!! Eta direction 5th downwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Deta_wenochar5th_do_numvar(df,f_loc,j)
!  implicit none
!  integer( kind = ik ), intent(in):: j
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(j == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(5:1:-1,iVar))
!    enddo
!  elseif(j == Ny - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(4, iVar) &
!              & + BCchar1(2) * f_loc(3, iVar) &
!              & + BCchar1(3) * f_loc(2, iVar)
!    enddo
!  elseif(j == Ny - 2) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(5:1:-1,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(5:1:-1,iVar))
!    enddo
!  endif 
!end subroutine Cal_Deri_Deta_wenochar5th_do_numvar
!
!! Zeta direction 5th upwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Dzeta_wenochar5th_up_numvar(df,f_loc,k)
!  implicit none
!  integer( kind = ik ), intent(in):: k
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(npz == 0 .and. k == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(3, iVar) &
!              & + BCchar1(2) * f_loc(4, iVar) &
!              & + BCchar1(3) * f_loc(5, iVar)
!    enddo
!  elseif(npz == 0 .and. k == 2)then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(:,iVar))
!    enddo
!  elseif(npz == npz0-1 .and. k == nz_local - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(:,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(:,iVar))
!    enddo
!  endif    
!end subroutine Cal_Deri_Dzeta_wenochar5th_up_numvar
!
!! Zeta direction 5th downwind WENO characteristic derivatives f(NumTemplate,NumVar)
!subroutine Cal_Deri_Dzeta_wenochar5th_do_numvar(df,f_loc,k)
!  implicit none
!  integer( kind = ik ), intent(in):: k
!  real( kind = rk ), intent(in):: f_loc(5,1:NumVar)
!  real( kind = rk ), intent(out):: df(1:NumVar)
!  integer( kind = ik ):: iVar 
!
!  df = 0.0_rk
!  ! Boundarys if
!  if(npz == 0 .and. k == 1) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_left(f_loc(5:1:-1,iVar))
!    enddo
!  elseif(npz == npz0-1 .and. k == nz_local - 1) then
!    do iVar = 1, NumVar
!      df(iVar) = BCchar1(1) * f_loc(4, iVar) &
!              & + BCchar1(2) * f_loc(3, iVar) &
!              & + BCchar1(3) * f_loc(2, iVar)
!    enddo
!  elseif(npz == npz0-1 .and. k == nz_local - 2) then
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue_del_right(f_loc(5:1:-1,iVar))
!    enddo
!  else
!    ! General WENO scheme
!    do iVar = 1, NumVar
!      df(iVar) = weno_halfvalue(f_loc(5:1:-1,iVar))
!    enddo
!  endif 
!end subroutine Cal_Deri_Dzeta_wenochar5th_do_numvar
!
!END MODULE FD5_Order_Weno
