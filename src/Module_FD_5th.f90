MODULE FD5_Order
  ! This is the module that contains all the subs and dates
  ! used for 5th order upwind/downwind and 6th order finite
  ! difference methods
  use SF_Constant,  only: ik, rk, NumVar, overLap
  use SF_CFD_Global,only: nx_local,Ny,nz_local,Nz,ModelType
  use MPI_Global,   only: if_parallel, MyId, NumProcess, MPI_COMM_WORLD, ierr,&
                        & npx,npx0,npz,npz0

  implicit none
  ! in order to improve the calculation speed, all the coefficients for finite
  ! difference calculation are used as constants
  
  real( kind = rk ), parameter::BC1(1:5) =(/ -25.d0/12.d0,      4.d0,  -3.d0, 4.d0/3.d0,    -0.25d0 /) 
  real( kind = rk ), parameter::BC2(1:5) =(/      -0.25d0,-5.d0/6.d0,  1.5d0,    -0.5d0, 1.d0/12.d0 /) 
  real( kind = rk ), parameter::BC3(1:5) =(/   1.d0/12.d0,-2.d0/3.d0,   0.d0, 2.d0/3.d0,-1.d0/12.d0 /)
  real( kind = rk ), parameter::BI = 60.d0;
  real( kind = rk ), parameter::alpha_up = -6.d0 
  real( kind = rk ), parameter::CO_FD_Upwind(-3:3) = (/( -1.d0 + 1.d0 / 12.d0 * alpha_up) / BI, &
                                                    &  (  9.d0 - 0.5d0 * alpha_up) / BI, &
                                                    &  (-45.d0 + 1.25d0 * alpha_up)/ BI, & 
                                                    &  ( -5.D0 / 3.D0 * alpha_up)/     BI, &
                                                    &  ( 45.d0 + 1.25d0 * alpha_up)/ BI, & 
                                                    &  ( -9.D0 - 0.5D0 * alpha_up) / BI, &
                                                    &  (  1.D0 + 1.D0 / 12.D0 * alpha_up) / BI/)

  real( kind = rk ), parameter::alpha_ce =  0.d0
  real( kind = rk ), parameter::CO_FD_Center(-3:3) = (/( -1.d0 + 1.d0 / 12.d0 * alpha_ce) / BI, &
                                                    &  (  9.d0 - 0.5d0 * alpha_ce) / BI, &
                                                    &  (-45.d0 + 1.25d0 * alpha_ce)/ BI, & 
                                                    &  ( -5.D0 / 3.D0 * alpha_ce)/     BI, &
                                                    &  ( 45.d0 + 1.25d0 * alpha_ce)/ BI, & 
                                                    &  ( -9.D0 - 0.5D0 * alpha_ce) / BI, &
                                                    &  (  1.D0 + 1.D0 / 12.D0 * alpha_ce) / BI/)

  real( kind = rk ), parameter::alpha_do =  6.d0 
  real( kind = rk ), parameter::CO_FD_Dowind(-3:3) = (/( -1.d0 + 1.d0 / 12.d0 * alpha_do) / BI, &
                                                    &  (  9.d0 - 0.5d0 * alpha_do) / BI, &
                                                    &  (-45.d0 + 1.25d0 * alpha_do)/ BI, & 
                                                    &  ( -5.D0 / 3.D0 * alpha_do)/BI     , &
                                                    &  ( 45.d0 + 1.25d0 * alpha_do)/ BI, & 
                                                    &  ( -9.D0 - 0.5D0 * alpha_do) / BI, &
                                                    &  (  1.D0 + 1.D0 / 12.D0 * alpha_do) / BI/) 

  contains

! Xi direction normal derivatives f(:,:,:,:)
subroutine Cal_Deri_Dxi_5th_up_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(iVar,ic - 3,jc,kc) &
                         &+ CO_FD_Upwind(-2) * f_loc(iVar,ic - 2,jc,kc) &
                         &+ CO_FD_Upwind(-1) * f_loc(iVar,ic - 1,jc,kc) &
                         &+ CO_FD_Upwind( 0) * f_loc(iVar,ic    ,jc,kc) &
                         &+ CO_FD_Upwind( 1) * f_loc(iVar,ic + 1,jc,kc) &
                         &+ CO_FD_Upwind( 2) * f_loc(iVar,ic + 2,jc,kc) &
                         &+ CO_FD_Upwind( 3) * f_loc(iVar,ic + 3,jc,kc)
        enddo
       enddo
     enddo
   enddo

   if( ModelType == 1)then 
   if(npx == 0) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = 1
      do iVar = 1, NumVar            
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC1(2) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC1(3) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC1(4) * f_loc(iVar,ic+3,jc,kc) &
                       &+ BC1(5) * f_loc(iVar,ic+4,jc,kc)
      enddo
                  
      ic = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC2(3) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC2(4) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC2(5) * f_loc(iVar,ic+3,jc,kc)
      enddo
    
      ic = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic-2,jc,kc) &
                       &+ BC3(2) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC3(4) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC3(5) * f_loc(iVar,ic+2,jc,kc)
      enddo
     enddo 
    enddo
   endif
   endif
   
   if(npx == npx0 - 1) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = nx_local
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC1(2) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC1(3) * f_loc(iVar,ic-2,jc,kc) &
                        &- BC1(4) * f_loc(iVar,ic-3,jc,kc) &
                        &- BC1(5) * f_loc(iVar,ic-4,jc,kc)
      enddo

      ic = nx_local-1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic+1,jc,kc) &
                        &- BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC2(3) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC2(4) * f_loc(iVar,ic-2,jc,kc) &
                        &- BC2(5) * f_loc(iVar,ic-3,jc,kc)
      enddo

      ic = nx_local-2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic+2,jc,kc) &
                        &- BC3(2) * f_loc(iVar,ic+1,jc,kc) &
                        &- BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC3(4) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC3(5) * f_loc(iVar,ic-2,jc,kc)
      enddo
     enddo 
    enddo
   endif
end subroutine Cal_Deri_Dxi_5th_up_numvar

subroutine Cal_Deri_Dxi_5th_ce_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Center(-3) * f_loc(iVar,ic - 3,jc,kc) &
                         &+ CO_FD_Center(-2) * f_loc(iVar,ic - 2,jc,kc) &
                         &+ CO_FD_Center(-1) * f_loc(iVar,ic - 1,jc,kc) &
                         &+ CO_FD_Center( 0) * f_loc(iVar,ic    ,jc,kc) &
                         &+ CO_FD_Center( 1) * f_loc(iVar,ic + 1,jc,kc) &
                         &+ CO_FD_Center( 2) * f_loc(iVar,ic + 2,jc,kc) &
                         &+ CO_FD_Center( 3) * f_loc(iVar,ic + 3,jc,kc)
        enddo
       enddo
     enddo
   enddo

   
   if( ModelType == 1)then 
   if(npx == 0) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = 1
      do iVar = 1, NumVar  
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC1(2) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC1(3) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC1(4) * f_loc(iVar,ic+3,jc,kc) &
                       &+ BC1(5) * f_loc(iVar,ic+4,jc,kc)
    
      enddo
         
      ic = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC2(3) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC2(4) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC2(5) * f_loc(iVar,ic+3,jc,kc)
      enddo
    
      ic = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic-2,jc,kc) &
                       &+ BC3(2) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC3(4) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC3(5) * f_loc(iVar,ic+2,jc,kc)
      enddo
     enddo 
    enddo
   endif
   endif
   
   if(npx == npx0 - 1) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = nx_local
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC1(2) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC1(3) * f_loc(iVar,ic-2,jc,kc) &
                        &- BC1(4) * f_loc(iVar,ic-3,jc,kc) &
                        &- BC1(5) * f_loc(iVar,ic-4,jc,kc)
      enddo

      ic = nx_local-1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic+1,jc,kc) &
                   &- BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                   &- BC2(3) * f_loc(iVar,ic-1,jc,kc) &
                   &- BC2(4) * f_loc(iVar,ic-2,jc,kc) &
                   &- BC2(5) * f_loc(iVar,ic-3,jc,kc)
      enddo

      ic = nx_local-2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic+2,jc,kc) &
                   &- BC3(2) * f_loc(iVar,ic+1,jc,kc) &
                   &- BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                   &- BC3(4) * f_loc(iVar,ic-1,jc,kc) &
                   &- BC3(5) * f_loc(iVar,ic-2,jc,kc)
      enddo
     enddo 
    enddo
   endif
   
   
end subroutine Cal_Deri_Dxi_5th_ce_numvar

subroutine Cal_Deri_Dxi_5th_do_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(iVar,ic - 3,jc,kc) &
                         &+ CO_FD_Dowind(-2) * f_loc(iVar,ic - 2,jc,kc) &
                         &+ CO_FD_Dowind(-1) * f_loc(iVar,ic - 1,jc,kc) &
                         &+ CO_FD_Dowind( 0) * f_loc(iVar,ic    ,jc,kc) &
                         &+ CO_FD_Dowind( 1) * f_loc(iVar,ic + 1,jc,kc) &
                         &+ CO_FD_Dowind( 2) * f_loc(iVar,ic + 2,jc,kc) &
                         &+ CO_FD_Dowind( 3) * f_loc(iVar,ic + 3,jc,kc)
        enddo
       enddo
     enddo
   enddo
   
   if( ModelType == 1)then  
   if(npx == 0) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC1(2) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC1(3) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC1(4) * f_loc(iVar,ic+3,jc,kc) &
                       &+ BC1(5) * f_loc(iVar,ic+4,jc,kc)    
      enddo
      ic = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC2(3) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC2(4) * f_loc(iVar,ic+2,jc,kc) &
                       &+ BC2(5) * f_loc(iVar,ic+3,jc,kc)
      enddo
    
      ic = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic-2,jc,kc) &
                       &+ BC3(2) * f_loc(iVar,ic-1,jc,kc) &
                       &+ BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                       &+ BC3(4) * f_loc(iVar,ic+1,jc,kc) &
                       &+ BC3(5) * f_loc(iVar,ic+2,jc,kc)
      enddo
     enddo 
    enddo
   endif
   endif
   
   if(npx == npx0 - 1) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = nx_local
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC1(2) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC1(3) * f_loc(iVar,ic-2,jc,kc) &
                        &- BC1(4) * f_loc(iVar,ic-3,jc,kc) &
                        &- BC1(5) * f_loc(iVar,ic-4,jc,kc)
      enddo

      ic = nx_local-1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic+1,jc,kc) &
                        &- BC2(2) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC2(3) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC2(4) * f_loc(iVar,ic-2,jc,kc) &
                        &- BC2(5) * f_loc(iVar,ic-3,jc,kc)
      enddo

      ic = nx_local-2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) =- BC3(1) * f_loc(iVar,ic+2,jc,kc) &
                        &- BC3(2) * f_loc(iVar,ic+1,jc,kc) &
                        &- BC3(3) * f_loc(iVar,ic  ,jc,kc) &
                        &- BC3(4) * f_loc(iVar,ic-1,jc,kc) &
                        &- BC3(5) * f_loc(iVar,ic-2,jc,kc)
      enddo
     enddo 
    enddo
   endif
  
   
end subroutine Cal_Deri_Dxi_5th_do_numvar

! Xi direction normal derivatives f(:,:,:) 
!subroutine Cal_Deri_Dxi_5th_up(df,f_loc)
!   implicit none
!   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   
!   integer( kind = ik )::ic,jc,kc 
!
!   df = 0.d0;
!   ! To ensure the memory consist 
!   do kc = 1, nz_local
!     do jc = 1, Ny
!       do ic = 1, nx_local
!        df(ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(ic - 3,jc,kc) &
!                    &+ CO_FD_Upwind(-2) * f_loc(ic - 2,jc,kc) &
!                    &+ CO_FD_Upwind(-1) * f_loc(ic - 1,jc,kc) &
!                    &+ CO_FD_Upwind( 0) * f_loc(ic    ,jc,kc) &
!                    &+ CO_FD_Upwind( 1) * f_loc(ic + 1,jc,kc) &
!                    &+ CO_FD_Upwind( 2) * f_loc(ic + 2,jc,kc) &
!                    &+ CO_FD_Upwind( 3) * f_loc(ic + 3,jc,kc)
!       enddo
!     enddo
!   enddo
!   
! if( ModelType == 1)then 
!   if(npx == 0) then
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      ic = 1
!      df(ic,jc,kc) = BC1(1) * f_loc(ic  ,jc,kc) &
!                  &+ BC1(2) * f_loc(ic+1,jc,kc) &
!                  &+ BC1(3) * f_loc(ic+2,jc,kc) &
!                  &+ BC1(4) * f_loc(ic+3,jc,kc) &
!                  &+ BC1(5) * f_loc(ic+4,jc,kc)
!      
!      ic = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic-1,jc,kc) &
!                  &+ BC2(2) * f_loc(ic  ,jc,kc) &
!                  &+ BC2(3) * f_loc(ic+1,jc,kc) &
!                  &+ BC2(4) * f_loc(ic+2,jc,kc) &
!                  &+ BC2(5) * f_loc(ic+3,jc,kc)
!    
!    
!      ic = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic-2,jc,kc) &
!                  &+ BC3(2) * f_loc(ic-1,jc,kc) &
!                  &+ BC3(3) * f_loc(ic  ,jc,kc) &
!                  &+ BC3(4) * f_loc(ic+1,jc,kc) &
!                  &+ BC3(5) * f_loc(ic+2,jc,kc)
!     enddo 
!    enddo
!   endif
!endif
!   if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      ic = nx_local
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic  ,jc,kc) &
!                   &- BC1(2) * f_loc(ic-1,jc,kc) &
!                   &- BC1(3) * f_loc(ic-2,jc,kc) &
!                   &- BC1(4) * f_loc(ic-3,jc,kc) &
!                   &- BC1(5) * f_loc(ic-4,jc,kc)
!
!      ic = nx_local-1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic+1,jc,kc) &
!                   &- BC2(2) * f_loc(ic  ,jc,kc) &
!                   &- BC2(3) * f_loc(ic-1,jc,kc) &
!                   &- BC2(4) * f_loc(ic-2,jc,kc) &
!                   &- BC2(5) * f_loc(ic-3,jc,kc)
!
!      ic = nx_local-2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic+2,jc,kc) &
!                   &- BC3(2) * f_loc(ic+1,jc,kc) &
!                   &- BC3(3) * f_loc(ic  ,jc,kc) &
!                   &- BC3(4) * f_loc(ic-1,jc,kc) &
!                   &- BC3(5) * f_loc(ic-2,jc,kc)
!     enddo 
!    enddo
!   endif
! 
!end subroutine Cal_Deri_Dxi_5th_up

subroutine Cal_Deri_Dxi_5th_ce(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Center(-3) * f_loc(ic - 3,jc,kc) &
                    &+ CO_FD_Center(-2) * f_loc(ic - 2,jc,kc) &
                    &+ CO_FD_Center(-1) * f_loc(ic - 1,jc,kc) &
                    &+ CO_FD_Center( 0) * f_loc(ic    ,jc,kc) &
                    &+ CO_FD_Center( 1) * f_loc(ic + 1,jc,kc) &
                    &+ CO_FD_Center( 2) * f_loc(ic + 2,jc,kc) &
                    &+ CO_FD_Center( 3) * f_loc(ic + 3,jc,kc)
       enddo
     enddo
   enddo
   
  if( ModelType == 1)then
   if(npx == 0) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = 1
      df(ic,jc,kc) = BC1(1) * f_loc(ic  ,jc,kc) &
                  &+ BC1(2) * f_loc(ic+1,jc,kc) &
                  &+ BC1(3) * f_loc(ic+2,jc,kc) &
                  &+ BC1(4) * f_loc(ic+3,jc,kc) &
                  &+ BC1(5) * f_loc(ic+4,jc,kc)
  
      ic = 2
      df(ic,jc,kc) = BC2(1) * f_loc(ic-1,jc,kc) &
                  &+ BC2(2) * f_loc(ic  ,jc,kc) &
                  &+ BC2(3) * f_loc(ic+1,jc,kc) &
                  &+ BC2(4) * f_loc(ic+2,jc,kc) &
                  &+ BC2(5) * f_loc(ic+3,jc,kc)
  
  
      ic = 3                                
      df(ic,jc,kc) = BC3(1) * f_loc(ic-2,jc,kc) &
                  &+ BC3(2) * f_loc(ic-1,jc,kc) &
                  &+ BC3(3) * f_loc(ic  ,jc,kc) &
                  &+ BC3(4) * f_loc(ic+1,jc,kc) &
                  &+ BC3(5) * f_loc(ic+2,jc,kc)
     enddo 
    enddo
   endif
  
  endif
  
   if(npx == npx0 - 1) then
    do kc = 1,nz_local
     do jc = 1,Ny
      ic = nx_local
      df(ic,jc,kc) =- BC1(1) * f_loc(ic  ,jc,kc) &
                   &- BC1(2) * f_loc(ic-1,jc,kc) &
                   &- BC1(3) * f_loc(ic-2,jc,kc) &
                   &- BC1(4) * f_loc(ic-3,jc,kc) &
                   &- BC1(5) * f_loc(ic-4,jc,kc)

      ic = nx_local-1
      df(ic,jc,kc) =- BC2(1) * f_loc(ic+1,jc,kc) &
                   &- BC2(2) * f_loc(ic  ,jc,kc) &
                   &- BC2(3) * f_loc(ic-1,jc,kc) &
                   &- BC2(4) * f_loc(ic-2,jc,kc) &
                   &- BC2(5) * f_loc(ic-3,jc,kc)

      ic = nx_local-2                                
      df(ic,jc,kc) = -BC3(1) * f_loc(ic+2,jc,kc) &
                   &- BC3(2) * f_loc(ic+1,jc,kc) &
                   &- BC3(3) * f_loc(ic  ,jc,kc) &
                   &- BC3(4) * f_loc(ic-1,jc,kc) &
                   &- BC3(5) * f_loc(ic-2,jc,kc)
     enddo 
    enddo
   endif
  
  
end subroutine Cal_Deri_Dxi_5th_ce

!subroutine Cal_Deri_Dxi_5th_do(df,f_loc)
!   implicit none
!   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   
!   integer( kind = ik )::ic,jc,kc 
!
!   df = 0.d0;
!   ! To ensure the memory consist 
!   do kc = 1, nz_local
!     do jc = 1, Ny
!       do ic = 1, nx_local
!        df(ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(ic - 3,jc,kc) &
!                    &+ CO_FD_Dowind(-2) * f_loc(ic - 2,jc,kc) &
!                    &+ CO_FD_Dowind(-1) * f_loc(ic - 1,jc,kc) &
!                    &+ CO_FD_Dowind( 0) * f_loc(ic    ,jc,kc) &
!                    &+ CO_FD_Dowind( 1) * f_loc(ic + 1,jc,kc) &
!                    &+ CO_FD_Dowind( 2) * f_loc(ic + 2,jc,kc) &
!                    &+ CO_FD_Dowind( 3) * f_loc(ic + 3,jc,kc)
!       enddo
!     enddo
!   enddo
!
!   
!   if( ModelType == 1)then  
!   if(npx == 0) then
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      ic = 1  
!      df(ic,jc,kc) = BC1(1) * f_loc(ic  ,jc,kc) &
!                  &+ BC1(2) * f_loc(ic+1,jc,kc) &
!                  &+ BC1(3) * f_loc(ic+2,jc,kc) &
!                  &+ BC1(4) * f_loc(ic+3,jc,kc) &
!                  &+ BC1(5) * f_loc(ic+4,jc,kc)
!      
!      ic = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic-1,jc,kc) &
!                  &+ BC2(2) * f_loc(ic  ,jc,kc) &
!                  &+ BC2(3) * f_loc(ic+1,jc,kc) &
!                  &+ BC2(4) * f_loc(ic+2,jc,kc) &
!                  &+ BC2(5) * f_loc(ic+3,jc,kc)
!
!
!      ic = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic-2,jc,kc) &
!                  &+ BC3(2) * f_loc(ic-1,jc,kc) &
!                  &+ BC3(3) * f_loc(ic  ,jc,kc) &
!                  &+ BC3(4) * f_loc(ic+1,jc,kc) &
!                  &+ BC3(5) * f_loc(ic+2,jc,kc)
!     enddo 
!    enddo
!   endif
!   endif
!   
!   if(npx == npx0 - 1) then
!    do kc = 1,nz_local
!     do jc = 1,Ny
!      ic = nx_local
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic  ,jc,kc) &
!                   &- BC1(2) * f_loc(ic-1,jc,kc) &
!                   &- BC1(3) * f_loc(ic-2,jc,kc) &
!                   &- BC1(4) * f_loc(ic-3,jc,kc) &
!                   &- BC1(5) * f_loc(ic-4,jc,kc)
!
!      ic = nx_local-1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic+1,jc,kc) &
!                   &- BC2(2) * f_loc(ic  ,jc,kc) &
!                   &- BC2(3) * f_loc(ic-1,jc,kc) &
!                   &- BC2(4) * f_loc(ic-2,jc,kc) &
!                   &- BC2(5) * f_loc(ic-3,jc,kc)
!
!      ic = nx_local-2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic+2,jc,kc) &
!                   &- BC3(2) * f_loc(ic+1,jc,kc) &
!                   &- BC3(3) * f_loc(ic  ,jc,kc) &
!                   &- BC3(4) * f_loc(ic-1,jc,kc) &
!                   &- BC3(5) * f_loc(ic-2,jc,kc)
!     enddo 
!    enddo
!   endif
!   
!   
!end subroutine Cal_Deri_Dxi_5th_do

! Eta direction normal derivatives f(:,:,:,:)
subroutine Cal_Deri_Deta_5th_up_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny - 3
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(iVar,ic,jc - 3,kc) &
                         &+ CO_FD_Upwind(-2) * f_loc(iVar,ic,jc - 2,kc) &
                         &+ CO_FD_Upwind(-1) * f_loc(iVar,ic,jc - 1,kc) &
                         &+ CO_FD_Upwind( 0) * f_loc(iVar,ic,jc    ,kc) &
                         &+ CO_FD_Upwind( 1) * f_loc(iVar,ic,jc + 1,kc) &
                         &+ CO_FD_Upwind( 2) * f_loc(iVar,ic,jc + 2,kc) &
                         &+ CO_FD_Upwind( 3) * f_loc(iVar,ic,jc + 3,kc)
        enddo
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC1(2) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC1(3) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC1(4) * f_loc(iVar,ic,jc+3,kc) &
                       &+ BC1(5) * f_loc(iVar,ic,jc+4,kc)
      enddo 
      jc = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC2(3) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC2(4) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC2(5) * f_loc(iVar,ic,jc+3,kc)
      enddo

      jc = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic,jc-2,kc) &
                       &+ BC3(2) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC3(4) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC3(5) * f_loc(iVar,ic,jc+2,kc)
      enddo
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC1(2) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC1(3) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC1(4) * f_loc(iVar,ic,jc-3,kc) &
                        &- BC1(5) * f_loc(iVar,ic,jc-4,kc)
      enddo

      jc = Ny - 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC2(3) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC2(4) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC2(5) * f_loc(iVar,ic,jc-3,kc)
      enddo

      jc = Ny - 2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic,jc+2,kc) &
                        &- BC3(2) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC3(4) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC3(5) * f_loc(iVar,ic,jc-2,kc)
      enddo
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_up_numvar

subroutine Cal_Deri_Deta_5th_ce_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny-3
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Center(-3) * f_loc(iVar,ic,jc - 3,kc) &
                         &+ CO_FD_Center(-2) * f_loc(iVar,ic,jc - 2,kc) &
                         &+ CO_FD_Center(-1) * f_loc(iVar,ic,jc - 1,kc) &
                         &+ CO_FD_Center( 0) * f_loc(iVar,ic,jc    ,kc) &
                         &+ CO_FD_Center( 1) * f_loc(iVar,ic,jc + 1,kc) &
                         &+ CO_FD_Center( 2) * f_loc(iVar,ic,jc + 2,kc) &
                         &+ CO_FD_Center( 3) * f_loc(iVar,ic,jc + 3,kc)
        enddo
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC1(2) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC1(3) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC1(4) * f_loc(iVar,ic,jc+3,kc) &
                       &+ BC1(5) * f_loc(iVar,ic,jc+4,kc)
      enddo 
      jc = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC2(3) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC2(4) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC2(5) * f_loc(iVar,ic,jc+3,kc)
      enddo

      jc = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic,jc-2,kc) &
                       &+ BC3(2) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC3(4) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC3(5) * f_loc(iVar,ic,jc+2,kc)
      enddo
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC1(2) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC1(3) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC1(4) * f_loc(iVar,ic,jc-3,kc) &
                        &- BC1(5) * f_loc(iVar,ic,jc-4,kc)
      enddo

      jc = Ny - 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC2(3) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC2(4) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC2(5) * f_loc(iVar,ic,jc-3,kc)
      enddo

      jc = Ny - 2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic,jc+2,kc) &
                        &- BC3(2) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC3(4) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC3(5) * f_loc(iVar,ic,jc-2,kc)
      enddo
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_ce_numvar

subroutine Cal_Deri_Deta_5th_do_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny-3
       do ic = 1, nx_local
        do iVar = 1, NumVar
        df(iVar,ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(iVar,ic,jc - 3,kc) &
                         &+ CO_FD_Dowind(-2) * f_loc(iVar,ic,jc - 2,kc) &
                         &+ CO_FD_Dowind(-1) * f_loc(iVar,ic,jc - 1,kc) &
                         &+ CO_FD_Dowind( 0) * f_loc(iVar,ic,jc    ,kc) &
                         &+ CO_FD_Dowind( 1) * f_loc(iVar,ic,jc + 1,kc) &
                         &+ CO_FD_Dowind( 2) * f_loc(iVar,ic,jc + 2,kc) &
                         &+ CO_FD_Dowind( 3) * f_loc(iVar,ic,jc + 3,kc)
        enddo
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC1(2) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC1(3) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC1(4) * f_loc(iVar,ic,jc+3,kc) &
                       &+ BC1(5) * f_loc(iVar,ic,jc+4,kc)
      enddo 
      jc = 2
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) = BC2(1) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC2(3) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC2(4) * f_loc(iVar,ic,jc+2,kc) &
                       &+ BC2(5) * f_loc(iVar,ic,jc+3,kc)
      enddo

      jc = 3
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = BC3(1) * f_loc(iVar,ic,jc-2,kc) &
                       &+ BC3(2) * f_loc(iVar,ic,jc-1,kc) &
                       &+ BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                       &+ BC3(4) * f_loc(iVar,ic,jc+1,kc) &
                       &+ BC3(5) * f_loc(iVar,ic,jc+2,kc)
      enddo
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC1(1) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC1(2) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC1(3) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC1(4) * f_loc(iVar,ic,jc-3,kc) &
                        &- BC1(5) * f_loc(iVar,ic,jc-4,kc)
      enddo

      jc = Ny - 1
      do iVar = 1, NumVar
      df(iVar,ic,jc,kc) =- BC2(1) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC2(2) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC2(3) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC2(4) * f_loc(iVar,ic,jc-2,kc) &
                        &- BC2(5) * f_loc(iVar,ic,jc-3,kc)
      enddo

      jc = Ny - 2
      do iVar = 1, NumVar                                
      df(iVar,ic,jc,kc) = -BC3(1) * f_loc(iVar,ic,jc+2,kc) &
                        &- BC3(2) * f_loc(iVar,ic,jc+1,kc) &
                        &- BC3(3) * f_loc(iVar,ic,jc  ,kc) &
                        &- BC3(4) * f_loc(iVar,ic,jc-1,kc) &
                        &- BC3(5) * f_loc(iVar,ic,jc-2,kc)
      enddo
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_do_numvar

! Eta direction normal derivatives f(:,:,:)
subroutine Cal_Deri_Deta_5th_up(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny-3
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(ic,jc - 3,kc) &
                    &+ CO_FD_Upwind(-2) * f_loc(ic,jc - 2,kc) &
                    &+ CO_FD_Upwind(-1) * f_loc(ic,jc - 1,kc) &
                    &+ CO_FD_Upwind( 0) * f_loc(ic,jc    ,kc) &
                    &+ CO_FD_Upwind( 1) * f_loc(ic,jc + 1,kc) &
                    &+ CO_FD_Upwind( 2) * f_loc(ic,jc + 2,kc) &
                    &+ CO_FD_Upwind( 3) * f_loc(ic,jc + 3,kc)
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc  ,kc) &
                  &+ BC1(2) * f_loc(ic,jc+1,kc) &
                  &+ BC1(3) * f_loc(ic,jc+2,kc) &
                  &+ BC1(4) * f_loc(ic,jc+3,kc) &
                  &+ BC1(5) * f_loc(ic,jc+4,kc)
                  
      jc = 2
      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc-1,kc) &
                  &+ BC2(2) * f_loc(ic,jc  ,kc) &
                  &+ BC2(3) * f_loc(ic,jc+1,kc) &
                  &+ BC2(4) * f_loc(ic,jc+2,kc) &
                  &+ BC2(5) * f_loc(ic,jc+3,kc)

      jc = 3                                
      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc-2,kc) &
                  &+ BC3(2) * f_loc(ic,jc-1,kc) &
                  &+ BC3(3) * f_loc(ic,jc  ,kc) &
                  &+ BC3(4) * f_loc(ic,jc+1,kc) &
                  &+ BC3(5) * f_loc(ic,jc+2,kc)
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc  ,kc) &
                   &- BC1(2) * f_loc(ic,jc-1,kc) &
                   &- BC1(3) * f_loc(ic,jc-2,kc) &
                   &- BC1(4) * f_loc(ic,jc-3,kc) &
                   &- BC1(5) * f_loc(ic,jc-4,kc)

      jc = Ny - 1
      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc+1,kc) &
                   &- BC2(2) * f_loc(ic,jc  ,kc) &
                   &- BC2(3) * f_loc(ic,jc-1,kc) &
                   &- BC2(4) * f_loc(ic,jc-2,kc) &
                   &- BC2(5) * f_loc(ic,jc-3,kc)

      jc = Ny - 2                                
      df(ic,jc,kc) =- BC3(1) * f_loc(ic,jc+2,kc) &
                   &- BC3(2) * f_loc(ic,jc+1,kc) &
                   &- BC3(3) * f_loc(ic,jc  ,kc) &
                   &- BC3(4) * f_loc(ic,jc-1,kc) &
                   &- BC3(5) * f_loc(ic,jc-2,kc)
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_up

subroutine Cal_Deri_Deta_5th_ce(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny-3
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Center(-3) * f_loc(ic,jc - 3,kc) &
                    &+ CO_FD_Center(-2) * f_loc(ic,jc - 2,kc) &
                    &+ CO_FD_Center(-1) * f_loc(ic,jc - 1,kc) &
                    &+ CO_FD_Center( 0) * f_loc(ic,jc    ,kc) &
                    &+ CO_FD_Center( 1) * f_loc(ic,jc + 1,kc) &
                    &+ CO_FD_Center( 2) * f_loc(ic,jc + 2,kc) &
                    &+ CO_FD_Center( 3) * f_loc(ic,jc + 3,kc)
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc  ,kc) &
                  &+ BC1(2) * f_loc(ic,jc+1,kc) &
                  &+ BC1(3) * f_loc(ic,jc+2,kc) &
                  &+ BC1(4) * f_loc(ic,jc+3,kc) &
                  &+ BC1(5) * f_loc(ic,jc+4,kc)
                  
      jc = 2
      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc-1,kc) &
                  &+ BC2(2) * f_loc(ic,jc  ,kc) &
                  &+ BC2(3) * f_loc(ic,jc+1,kc) &
                  &+ BC2(4) * f_loc(ic,jc+2,kc) &
                  &+ BC2(5) * f_loc(ic,jc+3,kc)

      jc = 3                                
      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc-2,kc) &
                  &+ BC3(2) * f_loc(ic,jc-1,kc) &
                  &+ BC3(3) * f_loc(ic,jc  ,kc) &
                  &+ BC3(4) * f_loc(ic,jc+1,kc) &
                  &+ BC3(5) * f_loc(ic,jc+2,kc)
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc  ,kc) &
                   &- BC1(2) * f_loc(ic,jc-1,kc) &
                   &- BC1(3) * f_loc(ic,jc-2,kc) &
                   &- BC1(4) * f_loc(ic,jc-3,kc) &
                   &- BC1(5) * f_loc(ic,jc-4,kc)

      jc = Ny - 1
      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc+1,kc) &
                   &- BC2(2) * f_loc(ic,jc  ,kc) &
                   &- BC2(3) * f_loc(ic,jc-1,kc) &
                   &- BC2(4) * f_loc(ic,jc-2,kc) &
                   &- BC2(5) * f_loc(ic,jc-3,kc)

      jc = Ny - 2                                
      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc+2,kc) &
                   &- BC3(2) * f_loc(ic,jc+1,kc) &
                   &- BC3(3) * f_loc(ic,jc  ,kc) &
                   &- BC3(4) * f_loc(ic,jc-1,kc) &
                   &- BC3(5) * f_loc(ic,jc-2,kc)
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_ce

subroutine Cal_Deri_Deta_5th_do(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 4, Ny-3
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(ic,jc - 3,kc) &
                    &+ CO_FD_Dowind(-2) * f_loc(ic,jc - 2,kc) &
                    &+ CO_FD_Dowind(-1) * f_loc(ic,jc - 1,kc) &
                    &+ CO_FD_Dowind( 0) * f_loc(ic,jc    ,kc) &
                    &+ CO_FD_Dowind( 1) * f_loc(ic,jc + 1,kc) &
                    &+ CO_FD_Dowind( 2) * f_loc(ic,jc + 2,kc) &
                    &+ CO_FD_Dowind( 3) * f_loc(ic,jc + 3,kc)
       enddo
     enddo
   enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = 1
      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc  ,kc) &
                  &+ BC1(2) * f_loc(ic,jc+1,kc) &
                  &+ BC1(3) * f_loc(ic,jc+2,kc) &
                  &+ BC1(4) * f_loc(ic,jc+3,kc) &
                  &+ BC1(5) * f_loc(ic,jc+4,kc)
                  
      jc = 2
      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc-1,kc) &
                  &+ BC2(2) * f_loc(ic,jc  ,kc) &
                  &+ BC2(3) * f_loc(ic,jc+1,kc) &
                  &+ BC2(4) * f_loc(ic,jc+2,kc) &
                  &+ BC2(5) * f_loc(ic,jc+3,kc)

      jc = 3                                
      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc-2,kc) &
                  &+ BC3(2) * f_loc(ic,jc-1,kc) &
                  &+ BC3(3) * f_loc(ic,jc  ,kc) &
                  &+ BC3(4) * f_loc(ic,jc+1,kc) &
                  &+ BC3(5) * f_loc(ic,jc+2,kc)
     enddo 
    enddo

    do kc = 1,nz_local
     do ic = 1,nx_local
      jc = Ny
      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc  ,kc) &
                   &- BC1(2) * f_loc(ic,jc-1,kc) &
                   &- BC1(3) * f_loc(ic,jc-2,kc) &
                   &- BC1(4) * f_loc(ic,jc-3,kc) &
                   &- BC1(5) * f_loc(ic,jc-4,kc)

      jc = Ny - 1
      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc+1,kc) &
                   &- BC2(2) * f_loc(ic,jc  ,kc) &
                   &- BC2(3) * f_loc(ic,jc-1,kc) &
                   &- BC2(4) * f_loc(ic,jc-2,kc) &
                   &- BC2(5) * f_loc(ic,jc-3,kc)

      jc = Ny - 2                                
      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc+2,kc) &
                   &- BC3(2) * f_loc(ic,jc+1,kc) &
                   &- BC3(3) * f_loc(ic,jc  ,kc) &
                   &- BC3(4) * f_loc(ic,jc-1,kc) &
                   &- BC3(5) * f_loc(ic,jc-2,kc)
     enddo 
    enddo
end subroutine Cal_Deri_Deta_5th_do

! Zeta direction normal derivatives f(:,:,:)
!subroutine Cal_Deri_Dzeta_5th_up(df,f_loc)
!   implicit none
!   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   
!   integer( kind = ik )::ic,jc,kc 
!
!   df = 0.d0;
!   ! To ensure the memory consist 
!   do kc = 1, nz_local
!     do jc = 1, Ny
!       do ic = 1, nx_local
!        df(ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(ic,jc,kc - 3) &
!                    &+ CO_FD_Upwind(-2) * f_loc(ic,jc,kc - 2) &
!                    &+ CO_FD_Upwind(-1) * f_loc(ic,jc,kc - 1) &
!                    &+ CO_FD_Upwind( 0) * f_loc(ic,jc,kc    ) &
!                    &+ CO_FD_Upwind( 1) * f_loc(ic,jc,kc + 1) &
!                    &+ CO_FD_Upwind( 2) * f_loc(ic,jc,kc + 2) &
!                    &+ CO_FD_Upwind( 3) * f_loc(ic,jc,kc + 3)
!       enddo
!     enddo
!   enddo
!
!   if(npz == 0) then
!    do jc = 1,Ny
!     do ic = 1,nx_local
!      kc = 1
!      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc,kc  ) &
!                  &+ BC1(2) * f_loc(ic,jc,kc+1) &
!                  &+ BC1(3) * f_loc(ic,jc,kc+2) &
!                  &+ BC1(4) * f_loc(ic,jc,kc+3) &
!                  &+ BC1(5) * f_loc(ic,jc,kc+4)
!                  
!      kc = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc,kc-1) &
!                  &+ BC2(2) * f_loc(ic,jc,kc  ) &
!                  &+ BC2(3) * f_loc(ic,jc,kc+1) &
!                  &+ BC2(4) * f_loc(ic,jc,kc+2) &
!                  &+ BC2(5) * f_loc(ic,jc,kc+3)
!
!      kc = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc,kc-2) &
!                  &+ BC3(2) * f_loc(ic,jc,kc-1) &
!                  &+ BC3(3) * f_loc(ic,jc,kc  ) &
!                  &+ BC3(4) * f_loc(ic,jc,kc+1) &
!                  &+ BC3(5) * f_loc(ic,jc,kc+2)
!     enddo 
!    enddo
!   endif
!
!   if(npz == npz0 - 1) then
!    do jc = 1,Ny
!     do ic = 1,nx_local
!      kc = nz_local
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc,kc  ) &
!                   &- BC1(2) * f_loc(ic,jc,kc-1) &
!                   &- BC1(3) * f_loc(ic,jc,kc-2) &
!                   &- BC1(4) * f_loc(ic,jc,kc-3) &
!                   &- BC1(5) * f_loc(ic,jc,kc-4)
!
!      kc = nz_local-1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc,kc+1) &
!                   &- BC2(2) * f_loc(ic,jc,kc  ) &
!                   &- BC2(3) * f_loc(ic,jc,kc-1) &
!                   &- BC2(4) * f_loc(ic,jc,kc-2) &
!                   &- BC2(5) * f_loc(ic,jc,kc-3)
!
!      kc = nz_local-2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc,kc+2) &
!                   &- BC3(2) * f_loc(ic,jc,kc+1) &
!                   &- BC3(3) * f_loc(ic,jc,kc  ) &
!                   &- BC3(4) * f_loc(ic,jc,kc-1) &
!                   &- BC3(5) * f_loc(ic,jc,kc-2)
!     enddo 
!    enddo
!   endif
!end subroutine Cal_Deri_Dzeta_5th_up

subroutine Cal_Deri_Dzeta_5th_ce(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Center(-3) * f_loc(ic,jc,kc - 3) &
                    &+ CO_FD_Center(-2) * f_loc(ic,jc,kc - 2) &
                    &+ CO_FD_Center(-1) * f_loc(ic,jc,kc - 1) &
                    &+ CO_FD_Center( 0) * f_loc(ic,jc,kc    ) &
                    &+ CO_FD_Center( 1) * f_loc(ic,jc,kc + 1) &
                    &+ CO_FD_Center( 2) * f_loc(ic,jc,kc + 2) &
                    &+ CO_FD_Center( 3) * f_loc(ic,jc,kc + 3)
       enddo
     enddo
   enddo

   if(npz == 0) then
    do jc = 1,Ny
     do ic = 1,nx_local
      kc = 1
      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc,kc  ) &
                  &+ BC1(2) * f_loc(ic,jc,kc+1) &
                  &+ BC1(3) * f_loc(ic,jc,kc+2) &
                  &+ BC1(4) * f_loc(ic,jc,kc+3) &
                  &+ BC1(5) * f_loc(ic,jc,kc+4)
                  
      kc = 2
      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc,kc-1) &
                  &+ BC2(2) * f_loc(ic,jc,kc  ) &
                  &+ BC2(3) * f_loc(ic,jc,kc+1) &
                  &+ BC2(4) * f_loc(ic,jc,kc+2) &
                  &+ BC2(5) * f_loc(ic,jc,kc+3)

      kc = 3                                
      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc,kc-2) &
                  &+ BC3(2) * f_loc(ic,jc,kc-1) &
                  &+ BC3(3) * f_loc(ic,jc,kc  ) &
                  &+ BC3(4) * f_loc(ic,jc,kc+1) &
                  &+ BC3(5) * f_loc(ic,jc,kc+2)
     enddo 
    enddo
   endif

   if(npz == npz0 - 1) then
    do jc = 1,Ny
     do ic = 1,nx_local
      kc = nz_local
      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc,kc  ) &
                   &- BC1(2) * f_loc(ic,jc,kc-1) &
                   &- BC1(3) * f_loc(ic,jc,kc-2) &
                   &- BC1(4) * f_loc(ic,jc,kc-3) &
                   &- BC1(5) * f_loc(ic,jc,kc-4)

      kc = nz_local-1
      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc,kc+1) &
                   &- BC2(2) * f_loc(ic,jc,kc  ) &
                   &- BC2(3) * f_loc(ic,jc,kc-1) &
                   &- BC2(4) * f_loc(ic,jc,kc-2) &
                   &- BC2(5) * f_loc(ic,jc,kc-3)

      kc = nz_local-2                                
      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc,kc+2) &
                   &- BC3(2) * f_loc(ic,jc,kc+1) &
                   &- BC3(3) * f_loc(ic,jc,kc  ) &
                   &- BC3(4) * f_loc(ic,jc,kc-1) &
                   &- BC3(5) * f_loc(ic,jc,kc-2)
     enddo 
    enddo
   endif
end subroutine Cal_Deri_Dzeta_5th_ce

!subroutine Cal_Deri_Dzeta_5th_do(df,f_loc)
!   implicit none
!   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
!   
!   integer( kind = ik )::ic,jc,kc 
!
!   df = 0.d0;
!   ! To ensure the memory consist 
!   do kc = 1, nz_local
!     do jc = 1, Ny
!       do ic = 1, nx_local
!        df(ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(ic,jc,kc - 3) &
!                    &+ CO_FD_Dowind(-2) * f_loc(ic,jc,kc - 2) &
!                    &+ CO_FD_Dowind(-1) * f_loc(ic,jc,kc - 1) &
!                    &+ CO_FD_Dowind( 0) * f_loc(ic,jc,kc    ) &
!                    &+ CO_FD_Dowind( 1) * f_loc(ic,jc,kc + 1) &
!                    &+ CO_FD_Dowind( 2) * f_loc(ic,jc,kc + 2) &
!                    &+ CO_FD_Dowind( 3) * f_loc(ic,jc,kc + 3)
!       enddo
!     enddo
!   enddo
!
!   if(npz == 0) then
!    do jc = 1,Ny
!     do ic = 1,nx_local
!      kc = 1
!      df(ic,jc,kc) = BC1(1) * f_loc(ic,jc,kc  ) &
!                  &+ BC1(2) * f_loc(ic,jc,kc+1) &
!                  &+ BC1(3) * f_loc(ic,jc,kc+2) &
!                  &+ BC1(4) * f_loc(ic,jc,kc+3) &
!                  &+ BC1(5) * f_loc(ic,jc,kc+4)
!                  
!      kc = 2
!      df(ic,jc,kc) = BC2(1) * f_loc(ic,jc,kc-1) &
!                  &+ BC2(2) * f_loc(ic,jc,kc  ) &
!                  &+ BC2(3) * f_loc(ic,jc,kc+1) &
!                  &+ BC2(4) * f_loc(ic,jc,kc+2) &
!                  &+ BC2(5) * f_loc(ic,jc,kc+3)
!
!      kc = 3                                
!      df(ic,jc,kc) = BC3(1) * f_loc(ic,jc,kc-2) &
!                  &+ BC3(2) * f_loc(ic,jc,kc-1) &
!                  &+ BC3(3) * f_loc(ic,jc,kc  ) &
!                  &+ BC3(4) * f_loc(ic,jc,kc+1) &
!                  &+ BC3(5) * f_loc(ic,jc,kc+2)
!     enddo 
!    enddo
!   endif
!
!   if(npz == npz0 - 1) then
!    do jc = 1,Ny
!     do ic = 1,nx_local
!      kc = nz_local
!      df(ic,jc,kc) =- BC1(1) * f_loc(ic,jc,kc  ) &
!                   &- BC1(2) * f_loc(ic,jc,kc-1) &
!                   &- BC1(3) * f_loc(ic,jc,kc-2) &
!                   &- BC1(4) * f_loc(ic,jc,kc-3) &
!                   &- BC1(5) * f_loc(ic,jc,kc-4)
!
!      kc = nz_local-1
!      df(ic,jc,kc) =- BC2(1) * f_loc(ic,jc,kc+1) &
!                   &- BC2(2) * f_loc(ic,jc,kc  ) &
!                   &- BC2(3) * f_loc(ic,jc,kc-1) &
!                   &- BC2(4) * f_loc(ic,jc,kc-2) &
!                   &- BC2(5) * f_loc(ic,jc,kc-3)
!
!      kc = nz_local-2                                
!      df(ic,jc,kc) = -BC3(1) * f_loc(ic,jc,kc+2) &
!                   &- BC3(2) * f_loc(ic,jc,kc+1) &
!                   &- BC3(3) * f_loc(ic,jc,kc  ) &
!                   &- BC3(4) * f_loc(ic,jc,kc-1) &
!                   &- BC3(5) * f_loc(ic,jc,kc-2)
!     enddo 
!    enddo
!   endif
!end subroutine Cal_Deri_Dzeta_5th_do

! Zeta direction periodic derivatives f(:,:,:,:)
subroutine Cal_Deri_Dzeta_5th_up_per_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1,NumVar
        df(iVar,ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(iVar,ic,jc,kc - 3) &
                         &+ CO_FD_Upwind(-2) * f_loc(iVar,ic,jc,kc - 2) &
                         &+ CO_FD_Upwind(-1) * f_loc(iVar,ic,jc,kc - 1) &
                         &+ CO_FD_Upwind( 0) * f_loc(iVar,ic,jc,kc    ) &
                         &+ CO_FD_Upwind( 1) * f_loc(iVar,ic,jc,kc + 1) &
                         &+ CO_FD_Upwind( 2) * f_loc(iVar,ic,jc,kc + 2) &
                         &+ CO_FD_Upwind( 3) * f_loc(iVar,ic,jc,kc + 3)
        enddo
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_up_per_numvar

subroutine Cal_Deri_Dzeta_5th_ce_per_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1,NumVar
        df(iVar,ic,jc,kc) = CO_FD_Center(-3) * f_loc(iVar,ic,jc,kc - 3) &
                         &+ CO_FD_Center(-2) * f_loc(iVar,ic,jc,kc - 2) &
                         &+ CO_FD_Center(-1) * f_loc(iVar,ic,jc,kc - 1) &
                         &+ CO_FD_Center( 0) * f_loc(iVar,ic,jc,kc    ) &
                         &+ CO_FD_Center( 1) * f_loc(iVar,ic,jc,kc + 1) &
                         &+ CO_FD_Center( 2) * f_loc(iVar,ic,jc,kc + 2) &
                         &+ CO_FD_Center( 3) * f_loc(iVar,ic,jc,kc + 3)
        enddo
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_ce_per_numvar

!subroutine cal_deri_dzeta_5th_ce_numvar(df,f_loc)
!   implicit none
!   real( kind = rk ), intent(in)::f_loc(1:numvar,1-overlap:nx_local+overlap,1:ny,1-overlap:nz_local+overlap)
!   real( kind = rk ), intent(out)::  df(1:numvar,1-overlap:nx_local+overlap,1:ny,1-overlap:nz_local+overlap)
!   
!   integer( kind = ik )::ic,jc,kc,ivar 
!
!   df = 0.d0;
!   ! to ensure the memory consist 
!   do kc = 1, nz_local
!     do jc = 1, ny
!       do ic = 1, nx_local
!        do ivar = 1,numvar
!        df(ivar,ic,jc,kc) = co_fd_center(-3) * f_loc(ivar,ic,jc,kc - 3) &
!                         &+ co_fd_center(-2) * f_loc(ivar,ic,jc,kc - 2) &
!                         &+ co_fd_center(-1) * f_loc(ivar,ic,jc,kc - 1) &
!                         &+ co_fd_center( 0) * f_loc(ivar,ic,jc,kc    ) &
!                         &+ co_fd_center( 1) * f_loc(ivar,ic,jc,kc + 1) &
!                         &+ co_fd_center( 2) * f_loc(ivar,ic,jc,kc + 2) &
!                         &+ co_fd_center( 3) * f_loc(ivar,ic,jc,kc + 3)
!        enddo
!       enddo
!     enddo
!   enddo
!   
!    if(npz == 0) then
!    do jc = 1,ny
!     do ic = 1,nx_local
!      kc = 1
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) = bc1(1) * f_loc(ivar,ic,jc,kc  ) &
!                       &+ bc1(2) * f_loc(ivar,ic,jc,kc+1) &
!                       &+ bc1(3) * f_loc(ivar,ic,jc,kc+2) &
!                       &+ bc1(4) * f_loc(ivar,ic,jc,kc+3) &
!                       &+ bc1(5) * f_loc(ivar,ic,jc,kc+4)
!      enddo  
!      
!      kc = 2
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) = bc2(1) * f_loc(ivar,ic,jc,kc-1) &
!                       &+ bc2(2) * f_loc(ivar,ic,jc,kc  ) &
!                       &+ bc2(3) * f_loc(ivar,ic,jc,kc+1) &
!                       &+ bc2(4) * f_loc(ivar,ic,jc,kc+2) &
!                       &+ bc2(5) * f_loc(ivar,ic,jc,kc+3)
!      enddo
!      kc = 3  
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) = bc3(1) * f_loc(ivar,ic,jc,kc-2) &
!                       &+ bc3(2) * f_loc(ivar,ic,jc,kc-1) &
!                       &+ bc3(3) * f_loc(ivar,ic,jc,kc  ) &
!                       &+ bc3(4) * f_loc(ivar,ic,jc,kc+1) &
!                       &+ bc3(5) * f_loc(ivar,ic,jc,kc+2)
!      enddo
!     enddo 
!    enddo
!   endif
!
!   if(npz == npz0 - 1) then
!    do jc = 1,ny
!     do ic = 1,nx_local
!      kc = nz_local
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) =- bc1(1) * f_loc(ivar,ic,jc,kc  ) &
!                        &- bc1(2) * f_loc(ivar,ic,jc,kc-1) &
!                        &- bc1(3) * f_loc(ivar,ic,jc,kc-2) &
!                        &- bc1(4) * f_loc(ivar,ic,jc,kc-3) &
!                        &- bc1(5) * f_loc(ivar,ic,jc,kc-4)
!      enddo
!      kc = nz_local-1
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) =- bc2(1) * f_loc(ivar,ic,jc,kc+1) &
!                        &- bc2(2) * f_loc(ivar,ic,jc,kc  ) &
!                        &- bc2(3) * f_loc(ivar,ic,jc,kc-1) &
!                        &- bc2(4) * f_loc(ivar,ic,jc,kc-2) &
!                        &- bc2(5) * f_loc(ivar,ic,jc,kc-3)
!      enddo
!      kc = nz_local-2 
!      do ivar = 1,numvar
!      df(ivar,ic,jc,kc) = -bc3(1) * f_loc(ivar,ic,jc,kc+2) &
!                        &- bc3(2) * f_loc(ivar,ic,jc,kc+1) &
!                        &- bc3(3) * f_loc(ivar,ic,jc,kc  ) &
!                        &- bc3(4) * f_loc(ivar,ic,jc,kc-1) &
!                        &- bc3(5) * f_loc(ivar,ic,jc,kc-2)
!      enddo
!     enddo 
!    enddo
!   endif
!   
!   
!   
!end subroutine cal_deri_dzeta_5th_ce_numvar


subroutine Cal_Deri_Dzeta_5th_do_per_numvar(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc,iVar 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        do iVar = 1,NumVar
        df(iVar,ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(iVar,ic,jc,kc - 3) &
                         &+ CO_FD_Dowind(-2) * f_loc(iVar,ic,jc,kc - 2) &
                         &+ CO_FD_Dowind(-1) * f_loc(iVar,ic,jc,kc - 1) &
                         &+ CO_FD_Dowind( 0) * f_loc(iVar,ic,jc,kc    ) &
                         &+ CO_FD_Dowind( 1) * f_loc(iVar,ic,jc,kc + 1) &
                         &+ CO_FD_Dowind( 2) * f_loc(iVar,ic,jc,kc + 2) &
                         &+ CO_FD_Dowind( 3) * f_loc(iVar,ic,jc,kc + 3)
        enddo
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_do_per_numvar

! Zeta direction periodic derivatives f(:,:,:)
subroutine Cal_Deri_Dzeta_5th_up_per(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Upwind(-3) * f_loc(ic,jc,kc - 3) &
                    &+ CO_FD_Upwind(-2) * f_loc(ic,jc,kc - 2) &
                    &+ CO_FD_Upwind(-1) * f_loc(ic,jc,kc - 1) &
                    &+ CO_FD_Upwind( 0) * f_loc(ic,jc,kc    ) &
                    &+ CO_FD_Upwind( 1) * f_loc(ic,jc,kc + 1) &
                    &+ CO_FD_Upwind( 2) * f_loc(ic,jc,kc + 2) &
                    &+ CO_FD_Upwind( 3) * f_loc(ic,jc,kc + 3)
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_up_per

subroutine Cal_Deri_Dzeta_5th_ce_per(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Center(-3) * f_loc(ic,jc,kc - 3) &
                    &+ CO_FD_Center(-2) * f_loc(ic,jc,kc - 2) &
                    &+ CO_FD_Center(-1) * f_loc(ic,jc,kc - 1) &
                    &+ CO_FD_Center( 0) * f_loc(ic,jc,kc    ) &
                    &+ CO_FD_Center( 1) * f_loc(ic,jc,kc + 1) &
                    &+ CO_FD_Center( 2) * f_loc(ic,jc,kc + 2) &
                    &+ CO_FD_Center( 3) * f_loc(ic,jc,kc + 3)
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_ce_per

subroutine Cal_Deri_Dzeta_5th_do_per(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,jc,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do jc = 1, Ny
       do ic = 1, nx_local
        df(ic,jc,kc) = CO_FD_Dowind(-3) * f_loc(ic,jc,kc - 3) &
                    &+ CO_FD_Dowind(-2) * f_loc(ic,jc,kc - 2) &
                    &+ CO_FD_Dowind(-1) * f_loc(ic,jc,kc - 1) &
                    &+ CO_FD_Dowind( 0) * f_loc(ic,jc,kc    ) &
                    &+ CO_FD_Dowind( 1) * f_loc(ic,jc,kc + 1) &
                    &+ CO_FD_Dowind( 2) * f_loc(ic,jc,kc + 2) &
                    &+ CO_FD_Dowind( 3) * f_loc(ic,jc,kc + 3)
       enddo
     enddo
   enddo
end subroutine Cal_Deri_Dzeta_5th_do_per

! 2D
! Xi direction normal derivatives f(:,:)
subroutine Cal_Deri_SurfDxi_5th_ce(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   do kc = 1, nz_local
     do ic = 1, nx_local
       df(ic,kc) = CO_FD_Center(-3) * f_loc(ic - 3,kc) &
                &+ CO_FD_Center(-2) * f_loc(ic - 2,kc) &
                &+ CO_FD_Center(-1) * f_loc(ic - 1,kc) &
                &+ CO_FD_Center( 0) * f_loc(ic    ,kc) &
                &+ CO_FD_Center( 1) * f_loc(ic + 1,kc) &
                &+ CO_FD_Center( 2) * f_loc(ic + 2,kc) &
                &+ CO_FD_Center( 3) * f_loc(ic + 3,kc)
     enddo
   enddo

   if( ModelType == 1)then
       
   if(npx == 0) then
    do kc = 1,nz_local
        
      ic = 1
      df(ic,kc) = BC1(1) * f_loc(ic  ,kc) &
               &+ BC1(2) * f_loc(ic+1,kc) &
               &+ BC1(3) * f_loc(ic+2,kc) &
               &+ BC1(4) * f_loc(ic+3,kc) &
               &+ BC1(5) * f_loc(ic+4,kc)
   
      ic = 2
      df(ic,kc) = BC2(1) * f_loc(ic-1,kc) &
               &+ BC2(2) * f_loc(ic  ,kc) &
               &+ BC2(3) * f_loc(ic+1,kc) &
               &+ BC2(4) * f_loc(ic+2,kc) &
               &+ BC2(5) * f_loc(ic+3,kc)
       
      ic = 3                                
      df(ic,kc) = BC3(1) * f_loc(ic-2,kc) &
               &+ BC3(2) * f_loc(ic-1,kc) &
               &+ BC3(3) * f_loc(ic  ,kc) &
               &+ BC3(4) * f_loc(ic+1,kc) &
               &+ BC3(5) * f_loc(ic+2,kc) 
    enddo
   endif
   endif
   !
   if(npx == npx0 - 1) then
    do kc = 1,nz_local
      ic = nx_local
      df(ic,kc) =- BC1(1) * f_loc(ic  ,kc) &
                &- BC1(2) * f_loc(ic-1,kc) &
                &- BC1(3) * f_loc(ic-2,kc) &
                &- BC1(4) * f_loc(ic-3,kc) &
                &- BC1(5) * f_loc(ic-4,kc)
      
      ic = nx_local-1
      df(ic,kc) =- BC2(1) * f_loc(ic+1,kc) &
                &- BC2(2) * f_loc(ic  ,kc) &
                &- BC2(3) * f_loc(ic-1,kc) &
                &- BC2(4) * f_loc(ic-2,kc) &
                &- BC2(5) * f_loc(ic-3,kc)
      
      ic = nx_local-2                                
      df(ic,kc) = -BC3(1) * f_loc(ic+2,kc) &
                &- BC3(2) * f_loc(ic+1,kc) &
                &- BC3(3) * f_loc(ic  ,kc) &
                &- BC3(4) * f_loc(ic-1,kc) &
                &- BC3(5) * f_loc(ic-2,kc)
    enddo
   endif
   
  
   
end subroutine Cal_Deri_SurfDxi_5th_ce

! Zeta direction normal derivatives f(:,:)
subroutine Cal_Deri_SurfDzeta_5th_ce(df,f_loc)
   implicit none
   real( kind = rk ), intent(in)::f_loc(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
   real( kind = rk ), intent(out)::  df(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
   
   integer( kind = ik )::ic,kc 

   df = 0.d0;
   ! To ensure the memory consist 
   
   do kc = 1, nz_local
       do ic = 1, nx_local
           
           df(ic,kc) = CO_FD_Center(-3) * f_loc(ic,kc - 3) &
                    &+ CO_FD_Center(-2) * f_loc(ic,kc - 2) &
                    &+ CO_FD_Center(-1) * f_loc(ic,kc - 1) &
                    &+ CO_FD_Center( 0) * f_loc(ic,kc    ) &
                    &+ CO_FD_Center( 1) * f_loc(ic,kc + 1) &
                    &+ CO_FD_Center( 2) * f_loc(ic,kc + 2) &
                    &+ CO_FD_Center( 3) * f_loc(ic,kc + 3)
           
          
       enddo
   enddo
   
if (ModelType == 1)then
   if(npz == 0) then
     do ic = 1,nx_local
      kc = 1
      df(ic,kc) = BC1(1) * f_loc(ic,kc  ) &
               &+ BC1(2) * f_loc(ic,kc+1) &
               &+ BC1(3) * f_loc(ic,kc+2) &
               &+ BC1(4) * f_loc(ic,kc+3) &
               &+ BC1(5) * f_loc(ic,kc+4)
                  
      kc = 2
      df(ic,kc) = BC2(1) * f_loc(ic,kc-1) &
               &+ BC2(2) * f_loc(ic,kc  ) &
               &+ BC2(3) * f_loc(ic,kc+1) &
               &+ BC2(4) * f_loc(ic,kc+2) &
               &+ BC2(5) * f_loc(ic,kc+3)

      kc = 3                                
      df(ic,kc) = BC3(1) * f_loc(ic,kc-2) &
               &+ BC3(2) * f_loc(ic,kc-1) &
               &+ BC3(3) * f_loc(ic,kc  ) &
               &+ BC3(4) * f_loc(ic,kc+1) &
               &+ BC3(5) * f_loc(ic,kc+2) 
    enddo
   endif

   if(npz == npz0 - 1) then
     do ic = 1,nx_local
      kc = nz_local
      df(ic,kc) =- BC1(1) * f_loc(ic,kc  ) &
                &- BC1(2) * f_loc(ic,kc-1) &
                &- BC1(3) * f_loc(ic,kc-2) &
                &- BC1(4) * f_loc(ic,kc-3) &
                &- BC1(5) * f_loc(ic,kc-4)

      kc = nz_local-1
      df(ic,kc) =- BC2(1) * f_loc(ic,kc+1) &
                &- BC2(2) * f_loc(ic,kc  ) &
                &- BC2(3) * f_loc(ic,kc-1) &
                &- BC2(4) * f_loc(ic,kc-2) &
                &- BC2(5) * f_loc(ic,kc-3)

      kc = nz_local-2                                
      df(ic,kc) = -BC3(1) * f_loc(ic,kc+2) &
                &- BC3(2) * f_loc(ic,kc+1) &
                &- BC3(3) * f_loc(ic,kc  ) &!
                &- BC3(4) * f_loc(ic,kc-1) &
                &- BC3(5) * f_loc(ic,kc-2)
     enddo 
   endif
   
endif

end subroutine Cal_Deri_SurfDzeta_5th_ce



END MODULE FD5_Order