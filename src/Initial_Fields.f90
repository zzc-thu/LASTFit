subroutine Initial_Fields
  ! Module dependencies
  use SF_Constant,   only: ik,rk,NumVar,FilesForContinue,pi,Unsteady_NS_Analysis
  use SF_CFD_Global, only: IF_Continue_Calculate,nx_local,nz_local,&
                           Ny,Sxi_Surface,TWall,Re_Ref,CSTHLND_REF,&
                           P_INF,T0_INF,Gamma,Mach_Ref,Mach_Normal
  use MPI_GLOBAL,    only: MyId,ierr,MPI_COMM_WORLD

  ! some arrays and vectors
  use SF_CFD_Global, only: U,V,W,P,T,Rho,Ucons0,&
                           dxidx,dxidy,dxidz,dxidt,&
                           detadx,detady,detadz,detadt,&
                           dzetadx,dzetady,dzetadz,dzetadt,&
                           BigU_Xi,BigU_Eta,BigU_Zeta,Cs,Mu,&
                           SHockH,ShockV,x_grid,y_grid,z_grid,&
                           Ucons0_steady,x_grid_steady,y_grid_steady,z_grid_steady,&
                           ShockH_steady,ShockV_steady,AnalysisType,Iteration

  ! some subroutines and functions
  use MPI_GLOBAL,     only: Parallel_Exchange,MPI_Barrier,&
                          & Parallel_Exchange_NumVar,Parallel_Exchange_Surface
  use OutputParaView, only: output_results
  use SFitting,       only: ShockRelation3D

  implicit none
  ! Local variables
  logical:: dir_exists
  character(len=100)::FILENAME_FLOW
  character(len=100)::FILENAME_SHK
  integer(kind=ik ):: ic,jc,kc 
  real( kind = rk ):: tempNym,MachNormalShock
  real( kind = rk ):: Stag_P,Stag_T,Tstar,Pwell
  real( kind = rk ):: surf_p(1:nx_local,1:nz_local)
  real( kind = rk ):: surf_T(1:nx_local,1:nz_local)
  real( kind = rk ):: DEG,SMS,DUM
  real( kind = rk ):: tempP,tempT,indexJ


 If(IF_Continue_Calculate == 1) then
    ! Here we read the initial fields from previous simulation
    ! in order to increase the efficient of the simulation
    ! all the processes read the sepearated files simultaneously
    if(MyId == 0) then
      write(*,'(A)')"We are reading the initial fields from previous simulation"
      inquire(file="./RESU/SF_Results00000.flowsfg", exist=dir_exists)
       if (.not. dir_exists) then
         write(*,'(A)')"The CONT file does not exist, you should start a new simulation"
          stop
       else
         write(*,'(A)')"The CONT file exists, we will read the initial fields from it"
       end if     
    endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
      
      ! Read the previous results
      write(filename_flow,'(A,A,I5.5,A)')'RESU/',trim(FilesForContinue),MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow"
      open(112,file=filename_flow,form='unformatted',status='old')
       read(112)Ucons0
       read(112)x_grid
       read(112)y_grid
       read(112)z_grid
      close(112)
      ! Transfer the conservative variables to the basic variables
      do kc = 1, nz_local
       do jc = 1,Ny
        do ic = 1, nx_local
         call CVtoU_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
         ! Calculate the temperature
         T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc)
        enddo
       enddo
      enddo

      ! Read the previous shock informations
      write(filename_shk,'(A,A,I5.5,A)')'RESU/',trim(FilesForContinue),MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Reading ShockH and ShockV"
     open(113,file=filename_shk,form='unformatted',status='old')
      read(113)ShockH
      read(113)ShockV
     close(113)

      ! Update the Jacobian
      call Calculate_Jaco
      
    if(MyId == 0) then
      write(*,'(A)')"We have read the initial fields from previous simulation"
    end if


 else
  if(MyId == 0) then
    write(*,'(A)')"We are starting a new simulation"
  endif

  ! Get the Values for jc = Ny
   call ShockRelation3D

!#ifdef DEBUG
  !if(MyId == 0) then
  ! write(*,*)"Output the shock surface informations for Check"
  !endif
  ! !Check Points: 
  ! open(123,file="CheckFiles/ShockVariables.dat",status="replace")
  ! write(123,*)"VARIABLES=ic Rho, U, V, W, P, T"
  ! do ic = 1,nx_local
  !  do kc = 1,nz_local 
  !    write(123,*)myid,ic,kc,Rho(ic,Ny,kc),U(ic,Ny,kc),V(ic,Ny,kc),&
  !                           W(ic,Ny,kc),P(ic,Ny,kc),T(ic,Ny,kc)
  !  enddo
  ! enddo 
  ! close(123) 

!#endif

  ! Here we have define two emperical parameters
  Tstar = pi * 50.d0 / 180.d0;
  Pwell = 2.22d0;
  ! Define the normal Mach numbers immediately
  ! after the shock at the stag points
  MachNormalShock  = sqrt((2.d0+(Gamma - 1.d0) * Mach_Normal**2.d0 )/ &
                          (2.d0*Gamma*Mach_Normal**2.d0-Gamma+1.d0))

  ! First, we calculate the pressure after the shock
  Stag_P = P_INF * (2.d0 * Gamma * Mach_Ref * Mach_Ref / (Gamma + 1.d0) - (Gamma - 1.d0)/(Gamma + 1.d0)); 
  ! Second, we calculate the pressure at the wall surfaces
  Stag_P = Stag_P * (1.D0+0.5D0*(Gamma-1.D0)*MachNormalShock*MachNormalShock)**(Gamma/(Gamma - 1.D0))
  
  Stag_T = T0_INF;
  
  do kc = 1,nz_local
   do ic = 1,nx_local
     DEG = MIN(abs(Sxi_Surface(ic,kc)),2.0_rk);
     SMS = (DEG / Tstar)**Pwell;
     DUM = 1.d0 + 0.5d0 * (Gamma - 1.d0) * SMS;
     if(abs(DUM)<1e-11) then
       surf_p(ic,kc) = Stag_P;
       surf_T(ic,kc) = Stag_T;
     else
       surf_p(ic,kc) = Stag_P / DUM**(Gamma / (Gamma - 1.d0));
       surf_T(ic,kc) = Stag_T / DUM;
     endif
   enddo
   if( TWall > 0.d0 ) then
     do ic = 1,nx_local
       surf_T(ic,kc) = TWall;
     enddo
   endif
  enddo

  ! Form the whole flow field(Basic Variables)
  ! Linear, Square of sqrt distribution are used.
  do kc = 1,nz_local
   do ic = 1,nx_local
     tempP = P(ic,Ny,kc) - surf_p(ic,kc);
     tempT = T(ic,Ny,kc) - surf_T(ic,kc);
   	! Initial field--Inner layers
     do jc = 1,Ny-1
          indexJ = REAL((jc - 1.d0) / (Ny - 1.d0),rk);
         U(ic,jc,kc) = indexJ * indexJ * U(ic,Ny,kc);
         V(ic,jc,kc) = indexJ * indexJ * V(ic,Ny,kc);
         W(ic,jc,kc) = indexJ * indexJ * W(ic,Ny,kc);
         P(ic,jc,kc) = indexJ * indexJ * tempP + surf_p(ic,kc);
         T(ic,jc,kc) = sqrt(indexJ) * tempT + surf_T(ic,kc);
       Rho(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / T(ic,jc,kc); 
     enddo
   enddo
  enddo
   ! For initial conservative variables
   do kc = 1,nz_local
    do jc = 1,Ny
     do ic = 1,nx_local
       Ucons0(1,ic,jc,kc) = Rho(ic,jc,kc);
       Ucons0(2,ic,jc,kc) = Rho(ic,jc,kc) * U(ic,jc,kc);
       Ucons0(3,ic,jc,kc) = Rho(ic,jc,kc) * V(ic,jc,kc); 
       Ucons0(4,ic,jc,kc) = Rho(ic,jc,kc) * W(ic,jc,kc); 
       Ucons0(5,ic,jc,kc) =   P(ic,jc,kc) / (Gamma - 1.d0) + 0.5d0 * Rho(ic,jc,kc) * &
                            ( U(ic,jc,kc)**2.d0 + V(ic,jc,kc)**2.d0 + W(ic,jc,kc)**2.d0 );
     enddo
    enddo
   enddo
    
 endif
 
  ! Here we need to define the initial BigXi, BigEta and BigZeta as well as the Mu
  DO kc = 1,nz_local
   DO jc = 1,Ny 
    DO ic = 1,nx_local
     ! Define the speed of sound
     Cs(ic,jc,kc) = sqrt(Gamma * P(ic,jc,kc) / Rho(ic,jc,kc))
     ! Define the velocity along xi, eta, zeta
     BigU_Xi(ic,jc,kc)= U(ic,jc,kc) *  dxidx(ic,jc,kc) &
                     &+ V(ic,jc,kc) *  dxidy(ic,jc,kc) &
                     &+ W(ic,jc,kc) *  dxidz(ic,jc,kc) &
                     &+ dxidt(ic,jc,kc)

    BigU_Eta(ic,jc,kc)= U(ic,jc,kc) * detadx(ic,jc,kc) &
                     &+ V(ic,jc,kc) * detady(ic,jc,kc) &
                     &+ W(ic,jc,kc) * detadz(ic,jc,kc) &
                     &+ detadt(ic,jc,kc)
    
    BigU_Zeta(ic,jc,kc)= U(ic,jc,kc) * dzetadx(ic,jc,kc) &
                      &+ V(ic,jc,kc) * dzetady(ic,jc,kc) &
                      &+ W(ic,jc,kc) * dzetadz(ic,jc,kc) &
                      &+ dzetadt(ic,jc,kc)

     ! Define the viscousity coefficient
     Mu(ic,jc,kc) = (T(ic,jc,kc)**1.5d0) * (1.d0 + Csthlnd_Ref)/(T(ic,jc,kc) + Csthlnd_Ref)/Re_Ref;
    END DO
   END DO
  END DO
  
  ! Output the initial fields
    ! Exchange flow related variables
    call Parallel_Exchange(Rho)
    call Parallel_Exchange(U)
    call Parallel_Exchange(V)
    call Parallel_Exchange(W)
    call Parallel_Exchange(P)
    call Parallel_Exchange(T)
    call Parallel_Exchange_NumVar(Ucons0)
    ! Exchange the shock related variables
    call Parallel_Exchange_Surface(ShockH)
    call Parallel_Exchange_Surface(ShockV)

  call output_results  

end subroutine Initial_Fields