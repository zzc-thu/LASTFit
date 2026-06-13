subroutine Singular_Initial_Fields
  ! This subroutine form the initial grids used for cone or cone like geometry
  ! Ths singular axis is located along the xi=0 line(plane)

  ! Module dependencies
  use SF_Constant,   only: ik,rk,NumVar,FilesForContinue,pi
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
                           ShockH,ShockV,ModelType,BigU_Xi_total,BigU_Eta_total,BigU_Zeta_total,&
                           Rho_total,U_total,V_total,W_total,P_total,T_total,shockH_total,shockV_total,Ucons0_total

  ! some subroutines and functions
  use MPI_GLOBAL,     only: Parallel_Exchange,MPI_Barrier,&
                          & Parallel_Exchange_NumVar,Parallel_Exchange_Surface
  use OutputParaView, only: output_results,output_results_singularity
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
 
    call Singular_Calculate_Jaco

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

    call Parallel_Exchange_singularity_negative(BigU_Xi,BigU_Xi_total)
    call Parallel_Exchange_singularity(BigU_Eta,BigU_Eta_total)
    call Parallel_Exchange_singularity(BigU_Zeta,BigU_Zeta_total)

    
    call Parallel_Exchange(BigU_Xi)
    call Parallel_Exchange(BigU_Eta)
    call Parallel_Exchange(BigU_Zeta)

    ! Exchange flow related variables

    CALL Parallel_Exchange_singularity(Rho,Rho_total)
    CALL Parallel_Exchange_singularity(U,U_total)
    CALL Parallel_Exchange_singularity(V,V_total)
    CALL Parallel_Exchange_singularity(W,W_total)
    CALL Parallel_Exchange_singularity(P,P_total)
    CALL Parallel_Exchange_singularity(T,T_total)
    CALL Parallel_Exchange_singularity_NumVar(Ucons0,Ucons0_total)
    ! Exchange the shock related variables
    CALL Parallel_Exchange_singularity_Surface(ShockH,ShockH_total)
    CALL Parallel_Exchange_singularity_Surface(ShockV,ShockV_total)
  
    
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

     if (ModelType == 1)then
         call output_results
     else
         call output_results_singularity
     endif

end subroutine Singular_Initial_Fields
    
subroutine Elliptical_Initial_Fields
  ! This subroutine forms the initial flow fields for elliptical cone geometry
  ! The initial field is designed to be physically consistent with elliptical geometry

  ! Module dependencies
  use SF_Constant,   only: ik, rk, NumVar, FilesForContinue, pi,overLap
  use SF_CFD_Global, only: IF_Continue_Calculate, nx_local, nz_local, &
                           Ny, Sxi_Surface, TWall, Re_Ref, CSTHLND_REF, &
                           P_INF, T0_INF, Gamma, Mach_Ref, Mach_Normal, &
                           Geometry_Type, R_head, theta_p_rad
  use MPI_GLOBAL,    only: MyId, ierr, MPI_COMM_WORLD

  ! Arrays and vectors
  use SF_CFD_Global, only: U, V, W, P, T, Rho, Ucons0, &
                           dxidx, dxidy, dxidz, dxidt, &
                           detadx, detady, detadz, detadt, &
                           dzetadx, dzetady, dzetadz, dzetadt, &
                           BigU_Xi, BigU_Eta, BigU_Zeta, Cs, Mu, &
                           ShockH, ShockV, ModelType, BigU_Xi_total, BigU_Eta_total, BigU_Zeta_total, &
                           Rho_total, U_total, V_total, W_total, P_total, T_total, &
                           shockH_total, shockV_total, Ucons0_total, &
                           WallNormalX, WallNormalY, WallNormalZ, &
                           WallSX, WallSY, WallSZ

  ! Subroutines and functions
  use MPI_GLOBAL,     only: Parallel_Exchange, MPI_Barrier, &
                          & Parallel_Exchange_NumVar, Parallel_Exchange_Surface
  use OutputParaView, only: output_results, output_results_singularity
  use SFitting,       only: ShockRelation3D

  implicit none
  
  ! Local variables
  logical :: dir_exists
  character(len=100) :: FILENAME_FLOW
  character(len=100) :: FILENAME_SHK
  integer(kind=ik) :: ic, jc, kc 
  real(kind=rk) :: tempNym, MachNormalShock
  real(kind=rk) :: Stag_P, Stag_T, Tstar, Pwell
  real(kind=rk) :: surf_p(1-overLAP:nx_local+overLAP, 1-overLAP:nz_local+overLAP)
  real(kind=rk) :: surf_T(1-overLAP:nx_local+overLAP, 1-overLAP:nz_local+overLAP)
  real(kind=rk) :: DEG, SMS, DUM
  real(kind=rk) :: tempP, tempT, indexJ
  real(kind=rk) :: norm_x, norm_y, norm_z
  real(kind=rk) :: wall_angle_x, wall_angle_y, wall_angle_z
  real(kind=rk) :: total_wall_angle, aspect_ratio
  real(kind=rk) :: local_mach, pressure_coef, temperature_coef
  real(kind=rk) :: curvature_effect, circumferential_angle
  real(kind=rk) :: effective_radius, local_radius

  ! Aspect ratio for elliptical cross-section (consistent with grid generation)
  aspect_ratio = 2.0_rk  

  ! ============================================================================
  ! CASE 1: Continue from previous simulation
  ! ============================================================================
  If (IF_Continue_Calculate == 1) then
    if (MyId == 0) then
      write(*, '(A)') "Reading initial fields from previous simulation"
      inquire(file="./RESU/SF_Results00000.flowsfg", exist=dir_exists)
      if (.not. dir_exists) then
        write(*, '(A)') "CONT file does not exist - starting new simulation"
        stop
      else
        write(*, '(A)') "CONT file exists - reading initial fields"
      endif     
    endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      
    ! Read previous flow results
    write(filename_flow, '(A,A,I5.5,A)') 'RESU/', trim(FilesForContinue), MyId, '.flowsfg'
    if (MyId == 0) write(*, *) "Reading flow field"
    open(112, file=filename_flow, form='unformatted', status='old')
    read(112) Ucons0
    close(112)
    
    ! Convert conservative to primitive variables
    do kc = 1-overLAP, nz_local+overLAP
      do jc = 1, Ny
        do ic = 1-overLAP, nx_local+overLAP
          call CVtoU_CPG(Rho(ic,jc,kc), U(ic,jc,kc), V(ic,jc,kc), W(ic,jc,kc), &
                        P(ic,jc,kc), Ucons0(1:NumVar,ic,jc,kc))
          T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc)
        enddo
      enddo
    enddo

    ! Read previous shock information
    write(filename_shk, '(A,A,I5.5,A)') 'RESU/', trim(FilesForContinue), MyId, '.shksfg'
    if (MyId == 0) write(*, *) "Reading shock height and velocity"
    open(113, file=filename_shk, form='unformatted', status='old')
    read(113) ShockH
    read(113) ShockV
    close(113)

    ! Update Jacobian
    call Singular_Calculate_Jaco

    if (MyId == 0) then
      write(*, '(A)') "Successfully read initial fields from previous simulation"
    endif

  ! ============================================================================
  ! CASE 2: Start new simulation with elliptical cone initial field
  ! ============================================================================
  else
    if (MyId == 0) then
      write(*, '(A)') "Starting new simulation with elliptical cone initial field"
    endif

    ! Get post-shock values at jc = Ny (shock surface)
    call ShockRelation3D

    ! ========================================================================
    ! ELLIPTICAL CONE SPECIFIC WALL PRESSURE AND TEMPERATURE DISTRIBUTION
    ! ========================================================================
    
    ! Parameters optimized for elliptical geometry
    Tstar = 45.0_rk * pi / 180.0_rk  ! Reference angle for elliptical cone
    Pwell = 2.5_rk                   ! Pressure distribution exponent
    
    ! Calculate post-shock stagnation conditions
    MachNormalShock = sqrt((2.0_rk + (Gamma - 1.0_rk) * Mach_Normal**2) / &
                          (2.0_rk * Gamma * Mach_Normal**2 - Gamma + 1.0_rk))

    Stag_P = P_INF * ((2.0_rk * Gamma * Mach_Ref * Mach_Ref / (Gamma + 1.0_rk)) - &
                     ((Gamma - 1.0_rk) / (Gamma + 1.0_rk)))
    Stag_P = Stag_P * (1.0_rk + 0.5_rk * (Gamma - 1.0_rk) * MachNormalShock**2)**(Gamma / (Gamma - 1.0_rk))
    Stag_T = T0_INF

    ! Calculate surface pressure and temperature distribution for elliptical cone
    do kc = 1-overLAP, nz_local+overLAP
      do ic = 1-overLAP, nx_local+overLAP
        ! Get local wall coordinates and normal vector
        norm_x = WallNormalX(ic, kc)
        norm_y = WallNormalY(ic, kc) 
        norm_z = WallNormalZ(ic, kc)
        
        ! Calculate circumferential position for elliptical correction
        circumferential_angle = atan2(WallSZ(ic,kc) * aspect_ratio, WallSY(ic,kc))
        
        ! Calculate local radius considering elliptical cross-section
        local_radius = sqrt((WallSY(ic,kc)**2) + (WallSZ(ic,kc) * aspect_ratio)**2)
        
        ! Effective radius for curvature effects (smaller radius = higher curvature)
        effective_radius = R_head * (1.0_rk + 0.3_rk * cos(2.0_rk * circumferential_angle))
        
        ! Calculate wall angles considering elliptical geometry
        wall_angle_x = acos(abs(norm_x))
        
        ! For elliptical cones, use combined angle that accounts for both meridional
        ! and circumferential curvature variations
        if (local_radius > 1.0e-10_rk) then
          curvature_effect = R_head / effective_radius
        else
          curvature_effect = 1.0_rk
        endif
        
        ! Enhanced wall angle calculation for elliptical geometry
        DEG = wall_angle_x * (1.0_rk + 0.2_rk * curvature_effect)
        
        ! Limit maximum angle for stability
        DEG = min(DEG, 0.8_rk * pi/2.0_rk)
        
        ! Modified pressure distribution for elliptical geometry
        ! Includes circumferential variation for elliptical cross-section
        SMS = (DEG / Tstar)**Pwell * (1.0_rk + 0.1_rk * cos(2.0_rk * circumferential_angle))
        DUM = 1.0_rk + 0.5_rk * (Gamma - 1.0_rk) * SMS
        
        if (abs(DUM) < 1.0e-11_rk) then
          surf_p(ic, kc) = Stag_P
          surf_T(ic, kc) = Stag_T
        else
          ! Elliptical cone pressure distribution
          pressure_coef = 1.0_rk / DUM**(Gamma / (Gamma - 1.0_rk))
          temperature_coef = 1.0_rk / DUM
          
          surf_p(ic, kc) = P_INF * pressure_coef
          surf_T(ic, kc) = T0_INF * temperature_coef
        endif
        
        ! Additional adjustment for nose region
        if (WallSX(ic,kc) < R_head * 0.5_rk) then
          ! Enhance pressure in nose region for elliptical cone
          surf_p(ic, kc) = surf_p(ic, kc) * (1.0_rk + 0.3_rk * exp(-WallSX(ic,kc)/(R_head*0.1_rk)))
        endif
      enddo
    enddo

    ! Apply wall temperature boundary condition if specified
    if (TWall > 0.0_rk) then
      do kc = 1-overLAP, nz_local+overLAP
        do ic = 1-overLAP, nx_local+overLAP
          surf_T(ic, kc) = TWall
        enddo
      enddo
    endif

    ! ========================================================================
    ! IMPROVED INITIAL FLOW FIELD DISTRIBUTION FOR ELLIPTICAL CONE
    ! ========================================================================
    do kc = 1-overLAP, nz_local+overLAP
      do ic = 1-overLAP, nx_local+overLAP
        tempP = P(ic, Ny, kc) - surf_p(ic, kc)
        tempT = T(ic, Ny, kc) - surf_T(ic, kc)
        
        ! Calculate circumferential variation for velocity distribution
        circumferential_angle = atan2(WallSZ(ic,kc) * aspect_ratio, WallSY(ic,kc))
        local_radius = sqrt((WallSY(ic,kc)**2) + (WallSZ(ic,kc) * aspect_ratio)**2)
        
        ! Initial field - Optimized for elliptical cone geometry
        do jc = 1, Ny-1
          indexJ = real((jc - 1.0_rk) / (Ny - 1.0_rk), rk)
          
          ! Velocity distribution with elliptical correction
          ! Reduced exponent for smoother transition
          U(ic, jc, kc) = indexJ**1.6_rk * U(ic, Ny, kc)
          V(ic, jc, kc) = indexJ**1.6_rk * V(ic, Ny, kc) 
          W(ic, jc, kc) = indexJ**1.6_rk * W(ic, Ny, kc)
          
          ! Pressure distribution - nearly linear for better stability
          P(ic, jc, kc) = indexJ * tempP + surf_p(ic, kc)
          
          ! Temperature distribution - smoother profile
          T(ic, jc, kc) = indexJ**0.7_rk * tempT + surf_T(ic, kc)
          
          ! Density from equation of state
          Rho(ic, jc, kc) = P(ic, jc, kc) * Gamma * Mach_Ref * Mach_Ref / T(ic, jc, kc)
        enddo
        
        ! Set boundary values explicitly
        jc = Ny
        P(ic, jc, kc) = P(ic, jc, kc)  ! Keep shock value
        T(ic, jc, kc) = T(ic, jc, kc)  ! Keep shock value
        Rho(ic, jc, kc) = Rho(ic, jc, kc)  ! Keep shock value
        
        jc = 1
        P(ic, jc, kc) = surf_p(ic, kc)
        T(ic, jc, kc) = surf_T(ic, kc)
        Rho(ic, jc, kc) = P(ic, jc, kc) * Gamma * Mach_Ref * Mach_Ref / T(ic, jc, kc)
      enddo
    enddo

    ! Initialize conservative variables
    do kc = 1-overLAP, nz_local+overLAP
      do jc = 1, Ny
        do ic = 1-overLAP, nx_local+overLAP
          Ucons0(1, ic, jc, kc) = Rho(ic, jc, kc)
          Ucons0(2, ic, jc, kc) = Rho(ic, jc, kc) * U(ic, jc, kc)
          Ucons0(3, ic, jc, kc) = Rho(ic, jc, kc) * V(ic, jc, kc)
          Ucons0(4, ic, jc, kc) = Rho(ic, jc, kc) * W(ic, jc, kc)
          Ucons0(5, ic, jc, kc) = P(ic, jc, kc) / (Gamma - 1.0_rk) + &
                                  0.5_rk * Rho(ic, jc, kc) * &
                                  (U(ic, jc, kc)**2 + V(ic, jc, kc)**2 + W(ic, jc, kc)**2)
        enddo
      enddo
    enddo
    
  endif  ! End of new simulation case

  ! ============================================================================
  ! INITIALIZE DERIVED QUANTITIES
  ! ============================================================================
  do kc = 1-overLAP, nz_local+overLAP
    do jc = 1, Ny
      do ic = 1-overLAP, nx_local+overLAP
        ! Speed of sound
        Cs(ic, jc, kc) = sqrt(Gamma * P(ic, jc, kc) / Rho(ic, jc, kc))
        
        ! Contravariant velocities
        BigU_Xi(ic, jc, kc) = U(ic, jc, kc) * dxidx(ic, jc, kc) + &
                              V(ic, jc, kc) * dxidy(ic, jc, kc) + &
                              W(ic, jc, kc) * dxidz(ic, jc, kc) + &
                              dxidt(ic, jc, kc)

        BigU_Eta(ic, jc, kc) = U(ic, jc, kc) * detadx(ic, jc, kc) + &
                               V(ic, jc, kc) * detady(ic, jc, kc) + &
                               W(ic, jc, kc) * detadz(ic, jc, kc) + &
                               detadt(ic, jc, kc)
        
        BigU_Zeta(ic, jc, kc) = U(ic, jc, kc) * dzetadx(ic, jc, kc) + &
                                V(ic, jc, kc) * dzetady(ic, jc, kc) + &
                                W(ic, jc, kc) * dzetadz(ic, jc, kc) + &
                                dzetadt(ic, jc, kc)

        ! Viscosity coefficient (Sutherland's law form)
        Mu(ic, jc, kc) = (T(ic, jc, kc)**1.5_rk) * (1.0_rk + CSTHLND_REF) / &
                         (T(ic, jc, kc) + CSTHLND_REF) / Re_Ref
      enddo
    enddo
  enddo

  ! ============================================================================
  ! PARALLEL DATA EXCHANGE
  ! ============================================================================
  
  ! Exchange contravariant velocities
  call Parallel_Exchange_singularity_negative(BigU_Xi, BigU_Xi_total)
  call Parallel_Exchange_singularity(BigU_Eta, BigU_Eta_total)
  call Parallel_Exchange_singularity(BigU_Zeta, BigU_Zeta_total)
  
  call Parallel_Exchange(BigU_Xi)
  call Parallel_Exchange(BigU_Eta)
  call Parallel_Exchange(BigU_Zeta)

  ! Exchange flow field variables
  call Parallel_Exchange_singularity(Rho, Rho_total)
  call Parallel_Exchange_singularity(U, U_total)
  call Parallel_Exchange_singularity(V, V_total)
  call Parallel_Exchange_singularity(W, W_total)
  call Parallel_Exchange_singularity(P, P_total)
  call Parallel_Exchange_singularity(T, T_total)
  call Parallel_Exchange_singularity_NumVar(Ucons0, Ucons0_total)
  
  ! Exchange shock variables
  call Parallel_Exchange_singularity_Surface(ShockH, ShockH_total)
  call Parallel_Exchange_singularity_Surface(ShockV, ShockV_total)
  
  ! Final exchanges for synchronization
  call Parallel_Exchange(Rho)
  call Parallel_Exchange(U)
  call Parallel_Exchange(V)
  call Parallel_Exchange(W)
  call Parallel_Exchange(P)
  call Parallel_Exchange(T)
  call Parallel_Exchange_NumVar(Ucons0)
  call Parallel_Exchange_Surface(ShockH)
  call Parallel_Exchange_Surface(ShockV)

  ! ============================================================================
  ! OUTPUT INITIAL FIELD
  ! ============================================================================
  if (ModelType == 1) then
    call output_results
  else
    call output_results_singularity
  endif

  if (MyId == 0) then
    write(*, '(A)') "Elliptical cone initial field setup completed"
  endif

end subroutine Elliptical_Initial_Fields