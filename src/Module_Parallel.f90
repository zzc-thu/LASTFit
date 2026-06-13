MODULE MPI_GLOBAL
  USE SF_Constant,  only:ik,rk,overLAP,NumVar
  !use SF_CFD_Global,only: nx_local,Ny,nz_local,Nz
  ! mpi module is provided by the compiler
  USE mpi
  implicit none
  ! Whether Parallel or Series
  logical :: If_parallel
  ! Define MPI Buffer Size
  ! Define the buffer layer for MPI sendrecv
  INTEGER( kind = ik ), PARAMETER:: IBUFFER_SIZE = 100000000;
  ! Define the MPI Communicator
  INTEGER( kind = ik ):: MyId           ! My process ID           
  INTEGER( kind = ik ):: NumProcess     ! Number of processes
  ! define the process along I directions
  ! ID_XM1 I- direction, ID_XP1 I+ direction
  INTEGER( kind = ik ):: ID_XM1, ID_XP1
  ! npx is the present process number along I directions
  ! npx0 is the total process number along I directions
  INTEGER( kind = ik ):: npx, npx0
  ! ID_ZM1 K- direction, ID_ZP1 K+ direction
  INTEGER( kind = ik ):: ID_ZM1, ID_ZP1
  ! npz is the present process number along K directions
  ! npz0 is the total process number along K directions
  INTEGER( kind = ik ):: npz, npz0
  
  ! user defined MPI DATA TYPE (for Send & Recv)
  INTEGER( kind = ik ):: TYPE_LAPX1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPZ2

  ! MPI error information
  INTEGER( kind = ik ):: ierr
  integer :: color, MPI_comm_npx0 

  ! Define some MPI arrays
  ! I_offset is the start grid point of the npx process
  ! I_nn is the number of grid points at the npx process
  INTEGER( kind = ik )::I_offset(0:1024),I_nn(0:1024)
  ! K_offset is the start grid point of the npz process
  ! K_nn is the number of grid points at the npz process
  INTEGER( kind = ik )::K_offset(0:1024),K_nn(0:1024)

  ! Define MPI Buffer Zone
  REAL( kind = rk ):: BUFFER_MPI(IBUFFER_SIZE)

  ! allocate some arrays for multi-variables transfer
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_2D_xi(:),Temp_Send2_2D_xi(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_2D_xi(:),Temp_Recv2_2D_xi(:)
  
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_2D_zeta(:),Temp_Send2_2D_zeta(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_2D_zeta(:),Temp_Recv2_2D_zeta(:)

  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_2D_xi_numvar(:),Temp_Send2_2D_xi_numvar(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_2D_xi_numvar(:),Temp_Recv2_2D_xi_numvar(:)

  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_2D_zeta_numvar(:),Temp_Send2_2D_zeta_numvar(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_2D_zeta_numvar(:),Temp_Recv2_2D_zeta_numvar(:)

  ! allocate some arrays for data transfer in 3D along xi and zeta directions
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_xi(:),Temp_Send2_3D_xi(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_xi(:),Temp_Recv2_3D_xi(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_zeta(:),Temp_Send2_3D_zeta(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_zeta(:),Temp_Recv2_3D_zeta(:)

  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_xi_interp(:),Temp_Send2_3D_xi_interp(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_xi_interp(:),Temp_Recv2_3D_xi_interp(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_zeta_interp(:),Temp_Send2_3D_zeta_interp(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_zeta_interp(:),Temp_Recv2_3D_zeta_interp(:)
  ! allocate some arrays for data transfer in 3D along xi and zeta directions with NumVar 
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_xi_numvar(:),Temp_Send2_3D_xi_numvar(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_xi_numvar(:),Temp_Recv2_3D_xi_numvar(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Send1_3D_zeta_numvar(:),Temp_Send2_3D_zeta_numvar(:)
  REAL( kind = rk ),ALLOCATABLE:: Temp_Recv1_3D_zeta_numvar(:),Temp_Recv2_3D_zeta_numvar(:)

  contains

! This subroutine is used to initialize the MPI environment
SUBROUTINE Parallel_Initilize
  Implicit none

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MyId, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcess, ierr)
  CALL MPI_BUFFER_ATTACH(BUFFER_MPI,8*IBUFFER_SIZE, ierr)

  
    If_parallel = .TRUE.
  
  if( NumProcess == 1) then
    
    If_parallel = .FALSE.

  endif

  if( MyId == 0 ) then                       
   write(*,"(A)")"|===================================================================================================|"
   write(*,"(A)")"|                                                                                                   |" 
   write(*,"(A)")"|                                                                                                   |"   
   write(*,"(A)")"|   _____ _                _   ______ _ _   _   _             _____       _                         |"
   write(*,"(A)")"|  /  ___| |              | |  |  ___(_) | | | (_)           /  ___|     | |                        |"
   write(*,"(A)")"|  \ \--\| |__   ___   ___| | _| |_   _| |_| |_ _ _ __   __ _\ \--\  ___ | |_   _____ _ __          |"
   write(*,"(A)")"|   \--\ \ |_ \ / _ \ / __| |/ /  _| | | __| __| | '_ \ / _` |\--\ \/ _ \| \ \ / / _ \ '__|         |"
   write(*,"(A)")"|  /\__/ / | | | (_) | (__|   <| |   | | |_| |_| | | | | (_| /\__/ / (_) | |\ V /  __/ |            |"
   write(*,"(A)")"|  \____/|_| |_|\___/ \___|_|\_\_|   |_|\__|\__|_|_| |_|\__  \____/ \___/|_| \_/ \___|_|            |"
   write(*,"(A)")"|                                                        __/ |                                      |"
   write(*,"(A)")"|                                                       |___/                                       |"
   write(*,"(A)")"|                                                                            Version 3.23.7         |"
   write(*,"(A)")"|                                                                                                   |"
   write(*,"(A)")"|                     The fitting solver is first developed at LAST Group @Tsinghua University      |"
   write(*,"(A)")"|                                                                                                   |"
   write(*,"(A)")"|                                                                         Developed By Youcheng Xi  |"
   write(*,"(A)")"|                                                                                                   |"
   write(*,"(A)")"|                                                                                                   |"
   write(*,"(A)")"|                                                                                                   |"
   write(*,"(A)")"|                                                                                                   |"  
   write(*,"(A)")"|===================================================================================================|"
  endif

END SUBROUTINE Parallel_Initilize

! This subroutine is used to finalize the MPI environment
! and deallocate the MPI buffer
SUBROUTINE Parallel_Finalize
   use mpi
   IMPLICIT NONE
    
    DEALLOCATE(Temp_Send1_2D_xi,Temp_Send2_2D_xi)
    DEALLOCATE(Temp_Recv1_2D_xi,Temp_Recv2_2D_xi)

    DEALLOCATE(Temp_Send1_2D_zeta,Temp_Send2_2D_zeta)
    DEALLOCATE(Temp_Recv1_2D_zeta,Temp_Recv2_2D_zeta)

    DEALLOCATE(Temp_Recv1_2D_xi_numvar,Temp_Recv2_2D_xi_numvar)
    DEALLOCATE(Temp_Send1_2D_xi_numvar,Temp_Send2_2D_xi_numvar)

    DEALLOCATE(Temp_Recv1_2D_zeta_numvar,Temp_Recv2_2D_zeta_numvar)
    DEALLOCATE(Temp_Send1_2D_zeta_numvar,Temp_Send2_2D_zeta_numvar)

    DEALLOCATE(Temp_Send1_3D_xi,Temp_Send2_3D_xi)
    DEALLOCATE(Temp_Recv1_3D_xi,Temp_Recv2_3D_xi)

    DEALLOCATE(Temp_Send1_3D_zeta,Temp_Send2_3D_zeta)
    DEALLOCATE(Temp_Recv1_3D_zeta,Temp_Recv2_3D_zeta)

    DEALLOCATE(Temp_Send1_3D_xi_interp,Temp_Send2_3D_xi_interp)
    DEALLOCATE(Temp_Recv1_3D_xi_interp,Temp_Recv2_3D_xi_interp)

    DEALLOCATE(Temp_Send1_3D_zeta_interp,Temp_Send2_3D_zeta_interp)
    DEALLOCATE(Temp_Recv1_3D_zeta_interp,Temp_Recv2_3D_zeta_interp)
    
    DEALLOCATE(Temp_Send1_3D_xi_numvar)
    DEALLOCATE(Temp_Send2_3D_xi_numvar)
    DEALLOCATE(Temp_Recv1_3D_xi_numvar)
    DEALLOCATE(Temp_Recv2_3D_xi_numvar)

    DEALLOCATE(Temp_Send1_3D_zeta_numvar)
    DEALLOCATE(Temp_Send2_3D_zeta_numvar)
    DEALLOCATE(Temp_Recv1_3D_zeta_numvar)
    DEALLOCATE(Temp_Recv2_3D_zeta_numvar)
    
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   CALL MPI_FINALIZE(ierr)

END SUBROUTINE Parallel_Finalize

! This subroutine is used to allocate the MPI transfer arrays
SUBROUTINE Parallel_Allocate
  USE SF_Constant,  Only: overLAP,NumVar
  USE SF_CFD_Global,Only: Ny,nx_local,nz_local
  implicit none

    ALLOCATE(Temp_Send1_2D_xi(overLAP * nz_local),Temp_Send2_2D_xi(overLAP * nz_local))
    ALLOCATE(Temp_Recv1_2D_xi(overLAP * nz_local),Temp_Recv2_2D_xi(overLAP * nz_local))
     Temp_Send1_2D_xi = 0.0_rk
     Temp_Recv1_2D_xi = 0.0_rk
     Temp_Send2_2D_xi = 0.0_rk
     Temp_Recv2_2D_xi = 0.0_rk

    ALLOCATE(Temp_Send1_2D_zeta(overLAP * nx_local),Temp_Send2_2D_zeta(overLAP * nx_local))
    ALLOCATE(Temp_Recv1_2D_zeta(overLAP * nx_local),Temp_Recv2_2D_zeta(overLAP * nx_local))
     Temp_Send1_2D_zeta = 0.0_rk
     Temp_Recv1_2D_zeta = 0.0_rk
     Temp_Send2_2D_zeta = 0.0_rk
     Temp_Recv2_2D_zeta = 0.0_rk

    ALLOCATE(Temp_Recv1_2D_xi_numvar(overLAP * nz_local * NumVar),Temp_Recv2_2D_xi_numvar(overLAP * nz_local * NumVar))
    ALLOCATE(Temp_Send1_2D_xi_numvar(overLAP * nz_local * NumVar),Temp_Send2_2D_xi_numvar(overLAP * nz_local * NumVar))
     Temp_Send1_2D_xi_numvar = 0.0_rk
     Temp_Recv1_2D_xi_numvar = 0.0_rk
     Temp_Send2_2D_xi_numvar = 0.0_rk
     Temp_Recv2_2D_xi_numvar = 0.0_rk

    ALLOCATE(Temp_Recv1_2D_zeta_numvar(overLAP * nx_local * NumVar),Temp_Recv2_2D_zeta_numvar(overLAP * nx_local * NumVar))
    ALLOCATE(Temp_Send1_2D_zeta_numvar(overLAP * nx_local * NumVar),Temp_Send2_2D_zeta_numvar(overLAP * nx_local * NumVar))
     Temp_Send1_2D_zeta_numvar = 0.0_rk
     Temp_Recv1_2D_zeta_numvar = 0.0_rk
     Temp_Send2_2D_zeta_numvar = 0.0_rk
     Temp_Recv2_2D_zeta_numvar = 0.0_rk

    ALLOCATE(Temp_Send1_3D_xi(overLAP * Ny * nz_local),Temp_Send2_3D_xi(overLAP * Ny * nz_local))
    ALLOCATE(Temp_Recv1_3D_xi(overLAP * Ny * nz_local),Temp_Recv2_3D_xi(overLAP * Ny * nz_local))
     Temp_Send1_3D_xi = 0.0_rk
     Temp_Recv1_3D_xi = 0.0_rk
     Temp_Send2_3D_xi = 0.0_rk
     Temp_Recv2_3D_xi = 0.0_rk
    
    ALLOCATE(Temp_Send1_3D_zeta(overLAP * Ny * (nx_local+2_ik*overLAP)),Temp_Send2_3D_zeta(overLAP * Ny * (nx_local+2_ik*overLAP)))
    ALLOCATE(Temp_Recv1_3D_zeta(overLAP * Ny * (nx_local+2_ik*overLAP)),Temp_Recv2_3D_zeta(overLAP * Ny * (nx_local+2_ik*overLAP)))
     Temp_Send1_3D_zeta = 0.0_rk
     Temp_Recv1_3D_zeta = 0.0_rk
     Temp_Send2_3D_zeta = 0.0_rk
     Temp_Recv2_3D_zeta = 0.0_rk

    ALLOCATE(Temp_Send1_3D_xi_interp(overLAP * (2*Ny-1) * nz_local),Temp_Send2_3D_xi_interp(overLAP * (2*Ny-1) * nz_local))
    ALLOCATE(Temp_Recv1_3D_xi_interp(overLAP * (2*Ny-1) * nz_local),Temp_Recv2_3D_xi_interp(overLAP * (2*Ny-1) * nz_local))
     Temp_Send1_3D_xi = 0.0_rk
     Temp_Recv1_3D_xi = 0.0_rk
     Temp_Send2_3D_xi = 0.0_rk
     Temp_Recv2_3D_xi = 0.0_rk
    
    ALLOCATE(Temp_Send1_3D_zeta_interp(overLAP * (2*Ny-1) * (nx_local+2_ik*overLAP)),Temp_Send2_3D_zeta_interp(overLAP * (2*Ny-1) * (nx_local+2_ik*overLAP)))
    ALLOCATE(Temp_Recv1_3D_zeta_interp(overLAP * (2*Ny-1) * (nx_local+2_ik*overLAP)),Temp_Recv2_3D_zeta_interp(overLAP * (2*Ny-1) * (nx_local+2_ik*overLAP)))
     Temp_Send1_3D_zeta = 0.0_rk
     Temp_Recv1_3D_zeta = 0.0_rk
     Temp_Send2_3D_zeta = 0.0_rk
     Temp_Recv2_3D_zeta = 0.0_rk 
     
    ALLOCATE(Temp_Send1_3D_xi_numvar(overLAP * Ny * nz_local * NumVar))
    ALLOCATE(Temp_Send2_3D_xi_numvar(overLAP * Ny * nz_local * NumVar))
    ALLOCATE(Temp_Recv1_3D_xi_numvar(overLAP * Ny * nz_local * NumVar))
    ALLOCATE(Temp_Recv2_3D_xi_numvar(overLAP * Ny * nz_local * NumVar))

    Temp_Send1_3D_xi_numvar = 0.d0;
    Temp_Send2_3D_xi_numvar = 0.d0;
    Temp_Recv1_3D_xi_numvar = 0.d0;
    Temp_Recv2_3D_xi_numvar = 0.d0;  
    
    ! 考虑了对应的角点处的数据传输
    ALLOCATE(Temp_Send1_3D_zeta_numvar(overLAP * Ny * (nx_local+2_ik*overLAP) * NumVar))
    ALLOCATE(Temp_Send2_3D_zeta_numvar(overLAP * Ny * (nx_local+2_ik*overLAP) * NumVar))
    ALLOCATE(Temp_Recv1_3D_zeta_numvar(overLAP * Ny * (nx_local+2_ik*overLAP) * NumVar))
    ALLOCATE(Temp_Recv2_3D_zeta_numvar(overLAP * Ny * (nx_local+2_ik*overLAP) * NumVar))
    
    Temp_Send1_3D_zeta_numvar = 0.d0;
    Temp_Send2_3D_zeta_numvar = 0.d0;
    Temp_Recv1_3D_zeta_numvar = 0.d0;
    Temp_Recv2_3D_zeta_numvar = 0.d0;


END SUBROUTINE Parallel_Allocate

! This subroutine is used to split the domain along I and K directions
! We obtain the npx and npz id for each process
! We also obtain the I_offset and I_nn for each process as well as K_offset and K_nn
SUBROUTINE Parallel_SplitXZ
    USE SF_Constant, Only: ik
    USE SF_CFD_Global, Only: Nx,Ny,Nz,nx_local,nz_local
    implicit none
    integer( kind = ik ):: k,ka
    integer( kind = ik ):: npx1,npx2
    integer( kind = ik ):: npz1,npz2

    
  ! If the corresponding number of parallel cores does 
  ! not match the specified number, the corresponding program exits
  if(Numprocess .NE. npx0*npz0) then
   if(myid.EQ.0) then
     write(*,*)"npx0*npz0!!! MPI numbers Wrong" 
   endif
   CALL MPI_FINALIZE(ierr)
   STOP
  endif

  npx = mod(MyId,npx0);
  npz = MyId/npx0;

  nx_local = Nx/npx0;
  nz_local = Nz/npz0;
  if(npx .LT. mod(Nx,npx0)) nx_local = nx_local + 1
  if(npz .LT. mod(Nz,npz0)) nz_local = nz_local + 1

  do k = 0, npx0-1
    ka = min(k,mod(Nx,npx0))
    I_offset(k) = int(Nx/npx0)*k + ka + 1;
    I_nn(k) = Nx / npx0;
    if(k .LT. mod(Nx,npx0)) I_nn(k) = I_nn(k) + 1
  end do

  do k = 0, npz0-1
    ka = min(k,mod(Nz,npz0))
    K_offset(k) = int(Nz/npz0)*k + ka + 1;
    K_nn(k) = Nz / npz0;
    if(k .LT. mod(Nz,npz0)) K_nn(k) = K_nn(k) + 1
  end do
 
 !! we define a new MPI datatype for data transfer
 ! call New_MPI_datatype

 !================= I direction ========================  
  npx1 = my_mod1(npx-1,npx0); 
  npx2 = my_mod1(npx+1,npx0); 

  ID_XM1 = npz*npx0+npx1
  ID_XP1 = npz*npx0+npx2
  
  ! Start and End Points should be treated differently 
  if(npx.eq.0     ) ID_XM1 = MPI_PROC_NULL
  if(npx.eq.npx0-1) ID_XP1 = MPI_PROC_NULL

 !================= z direction ========================
 ! periodic boundary conditions 
  npz1 = my_mod1(npz-1,npz0)
  npz2 = my_mod1(npz+1,npz0)

  ID_ZM1 = npz1*npx0+npx
  ID_ZP1 = npz2*npx0+npx
  
  ! This part is only used for non-periodic boundary conditions
  ! Start and End Points should be treated differently 
  ! if(npz.eq.0     ) ID_ZM1 = MPI_PROC_NULL
  ! if(npz.eq.npz0-1) ID_ZP1 = MPI_PROC_NULL 

 !=======================================================
  
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  color = merge(0, MPI_UNDEFINED, npx == 0)  ! npx==0的进程color=0，否则为MPI_UNDEFINED 
  CALL MPI_Comm_split(MPI_COMM_WORLD, color, MyId, MPI_comm_npx0, ierr) 
  
  ! Output the Grid Information for each processes
    WRITE(*,1123)MyId+1, NumProcess, nx_local, Ny, nz_local
1123 FORMAT(' Process-',I5,'  /',I5, ' MPI Grid is activated. Allocated mesh size: ', I5, ' X', I5, ' X', I5, '.')
    
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)  
END SUBROUTINE Parallel_SplitXZ

subroutine New_MPI_datatype
   use SF_CFD_Global, only:nx_local,Ny,nz_local
   implicit none
   integer( kind = ik)::TYPE_tmp
  
   call MPI_TYPE_Vector(     Ny, overLAP,                 nx_local+2*overLAP,MPI_DOUBLE_PRECISION,TYPE_LAPX1,ierr)
  !call MPI_TYPE_Vector(overLAP,nx_local,                 nx_local+2*overLAP,MPI_DOUBLE_PRECISION,TYPE_LAPY1,ierr)
   call MPI_TYPE_Vector(overLAP,nx_local,(nx_local+2*overLAP)*(Ny+2*overLAP),MPI_DOUBLE_PRECISION,TYPE_LAPZ1,ierr)

   call MPI_TYPE_Vector(     Ny, nx_local,                nx_local+2*overLAP,MPI_DOUBLE_PRECISION,TYPE_tmp,ierr)

   call MPI_TYPE_HVector( nz_local,1,(nx_local+2*overLAP)*(Ny+2*overLAP)*rk, TYPE_LAPX1,TYPE_LAPX2,ierr)
  !call MPI_TYPE_HVector( nz_local,1,(nx_local+2*overLAP)*(Ny+2*overLAP)*rk, TYPE_LAPY1,TYPE_LAPY2,ierr)
   call MPI_TYPE_HVector(  overLAP,1,(nx_local+2*overLAP)*(Ny+2*overLAP)*rk, TYPE_tmp,TYPE_LAPZ2,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX1,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ1,ierr)

   call MPI_TYPE_COMMIT(TYPE_LAPX2,ierr)
   call MPI_TYPE_COMMIT(TYPE_LAPZ2,ierr)

   call MPI_barrier(MPI_COMM_WORLD,ierr)
end subroutine New_MPI_datatype

SUBROUTINE Parallel_Exchange_NumVar(f)
  USE SF_Constant,   only:NumVar
  USE SF_CFD_Global, only:nx_local,Ny,nz_local
  implicit none
  real( kind = rk )::f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)


  call exchange_boundary_x_numvar(f)

  call exchange_boundary_z_numvar(f)

END SUBROUTINE Parallel_Exchange_NumVar

SUBROUTINE Parallel_Exchange(f)
  use SF_CFD_Global, only:nx_local,Ny,nz_local
  implicit none
  real( kind = rk )::f(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)

    call exchange_boundary_x(f)

    call exchange_boundary_z(f)

  !return 
END SUBROUTINE Parallel_Exchange



SUBROUTINE Parallel_Exchange_interp(f)
  use SF_CFD_Global, only:nx_local,Ny,nz_local
  implicit none
  real( kind = rk )::f(1-overLap:nx_local+overLap,1:2*Ny-1,1-overLap:nz_local+overLap)

    call exchange_boundary_x_interp(f)

    call exchange_boundary_z_interp(f)

  !return 
END SUBROUTINE Parallel_Exchange_interp


SUBROUTINE Parallel_Exchange_NumVar_Surface(f)
  USE SF_Constant,   only:NumVar
  use SF_CFD_Global, only:nx_local,nz_local
  implicit none
  real( kind = rk )::f(1:NumVar,1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)

    call exchange_boundary_x_surface_numvar(f)
    call exchange_boundary_z_surface_numvar(f)

END SUBROUTINE Parallel_Exchange_NumVar_Surface

SUBROUTINE exchange_boundary_x_surface_numvar(f)
  use SF_Constant,       only:NumVar
  use SF_CFD_Global,     only:nx_local,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1:NumVar,1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  integer:: ic,kc,iVar,k1
  integer:: NumInt_MPI

  ! Here we perform the data transfer along I direction only
  do kc = 1,nz_local
      do ic = 1,overLAP
        do iVar = 1,NumVar
        k1 = (kc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
        Temp_Send1_2D_xi_numvar(k1) = f(iVar,ic,kc)                  ! left boundary of f
        Temp_Send2_2D_xi_numvar(k1) = f(iVar,nx_local-overLAP+ic,kc) ! right boundary of f
        enddo
      enddo
  enddo
   
  ! Here we count the total number of data used in MPI_SendRecv 
  NumInt_MPI = overLAP*nz_local*NumVar

  call MPI_Sendrecv(  Temp_Send1_2D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XM1, 9001,   &
                      Temp_Recv2_2D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_2D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XP1, 8001,   &           
                      Temp_Recv1_2D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XM1, 8001,MPI_COMM_WORLD,Status,ierr)   
  
    if( ID_XM1 .NE. MPI_PROC_NULL) then
      do kc = 1, nz_local
        do ic = 1, overLAP
          do iVar = 1, NumVar
            k1 = (kc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
            f(iVar,ic-overLAP,kc) = Temp_Recv1_2D_xi_numvar(k1)
          enddo
        enddo
      enddo
    endif

    if( ID_XP1 .NE. MPI_PROC_NULL) then
      do kc = 1, nz_local
        do ic = 1, overLAP
          do iVar = 1, NumVar
            k1 = (kc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
            f(iVar,nx_local+ic,kc) = Temp_Recv2_2D_xi_numvar(k1)
          enddo
        enddo
      enddo
    endif
END SUBROUTINE exchange_boundary_x_surface_numvar

SUBROUTINE exchange_boundary_z_surface_numvar(f)
  use SF_Constant,       only:NumVar
  use SF_CFD_Global,     only:nx_local,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1:NumVar,1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  integer:: ic,kc,iVar,k1
  integer:: NumInt_MPI

  ! Here we perform the data transfer along K direction only
  do kc = 1,overLAP
      do ic = 1,nx_local
        do iVar = 1,NumVar
        k1 = (kc-1)*nx_local*NumVar + (ic-1)*NumVar + iVar
        Temp_Send1_2D_zeta_numvar(k1) = f(iVar,ic,kc)                  ! left boundary of f
        Temp_Send2_2D_zeta_numvar(k1) = f(iVar,ic,kc+nz_local-overLAP) ! right boundary of f
        enddo
      enddo
  enddo
   
  ! Here we count the total number of data used in MPI_SendRecv 
  NumInt_MPI = overLAP*nx_local*NumVar

  call MPI_Sendrecv(  Temp_Send1_2D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_ZM1, 9001,   &
                      Temp_Recv2_2D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_ZP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_2D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_ZP1, 8001,   &           
                      Temp_Recv1_2D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_ZM1, 8001,MPI_COMM_WORLD,Status,ierr)   
  
    if( ID_ZM1 .NE. MPI_PROC_NULL) then
      do kc = 1, overLAP
        do ic = 1, nx_local
          do iVar = 1, NumVar
            k1 = (kc-1)*nx_local*NumVar + (ic-1)*NumVar + iVar
            f(iVar,ic,kc-overLAP) = Temp_Recv1_2D_zeta_numvar(k1)
          enddo
        enddo
      enddo
    endif

    if( ID_ZP1 .NE. MPI_PROC_NULL) then
      do kc = 1, overLAP
        do ic = 1, nx_local
          do iVar = 1, NumVar
            k1 = (kc-1)*nx_local*NumVar + (ic-1)*NumVar + iVar
            f(iVar,ic,nz_local+kc) = Temp_Recv2_2D_zeta_numvar(k1)
          enddo
        enddo
      enddo
    endif  
END SUBROUTINE exchange_boundary_z_surface_numvar

SUBROUTINE Parallel_Exchange_Surface(f)
  use SF_CFD_Global, only:nx_local,nz_local
  implicit none
  real( kind = rk )::f(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
    call exchange_boundary_x_surface(f)
    call exchange_boundary_z_surface(f)

END SUBROUTINE Parallel_Exchange_Surface

SUBROUTINE exchange_boundary_x_surface(f)
  use SF_CFD_Global, only:nx_local,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real( kind = rk )::f(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  integer:: ic,kc,k1

  ! Here we perform the data transfer along I direction only
  do kc = 1,nz_local
      do ic = 1,overLAP
        k1 = (kc-1)*overLAP + ic
        Temp_Send1_2D_xi(k1) = f(ic,kc)                  ! left boundary of f
        Temp_Send2_2D_xi(k1) = f(nx_local-overLAP+ic,kc) ! right boundary of f
      enddo
  enddo

  call MPI_Sendrecv(  Temp_Send1_2D_xi,overLAP*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 9001,   &
                      Temp_Recv2_2D_xi,overLAP*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_2D_xi,overLAP*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 8001,   &           
                      Temp_Recv1_2D_xi,overLAP*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 8001,MPI_COMM_WORLD,Status,ierr)   

   if( ID_XM1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local
      do ic = 1,overLAP
       k1 = (kc-1)*overLAP + ic
       f(ic-overLAP,kc) = Temp_Recv1_2D_xi(k1);
      enddo  
    enddo
   endif
  
   if( ID_XP1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local  
      do ic = 1,overLAP
       k1 = (kc-1)*overLAP + ic
       f(nx_local+ic,kc) = Temp_Recv2_2D_xi(k1);
      enddo  
    enddo
   endif
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE exchange_boundary_x_surface

SUBROUTINE exchange_boundary_z_surface(f)
  use SF_CFD_Global, only:nx_local,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real( kind = rk )::f(1-overLap:nx_local+overLap,1-overLap:nz_local+overLap)
  integer:: ic,kc,k1

  ! Here we perform the data transfer along K direction only
  do kc = 1,overLAP
      do ic = 1,nx_local
        k1 = (kc-1)*nx_local + ic
        Temp_Send1_2D_zeta(k1) = f(ic,kc)                  ! left boundary of f
        Temp_Send2_2D_zeta(k1) = f(ic,nz_local-overLAP+kc) ! right boundary of f
      enddo
  enddo

  call MPI_Sendrecv(  Temp_Send1_2D_zeta,nx_local*overLAP,MPI_DOUBLE_PRECISION, ID_ZM1, 9001,   &
                      Temp_Recv2_2D_zeta,nx_local*overLAP,MPI_DOUBLE_PRECISION, ID_ZP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_2D_zeta,nx_local*overLAP,MPI_DOUBLE_PRECISION, ID_ZP1, 8001,   &           
                      Temp_Recv1_2D_zeta,nx_local*overLAP,MPI_DOUBLE_PRECISION, ID_ZM1, 8001,MPI_COMM_WORLD,Status,ierr)   

   if( ID_ZM1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP
      do ic = 1,nx_local
        k1 = (kc-1)*nx_local + ic
       f(ic,kc-overLAP) = Temp_Recv1_2D_zeta(k1);
      enddo  
    enddo
   endif
  
   if( ID_ZP1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP
      do ic = 1,nx_local
        k1 = (kc-1)*nx_local + ic
       f(ic,nz_local+kc) = Temp_Recv2_2D_zeta(k1);
      enddo  
    enddo
   endif
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE exchange_boundary_z_surface

! mpi message send and recv, using user defined data type
! which can greatly simplify the data transfer
SUBROUTINE exchange_boundary_x(f)
  use SF_CFD_Global, only:nx_local,Ny,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
  integer:: ic,jc,kc,k1

  ! Here we perform the data transfer along I direction only
  do kc = 1,nz_local
    do jc = 1,Ny
      do ic = 1,overLAP
        k1 = (kc-1)*Ny*overLAP + (jc-1)*overLAP + ic
        Temp_Send1_3D_xi(k1) = f(ic,jc,kc)                  ! left boundary of f
        Temp_Send2_3D_xi(k1) = f(nx_local-overLAP+ic,jc,kc) ! right boundary of f
      enddo
    enddo
  enddo

  call MPI_Sendrecv(  Temp_Send1_3D_xi,overLAP*Ny*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 9001,   &
                      Temp_Recv2_3D_xi,overLAP*Ny*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_3D_xi,overLAP*Ny*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 8001,   &           
                      Temp_Recv1_3D_xi,overLAP*Ny*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 8001,MPI_COMM_WORLD,Status,ierr)   

   if( ID_XM1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local
     do jc = 1,Ny 
      do ic = 1,overLAP
       k1 = (kc-1)*Ny*overLAP + (jc-1)*overLAP + ic
       f(ic-overLAP,jc,kc) = Temp_Recv1_3D_xi(k1);
      enddo  
     enddo
    enddo
   endif
  
   if( ID_XP1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local  
     do jc = 1,Ny 
      do ic = 1,overLAP
       k1 = (kc-1)*Ny*overLAP + (jc-1)*overLAP + ic
       f(nx_local+ic,jc,kc) = Temp_Recv2_3D_xi(k1);
      enddo  
     enddo
    enddo
   endif
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine exchange_boundary_x

subroutine exchange_boundary_z(f)
   use SF_CFD_Global, only:nx_local,Ny,nz_local
   implicit none
   integer Status(MPI_status_Size)
   real(kind=rk):: f(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
   integer:: ic,jc,kc,k1
   ! this integer is used to represent the total number of grid points in the x direction
   integer:: temp_Int
  
   temp_Int = nx_local+2*overLAP
  ! Here we perform the data transfer along K direction only
  do kc = 1,overLAP
    do jc = 1,Ny
      do ic = 1-overLap,nx_local+overLap      
        k1 = (kc-1)*Ny*temp_Int + (jc-1)*temp_Int + ic + overLAP
        Temp_Send1_3D_zeta(k1) = f(ic,jc,kc)                  ! left boundary of f
        Temp_Send2_3D_zeta(k1) = f(ic,jc,nz_local-overLAP+kc) ! right boundary of f
      enddo
    enddo
  enddo

  call MPI_SendRecv(Temp_Send1_3D_zeta,overLAP*Ny*temp_Int,MPI_DOUBLE_PRECISION,ID_ZM1,9003,&
                    Temp_Recv2_3D_zeta,overLAP*Ny*temp_Int,MPI_DOUBLE_PRECISION,ID_ZP1,9003,MPI_COMM_WORLD,Status,ierr);

  call MPI_SendRecv(Temp_Send2_3D_zeta,overLAP*Ny*temp_Int,MPI_DOUBLE_PRECISION,ID_ZP1,8003,&
                    Temp_Recv1_3D_zeta,overLAP*Ny*temp_Int,MPI_DOUBLE_PRECISION,ID_ZM1,8003,MPI_COMM_WORLD,Status,ierr);


  !call MPI_SendRecv(Temp_Send1_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZM1,9003,&
  !                  Temp_Recv2_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZP1,9003,MPI_COMM_WORLD,Status,ierr);
  !
  !call MPI_SendRecv(Temp_Send2_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZP1,8003,&
  !                  Temp_Recv1_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZM1,8003,MPI_COMM_WORLD,Status,ierr);
  
  if(ID_ZM1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP
     do jc = 1,Ny 
      do ic = 1-overLAP,nx_local+overLAP
       k1 = (kc-1)*Ny*temp_Int + (jc-1)*temp_Int + ic + overLAP
       f(ic,jc,kc - overLAP) = Temp_Recv1_3D_zeta(k1);
      enddo  
     enddo
    enddo
  endif
  
  if(ID_ZP1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP  
     do jc = 1,Ny 
      do ic = 1-overLAP,nx_local+overLAP
       k1 = (kc-1)*Ny*temp_Int + (jc-1)*temp_Int + ic + overLAP
       f(ic,jc,kc+nz_local) = Temp_Recv2_3D_zeta(k1);
      enddo  
     enddo
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)                  
end subroutine exchange_boundary_z

SUBROUTINE exchange_boundary_x_interp(f)
  use SF_CFD_Global, only:nx_local,Ny,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1-overLAP:nx_local+overLAP,1:2*Ny-1,1-overLAP:nz_local+overLAP)
  integer:: ic,jc,kc,k1

  ! Here we perform the data transfer along I direction only
  do kc = 1,nz_local
    do jc = 1,2*Ny-1
      do ic = 1,overLAP
        k1 = (kc-1)*(2*Ny-1)*overLAP + (jc-1)*overLAP + ic
        Temp_Send1_3D_xi_interp(k1) = f(ic,jc,kc)                  ! left boundary of f
        Temp_Send2_3D_xi_interp(k1) = f(nx_local-overLAP+ic,jc,kc) ! right boundary of f
      enddo
    enddo
  enddo

  call MPI_Sendrecv(  Temp_Send1_3D_xi_interp,overLAP*(2*Ny-1)*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 9001,   &
                      Temp_Recv2_3D_xi_interp,overLAP*(2*Ny-1)*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_3D_xi_interp,overLAP*(2*Ny-1)*nz_local,MPI_DOUBLE_PRECISION, ID_XP1, 8001,   &           
                      Temp_Recv1_3D_xi_interp,overLAP*(2*Ny-1)*nz_local,MPI_DOUBLE_PRECISION, ID_XM1, 8001,MPI_COMM_WORLD,Status,ierr)   

   if( ID_XM1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local
     do jc = 1,2*Ny-1 
      do ic = 1,overLAP
       k1 = (kc-1)*(2*Ny-1)*overLAP + (jc-1)*overLAP + ic
       f(ic-overLAP,jc,kc) = Temp_Recv1_3D_xi_interp(k1);
      enddo  
     enddo
    enddo
   endif
  
   if( ID_XP1 .NE. MPI_PROC_NULL) then
    do kc = 1,nz_local  
     do jc = 1,2*Ny-1
      do ic = 1,overLAP
       k1 = (kc-1)*(2*Ny-1)*overLAP + (jc-1)*overLAP + ic
       f(nx_local+ic,jc,kc) = Temp_Recv2_3D_xi_interp(k1);
      enddo  
     enddo
    enddo
   endif
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine exchange_boundary_x_interp

subroutine exchange_boundary_z_interp(f)
   use SF_CFD_Global, only:nx_local,Ny,nz_local
   implicit none
   integer Status(MPI_status_Size)
   real(kind=rk):: f(1-overLAP:nx_local+overLAP,1:2*Ny-1,1-overLAP:nz_local+overLAP)
   integer:: ic,jc,kc,k1
   ! this integer is used to represent the total number of grid points in the x direction
   integer:: temp_Int
  
   temp_Int = nx_local+2*overLAP
  ! Here we perform the data transfer along K direction only
  do kc = 1,overLAP
    do jc = 1,2*Ny-1
      do ic = 1-overLap,nx_local+overLap      
        k1 = (kc-1)*(2*Ny-1)*temp_Int + (jc-1)*temp_Int + ic + overLAP
        Temp_Send1_3D_zeta_interp(k1) = f(ic,jc,kc)                  ! left boundary of f
        Temp_Send2_3D_zeta_interp(k1) = f(ic,jc,nz_local-overLAP+kc) ! right boundary of f
      enddo
    enddo
  enddo

  call MPI_SendRecv(Temp_Send1_3D_zeta_interp,overLAP*(2*Ny-1)*temp_Int,MPI_DOUBLE_PRECISION,ID_ZM1,9003,&
                    Temp_Recv2_3D_zeta_interp,overLAP*(2*Ny-1)*temp_Int,MPI_DOUBLE_PRECISION,ID_ZP1,9003,MPI_COMM_WORLD,Status,ierr);

  call MPI_SendRecv(Temp_Send2_3D_zeta_interp,overLAP*(2*Ny-1)*temp_Int,MPI_DOUBLE_PRECISION,ID_ZP1,8003,&
                    Temp_Recv1_3D_zeta_interp,overLAP*(2*Ny-1)*temp_Int,MPI_DOUBLE_PRECISION,ID_ZM1,8003,MPI_COMM_WORLD,Status,ierr);


  !call MPI_SendRecv(Temp_Send1_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZM1,9003,&
  !                  Temp_Recv2_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZP1,9003,MPI_COMM_WORLD,Status,ierr);
  !
  !call MPI_SendRecv(Temp_Send2_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZP1,8003,&
  !                  Temp_Recv1_3D_zeta,overLAP*Ny*nx_local,MPI_DOUBLE_PRECISION,ID_ZM1,8003,MPI_COMM_WORLD,Status,ierr);
  
  if(ID_ZM1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP
     do jc = 1,2*Ny-1 
      do ic = 1-overLAP,nx_local+overLAP
       k1 = (kc-1)*(2*Ny-1)*temp_Int + (jc-1)*temp_Int + ic + overLAP
       f(ic,jc,kc - overLAP) = Temp_Recv1_3D_zeta_interp(k1);
      enddo  
     enddo
    enddo
  endif
  
  if(ID_ZP1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP  
     do jc = 1,2*Ny-1 
      do ic = 1-overLAP,nx_local+overLAP
       k1 = (kc-1)*(2*Ny-1)*temp_Int + (jc-1)*temp_Int + ic + overLAP
       f(ic,jc,kc+nz_local) = Temp_Recv2_3D_zeta_interp(k1);
      enddo  
     enddo
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)                  
end subroutine exchange_boundary_z_interp


subroutine exchange_boundary_x_numvar(f)
  use SF_Constant,       only:NumVar
  use SF_CFD_Global,     only:nx_local,Ny,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  integer:: ic,jc,kc,iVar,k1
  integer:: NumInt_MPI

  ! Here we perform the data transfer along I direction only
  do kc = 1,nz_local
    do jc = 1,Ny
      do ic = 1,overLAP
        do iVar = 1,NumVar
        k1 = (kc-1)*Ny*overLAP*NumVar + (jc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
        Temp_Send1_3D_xi_numvar(k1) = f(iVar,ic,jc,kc)                  ! left boundary of f
        Temp_Send2_3D_xi_numvar(k1) = f(iVar,nx_local-overLAP+ic,jc,kc) ! right boundary of f
        enddo
      enddo
    enddo
  enddo
   
  ! Here we count the total number of data used in MPI_SendRecv 
  NumInt_MPI = overLAP*Ny*nz_local*NumVar

  call MPI_Sendrecv(  Temp_Send1_3D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XM1, 9001,   &
                      Temp_Recv2_3D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XP1, 9001,MPI_COMM_WORLD,Status,ierr)   
  call MPI_Sendrecv(  Temp_Send2_3D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XP1, 8001,   &           
                      Temp_Recv1_3D_xi_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION, ID_XM1, 8001,MPI_COMM_WORLD,Status,ierr)   
  
    if( ID_XM1 .NE. MPI_PROC_NULL) then
      do kc = 1, nz_local
        do jc = 1, Ny
          do ic = 1, overLAP
            do iVar = 1, NumVar
              k1 = (kc-1)*Ny*overLAP*NumVar + (jc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
              f(iVar,ic-overLAP,jc,kc) = Temp_Recv1_3D_xi_numvar(k1)
            enddo
          enddo
        enddo
      enddo
    endif

    if( ID_XP1 .NE. MPI_PROC_NULL) then
      do kc = 1, nz_local
        do jc = 1, Ny
          do ic = 1, overLAP
            do iVar = 1, NumVar
              k1 = (kc-1)*Ny*overLAP*NumVar + (jc-1)*overLAP*NumVar + (ic-1)*NumVar + iVar
              f(iVar,nx_local+ic,jc,kc) = Temp_Recv2_3D_xi_numvar(k1)
            enddo
          enddo
        enddo
      enddo
    endif

end subroutine exchange_boundary_x_numvar

subroutine exchange_boundary_z_numvar(f)
  use SF_Constant,       only:NumVar
  use SF_CFD_Global,     only:nx_local,Ny,nz_local
  implicit none
  integer Status(MPI_status_Size)
  real(kind=rk):: f(1:NumVar,1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap)
  integer:: ic,jc,kc,iVar,k1
  integer:: NumInt_MPI

  ! Here we perform the data transfer along K direction only
  do kc = 1,overLAP
    do jc = 1,Ny
      do ic = 1-overLap,nx_local+overLap      
        do iVar = 1,NumVar
          k1 = (kc-1)*Ny*(nx_local+2*overLap)*NumVar &
            &+ (jc-1)*(nx_local+2*overLap)*NumVar &
            &+ (ic-1+overLAP)*NumVar &
            &+ iVar
          Temp_Send1_3D_zeta_numvar(k1) = f(iVar,ic,jc,kc)                  ! left boundary of f
          Temp_Send2_3D_zeta_numvar(k1) = f(iVar,ic,jc,nz_local-overLAP+kc) ! right boundary of f
        enddo
      enddo
    enddo
  enddo

  ! Here we count the total number of data used in MPI_SendRecv 
  NumInt_MPI = overLAP*Ny*(nx_local+2*overLap)*NumVar

  call MPI_SendRecv(Temp_Send1_3D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION,ID_ZM1,9003,&
                    Temp_Recv2_3D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION,ID_ZP1,9003,MPI_COMM_WORLD,Status,ierr);

  call MPI_SendRecv(Temp_Send2_3D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION,ID_ZP1,8003,&
                    Temp_Recv1_3D_zeta_numvar,NumInt_MPI,MPI_DOUBLE_PRECISION,ID_ZM1,8003,MPI_COMM_WORLD,Status,ierr);

  if(ID_ZM1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP
     do jc = 1,Ny 
      do ic = 1-overLAP,nx_local+overLAP
       do iVar = 1,NumVar
         k1 = (kc-1)*Ny*(nx_local+2*overLap)*NumVar &
              &+ (jc-1)*(nx_local+2*overLap)*NumVar &
              &+ (ic-1+overLAP)*NumVar &
              &+ iVar
         f(iVar,ic,jc,kc - overLAP) = Temp_Recv1_3D_zeta_numvar(k1);
       enddo
      enddo  
     enddo
    enddo
  endif
  
  if(ID_ZP1 .NE. MPI_PROC_NULL) then
    do kc = 1,overLAP  
     do jc = 1,Ny 
      do ic = 1-overLAP,nx_local+overLAP
       do iVar = 1,NumVar
        k1 =   (kc-1)*Ny*(nx_local+2*overLap)*NumVar &
            &+ (jc-1)*(nx_local+2*overLap)*NumVar &
            &+ (ic-1+overLAP)*NumVar &
            &+ iVar
         f(iVar,ic,jc,kc+nz_local) = Temp_Recv2_3D_zeta_numvar(k1);
       enddo
      enddo  
     enddo
    enddo
  endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine exchange_boundary_z_numvar


! Internal function to handle operation
integer function my_mod1(i,n)
  implicit none
  integer,intent(in)::i,n 
  if(i.lt.0) then
   my_mod1 = i + n
  else if(i.gt.n-1) then
   my_mod1 = i - n
  else
   my_mod1 = i
  endif
end function my_mod1


END MODULE MPI_GLOBAL