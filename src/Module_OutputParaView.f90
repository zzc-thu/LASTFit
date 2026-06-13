! Here, for convenience, we put all the variables output by ParaView in a module
! So when we call it elsewhere, we only need to call this module
! The advantage of this is that we can put all the output variables in one place for easy management
Module OutputParaView
  ! define the precisions
  use SF_Constant, only: ik,rk 

  ! define the output dimensional variables from global variables
  use SF_CFD_Global,only: Nx,Ny,Nz,nx_local,nz_local,OutputFileNo
  ! define the output array variables from global variables
  use SF_CFD_Global,only: X_grid,Y_grid,Z_grid,x_grid_steady,y_grid_steady,z_grid_steady,&
                          x_grid_interp,y_grid_interp,z_grid_interp
  use SF_CFD_Global,only: U,V,W,P,T,Rho,U_pert,V_pert,W_pert,P_pert,T_pert,Rho_pert,&
                          U_inte,V_inte,W_inte,P_inte,T_inte,Rho_inte              
  use SF_CFD_Global,only: dxdxi,dydxi,dzdxi,dxdeta,dydeta,dzdeta,dxdzeta,dydzeta,dzdzeta,Jaco,&
                          dxidx,dxidy,dxidz,detadx,detady,detadz,dzetadx,dzetady,dzetadz,&
                          dUcons

  ! define the MPI related variables
  use MPI_Global,   only: MyId,npx0,npz0,npx,npz,&
                          I_nn,K_nn,I_offset,K_offset,If_parallel
  implicit none
  
  integer( kind = ik ),parameter:: npy0 = 1
  integer( kind = ik ),parameter:: npy = 0

  contains 

  subroutine output_initial_grid
      implicit none
      character(128) :: filename_output_vts,filename_pvts
      logical :: dir_exists
      integer(kind = ik)::global_Istart,global_Iend
      integer(kind = ik)::global_Jstart,global_Jend
      integer(kind = ik)::global_Kstart,global_Kend
      integer(kind = ik)::i_local_start,i_local_end
      integer(kind = ik)::j_local_start,j_local_end
      integer(kind = ik)::k_local_start,k_local_end
  
      ! output the initial grid .pvts file
      if(myid == 0) then
       !inquire(file="INIT", exist=dir_exists)
       !if (.not. dir_exists) then
       !  call system("mkdir -p INIT")
       !endif
  
       filename_pvts = "INIT/Initial_grid.pvts"
       call output_global_pvts_initial_grid(filename_pvts)
       !    
      endif
      
      if(myid == 0) then
       write(*,*)"Output the initial grid .vts files"
      endif
  
      global_Istart = i_offset(npx);
      global_Iend   = i_offset(npx) + i_nn(npx) - 1;
    
      global_Jstart = 1_ik;
      global_Jend   = Ny;
    
      global_Kstart = k_offset(npz);
      global_Kend   = k_offset(npz) + k_nn(npz) - 1;
    
    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Iend   = global_Iend + 1;   
       i_local_start = 1; 
       i_local_end = nx_local + 1;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = Ny;  
   
     if(npz == 0) then
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       k_local_start = 0; 
       k_local_end = nz_local;
     else
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;  
     endif
    else
      i_local_start = 1; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = Ny; 
      k_local_start = 1; 
      k_local_end = nz_local;
    endif
   
    write(filename_output_vts,'("INIT/INIT_GridPart",I5.5,".vts")')myid
    
    CALL output_vts_files_initial_grid_local_core(filename_output_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  end subroutine output_initial_grid

  subroutine output_global_pvts_initial_grid(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float

     !vtk_float = "Float64"
     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"

     open(unit=123, file=filename_pvts, status="replace")
     
     xml_part = ' &
     & <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <PStructuredGrid WholeExtent="1 '//int2str(Nx)//' 1 '//int2str(Ny)//' 1 '//int2str(Nz)//'" GhostLevel="1"> &
     &   <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &   </PPoints> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;

     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;
            
            Piece_SourceFile_vts = "";
            write(Piece_SourceFile_vts,'("INIT_GridPart",I5.5,".vts")')local_proc_num

            global_istart = i_offset(ic);
            global_iend   = i_offset(ic) + i_nn(ic) - 1;
  
            !global_jstart = j_offset(jc);
            !global_jend   = j_offset(jc) + j_nn(jc) - 1;

            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = Ny;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;
           
          if(if_parallel) then  
           if(ic == 0 ) then
             global_Iend   = global_Iend + 1;  
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend + 1;  
           endif

           !if(jc == 0 ) then
           !  global_Jend   = global_Jend; 
           !elseif(jc == npy0 - 1) then
           !  global_Jstart = global_Jstart - 1; 
           !else
           !  global_Jstart = global_Jstart - 1;
           !  global_Jend   = global_Jend + 1; 
           !endif
           global_jstart = 1;
           global_jend   = Ny;

           if(kc == 0 ) then
             global_Kend   = global_Kend + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
           else
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend + 1; 
           endif
          else
           ! global_Istart = global_Istart;
           ! global_Iend   = global_Iend;
           ! global_Jstart = global_Jstart;
           ! global_Jend   = global_Jend;
           ! global_Kstart = global_Kstart;
           ! global_Kend   = global_Kend;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '

         enddo
     enddo

     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)

  end subroutine output_global_pvts_initial_grid

  ! output the local .vts file
  subroutine output_vts_files_initial_grid_local_core(filename_output_vts,i_local_start,i_local_end,&
                                                 j_local_start,j_local_end,k_local_start,k_local_end,&
                                                 istart,iend,jstart,jend,kstart,kend)
    use mpi, only:MPI_OFFSET_KIND
    implicit none
    integer::i,j,k
    integer::icounts
    integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
    integer::loc_i_length,loc_j_length,loc_k_length 
    integer::istart,iend,jstart,jend,kstart,kend
    integer, parameter :: int64_kind = selected_int_kind(2*range(1))
    integer(int64_kind) :: gridsize_64
    integer::size_real
    integer (KIND=MPI_OFFSET_KIND) ::offset_x,delta_offset_w
    
    character(len=65536) :: xml_part
    character(len=7) :: vtk_float
    character(128)::filename_output_vts
 
    real(4),allocatable::pointData2(:,:)
 
    loc_i_length = i_local_end - i_local_start + 1;
    loc_j_length = j_local_end - j_local_start + 1;
    loc_k_length = k_local_end - k_local_start + 1;
     
    allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
 
    pointData2 = 0.0;
 
    ! single is enough
    size_real = 4;
    vtk_float = "Float32";
 
    !print*,size_real,int64_kind
 
    gridsize_64 = int(size_real,int64_kind) * int(nx,int64_kind) * int(ny,int64_kind) * int(nz,int64_kind)
    delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
    offset_x = 0
 
    open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
    ! forming the xml_part files
    xml_part = ' &
    & <?xml version="1.0"?> &
    & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
    &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
    &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
    &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
    &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
    &    <Points> &
    &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
    &    </Points> '
 
    xml_part = trim(adjustl(xml_part)) // ' &
    &    </Piece> &
    &   </StructuredGrid> &
    &   <AppendedData encoding="raw"> '
 
    ! the xml_part is finished
    !write(*,*)"write head"
      write(1231) trim(adjustl(xml_part))
    !write(*,*)"write xx,yy,zz"
      write(1231) "_"
      write(1231) 3*size_real*int(gridsize_64,int64_kind)
      icounts = 1
      do k = k_local_start,k_local_end
        do j = j_local_start,j_local_end
          do i = i_local_start,i_local_end
            pointData2(1,icounts) = real(X_grid(i,j,k),kind=4)
            pointData2(2,icounts) = real(Y_grid(i,j,k),kind=4)
            pointData2(3,icounts) = real(Z_grid(i,j,k),kind=4)
            icounts = icounts + 1     
          enddo
        enddo
      enddo
       ! Write the all grid points data
      write(1231)pointData2
 
    write(1231) ' &
     &    </AppendedData> &
     &  </VTKFile>' 
    close(1231)
 
    deallocate(pointData2)
   
  end subroutine output_vts_files_initial_grid_local_core  
  
  subroutine output_results
    implicit none
      character(128) :: filename_output_vts,filename_pvts
      logical :: dir_exists
      integer(kind = ik)::global_Istart,global_Iend
      integer(kind = ik)::global_Jstart,global_Jend
      integer(kind = ik)::global_Kstart,global_Kend
      integer(kind = ik)::i_local_start,i_local_end
      integer(kind = ik)::j_local_start,j_local_end
      integer(kind = ik)::k_local_start,k_local_end
      integer(kind = ik)::cmd_status
  
      ! output the initial grid .pvts file
      if(myid == 0) then
       inquire(file="RESU", exist=dir_exists)
       if (.not. dir_exists) then
         call system("mkdir -p RESU")
       endif

       filename_pvts = "";
       write(filename_pvts,'("RESU/Result_",I8.8,".pvts")')OutputFileNo
       call output_global_pvts_resu(filename_pvts)    
      endif
      
      if(myid == 0) then
       write(*,*)"Output the result .vts files"
      endif
  
      global_Istart = i_offset(npx);
      global_Iend   = i_offset(npx) + i_nn(npx) - 1;
    
      global_Jstart = 1_ik;
      global_Jend   = Ny;
    
      global_Kstart = k_offset(npz);
      global_Kend   = k_offset(npz) + k_nn(npz) - 1;
    
    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Iend   = global_Iend + 1;
       i_local_start = 1; 
       i_local_end = nx_local + 1;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = Ny;  
   
     if(npz == 0) then
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       k_local_start = 0; 
       k_local_end = nz_local;
     else
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;  
     endif
    else
      i_local_start = 1; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = Ny; 
      k_local_start = 1; 
      k_local_end = nz_local+1;
    endif
   
    write(filename_output_vts,'("RESU/RESU_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,myid

    CALL output_vts_files_local_core(filename_output_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  end subroutine output_results

   subroutine output_results_singularity
      implicit none
      character(128) :: filename_output_vts,filename_pvts
      logical :: dir_exists
      integer(kind = ik)::global_Istart,global_Iend
      integer(kind = ik)::global_Jstart,global_Jend
      integer(kind = ik)::global_Kstart,global_Kend
      integer(kind = ik)::i_local_start,i_local_end
      integer(kind = ik)::j_local_start,j_local_end
      integer(kind = ik)::k_local_start,k_local_end
      integer(kind = ik)::cmd_status
  
      ! output the initial grid .pvts file
      if(myid == 0) then
       inquire(file="RESU", exist=dir_exists)
       if (.not. dir_exists) then
         call system("mkdir -p RESU")
       endif
  
       filename_pvts = "";
       write(filename_pvts,'("RESU/Result_",I8.8,".pvts")')OutputFileNo
       call output_global_pvts_resu_singularity(filename_pvts)    
      endif
      
      if(myid == 0) then
       write(*,*)"Output the result of singularity .vts files"
      endif
  
      global_Istart = i_offset(npx);
      global_Iend   = i_offset(npx) + i_nn(npx) - 1;
      
      global_Jstart = 1_ik;
      global_Jend   = Ny;
    
      global_Kstart = k_offset(npz);
      global_Kend   = k_offset(npz) + k_nn(npz) - 1;
    
    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend;   
       i_local_start = 0; 
       i_local_end = nx_local ;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = Ny;  
   
     if(npz == 0) then
       global_Kstart = global_Kstart;
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;
     else
       global_Kstart = global_Kstart-1 ;
       global_Kend   = global_Kend+1 ;
       k_local_start = 0; 
       k_local_end = nz_local+1 ;  
     endif
    else
      global_Istart = global_Istart-1
      global_Iend   = global_Iend  
      global_Jstart = global_Jstart 
      global_Jend   = global_Jend
      global_Kstart = global_Kstart 
      global_Kend   = global_Kend + 1
      
      i_local_start = 0; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = Ny; 
      k_local_start = 1; 
      k_local_end = nz_local+1;
    endif
   
    write(filename_output_vts,'("RESU/RESU_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,myid
    
    CALL output_vts_files_local_core_singularity(filename_output_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  end subroutine output_results_singularity

  
  subroutine output_Pert_results
    implicit none
      character(128) :: filename_output_vts,filename_pvts
      logical :: dir_exists
      integer(kind = ik)::global_Istart,global_Iend
      integer(kind = ik)::global_Jstart,global_Jend
      integer(kind = ik)::global_Kstart,global_Kend
      integer(kind = ik)::i_local_start,i_local_end
      integer(kind = ik)::j_local_start,j_local_end
      integer(kind = ik)::k_local_start,k_local_end
      integer(kind = ik)::cmd_status
  
      ! output the initial grid .pvts file
      if(myid == 0) then
      ! inquire(file="Pert", exist=dir_exists)
      ! if (.not. dir_exists) then
      !   call system("mkdir -p Pert")
      ! endif

       filename_pvts = "";
       write(filename_pvts,'("Pert/Pert_Result_",I8.8,".pvts")')OutputFileNo
       call output_global_pvts_Pert_resu(filename_pvts)    
      endif
      
      !if(myid == 0) then
      ! write(*,*)"Output the Perturation result .vts files"
      !endif
  
      global_Istart = i_offset(npx);
      global_Iend   = i_offset(npx) + i_nn(npx) - 1;
    
      global_Jstart = 1_ik;
      global_Jend   = 2*Ny-1;
    
      global_Kstart = k_offset(npz);
      global_Kend   = k_offset(npz) + k_nn(npz) - 1;
    
    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Iend   = global_Iend + 1;   
       i_local_start = 1; 
       i_local_end = nx_local + 1;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = 2*Ny-1;  
   
     if(npz == 0) then
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       k_local_start = 0; 
       k_local_end = nz_local;
     else
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;  
     endif
    else
      i_local_start = 1; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = 2*Ny-1; 
      k_local_start = 1; 
      k_local_end = nz_local;
    endif
   
    write(filename_output_vts,'("Pert/Pert_Result_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,myid

    CALL output_vts_files_Pert_local_core(filename_output_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  end subroutine output_Pert_results
  
   subroutine output_LNSResults

   implicit none
      character(128) :: filename_output_vts,filename_pvts
      logical :: dir_exists
      integer(kind = ik)::global_Istart,global_Iend
      integer(kind = ik)::global_Jstart,global_Jend
      integer(kind = ik)::global_Kstart,global_Kend
      integer(kind = ik)::i_local_start,i_local_end
      integer(kind = ik)::j_local_start,j_local_end
      integer(kind = ik)::k_local_start,k_local_end
      integer(kind = ik)::cmd_status
  
      ! output the initial grid .pvts file
      if(myid == 0) then
       inquire(file="LNSResults", exist=dir_exists)

       filename_pvts = "";
       write(filename_pvts,'("LNSResults/LNSResults_",I8.8,".pvts")')OutputFileNo
       call output_global_pvts_LNSResults(filename_pvts)
      endif
      
      if(myid == 0) then
       write(*,*)"Output the LNS result .vts files"
      endif
  
      global_Istart = i_offset(npx);
      global_Iend   = i_offset(npx) + i_nn(npx) - 1;
    
      global_Jstart = 1_ik;
      global_Jend   = 2*Ny-1;
    
      global_Kstart = k_offset(npz);
      global_Kend   = k_offset(npz) + k_nn(npz) - 1;
    
    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Iend   = global_Iend + 1;   
       i_local_start = 1; 
       i_local_end = nx_local + 1;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = 2*Ny-1;  
   
     if(npz == 0) then
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       k_local_start = 0; 
       k_local_end = nz_local;
     else
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;  
     endif
    else
      i_local_start = 1; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = 2*Ny-1; 
      k_local_start = 1; 
      k_local_end = nz_local;
    endif
    
    write(filename_output_vts,'("LNSResults/LNSResults_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,myid

    CALL output_vts_files_LNSResults_local_core(filename_output_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  
  end subroutine output_LNSResults
  
  ! output the global .pvts file
  subroutine output_global_pvts_resu(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float

     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"

     open(unit=123, file=filename_pvts, status="replace")
     
     xml_part = ' &
     & <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <PStructuredGrid WholeExtent="1 '//int2str(Nx)//' 1 '//int2str(Ny)//' 1 '//int2str(Nz)//'" GhostLevel="1"> &
     &   <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &   </PPoints> &  
     &   <PPointData> &
     &    <PDataArray type="'//vtk_float//'" Name="Rho" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="U" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="V" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="W" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="P" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="T" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU1" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU2" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU3" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU4" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU5" format="appended"/> &     
     &   </PPointData> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;

     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;

            write(Piece_SourceFile_vts,'("RESU_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,local_proc_num
           
            global_istart = i_offset(ic);
            global_iend   = i_offset(ic) + i_nn(ic) - 1;
  
            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = Ny;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;

          if(If_parallel) then 
           if(ic == 0 ) then
             global_Iend   = global_Iend + 1;
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend + 1;  
           endif

           if(kc == 0 ) then
             global_Kend   = global_Kend + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
           else
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend + 1; 
           endif
          else
            global_Istart = global_Istart;
            global_Iend   = global_Iend;
            global_Jstart = global_Jstart;
            global_Jend   = global_Jend;
            global_Kstart = global_Kstart;
            global_Kend   = global_Kend;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '

         enddo
     enddo

     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)

  end subroutine output_global_pvts_resu

  subroutine output_global_pvts_resu_singularity(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
  
     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"
  
     open(unit=123, file=filename_pvts, status="replace")

     xml_part = ' &
     &    <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &    <PStructuredGrid WholeExtent="0 '//int2str(Nx)//' 1 '//int2str(Ny)//' 1 '//int2str(Nz+1)//'" GhostLevel="1"> &
     &    <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &    </PPoints> &  
     &    <PPointData> &
     &    <PDataArray type="'//vtk_float//'" Name="Rho" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="U" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="V" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="W" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="P" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="T" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU1" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU2" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU3" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU4" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dU5" format="appended"/> &     
     &    </PPointData> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;
  
     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;
  
            write(Piece_SourceFile_vts,'("RESU_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,local_proc_num
           
            global_istart = i_offset(ic);
            global_iend   = i_offset(ic) + i_nn(ic) - 1;

            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = Ny;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;
  
          if(If_parallel) then 
           if(ic == 0 ) then
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend  ;  
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
             global_Iend   = global_Iend;
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend   + 1;  
           endif
  
           if(kc == 0 ) then
             global_Kend   = global_Kend   + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend   + 1; 
           else
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend   + 1; 
           endif
          else
            global_Istart = global_Istart - 1;
            global_Iend   = global_Iend ;
            global_Jstart = global_Jstart;
            global_Jend   = global_Jend;
            global_Kstart = global_Kstart;
            global_Kend   = global_Kend + 1;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '
  
         enddo
     enddo
  
     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)

  end subroutine output_global_pvts_resu_singularity
  
  
  subroutine output_global_pvts_Pert_resu(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float

     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"

     open(unit=123, file=filename_pvts, status="replace")
     
     xml_part = ' &
     & <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <PStructuredGrid WholeExtent="1 '//int2str(Nx)//' 1 '//int2str(2*Ny-1)//' 1 '//int2str(Nz)//'" GhostLevel="1"> &
     &   <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &   </PPoints> &  
     &   <PPointData> &
     &    <PDataArray type="'//vtk_float//'" Name="Rho_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="U_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="V_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="W_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="P_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="T_pert" format="appended"/> &
     &   </PPointData> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;

     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;

            write(Piece_SourceFile_vts,'("Pert_Result_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,local_proc_num
           
            global_istart = i_offset(ic);
            global_iend   = i_offset(ic) + i_nn(ic) - 1;
  
            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = 2*Ny-1;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;

          if(If_parallel) then 
           if(ic == 0 ) then
             global_Iend   = global_Iend + 1;  
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend + 1;  
           endif

           if(kc == 0 ) then
             global_Kend   = global_Kend + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
           else
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend + 1; 
           endif
          else
            global_Istart = global_Istart;
            global_Iend   = global_Iend;
            global_Jstart = global_Jstart;
            global_Jend   = global_Jend;
            global_Kstart = global_Kstart;
            global_Kend   = global_Kend;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '

         enddo
     enddo

     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)

  end subroutine output_global_pvts_Pert_resu
  
  subroutine output_global_pvts_LNSResults(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float

     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"

     open(unit=123, file=filename_pvts, status="replace")
     
     xml_part = ' &
     & <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <PStructuredGrid WholeExtent="1 '//int2str(Nx)//' 1 '//int2str(2*Ny-1)//' 1 '//int2str(Nz)//'" GhostLevel="1"> &
     &   <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &   </PPoints> &  
     &   <PPointData> &
     &    <PDataArray type="'//vtk_float//'" Name="Rho_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="U_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="V_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="W_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="P_pert" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="T_pert" format="appended"/> &
     &   </PPointData> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;

     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;

            write(Piece_SourceFile_vts,'("LNSResults_Part_",I8.8,"_",I5.5,".vts")')OutputFileNo,local_proc_num
           
            global_istart = i_offset(ic);    
            global_iend   = i_offset(ic) + i_nn(ic) - 1;
  
            !global_jstart = j_offset(jc);
            !global_jend   = j_offset(jc) + j_nn(jc) - 1;

            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = 2*Ny-1;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;

          if(If_parallel) then 
           if(ic == 0 ) then
             global_Iend   = global_Iend + 1;  
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend + 1;  
           endif

           !if(jc == 0 ) then
           !  global_Jend   = global_Jend + 1; 
           !elseif(jc == npy0 - 1) then
           !  global_Jstart = global_Jstart - 1; 
           !else
           !  global_Jstart = global_Jstart - 1;
           !  global_Jend   = global_Jend + 1; 
           !endif

           if(kc == 0 ) then
             global_Kend   = global_Kend + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
           else
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend + 1; 
           endif
          else
            global_Istart = global_Istart;
            global_Iend   = global_Iend;
            global_Jstart = global_Jstart;
            global_Jend   = global_Jend;
            global_Kstart = global_Kstart;
            global_Kend   = global_Kend;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '

         enddo
     enddo

     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)

  end subroutine output_global_pvts_LNSResults
  
  
  ! output the local .vts file
  subroutine output_vts_files_local_core(filename_output_vts,i_local_start,i_local_end,&
                                     j_local_start,j_local_end,k_local_start,k_local_end,&
                                                     istart,iend,jstart,jend,kstart,kend)
     use mpi, only:MPI_OFFSET_KIND
     implicit none
     integer::i,j,k
     integer::icounts
     integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
     integer::loc_i_length,loc_j_length,loc_k_length 
     integer::istart,iend,jstart,jend,kstart,kend
     integer, parameter :: int64_kind = selected_int_kind(2*range(1))
     integer(int64_kind) :: gridsize_64
     integer::size_real
     integer (KIND=MPI_OFFSET_KIND) :: offset,offset_x,delta_offset_w
     
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
     character(128)::filename_output_vts
  
     real(4),allocatable::pointData2(:,:)
     real(4),allocatable::tempVar(:,:,:)
  
     loc_i_length = i_local_end - i_local_start + 1;
     loc_j_length = j_local_end - j_local_start + 1;
     loc_k_length = k_local_end - k_local_start + 1;
      
     allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
     allocate(tempVar(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end))
  
     pointData2 = 0.0;
     tempVar = 0.0;
  
     ! single is enough
     size_real = 4;
     vtk_float = "Float32";
  
     !print*,size_real,int64_kind
  
     !gridsize_64 = int(size_real,int64_kind) * int(nx,int64_kind) * int(ny,int64_kind) * int(nz,int64_kind)

     gridsize_64 = int(size_real,int64_kind) * int(loc_i_length,int64_kind) &
                                          &  * int(loc_j_length,int64_kind) &
                                          &  * int(loc_k_length,int64_kind)

     delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
     offset_x = 0
  
     open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
     ! forming the xml_part files
     xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
     &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
     &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &    <Points> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
     &    </Points> &
     &   <PointData> '
  
     offset = offset_x + 3 * gridsize_64 + storage_size(gridsize_64)/8 
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Rho" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="U" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="V" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="W" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="P" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="T" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU1" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU2" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU3" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU4" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU5" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
     xml_part = trim(adjustl(xml_part)) // ' &
     &       </PointData> &
     &     </Piece> &
     &   </StructuredGrid> &
     &   <AppendedData encoding="raw"> '
  
     ! the xml_part is finished
     !write(*,*)"write head"
       write(1231) trim(adjustl(xml_part))
     !write(*,*)"write xx,yy,zz"
       write(1231) "_"
       write(1231) 3*size_real*int(gridsize_64,int64_kind)
       icounts = 1
       do k = k_local_start,k_local_end
         do j = j_local_start,j_local_end
           do i = i_local_start,i_local_end
             pointData2(1,icounts) = real(X_grid(i,j,k),kind=4)
             pointData2(2,icounts) = real(Y_grid(i,j,k),kind=4)
             pointData2(3,icounts) = real(Z_grid(i,j,k),kind=4)
              icounts = icounts + 1     
           enddo
         enddo
       enddo
       !write(1231)pointData
       write(1231)pointData2
  
     !write(*,*)"write Rho"  
       write(1231)gridsize_64
       !write(1231)Rho
       tempVar = real(Rho(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write U"
       write(1231)gridsize_64
       !write(1231)U
       tempVar = real(U(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write V"
       write(1231)gridsize_64
       !write(1231)V
       tempVar = real(V(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write W"
       write(1231)gridsize_64
       !write(1231)W
       tempVar = real(W(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write P"
       write(1231)gridsize_64
       !write(1231)P
       tempVar = real(P(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write T"
       write(1231)gridsize_64
       !write(1231)T
       tempVar = real(T(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dU1"
       write(1231)gridsize_64
       tempVar = real(dUcons(1,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar       

     !write(*,*)"write dU2"
       write(1231)gridsize_64
       tempVar = real(dUcons(2,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dU3"
       write(1231)gridsize_64
       tempVar = real(dUcons(3,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dU4"
       write(1231)gridsize_64
       tempVar = real(dUcons(4,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dU5"
       write(1231)gridsize_64
       tempVar = real(dUcons(5,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     write(1231) ' &
      &    </AppendedData> &
      &  </VTKFile>' 
     close(1231)
     
     deallocate(pointData2)
     deallocate(tempVar)
  
   end subroutine output_vts_files_local_core
  
 subroutine output_vts_files_local_core_singularity(filename_output_vts,i_local_start,i_local_end,&
                                     j_local_start,j_local_end,k_local_start,k_local_end,&
                                                     istart,iend,jstart,jend,kstart,kend)
     use mpi, only:MPI_OFFSET_KIND
     implicit none
     integer::i,j,k
     integer::icounts
     integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
     integer::loc_i_length,loc_j_length,loc_k_length 
     integer::istart,iend,jstart,jend,kstart,kend
     integer, parameter :: int64_kind = selected_int_kind(2*range(1))
     integer(int64_kind) :: gridsize_64
     integer::size_real
     integer (KIND=MPI_OFFSET_KIND) :: offset,offset_x,delta_offset_w
     
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
     character(128)::filename_output_vts
  
     real(4),allocatable::pointData2(:,:)
     real(4),allocatable::tempVar(:,:,:)
  
     loc_i_length = i_local_end - i_local_start + 1;
     loc_j_length = j_local_end - j_local_start + 1;
     loc_k_length = k_local_end - k_local_start + 1;
      
     allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
     allocate(tempVar(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end))
  
     pointData2 = 0.0;
     tempVar = 0.0;
  
     ! single is enough
     size_real = 4;
     vtk_float = "Float32";
  
     !print*,size_real,int64_kind
  
     !gridsize_64 = int(size_real,int64_kind) * int(nx,int64_kind) * int(ny,int64_kind) * int(nz,int64_kind)
 
     gridsize_64 = int(size_real,int64_kind) * int(loc_i_length,int64_kind) &
                                          &  * int(loc_j_length,int64_kind) &
                                          &  * int(loc_k_length,int64_kind)
 
     delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
     offset_x = 0
  
     open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
     ! forming the xml_part files
     xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
     &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
     &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &    <Points> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
     &    </Points> &
     &   <PointData> '
  
     offset = offset_x + 3 * gridsize_64 + storage_size(gridsize_64)/8 
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Rho" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="U" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="V" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="W" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="P" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="T" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
 
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU1" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
 
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU2" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
 
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU3" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
 
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU4" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
 
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dU5" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
     xml_part = trim(adjustl(xml_part)) // ' &
     &       </PointData> &
     &     </Piece> &
     &   </StructuredGrid> &
     &   <AppendedData encoding="raw"> '
  
     ! the xml_part is finished
     !write(*,*)"write head"
       write(1231) trim(adjustl(xml_part))
     !write(*,*)"write xx,yy,zz"
       write(1231) "_"
       write(1231) 3*size_real*int(gridsize_64,int64_kind)
       icounts = 1
       do k = k_local_start,k_local_end
         do j = j_local_start,j_local_end
           do i = i_local_start,i_local_end
             pointData2(1,icounts) = real(X_grid(i,j,k),kind=4)
             pointData2(2,icounts) = real(Y_grid(i,j,k),kind=4)
             pointData2(3,icounts) = real(Z_grid(i,j,k),kind=4)
              icounts = icounts + 1     
           enddo
         enddo
       enddo
       !write(1231)pointData
       write(1231)pointData2
  
     !write(*,*)"write Rho"  
       write(1231)gridsize_64
       !write(1231)Rho
       tempVar = real(Rho(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write U"
       write(1231)gridsize_64
       !write(1231)U
       tempVar = real(U(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write V"
       write(1231)gridsize_64
       !write(1231)V
       tempVar = real(V(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write W"
       write(1231)gridsize_64
       !write(1231)W
       tempVar = real(W(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write P"
       write(1231)gridsize_64
       !write(1231)P
       tempVar = real(P(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write T"
       write(1231)gridsize_64
       !write(1231)T
       tempVar = real(T(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
 
     !write(*,*)"write dU1"
       write(1231)gridsize_64
       tempVar = real(dUcons(1,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar       
 
     !write(*,*)"write dU2"
       write(1231)gridsize_64
       tempVar = real(dUcons(2,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
 
     !write(*,*)"write dU3"
       write(1231)gridsize_64
       tempVar = real(dUcons(3,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
 
     !write(*,*)"write dU4"
       write(1231)gridsize_64
       tempVar = real(dUcons(4,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
 
     !write(*,*)"write dU5"
       write(1231)gridsize_64
       tempVar = real(dUcons(5,i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
 
     write(1231) ' &
      &    </AppendedData> &
      &  </VTKFile>' 
     close(1231)
     
     deallocate(pointData2)
     deallocate(tempVar)
  
  end subroutine output_vts_files_local_core_singularity
                                                     
   subroutine output_vts_files_Pert_local_core(filename_output_vts,i_local_start,i_local_end,&
                                     j_local_start,j_local_end,k_local_start,k_local_end,&
                                                     istart,iend,jstart,jend,kstart,kend)
     use mpi, only:MPI_OFFSET_KIND
     implicit none
     integer::i,j,k
     integer::icounts
     integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
     integer::loc_i_length,loc_j_length,loc_k_length 
     integer::istart,iend,jstart,jend,kstart,kend
     integer, parameter :: int64_kind = selected_int_kind(2*range(1))
     integer(int64_kind) :: gridsize_64
     integer::size_real
     integer (KIND=MPI_OFFSET_KIND) :: offset,offset_x,delta_offset_w
     
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
     character(128)::filename_output_vts
  
     real(4),allocatable::pointData2(:,:)
     real(4),allocatable::tempVar(:,:,:)
  
     loc_i_length = i_local_end - i_local_start + 1;
     loc_j_length = j_local_end - j_local_start + 1;
     loc_k_length = k_local_end - k_local_start + 1;
      
     allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
     allocate(tempVar(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end))
  
     pointData2 = 0.0;
     tempVar = 0.0;
  
     ! single is enough
     size_real = 4;
     vtk_float = "Float32";
  
     !print*,size_real,int64_kind
  
     !gridsize_64 = int(size_real,int64_kind) * int(nx,int64_kind) * int(ny,int64_kind) * int(nz,int64_kind)

     gridsize_64 = int(size_real,int64_kind) * int(loc_i_length,int64_kind) &
                                          &  * int(loc_j_length,int64_kind) &
                                          &  * int(loc_k_length,int64_kind)

     delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
     offset_x = 0
  
     open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
     ! forming the xml_part files
     xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
     &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
     &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &    <Points> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
     &    </Points> &
     &   <PointData> '
  
     offset = offset_x + 3 * gridsize_64 + storage_size(gridsize_64)/8 
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Rho_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="U_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="V_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="W_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="P_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="T_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
     xml_part = trim(adjustl(xml_part)) // ' &
     &       </PointData> &
     &     </Piece> &
     &   </StructuredGrid> &
     &   <AppendedData encoding="raw"> '
  
     ! the xml_part is finished
     !write(*,*)"write head"
       write(1231) trim(adjustl(xml_part))
     !write(*,*)"write xx,yy,zz"
       write(1231) "_"
       write(1231) 3*size_real*int(gridsize_64,int64_kind)
       icounts = 1
       do k = k_local_start,k_local_end
         do j = j_local_start,j_local_end
           do i = i_local_start,i_local_end
             pointData2(1,icounts) = real(X_grid_interp(i,j,k),kind=4)
             pointData2(2,icounts) = real(Y_grid_interp(i,j,k),kind=4)
             pointData2(3,icounts) = real(Z_grid_interp(i,j,k),kind=4)
              icounts = icounts + 1     
           enddo
         enddo
       enddo
       !write(1231)pointData
       write(1231)pointData2
  
     !write(*,*)"write Rho"  
       write(1231)gridsize_64
       !write(1231)Rho
       tempVar = real(Rho_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write U"
       write(1231)gridsize_64
       !write(1231)U
       tempVar = real(U_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write V"
       write(1231)gridsize_64
       !write(1231)V
       tempVar = real(V_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write W"
       write(1231)gridsize_64
       !write(1231)W
       tempVar = real(W_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write P"
       write(1231)gridsize_64
       !write(1231)P
       tempVar = real(P_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write T"
       write(1231)gridsize_64
       !write(1231)T
       tempVar = real(T_pert(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     write(1231) ' &
      &    </AppendedData> &
      &  </VTKFile>' 
     close(1231)
     
     deallocate(pointData2)
     deallocate(tempVar)
  
   end subroutine output_vts_files_Pert_local_core                                                  
                                                     
    subroutine output_vts_files_LNSResults_local_core(filename_output_vts,i_local_start,i_local_end,&
                                     j_local_start,j_local_end,k_local_start,k_local_end,&
                                                     istart,iend,jstart,jend,kstart,kend)
     use mpi, only:MPI_OFFSET_KIND
     implicit none
     integer::i,j,k
     integer::icounts
     integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
     integer::loc_i_length,loc_j_length,loc_k_length 
     integer::istart,iend,jstart,jend,kstart,kend
     integer, parameter :: int64_kind = selected_int_kind(2*range(1))
     integer(int64_kind) :: gridsize_64
     integer::size_real
     integer (KIND=MPI_OFFSET_KIND) :: offset,offset_x,delta_offset_w
     
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
     character(128)::filename_output_vts
  
     real(4),allocatable::pointData2(:,:)
     real(4),allocatable::tempVar(:,:,:)
  
     loc_i_length = i_local_end - i_local_start + 1;
     loc_j_length = j_local_end - j_local_start + 1;
     loc_k_length = k_local_end - k_local_start + 1;
      
     allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
     allocate(tempVar(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end))
  
     pointData2 = 0.0;
     tempVar = 0.0;
  
     ! single is enough
     size_real = 4;
     vtk_float = "Float32";
  
     !print*,size_real,int64_kind
  
     !gridsize_64 = int(size_real,int64_kind) * int(nx,int64_kind) * int(ny,int64_kind) * int(nz,int64_kind)

     gridsize_64 = int(size_real,int64_kind) * int(loc_i_length,int64_kind) &
                                          &  * int(loc_j_length,int64_kind) &
                                          &  * int(loc_k_length,int64_kind)

     delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
     offset_x = 0
  
     open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
     ! forming the xml_part files
     xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
     &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
     &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &    <Points> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
     &    </Points> &
     &   <PointData> '
  
     offset = offset_x + 3 * gridsize_64 + storage_size(gridsize_64)/8 
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Rho_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="U_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="V_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="W_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="P_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="T_pert" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
     xml_part = trim(adjustl(xml_part)) // ' &
     &       </PointData> &
     &     </Piece> &
     &   </StructuredGrid> &
     &   <AppendedData encoding="raw"> '
  
     ! the xml_part is finished
     !write(*,*)"write head"
       write(1231) trim(adjustl(xml_part))
     !write(*,*)"write xx,yy,zz"
       write(1231) "_"
       write(1231) 3*size_real*int(gridsize_64,int64_kind)
       icounts = 1
       do k = k_local_start,k_local_end
         do j = j_local_start,j_local_end
           do i = i_local_start,i_local_end
             pointData2(1,icounts) = real(X_grid_interp(i,j,k),kind=4)
             pointData2(2,icounts) = real(Y_grid_interp(i,j,k),kind=4)
             pointData2(3,icounts) = real(Z_grid_interp(i,j,k),kind=4)
              icounts = icounts + 1     
           enddo
         enddo
       enddo
       !write(1231)pointData
       write(1231)pointData2
  
     !write(*,*)"write Rho"  
       write(1231)gridsize_64
       !write(1231)Rho
       tempVar = real(Rho_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write U"
       write(1231)gridsize_64
       !write(1231)U
       tempVar = real(U_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write V"
       write(1231)gridsize_64
       !write(1231)V
       tempVar = real(V_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write W"
       write(1231)gridsize_64
       !write(1231)W
       tempVar = real(W_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write P"
       write(1231)gridsize_64
       !write(1231)P
       tempVar = real(P_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write T"
       write(1231)gridsize_64
       !write(1231)T
       tempVar = real(T_inte(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     write(1231) ' &
      &    </AppendedData> &
      &  </VTKFile>' 
     close(1231)
     
     deallocate(pointData2)
     deallocate(tempVar)
  
   end subroutine output_vts_files_LNSResults_local_core                                                                  
                                                     
  ! output the jacobian .pvts file and .vts files
  subroutine output_jacobian
    implicit none
    character(128) :: filename_pvts
    character(128) :: filename_output_jaco_vts    
    logical :: dir_exists
    integer(kind = ik)::global_Istart,global_Iend
    integer(kind = ik)::global_Jstart,global_Jend
    integer(kind = ik)::global_Kstart,global_Kend
    integer(kind = ik)::i_local_start,i_local_end
    integer(kind = ik)::j_local_start,j_local_end
    integer(kind = ik)::k_local_start,k_local_end
  
    ! output the initial grid .pvts file
    if(myid == 0) then
     !inquire(file="JACO", exist=dir_exists)
     !if(.not. dir_exists) then
     !  call system("mkdir -p JACO")
     !endif
  
     filename_pvts = "JACO/Jacobian.pvts"
     call output_global_pvts_jacobian(filename_pvts)    
    endif
    
    !if(myid == 0) then
    ! write(*,*)"Output the jacobian .vts files"
    !endif
  
    global_Istart = i_offset(npx);
    global_Iend   = i_offset(npx) + i_nn(npx) - 1;
  
    global_Jstart = 1_ik;
    global_Jend   = Ny;
  
    global_Kstart = k_offset(npz);
    global_Kend   = k_offset(npz) + k_nn(npz) - 1;

    if( if_parallel) then  
     !Output the VTS files, we have 9 cases, we need to handle all the things
     if(npx == 0) then
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend;   
       i_local_start = 0; 
       i_local_end = nx_local ;
     else if( npx == npx0 - 1) then
       global_Istart = global_Istart - 1;
       i_local_start = 0; 
       i_local_end = nx_local;
     else
       global_Istart = global_Istart - 1;
       global_Iend   = global_Iend   + 1;
       i_local_start = 0; 
       i_local_end = nx_local + 1;     
     endif 
   
       j_local_start = 1; 
       j_local_end = Ny;  
   
     if(npz == 0) then
       global_Kstart = global_Kstart;
       global_Kend   = global_Kend + 1; 
       k_local_start = 1; 
       k_local_end = nz_local + 1;    
     else if( npz == npz0 - 1) then
       global_Kstart = global_Kstart - 1;
       global_Kend   = global_Kend   + 1;
       k_local_start = 0; 
       k_local_end = nz_local + 1;
     else
       global_Kstart = global_Kstart-1 ;
       global_Kend   = global_Kend+1 ;
       k_local_start = 0; 
       k_local_end = nz_local+1 ;  
     endif
    else
      global_Istart = global_Istart-1
      global_Iend   = global_Iend  
      global_Jstart = global_Jstart 
      global_Jend   = global_Jend
      global_Kstart = global_Kstart 
      global_Kend   = global_Kend + 1
      
      i_local_start = 0; 
      i_local_end = nx_local; 
      j_local_start = 1; 
      j_local_end = Ny; 
      k_local_start = 1; 
      k_local_end = nz_local + 1;
    endif
    
    filename_output_jaco_vts ="";
    write(filename_output_jaco_vts,'("JACO/JACO_Part",I5.5,".vts")')MyId

    CALL output_vts_files_jacobian_local_core(filename_output_jaco_vts,i_local_start,i_local_end,&
                                            j_local_start,j_local_end,k_local_start,k_local_end,&
                                            global_Istart,global_Iend,global_Jstart,global_Jend,&
                                            global_Kstart,global_Kend)
  end subroutine output_jacobian
   
  subroutine output_global_pvts_jacobian(filename_pvts)
     implicit none
     integer::ic,global_istart,global_iend
     integer::jc,global_jstart,global_jend
     integer::kc,global_kstart,global_kend
     integer::local_proc_num
     character(128) :: filename_pvts
     character(128) :: Piece_SourceFile_vts
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float

     !vtk_float = "Float64"
     vtk_float = "Float32"
     
     !write(*,*)"Output the global pvts files"

     open(unit=123, file=filename_pvts, status="replace")
     
     xml_part = ' &
     & <VTKFile type="PStructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <PStructuredGrid WholeExtent="0 '//int2str(Nx)//' 1 '//int2str(Ny)//' 1 '//int2str(Nz+1)//'" GhostLevel="1"> &
     &   <PPoints> &
     &    <PDataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended"/> &
     &   </PPoints> &  
     &   <PPointData> &
     &    <PDataArray type="'//vtk_float//'" Name="dxidx" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dxidy" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dxidz" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="detadx" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="detady" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="detadz" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzetadx" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzetady" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzetadz" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dxdxi" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dxdeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dxdzeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dydxi" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dydeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dydzeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzdxi" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzdeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="dzdzeta" format="appended"/> &
     &    <PDataArray type="'//vtk_float//'" Name="Jaco" format="appended"/> &
     &   </PPointData> '
     ! Here we have only one mpi_parts along the eta direction
     ! therefore, the npy0 = 1, and the jc = 0
     jc = 0_ik;
     !npy0 = 1_ik;

     do ic = 0, npx0 - 1    ! \xi processes
         do kc = 0, npz0 - 1  ! \zeta processes
            ! Get the local processes number
            local_proc_num = kc*(npx0*npy0)+jc*npx0+ic;

            Piece_SourceFile_vts = "";
            write(Piece_SourceFile_vts,'("JACO_Part",I5.5,".vts")')local_proc_num
           
            global_istart = i_offset(ic);
            global_iend   = i_offset(ic) + i_nn(ic) - 1;
  
            ! Here we have only one mpi_parts along the eta direction
            global_jstart = 1;
            global_jend   = Ny;
  
            global_kstart = k_offset(kc);
            global_kend   = k_offset(kc) + k_nn(kc) - 1;

          if(If_parallel) then 
            if(ic == 0 ) then
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend ;  
           elseif(ic == npx0 - 1) then
             global_Istart = global_Istart - 1;  
             global_Iend   = global_Iend;
           else
             global_Istart = global_Istart - 1;
             global_Iend   = global_Iend   + 1;  
           endif
  
           if(kc == 0 ) then
             global_Kstart = global_Kstart ;
             global_Kend   = global_Kend   + 1; 
           elseif(kc == npz0 - 1) then
             global_Kstart = global_Kstart - 1;
             global_Kend   = global_Kend   + 1; 
           else
             global_Kstart = global_Kstart-1 ;
             global_Kend   = global_Kend+1   ; 
           endif
          else
            global_Istart = global_Istart - 1;
            global_Iend   = global_Iend ;
            global_Jstart = global_Jstart;
            global_Jend   = global_Jend;
            global_Kstart = global_Kstart;
            global_Kend   = global_Kend + 1;
          endif

          xml_part = trim(adjustl(xml_part)) // ' &
         & <Piece Extent= " '//int2str(global_istart)//' '//int2str(global_iend)//' &
         & '//int2str(global_jstart)//' '//int2str(global_jend)//' &
         & '//int2str(global_kstart)//' '//int2str(global_kend)//'" Source="'//trim(Piece_SourceFile_vts)//'"/> '

         enddo
     enddo

     xml_part = trim(adjustl(xml_part)) // ' &
     &  </PStructuredGrid> &
     & </VTKFile> ' 
      write(123,'(A)')trim(adjustl(xml_part))
     close(123)
  
  end subroutine output_global_pvts_jacobian

  subroutine output_vts_files_jacobian_local_core(filename_output_vts,i_local_start,i_local_end,&
                                     j_local_start,j_local_end,k_local_start,k_local_end,&
                                                     istart,iend,jstart,jend,kstart,kend)
     use mpi, only:MPI_OFFSET_KIND
     implicit none
     integer::i,j,k
     integer::icounts
     integer::i_local_start,i_local_end,j_local_start,j_local_end,k_local_start,k_local_end
     integer::loc_i_length,loc_j_length,loc_k_length 
     integer::istart,iend,jstart,jend,kstart,kend
     integer, parameter :: int64_kind = selected_int_kind(2*range(1))
     integer(int64_kind) :: gridsize_64
     integer::size_real
     integer (KIND=MPI_OFFSET_KIND) :: offset,offset_x,delta_offset_w
     
     !character(len=16) :: int2str
     !character(len=32) :: int2str_o
     character(len=65536) :: xml_part
     character(len=7) :: vtk_float
     character(128)::filename_output_vts
  
     real(4),allocatable::pointData2(:,:)
     real(4),allocatable::tempVar(:,:,:)
  
     loc_i_length = i_local_end - i_local_start + 1;
     loc_j_length = j_local_end - j_local_start + 1;
     loc_k_length = k_local_end - k_local_start + 1;
      
     allocate(pointData2(3,loc_i_length*loc_j_length*loc_k_length))
     allocate(tempVar(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end))
  
     pointData2 = 0.0;
     tempVar = 0.0;
  
     ! single is enough
     size_real = 4;
     vtk_float = "Float32";
  
     !print*,size_real,int64_kind
  
     gridsize_64 = int(size_real,int64_kind) * int(loc_i_length,int64_kind) &
                                            &* int(loc_j_length,int64_kind) &
                                            &* int(loc_k_length,int64_kind)
     delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 
     offset_x = 0
  
     open(unit=1231, file=trim(filename_output_vts), access="stream", form="unformatted", status="replace")
     ! forming the xml_part files
     xml_part = ' &
     & <?xml version="1.0"?> &
     & <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
     &  <StructuredGrid WholeExtent=" '//int2str(istart)//' '//int2str(iend)//' &
     &'//int2str(jstart)//' '//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &   <Piece Extent=" '//int2str(istart)//' '//int2str(iend)//' '//int2str(jstart)//' &
     &'//int2str(jend)//' '//int2str(kstart)//' '//int2str(kend)//'"> &
     &    <Points> &
     &     <DataArray type="'//vtk_float//'" NumberOfComponents="3" format="appended" offset="'//int2str_o(offset_x)//'"/> &
     &    </Points> &
     &   <PointData> '
  
     offset = offset_x + 3 * gridsize_64 + storage_size(gridsize_64)/8 
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxidx" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxidy" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxidz" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="detadx" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="detady" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="detadz" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzetadx" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzetady" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzetadz" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxdxi" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxdeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dxdzeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dydxi" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dydeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w
  
      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dydzeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzdxi" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzdeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="dzdzeta" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w      

      xml_part = trim(adjustl(xml_part)) // ' &
     & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Jaco" format="appended" &
     &  offset="'//int2str_o(offset)//'"/>'
      offset = offset + delta_offset_w 

     xml_part = trim(adjustl(xml_part)) // ' &
     &       </PointData> &
     &     </Piece> &
     &   </StructuredGrid> &
     &   <AppendedData encoding="raw"> '

     ! the xml_part is finished
     !write(*,*)"write head"
       write(1231) trim(adjustl(xml_part))
     !write(*,*)"write xx,yy,zz"
       write(1231) "_"
       write(1231) 3*size_real*int(gridsize_64,int64_kind)
       icounts = 1
       do k = k_local_start,k_local_end
         do j = j_local_start,j_local_end
           do i = i_local_start,i_local_end
             pointData2(1,icounts) = real(X_grid(i,j,k),kind=4)
             pointData2(2,icounts) = real(Y_grid(i,j,k),kind=4)
             pointData2(3,icounts) = real(Z_grid(i,j,k),kind=4)
              icounts = icounts + 1     
           enddo
         enddo
       enddo
       !write(1231)pointData
       write(1231)pointData2
  
       write(1231)gridsize_64
       tempVar = real(dxidx(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dxidy"
       write(1231)gridsize_64
       tempVar = real(dxidy(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write dxidz"
       write(1231)gridsize_64
       tempVar = real(dxidz(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write detadx"
       write(1231)gridsize_64
       tempVar = real(detadx(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write detady"
       write(1231)gridsize_64
       tempVar = real(detady(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write detadz"
       write(1231)gridsize_64
       tempVar = real(detadz(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dzetadx"
       write(1231)gridsize_64
       tempVar = real(dzetadx(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dzetady"
       write(1231)gridsize_64
       tempVar = real(dzetady(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dzetadz"
       write(1231)gridsize_64
       tempVar = real(dzetadz(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dxdxi"  
       write(1231)gridsize_64
       tempVar = real(dxdxi(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dxdeta"
       write(1231)gridsize_64
       tempVar = real(dxdeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar   
  
     !write(*,*)"write dxdzeta"
       write(1231)gridsize_64
       tempVar = real(dxdzeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dydxi"
       write(1231)gridsize_64
       tempVar = real(dydxi(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dydeta"
       write(1231)gridsize_64
       tempVar = real(dydeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dydzeta"
       write(1231)gridsize_64
       tempVar = real(dydzeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write dzdxi"
       write(1231)gridsize_64
       tempVar = real(dzdxi(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dzdeta"
       write(1231)gridsize_64
       tempVar = real(dzdeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar
  
     !write(*,*)"write dzdzeta"
       write(1231)gridsize_64
       tempVar = real(dzdzeta(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar

     !write(*,*)"write Jaco"
       write(1231)gridsize_64
       tempVar = real(Jaco(i_local_start:i_local_end,j_local_start:j_local_end,k_local_start:k_local_end),kind=4)
       write(1231)tempVar           

     write(1231) ' &
      &    </AppendedData> &
      &  </VTKFile>' 
     close(1231)
     
     deallocate(pointData2)
     deallocate(tempVar)
  
  end subroutine output_vts_files_jacobian_local_core

  function int2str(int_num)
    implicit none
    integer :: int_num
    character(len=16) :: int2str, ret_value
    write(ret_value, "(I0)") int_num
    int2str = ret_value
  end function int2str

  function int2str_o(int_num)
    use mpi, only:MPI_OFFSET_KIND
    implicit none
    integer(KIND=MPI_OFFSET_KIND) :: int_num
    character(len=32) :: int2str_o, ret_value
    write(ret_value, "(I0)") int_num
    int2str_o = ret_value
  end function int2str_o

END MODULE OutputParaView