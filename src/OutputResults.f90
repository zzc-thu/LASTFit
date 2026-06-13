subroutine OutputResults(Iter)
  use SF_Constant,       only:ik,FilesForContinue
  use SF_CFD_Global,     only:Ny,nx_local,nz_local,StepsWriteData,OutputFileNo,&
                              Nx,Nz,ModelType
  use MPI_GLOBAL,        only:MyId,NumProcess!Parallel_OutputFunction3D
  ! matrix and array
  use SF_CFD_Global,     only:x_grid,y_grid,z_grid,Rho,U,V,W,P,T,&
                              ShockH,ShockV,Ucons0
  ! subs and functions
  use OutputParaView,    only:output_results,output_results_singularity

  implicit none 
  integer( kind = ik ),intent(in):: Iter
  character(len=100)filename
  character(len=100)filename_shk,filename_flow
  character(len=10) :: timestep_str
  
  if(mod(Iter,StepsWriteData) == 0) then
    ! 这里输出分割符号来区别输出和计算输出残差
    if(MyId == 0) then
      write(*,"(A)")"========================================================================================="
    endif
  ! Output Continue files
   ! Output the shock info
     write(filename_shk,'(A,A,I5.5,A)')'RESU/',trim(FilesForContinue),MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Saving ShockH and ShockV"
     open(111,file=filename_shk,form='unformatted',status='replace')
     write(111)ShockH
     write(111)ShockV
     close(111)
     if(MyId == 0 ) write(*,*)"Saving finished"    

   ! Output the flow info
     write(filename_flow,'(A,A,I5.5,A)')'RESU/',trim(FilesForContinue),MyId,'.flowsfg'
     if(MyId == 0 ) write(*,*)"Saving flow DATA Ucons0"
     open(112,file=filename_flow,form='unformatted',status='replace')
     write(112)Ucons0
     write(112)x_grid
     write(112)y_grid
     write(112)z_grid
     close(112)
     if(MyId == 0 ) write(*,*)"Saving finished"
    
     !if(MyId == 0 ) write(*,*)"Saving Iteration Number"
     !open(113,file='RESU/IterationNumber.dat',form='unformatted',status='replace')
     !write(113)Iter + 1
     !close(113)
     
    !Output the Total Results
    !if(MyId == 0) then
    !write(filename_flow,'(A,A,"_",I5.5,A)')'RESU/',trim(FilesForContinue),OutputFileNo,'.dat'
    !open(55,file=filename_flow,form='unformatted',access = 'stream')
    !endif
    !call Parallel_OutputFunction(55,x_grid(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,y_grid(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,Rho(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,U(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,V(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,W(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,P(:,:,1),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,T(:,:,1),Nx,Ny,nx_local,Ny)    
    !call Parallel_OutputFunction(55,x_grid(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,y_grid(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,Rho(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,U(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,V(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,W(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,P(:,:,21),Nx,Ny,nx_local,Ny)
    !call Parallel_OutputFunction(55,T(:,:,21),Nx,Ny,nx_local,Ny)       
    !if(MyId == 0) then
    !close(55)
    !endif
    ! 
   ! Output the .pvts and .vts files
   ! Update the OutputFileNo in consequence
     OutputFileNo = OutputFileNo + 1; 
     
     if (ModelType == 1)then
         call output_results
     else
         call output_results_singularity
     endif
     
  endif
  
  end subroutine OutputResults
  
    
    
SUBROUTINE Parallel_OutputFunction(FileNumber,Uin,Nx,Ny,nx_local,ny_local)
  use SF_Constant
  use MPI_Global
  implicit none
  ! This function transfer all the variables
  ! into main process and output the whole file
  ! for calculations.
  INTEGER(kind = ik)::ki,ia,I,FileNumber
  INTEGER(kind = ik)::np_recv,i1,j1
  INTEGER(kind = ik)::Nx,Ny,nx_local,ny_local
  integer::Status(MPI_STATUS_SIZE)
  REAL(kind = rk)::Uin(1-overLAP:nx_local+overLAP,1:ny_local)
  REAL(kind = rk)::U(nx_local,ny_local)
  REAL(kind = rk)::buff1(nx_local*ny_local),buff2(nx_local,ny_local),F(Nx,Ny)
  !write(*,*)Ny,ny_local
  U = Uin(1:nx_local,1:ny_local)
  !write(*,*)myid,nx_local,Ny,ny_local,npx
  if(myid.EQ.0) then
    do ki=0,numprocess - 1	    
      np_recv=ki
     if(np_recv.eq.0) then
       do i1=1,nx_local
        do j1=1,ny_local
         F(i1,j1)=U(i1,j1)
        enddo
       enddo
     else
       
       call MPI_Recv(buff1,i_nn(ki)*Ny,MPI_DOUBLE_PRECISION, np_recv,9001,MPI_COMM_WORLD,Status,ierr)
       do j1=1,Ny
        do i1=1,i_nn(ki)
	       i=i_offset(ki)+i1-1
	       ia=(j1-1)*i_nn(ki)+i1
	       F(i,j1)=buff1(ia)
        enddo
       enddo
     endif
    enddo
    WRITE(FileNumber)F
  else    
   	do i1=1,nx_local
	   do j1=1,Ny
	    buff2(i1,j1)=U(i1,j1)
	   enddo
    enddo
    
    ! Attension this buff2 is used along the colum index first
    call MPI_Bsend(buff2(1,1),nx_local*ny_local,MPI_DOUBLE_PRECISION,0,9001,MPI_COMM_WORLD,ierr)
  endif
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE Parallel_OutputFunction        