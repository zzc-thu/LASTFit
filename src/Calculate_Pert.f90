subroutine Calculate_Pert(Iter)
  ! Module dependencies
  use SF_Constant,   only: ik,rk,NumVar,FilesForContinue,pi,OverLap,PertForContinue
  use SF_CFD_Global, only: IF_Continue_Calculate,nx_local,nz_local,&
                           Ny,Gamma,Mach_Ref
  use MPI_GLOBAL,    only: MyId,ierr,MPI_COMM_WORLD

  ! some arrays and vectors
  use SF_CFD_Global, only: U0,V0,W0,P0,T0,Rho0,Ucons0_steady,&
                           U,V,W,P,T,Rho,Ucons0,&
                           U_interp,V_interp,W_interp,P_interp,T_interp,Rho_interp,&
                           U_inf,V_inf,W_inf,P_inf,Rho_inf,T_inf,&     
                           U_free,V_free,W_free,P_free,Rho_free,T_free,&   
                           U_pert,V_pert,W_pert,P_pert,T_pert,Rho_pert,&
                           ShockH_steady,ShockV_steady,shockH,shockV,&
                           Heta,X_grid_steady,Y_grid_steady,Z_grid_steady,&
                           x_grid_interp,y_grid_interp,z_grid_interp,&
                           WallSx,WallSy,WallSz,WallNormalX,WallNormalY,WallNormalZ,&
                           x_grid,y_grid,z_grid,ShockH,ShockV,&
                           Hdst,Hdst_steady,StepsWriteData,OutputFileNo,Nx,Nz,&
                           etaSC,Hdst_interp,Rho0_interp,U0_interp,V0_interp,W0_interp,P0_interp,T0_interp
  ! some subroutines and functions
  use MPI_GLOBAL,     only: Parallel_Exchange,MPI_Barrier,&
                          & Parallel_Exchange_NumVar,Parallel_Exchange_Surface,Parallel_Exchange_interp
  use OutputParaView, only: output_Pert_results
  use SFitting,       only: ShockRelation3D
  use grid_distribution, only: GridDistributions
  
  implicit none
  ! Local variables
  integer(kind = ik),intent(in):: Iter
  logical:: dir_exists
  character(len=100)::FILENAME_FLOW
  character(len=100)::FILENAME_SHK
  integer(kind=ik ):: ic,jc,kc,loc_jc
  real( kind = rk ):: tempP,tempT,indexJ
  
    ! Here we read the steady state fields from previous simulation
    ! in order to increase the efficient of the simulation
    ! all the processes read the sepearated files simultaneously
  
    !if(MyId == 0) then
    !  write(*,'(A)')"We are reading the initial fields from previous simulation"
    !  inquire(file="./RESU_STEADY/SF_Results00000.flowsfg", exist=dir_exists)
    !   if (.not. dir_exists) then
    !     write(*,'(A)')"The steady file does not exist"
    !      stop
    !   else
    !     write(*,'(A)')"The steady file exists, we will read the steady fields from it"
    !   end if     
    !endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
      
      ! Read the previous results
      write(filename_flow,'(A,A,I5.5,A)')'RESU_STEADY/',trim(FilesForContinue),MyId,'.flowsfg'
      !if(MyId == 0 ) write(*,*)"Reading steady state flow"
      open(112,file=filename_flow,form='unformatted',status='old')
       read(112)Ucons0_steady
       read(112)x_grid_steady
       read(112)y_grid_steady
       read(112)z_grid_steady
      close(112)
      
       ! Read the previous shock informations
      write(filename_shk,'(A,A,I5.5,A)')'RESU_STEADY/',trim(FilesForContinue),MyId,'.shksfg'
     !if(MyId == 0 ) write(*,*)"Reading ShockH and ShockV"
     open(113,file=filename_shk,form='unformatted',status='old')
      read(113)ShockH_steady
      read(113)ShockV_steady
     close(113)
      
      ! Transfer the conservative variables to the basic variables
      do kc = 1, nz_local
       do jc = 1,Ny
        do ic = 1, nx_local
         call CVtoU_CPG(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),Ucons0_steady(1:NumVar,ic,jc,kc))
         ! Calculate the temperature
         T0(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P0(ic,jc,kc) / Rho0(ic,jc,kc)
        enddo
       enddo
      enddo
   !  generate the perturbed grid
       
   call GridDistributions(Ny,(Ny+1)/2,(Ny+1)/2,1.05d0,etaSC)  
      
   ! Calculate the Distance to the wall
   do kc = 1-overLap, nz_local+overLap
    do jc = 1,Ny
     do ic = 1-overLap, nx_local+overLap
         Hdst(ic,jc,kc)          = ShockH(ic,kc)* Heta(ic,jc,kc);
         Hdst_steady(ic,jc,kc)   = ShockH_steady(ic,kc)* Heta(ic,jc,kc);
     enddo
    enddo
   enddo
   
   do kc = 1-overLap, nz_local+overLap
    do jc = 1,2*Ny-1
     do ic = 1-overLap, nx_local+overLap
          Hdst_interp(ic,jc,kc)   = shockH_steady(ic,kc) * etaSC(jc)
          X_grid_interp(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * Hdst_interp(ic,jc,kc);
          Y_grid_interp(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * Hdst_interp(ic,jc,kc);
          Z_grid_interp(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * Hdst_interp(ic,jc,kc);
       enddo
    enddo
   enddo
 
    do kc = 1, nz_local
      do ic = 1,nx_local
          do jc = 1,2*Ny-1
          	 if (Hdst_interp(ic,jc,kc) .GE. ShockH_steady(ic,kc))then
              loc_jc = jc
             exit 
             endif
          end do
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),Rho0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),Rho0_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),U0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),U0_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),V0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),V0_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),W0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),W0_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),P0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),P0_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst_steady(ic,1:Ny,kc),T0(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),T0_interp(ic,1:loc_jc,kc),loc_jc)
 
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),Rho0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,Rho0_interp(ic,1:loc_jc,kc))
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),U0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,U0_interp(ic,1:loc_jc,kc))
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),V0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,V0_interp(ic,1:loc_jc,kc))
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),W0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,W0_interp(ic,1:loc_jc,kc))
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),P0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,P0_interp(ic,1:loc_jc,kc))
    !call NewtonInterpHighOrder(Hdst_steady(ic,1:Ny,kc),T0(ic,1:Ny,kc),ny,Hdst_interp(ic,1:loc_jc,kc),loc_jc,T0_interp(ic,1:loc_jc,kc))          
          
          Rho0_interp(ic,loc_jc:2*Ny-1,kc)=Rho_inf
          U0_interp(ic,loc_jc:2*Ny-1,kc)=U_inf
          V0_interp(ic,loc_jc:2*Ny-1,kc)=V_inf
          W0_interp(ic,loc_jc:2*Ny-1,kc)=W_inf
          P0_interp(ic,loc_jc:2*Ny-1,kc)=P_inf
          T0_interp(ic,loc_jc:2*Ny-1,kc)=T_inf
      enddo 
    enddo

 ! 打开文件
!open(unit=10, file='output.dat',  status='replace', action='write')
!
!! 假设数组Rho0(1,:,1)有n个元素，我们使用循环
!do ic = 1, size(Rho0, 2)   ! 第二个维度的大小
!    write(10, '(F10.5)') Rho0(1, ic, 1)
!end do
!
!do ic = 1, size(Rho0_interp, 2)   ! 第二个维度的大小
!    write(10, '(F10.5)') Rho0_interp(1, ic, 1)
!end do
!! 关闭文件
!close(10)   
    
    do kc = 1, nz_local
      do ic = 1,nx_local
          do jc = 1,2*Ny-1
          	 if (Hdst_interp(ic,jc,kc) .GE. ShockH(ic,kc))then
              loc_jc = jc
             exit 
             endif
          end do
          
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),Rho(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),Rho_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),U(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),U_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),V(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),V_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),W(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),W_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),P(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),P_interp(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),T(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),T_interp(ic,1:loc_jc,kc),loc_jc)
             
          Rho_interp(ic,loc_jc:2*Ny-1,kc)=Rho_free(ic,Ny,kc)
          U_interp(ic,loc_jc:2*Ny-1,kc)=U_free(ic,Ny,kc)
          V_interp(ic,loc_jc:2*Ny-1,kc)=V_free(ic,Ny,kc)
          W_interp(ic,loc_jc:2*Ny-1,kc)=W_free(ic,Ny,kc)
          P_interp(ic,loc_jc:2*Ny-1,kc)=P_free(ic,Ny,kc)
          T_interp(ic,loc_jc:2*Ny-1,kc)=T_free(ic,Ny,kc)
      enddo 
    enddo
    
    !write(*,*)"loc_jc:",loc_jc
    !write(*,*)"Rho_interp：",Rho_interp(1,:,1)
    !calculate the perturbed fields
    do kc = 1, nz_local
     do jc = 1,2*Ny-1
      do ic = 1, nx_local
        Rho_pert(ic,jc,kc) = Rho_interp(ic,jc,kc) - Rho0_interp(ic,jc,kc)
          U_pert(ic,jc,kc) = U_interp(ic,jc,kc) - U0_interp(ic,jc,kc)
          V_pert(ic,jc,kc) = V_interp(ic,jc,kc) - V0_interp(ic,jc,kc)
          W_pert(ic,jc,kc) = W_interp(ic,jc,kc) - W0_interp(ic,jc,kc)
          P_pert(ic,jc,kc) = P_interp(ic,jc,kc) - P0_interp(ic,jc,kc)
          T_pert(ic,jc,kc) = T_interp(ic,jc,kc) - T0_interp(ic,jc,kc)
      enddo
     enddo
    enddo
    
!write(*,*)"Rho_pert：",Rho_pert(1,:,1)
  
    call Parallel_Exchange(X_grid_steady)
    call Parallel_Exchange(Y_grid_steady)
    call Parallel_Exchange(Z_grid_steady)
    
    call Parallel_Exchange_interp(X_grid_interp)
    call Parallel_Exchange_interp(Y_grid_interp)
    call Parallel_Exchange_interp(Z_grid_interp)    
    
    call Parallel_Exchange(Rho0)
    call Parallel_Exchange(U0)
    call Parallel_Exchange(V0)
    call Parallel_Exchange(W0)
    call Parallel_Exchange(P0)
    call Parallel_Exchange(T0)
    
    call Parallel_Exchange_interp(Rho0_interp)
    call Parallel_Exchange_interp(U0_interp)
    call Parallel_Exchange_interp(V0_interp)
    call Parallel_Exchange_interp(W0_interp)
    call Parallel_Exchange_interp(P0_interp)
    call Parallel_Exchange_interp(T0_interp)
    
    call Parallel_Exchange_interp(Rho_interp)
    call Parallel_Exchange_interp(U_interp)
    call Parallel_Exchange_interp(V_interp)
    call Parallel_Exchange_interp(W_interp)
    call Parallel_Exchange_interp(P_interp)
    call Parallel_Exchange_interp(T_interp)
    
    call Parallel_Exchange_interp(Rho_pert)
    call Parallel_Exchange_interp(U_pert)
    call Parallel_Exchange_interp(V_pert)
    call Parallel_Exchange_interp(W_pert)
    call Parallel_Exchange_interp(P_pert)
    call Parallel_Exchange_interp(T_pert)
    
    call Parallel_Exchange_NumVar(Ucons0_steady)
    call Parallel_Exchange_Surface(ShockH_steady)
    call Parallel_Exchange_Surface(ShockV_steady)
    
    
!if(mod(Iter,StepsWriteData) == 0) then
  ! write(*,*),size(Rho_pert,1),size(Rho_pert,2),size(Rho_pert,3)
  !if(MyId == 0 ) write(*,*)"Saving Perturbed fields finished" 
   ! Output the flow info
     !write(filename_flow,'(A,A,I5.5,"_",I5.5,A)')'Pert/',trim(PertForContinue),MyId,OutputFileNo,'.dat'
     !open(233,file=filename_flow,form='unformatted',access = 'stream')
     !write(233)x_grid_interp(1:nx_local,:,1)
     !write(233)y_grid_interp(1:nx_local,:,1)
     !write(233)Rho_pert(1:nx_local,:,1)
     !write(233)U_pert(1:nx_local,:,1)
     !write(233)V_pert(1:nx_local,:,1)
     !write(233)W_pert(1:nx_local,:,1)
     !write(233)P_pert(1:nx_local,:,1)
     !write(233)T_pert(1:nx_local,:,1)
     !
     !close(233)
!100 FORMAT(F33.16)    
   
   call output_Pert_results
   
    !Output the Total Results
    !if(MyId == 0) then
    !write(filename_flow,'(A,A,"_",I5.5,A)')'Pert/',trim(PertForContinue),OutputFileNo,'.dat'
    !open(55,file=filename_flow,form='unformatted',access = 'stream')
    !endif
    !call Parallel_OutputFunction(55,x_grid_interp(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,y_grid_interp(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,Rho_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,U_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,V_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,W_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,P_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)
    !call Parallel_OutputFunction(55,T_pert(:,:,4),Nx,2*Ny-1,nx_local,2*Ny-1)    
    !
    !if(MyId == 0) then
    !close(55)
    !endif
    
OutputFileNo = OutputFileNo + 1;  

 end subroutine Calculate_Pert
 
 
!subroutine InterpToSteadyGrid(h_unsteady, phi_unsteady, h_steady, phi_interp)
!    !--- 输入输出参数声明 ---
!    use SF_Constant,  only: ik,rk,OverLap
!    use SF_CFD_Global,only: nx_local,nz_local,Ny
!    implicit none
!    !real( kind = rk ), intent(in)  :: h_unsteady  (1-OverLap:nx_local+OverLap, 1:ny, 1-OverLap:nz_local+OverLap)   ! 非定常网格高度(3D)
!    !real( kind = rk ), intent(in)  :: phi_unsteady(1-OverLap:nx_local+OverLap, 1:ny, 1-OverLap:nz_local+OverLap)  ! 非定常物理量(3D)
!    !real( kind = rk ), intent(in)  :: h_steady    (1-OverLap:nx_local+OverLap, 1:2*ny-1, 1-OverLap:nz_local+OverLap)     ! 定常网格高度(3D)
!    !real( kind = rk ), intent(out) :: phi_interp  (1-OverLap:nx_local+OverLap, 1:2*ny-1, 1-OverLap:nz_local+OverLap)  ! 输出物理量(3D)
!    real( kind = rk ), intent(in)  :: h_unsteady  (:, :, :)   ! 非定常网格高度(3D)
!    real( kind = rk ), intent(in)  :: phi_unsteady(:, :, :)  ! 非定常物理量(3D)
!    real( kind = rk ), intent(in)  :: h_steady    (:, :, :)   ! 定常网格高度(3D)
!    real( kind = rk ), intent(out) :: phi_interp  (:, :, :)  ! 输出物理量(3D)    
!    
!    !--- 局部变量 ---
!    integer( kind = ik ) :: ic, jc, kc
!    
!    !--- 三维网格遍历 ---
!        do kc = 1, nz_local
!            do ic = 1, nx_local
!                
!                !call LinearInterp1D(h_unsteady(ic,:,kc), phi_unsteady(ic,:,kc), h_steady(ic,:,kc), phi_interp(ic,:,kc))
!                call CubicSplineInterp1D(h_unsteady(ic,:,kc), phi_unsteady(ic,:,kc), h_steady(ic,:,kc), phi_interp(ic,:,kc),)
!            end do
!        end do
!
!    end subroutine InterpToSteadyGrid
    
    !subroutine LinearInterp1D(x, y, xi, yi) 
    !!--- 输入输出参数声明 --- 
    !use SF_Constant,  only: ik,rk 
    !use SF_CFD_Global,only: ny 
    !
    !implicit none 
    !real(kind = rk), intent(in) :: x(ny)    ! 原始网格点 
    !real(kind = rk), intent(in) :: y(ny)    ! 原始物理量 
    !real(kind = rk), intent(in) :: xi(ny)   ! 插值网格点 
    !real(kind = rk), intent(out) :: yi(ny)  ! 插值后的物理量 
    !
    !
    !!--- 局部变量 --- 
    !integer(kind = ik) :: i, j 
    !real(kind = rk) :: slope 
    !
    !!--- 遍历插值点 --- 
    !do i = 1, ny 
    !    ! 寻找插值点所在的区间 
    !    do j = 1, ny - 1 
    !        if (xi(i) >= x(j) .and. xi(i) <= x(j+1)) then 
    !            ! 计算线性插值的斜率 
    !            slope = (y(j+1) - y(j)) / (x(j+1) - x(j)) 
    !            ! 进行线性插值 
    !            yi(i) = y(j) + slope * (xi(i) - x(j)) 
    !            exit 
    !        end if 
    !    end do 
    !end do 
    !end subroutine LinearInterp1D 
    
    ! 导入必要的模块 
subroutine CubicSplineInterp1D(x, y, xi, yi,local_ny)
    !--- 输入输出参数声明 --- 
    use SF_Constant,  only: ik,rk 
    use SF_CFD_Global,only: ny 
    implicit none 
    integer(kind = ik), intent(in) :: local_ny
    real(kind = rk), intent(in) :: x(1:ny)    ! 原始网格点 
    real(kind = rk), intent(in) :: y(1:ny)    ! 原始物理量 
    real(kind = rk), intent(in) :: xi(1:local_ny)   ! 插值网格点
    real(kind = rk), intent(out) :: yi(1:local_ny)  ! 插值后的物理量
 
    !--- 局部变量 --- 
    integer(kind = ik) :: i, j 
    real(kind = rk) :: h(ny-1), alpha(ny-1), l(ny), mu(ny), z(ny) 
    real(kind = rk) :: a(ny), b(ny), c(ny), d(ny) 
    real(kind = rk) :: dx, dy 
 
    ! 计算步长 
    do i = 1, ny-1 
        h(i) = x(i+1) - x(i) 
    end do 
 
    ! 计算 alpha 
    do i = 2, ny-1 
        alpha(i) = 3.0_rk / h(i) * (y(i+1) - y(i)) - 3.0_rk / h(i-1) * (y(i) - y(i-1)) 
    end do 
 
    ! 求解三对角线性方程组 
    l(1) = 1.0_rk 
    mu(1) = 0.0_rk 
    z(1) = 0.0_rk 
    do i = 2, ny-1 
        l(i) = 2.0_rk * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1) 
        mu(i) = h(i) / l(i) 
        z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i) 
    end do 
    l(ny) = 1.0_rk 
    z(ny) = 0.0_rk 
    c(ny) = 0.0_rk 
    do j = ny-1, 1, -1 
        c(j) = z(j) - mu(j) * c(j+1) 
        b(j) = (y(j+1) - y(j)) / h(j) - h(j) * (c(j+1) + 2.0_rk * c(j)) / 3.0_rk 
        d(j) = (c(j+1) - c(j)) / (3.0_rk * h(j)) 
        a(j) = y(j) 
    end do 
 
    ! 遍历插值点 
    do i = 1, local_ny
        ! 寻找插值点所在的区间 
        do j = 1, ny - 1 
            if (xi(i) >= x(j) .and. xi(i) <= x(j+1)) then 
                dx = xi(i) - x(j) 
                ! 进行三次样条插值 
                yi(i) = a(j) + b(j) * dx + c(j) * dx**2 + d(j) * dx**3 
                exit 
            end if 
        end do 
    end do 
end subroutine CubicSplineInterp1D 

SUBROUTINE NewtonInterpHighOrder(x, y, n, xi, m, yi)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: n, m 
    REAL*8, INTENT(IN) :: x(n), y(n), xi(m)
    REAL*8, INTENT(OUT) :: yi(m)
    REAL*8, ALLOCATABLE :: F(:,:)  ! 动态分配差商表
    INTEGER :: i, j, k 
    REAL*8 :: term, prod_term 
 
    ! 动态分配内存以避免堆栈溢出[4]()
    ALLOCATE(F(n, n))
    
    ! 初始化差商表
    F(:, 1) = y(:)
    
    ! 构建差商表（至n-1阶）
    DO j = 2, n 
        DO i = j, n 
            F(i, j) = (F(i, j-1) - F(i-1, j-1)) / (x(i) - x(i-j+1))
        END DO 
    END DO 
 
    ! 牛顿插值公式 
    DO k = 1, m 
        term = F(1, 1)  ! 常数项 
        prod_term = 1.0D0 
        DO j = 2, n 
            prod_term = prod_term * (xi(k) - x(j-1))
            term = term + F(j, j) * prod_term 
        END DO 
        yi(k) = term 
    END DO 
 
    ! 释放内存 
    DEALLOCATE(F)
END SUBROUTINE NewtonInterpHighOrder 
    
!    SUBROUTINE Parallel_OutputFunction3D(FileNumber,Uin,Nx,Ny,Nz,nx_local,ny_local,nz_local) 
!    use SF_Constant 
!    use MPI_Global
!    implicit none 
! 
!    ! 三维输出函数说明： 
!    ! - z方向包含overlap（类似原x方向） 
!    ! - y方向无overlap 
!    ! - 数据按x-y-z维度存储 
! 
!    INTEGER(kind = ik)::ki, ia, i, j, k, FileNumber 
!    INTEGER(kind = ik)::np_recv, i1, j1, k1 
!    INTEGER(kind = ik)::Nx, Ny, Nz, nx_local, ny_local, nz_local 
!    INTEGER Status(MPI_STATUS_SIZE) 
!    !INTEGER(kind = ik)::ierr  ! 定义MPI错误变量 
!    REAL(kind = rk)::Uin(1-overLAP:nx_local+overLAP, 1:ny_local, 1-overLAP:nz_local+overLAP) 
!    REAL(kind = rk)::U(nx_local, ny_local, nz_local) 
!    REAL(kind = rk)::buff1(nx_local*ny_local*nz_local), buff2(nx_local, ny_local, nz_local) 
!    REAL(kind = rk)::F(Nx, Ny, Nz) 
! 
!    ! 初始化本地数据（去除overlap） 
!    U = Uin(1:nx_local, 1:ny_local, 1:nz_local) 
! 
!    ! 主进程处理 
!    if(myid == 0) then 
!        !allocate(buff1(nx_local*ny_local*nz_local), stat=ierr) 
!        !if (ierr /= 0) then 
!        !    write(*,*) 'Memory allocation for buff1 failed!' 
!        !    stop 
!        !end if 
!        F = 0.0_rk 
! 
!        do ki = 0, numprocess-1 
!            np_recv = ki 
! 
!            if(np_recv == 0) then 
!                ! 主进程直接写入本地数据 
!                do i1 = 1, nx_local 
!                    do j1 = 1, ny_local 
!                        do k1 = 1, nz_local 
!                            F(i1,j1,k1) = U(i1,j1,k1) 
!                        end do 
!                    end do 
!                end do 
!            else 
!                ! 接收其他进程数据 
!                call MPI_Recv(buff1, i_nn(ki)*ny*k_nn(ki), MPI_DOUBLE_PRECISION, & 
!                              np_recv, 90001, MPI_COMM_WORLD, Status, ierr) 
!                !if (ierr /= 0) then 
!                !    write(*,*) 'MPI_Recv failed!' 
!                !    stop 
!                !end if 
! 
!                ! 重组三维数据 
!                do k1 = 1, k_nn(ki) 
!                    do j1 = 1, ny_local 
!                        do i1 = 1, i_nn(ki) 
!                            i = i_offset(ki) + i1 - 1 
!                            k = k_offset(ki) + k1 - 1 
!                            ia = (k1-1)*i_nn(ki)*ny_local + (j1-1)*i_nn(ki) + i1 
!                            F(i,j1,k) = buff1(ia) 
!                        end do 
!                    end do 
!                end do 
!            end if 
!        end do 
! 
!        ! 写入全局文件 
!        WRITE(FileNumber) F 
!        !deallocate(buff1) 
!    else 
!        ! 非主进程发送数据 
!        ! allocate(buff2(nx_local, ny_local, nz_local), stat=ierr) 
!        !if (ierr /= 0) then 
!        !    write(*,*) 'Memory allocation for buff2 failed!' 
!        !    stop 
!        !end if 
!        buff2 = U 
! 
!        ! 展平三维数组为一维 
!        call MPI_Bsend(buff2(1,1,1), nx_local*ny_local*nz_local, MPI_DOUBLE_PRECISION, & 
!                       0, 90001, MPI_COMM_WORLD, ierr) 
!        !if (ierr /= 0) then 
!        !    write(*,*) 'MPI_Bsend failed!' 
!        !    stop 
!        !end if 
!        !deallocate(buff2) 
!    end if 
! 
!    ! 同步所有进程 
!    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!    !if (ierr /= 0) then 
!    !    write(*,*) 'MPI_BARRIER failed!' 
!    !    stop 
!    !end if 
! 
! 
!END SUBROUTINE Parallel_OutputFunction3D 
!    
