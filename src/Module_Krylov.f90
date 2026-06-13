Module Krylov_SubSpace
! This module contains the major processes used for Krylov subspace methods
! The module is designed to be used independently of the main program with 
! the user providing the necessary Matrix-Vector product and preconditioner(optional) 
  use SF_Constant,    only: ik,rk,NumVar,overLAP,Fre_eps,&
                          & InnerGmresTol,MaxKrylovSubSpace,MaxRestartedNumber
  use SF_CFD_Global,  only: nx_local,Ny,nz_local,DT0,Nsize,Nx,Nz,Select_Kind_Precon_Solver
  use mpi    ! MPI module, the Krylov subspace methods are parallelized using MPI

  use SF_CFD_Global,  only: dUcons0_old,ShockH_old,ShockV_old,ShockAc_old,&
                          & Ucons0,ShockH,ShockV,ShockAc,UconsOld,dUcons0,&
                          & dUcons,dShockH,dShockV,RHS_Krylov,DeltaSolu

  ! some functions and subroutines
  use SFitting,       only: Calculate_ShockAC3D
  implicit none
  real( kind = rk ):: Fre_epsUse
  real( kind = rk ), external :: ddot  ! blas 1 function
    
contains
!=================== 这些部分是独立于CFD计算的 ======================
! 理论上适用于任何的矩阵向量乘法和预处理器
 ! FUNCTION: norm_vec
  function norm_vec(v_in)
    
    implicit none
    ! Define some variables
    real( kind = rk )                            :: norm_vec, norm_vec_local
    integer( kind = ik )                         :: Nlength,ic
    integer( kind = ik )                         :: ierr
    ! Define some vectors and arrays
    real( kind = rk ), dimension(:), intent(in)  :: v_in
    
    norm_vec       = 0.0_rk;
    norm_vec_local = 0.0_rk;

    Nlength = size(v_in);
    
    ! 为了程序性能调优以及便于程序理解，我们把naive的实现方式和Blas/Lapack
    ! 的实现方式都标记在注释之中
    !! 这里实际上就是需要计算一个对应的向量的点乘过程，naive的实现方式如下
    !do ic = 1, Nlength
    !    norm_vec_local = norm_vec_local + v_in(ic)*v_in(ic);
    !end do
    ! 调用Blas 1的函数
    norm_vec_local = ddot(Nlength,v_in,1,v_in,1);
    
    ! 这里是一个全局的归约操作，将所有的local的norm_vec_local进行求和
    call MPI_ALLReduce(norm_vec_local,norm_vec,1,MPI_DOUBLE_PRECISION,&
                      &MPI_SUM,MPI_COMM_WORLD,ierr);
    
    ! 这里是一个开方操作
    norm_vec = sqrt(norm_vec);

  end function norm_vec
 
   ! FUNCTION: DOR_PROD
  function dot_prod(v_in,w_in)
          
    implicit none
  
    real(kind=rk)                           :: dot_prod, dot_prod_local
    real(kind=rk), dimension(:), intent(in) :: v_in, w_in
    integer(kind=ik)                        :: Nlength,ic
    integer(kind=ik)                        :: ierr
          
    dot_prod = 0.d0;
    dot_prod_local = 0.d0;
  
    Nlength = size(v_in);
          
    ! Naive implementation
    !do i = 1,Nlength
    !  dot_prod_local = dot_prod_local + v_in(i) * w_in(i);
    !enddo
    ! Blas 1 implementation
    dot_prod_local = ddot(Nlength,v_in,1_ik,w_in,1_ik);
          
    ! We should use MPI_ALLReduce
  
    call MPI_AllReduce(dot_prod_local,dot_prod,1,MPI_DOUBLE_PRECISION,&
                                                 MPI_SUM,MPI_COMM_WORLD,ierr)
  
    !  每一个进程都得到了对应的dot_prod值
  
  end function dot_prod
 
 ! SUBROUTINE Arnodi processes
  subroutine arnoldi(Qmat, Hmat, kstep, in, im)
    implicit none
    integer,                            intent(in)    :: kstep, in, im
    real(kind = rk),dimension(in,im+1), intent(inout) :: Qmat
    real(kind = rk),dimension(im+1,im), intent(inout) :: Hmat
    ! Define some local arguments
    real(kind = rk),dimension(in)       :: x_temp
    real(kind = rk)                     :: dotsum,dotsum_local
    integer(kind = ik)                  :: ic,index_dot
    integer(kind = ik)                  :: ierr

    ! subroutine content
    ! The linear operator is defined by user, 
    ! kstep loop is the main loop for arnoldi process
    ! w_j = A * q_j, w_j = Qmat(1:in,kstep+1)
     call LinearOperatorAX(in,Qmat(1:in,kstep),Qmat(1:in,kstep+1))

     do ic = 1, kstep

        dotsum_local = 0.d0;
        ! Hij = (A*q_j,q_i) = (w_j, q_i)
          ! naive implementation
          !do index_dot = 1, in
          !  dotsum_local = dotsum_local + Qmat(index_dot,ic)*Qmat(index_dot,kstep+1);
          !end do
          ! Blas 1 implementation
          dotsum_local = ddot(in,Qmat(1:in,ic),1,Qmat(1:in,kstep+1),1);
        
          call MPI_AllReduce(dotsum_local,dotsum,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr);
          ! j = kstep  
          ! Hji = (q_i, w_j)               
          Hmat(ic,kstep) = dotsum;
         
        ! w_j = w_j - Hij * v_i
        Qmat(1:in,kstep+1) = Qmat(1:in,kstep+1) - Hmat(ic,kstep)*Qmat(1:in,ic);

     enddo
        ! H(j+1,j) = ||w_j||
        Hmat(kstep+1,kstep) = norm_vec(Qmat(1:in,kstep+1));
        ! q_(j+1) = w_j / H(j+1,j) renormalize
        Qmat(1:in,kstep+1) = Qmat(1:in,kstep+1)/Hmat(kstep+1,kstep);

  end subroutine arnoldi

  subroutine apply_givens_rotation(H, cs, sn, k)
        implicit none
        real(kind=rk),    dimension(:,:), intent(inout)   :: H
        integer,                          intent(in)      :: k
        real(kind=rk),    dimension(:),   intent(inout)   :: cs, sn
      ! Local arguments
        real(kind=rk)    :: temp
        integer          :: i

        do i=1,k-1
            temp     =  cs(i)*H(i,k) + sn(i)*H(i+1,k)
            H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
            H(i,k)   =  temp
        end do

        if (H(k,k)==0) then
            cs(k) = 0.0D0
            sn(k) = 1.0D0
        else
            temp  = sqrt(H(k,k)**2 + H(k+1,k)**2)
            cs(k) = H(k,k) / temp
            sn(k) = cs(k) * H(k+1,k) / H(k,k)
        end if
        H(k,k)   = cs(k)*H(k,k) + sn(k)*H(k+1,k)
        H(k+1,k) = 0.0D0

  end subroutine apply_givens_rotation

  subroutine back_substitute(H, beta)
      implicit none
      real( kind = rk ), dimension(:,:), intent(in)    :: H
      real( kind = rk ), dimension(:),   intent(inout) :: beta
      real( kind = rk )                                :: dotsum
      integer :: i, k, index_dot  
      k = size(beta)
      beta(k) = beta(k)/H(k,k)
      do i=k-1,1,-1
          dotsum = 0.d0;
            ! naive implementation
            !do index_dot = i+1,k,1
            !  dotsum = dotsum + H(i,index_dot)*beta(index_dot)
            !enddo
            ! Blas 1 implementation
            dotsum = ddot(k-i,H(i,i+1:k),1,beta(i+1:k),1) 

          beta(i) = (beta(i) - dotsum)/H(i,i)
      end do

  end subroutine back_substitute

  subroutine gmres(n, b, x, m, tol, mrn, verbose)
    ! GMRES (Generalized Minimal Residual Method)迭代求解线性方程组的子程序
    ! 参数说明：
    ! n: 方程组的维度
    ! b: 线性方程组的右端向量
    ! x: 输入/输出的解向量
    ! m: Krylov子空间的维度
    ! tol: 收敛容差
    ! verbose: 是否输出详细迭代信息的标志
  
    implicit none
  
    ! 输入/输出参数声明 -----------------------------------------------------------
    integer,                              intent(in)    :: n    ! 方程维度
    integer,                              intent(in)    :: m    ! Krylov子空间维度
    integer,                              intent(in)    :: mrn  ! max_restarted_number  ! 最大重启次数
    integer,                              intent(in)    :: verbose  ! 输出详细信息标志
    real( kind = rk ),    dimension(n),   intent(in)    :: b    ! 方程右端向量
    real( kind = rk ),    dimension(n),   intent(inout) :: x    ! 解向量
    real( kind = rk ),                    intent(in)    :: tol  ! 收敛容差

  
    ! 局部变量声明 ----------------------------------------------------------------
    integer                             :: k,k_out,i,j  ! 循环和临时变量
    integer                             :: this_proc, ierr  ! MPI进程相关
    integer                             :: rest_num  ! 重启次数计数器
    real( kind = rk )                   :: b_norm    ! 右端向量的范数
    real( kind = rk )                   :: res_norm  ! 残差范数
    real( kind = rk )                   :: error,error_old  ! 误差值
    real( kind = rk )                   :: tempVy    ! 临时变量
    real( kind = rk ), dimension(n)     :: res, x_temp  ! 残差和临时解向量
    real( kind = rk ), dimension(m+1)   :: sn, cs    ! Givens旋转的正弦和余弦
    real( kind = rk ), dimension(n,m+1) :: Q         ! Q矩阵（正交基）
    real( kind = rk ), dimension(m+1,m) :: H         ! Hessenberg矩阵
    real( kind = rk ), dimension(m+1)   :: beta      ! beta向量
  
    ! 子程序主体 -----------------------------------------------------------------
      
    ! 获取当前MPI进程的秩
    call mpi_comm_rank(mpi_comm_world, this_proc, ierr)
  
    ! 初始解向量置零
    x = 0.d0
  
    ! 外层重启循环
    do rest_num = 1,mrn
        ! 重置所有临时变量,每次重启的时候都要重置
        res = 0.d0
        sn = 0.d0
        cs = 0.d0
        Q = 0.d0
        H = 0.d0
        beta = 0.d0

        
        ! 初始解向量置零
        if( rest_num == 1 ) then
          x_temp = 0.d0;
        else
          call LinearOperatorAX(n,x,x_temp)
        endif
  
        ! 计算右端向量范数
        b_norm = norm_vec(b)
  
        ! 计算初始残差
        res = b - x_temp                 ! r0 = b - A x0
        res_norm = norm_vec(res)         ! |r0|
  
        ! 计算相对误差
        error   = res_norm / b_norm      ! error = 1 if x0 = 0
        error_old = error                ! 
  
        ! 初始化beta向量和Q矩阵的第一列
        beta(1) = 1.d0                   ! beta0 = |r0|
        Q(1:n,1)  = res(1:n) / res_norm  ! q1 = r0 / |r0|
          
        ! Arnoldi迭代过程
        do k=1,m
           ! 执行Arnoldi迭代生成正交基，Hessenberg矩阵
           ! 实际上就是利用待定系数法配合Gram-Schmidt正交化过程
           call arnoldi(Q, H, k, n, m)
           
           ! 应用Givens旋转进行QR分解
           call apply_givens_rotation(H, cs, sn, k)
           
           ! 更新beta向量
           beta(k+1) = -sn(k)*beta(k)
           beta(k)   =  cs(k)*beta(k)
           
           ! 计算误差
           error = abs(beta(k+1))/b_norm
           
           ! 检查数值异常
           if(ISNAN(error)) then
             WRITE(*,*)'GMRES IS NAN. WRONG! WRONG!'
             stop
           endif
             
           ! 输出详细迭代信息
           if (this_proc==0 .and. verbose==1) then
               write(*,"(A,I1,A,I3,A,E14.6)")"Gmres Rest:",rest_num,",Subspace:",k,",Error:",error
           end if
             
           k_out = k
           
           ! 判断是否收敛
           if(error<tol) exit
        end do
          
        ! 回代求解最小二乘问题
        call back_substitute(H(1:k_out,1:k_out), beta(1:k_out))
          
        ! 更新解向量
        do i = 1,n
           tempVy = 0.d0 
           do j = 1,k
             tempVy = tempVy + Q(i,j)*beta(j) 
           enddo
           x(i) = x(i) + tempVy
        enddo
          
        ! 判断是否收敛并退出
        if(error<tol) exit
    enddo
  end subroutine gmres

  subroutine LinearOperatorAX(n,v_in,v_out)
    ! 这个子程序计算了对应的矩阵向量乘:
    !   n: 对应的是向量的维度                   
    !   v_in(n):   对应的输入向量
    !   v_out(n):  对应的输出向量
    !   每一个block上向量的维度为（三部分组成
    !   守恒变量（不包含边界），激波高度和激波速度）
    !   物面边界上的密度和压力会发生变化
    use SF_Constant,           only: ik,rk,NumVar,overLAP
    use SF_CFD_Global,         only: nx_local,Ny,nz_local
    implicit none
    integer( kind = ik ),                intent(in)    :: n
    real( kind = rk ), dimension(n),     intent(in)    :: v_in
    real( kind = rk ), dimension(n),     intent(out)   :: v_out
    ! local variables
    integer( kind = ik )         :: ic, jc, kc, iVar
    real( kind = rk )            :: NormOfv_in
    integer( kind = ik )         :: index_vector

   ! Step 1: 
     NormOfv_in = norm_vec(v_in);
     Fre_epsUse = Fre_eps / NormOfv_in;
    
    !Fre_epsUse = Fre_eps / sqrt(real(n,rk));

   ! Step 2: ===============================更新对应的求解变量==========================================
    !                                   \varphi + Fre_eps * v_in   
    ! 计算对应的通量函数R_snew，关于\varphi + Fre_eps * v_in
    ! Ucons0, ShockH, ShockV是参与通量计算的变量
    index_vector = 0;
    ! Flux Update
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            index_vector = index_vector + 1;
            Ucons0(iVar,ic,jc,kc) = UconsOld(iVar,ic,jc,kc) + Fre_epsUse * v_in(index_vector);
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! ShockH and ShockV
    DO kc = 1, nz_local
      DO ic = 1, nx_local
        index_vector = index_vector + 1;
        ShockH(ic,kc) = ShockH_old(ic,kc) + Fre_epsUse * v_in(index_vector);
        index_vector = index_vector + 1;
        ShockV(ic,kc) = ShockV_old(ic,kc) + Fre_epsUse * v_in(index_vector);
      ENDDO
    ENDDO
   ! Step 3: ===============================更新新的通量函数==========================================
    ! 更新对应的通量函数,注意我们还需要更新对应的网格系数矩阵
   CALL Calculate_Jaco
   CALL Update_Variables
   CALL Calculate_Flux
   CALL Calculate_ShockAC3D

   ! Step 4: ===============================计算对应的矩阵向量乘法=====================================
   index_vector = 0;
    ! Flux Update
    DO kc = 1, nz_local
     DO jc = 1, Ny-1
      DO ic = 1, nx_local
       DO iVar = 1, NumVar
        index_vector = index_vector + 1;
        v_out(index_vector) = v_in(index_vector)/DT0 - &
                              & (dUcons0(iVar,ic,jc,kc) - dUcons0_old(iVar,ic,jc,kc))/Fre_epsUse;
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! ShockH and ShockV
    DO kc = 1, nz_local
     DO ic = 1, nx_local
      index_vector = index_vector + 1;
      !v_out(index_vector) = v_in(index_vector)/DT0 - &
      !                      & (ShockV(ic,kc) - ShockV_old(ic,kc))/Fre_epsUse;
      v_out(index_vector) = v_in(index_vector)/DT0 - v_in(index_vector+1);
      index_vector = index_vector + 1;
      v_out(index_vector) = v_in(index_vector)/DT0 - &
                            & (ShockAc(ic,kc) -ShockAc_old(ic,kc))/Fre_epsUse;
     ENDDO
    ENDDO

  end subroutine LinearOperatorAX

  subroutine FormRHS(n,vec_b)
    ! 这个子程序用于计算右端向量
    !   vec_b(n):  对应的右端向量
    !   每一个block上向量的维度为（四部分组成，边界条件
    !   网格中的守恒变量（不包含边界），激波高度和激波速度）
    !   物面边界上的密度和压力会发生变化
    use SF_Constant,           only: ik,rk,NumVar,overLAP
    use SF_CFD_Global,         only: nx_local,Ny,nz_local
    implicit none
    integer( kind = ik ),   intent(in)   :: n
    real( kind = rk ),      dimension(n) :: vec_b
    ! local variables
    integer( kind = ik )         :: ic, jc, kc, iVar, index_vector

    ! Ucons0, ShockH, ShockV是参与通量计算的变量
    index_vector = 0;
    ! Flux Update
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            index_vector = index_vector + 1;
            vec_b(index_vector) = dUcons0(iVar,ic,jc,kc);
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    dUcons0_old = dUcons0;
    ! ShockH and ShockV
    DO kc = 1, nz_local
      DO ic = 1, nx_local
        index_vector = index_vector + 1;
        vec_b(index_vector) = ShockV(ic,kc);
        index_vector = index_vector + 1;
        vec_b(index_vector) = ShockAc(ic,kc);
      ENDDO
    ENDDO
    ShockV_old = ShockV;
    ShockAc_old = ShockAc;

  end subroutine FormRHS 

  subroutine GMRES_iteration
    ! 这里相当于是GMRES求解器和CFD求解器的打包，
    ! 我们在这里把所有需要的过程进行封装
    implicit none
    integer( kind = ik ):: index_vector
    integer( kind = ik ):: ic,jc,kc,iVar
    integer( kind = ik ):: output_gmres_debug = 0
    
    ! The present state of the solution
    UconsOld = Ucons0;
    ShockH_old= ShockH;
    ShockV_old= ShockV;

    call Calculate_Flux
    call Calculate_ShockAC3D

    ! save the present state of the residual
    dUcons0_old = dUcons0;  ! dU0/dt = R(U0)
    !ShockV_old = ShockV;   ! dH /dt = ShockV
    ShockAc_old = ShockAc;  ! dV /dt = ShockAc
    
    ! Form the right hand side
    call FormRHS(NSize, RHS_Krylov)

    ! Without the preconditioner
    call gmres(NSize, RHS_Krylov, DeltaSolu, MaxKrylovSubSpace, &
              &InnerGmresTol, MaxRestartedNumber, output_gmres_debug)
    
    ! Update the solution
    ! 为了形式上和其他Implicit方法保持一致，我们需要把对应的DeltaX转换为相应的增量 
    index_vector = 0;
    ! Flux Update
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            index_vector = index_vector + 1;
            dUcons(iVar,ic,jc,kc) = DeltaSolu(index_vector);
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! ShockH and ShockV
    DO kc = 1, nz_local
      DO ic = 1, nx_local
        index_vector = index_vector + 1;
        dShockH(ic,kc) = DeltaSolu(index_vector);;
        index_vector = index_vector + 1;
        dShockV(ic,kc) = DeltaSolu(index_vector);;
      ENDDO
    ENDDO     

    
      Ucons0 = UconsOld   + dUcons;
      ShockH = ShockH_old + dShockH;
      ShockV = ShockV_old + dShockV;

    ! 
    call Calculate_Jaco
    call Update_Variables
    

  end subroutine GMRES_iteration
  
  ! 这个是优化后的GMRES Solver，其通用性会大大下降，但是和当前的CFD求解器
  ! 结合的更好
  ! 为了求解效率，我们需要把Gmres过程直接和
  ! CFD求解过程进行打包，这样可以减少中间变量的存储和计算
  subroutine Gmres_SF_Solver
    use SF_CFD_Global,  only: BSS_A,BSS_U,BSS_V,HessenbergMat,Givens,Ve
    use mpi    ! MPI module, the Krylov subspace methods are parallelized using MPI
    IMPLICIT NONE ! 
    INTEGER( kind = ik ):: ic,jc,kc,iVar
    INTEGER( kind = ik ):: IM,JM,M_DIM,Lsub
    INTEGER::ierr
    REAL( kind = rk ):: Beta_local, Beta_global
    REAL( kind = rk ):: Prod_local, Prod_global
    REAL( kind = rk ):: H_sqrt, HIJ1, HIJ2, QJ1, QJ2
    REAL( kind = rk ):: Norm_cri
    REAL( kind = rk ):: eps_fre

    ! The present state of the solution
    UconsOld = Ucons0;
    ShockH_old= ShockH;
    ShockV_old= ShockV;

    ! Calculate the right hand side of the present state
    call Calculate_Flux
    call Calculate_ShockAC3D

    ! save the present state of the residual
    dUcons0_old = dUcons0;  ! dU0/dt = R(U0)
    !ShockV_old = ShockV;   ! dH /dt = ShockV
    ShockAc_old = ShockAc;  ! dV /dt = ShockAc

    ! define the max krylov subspace
    M_DIM = MaxKrylovSubSpace
    ! Norm Criterion for GMRES
    Norm_cri = sqrt(REAL(Nx * Nz * (NumVar * (Ny - 1) + 2),rk))
    ! GMRES solver
    ! r = b - A*x0, x0 = 0, r = b
    ! beta = norm(r0)

    Beta_local = 0.0_rk;
    DO kc = 1, nz_local
     DO ic = 1, nx_local
         Beta_local = Beta_local +  ShockV(ic,kc) *  ShockV(ic,kc);
         Beta_local = Beta_local + ShockAc(ic,kc) * ShockAc(ic,kc);
     ENDDO
    ENDDO
    DO kc = 1, nz_local
     DO jc = 1, Ny-1
      DO ic = 1, nx_local
       DO iVar = 1, NumVar
         Beta_local = Beta_local + dUcons0(iVar,ic,jc,kc) * dUcons0(iVar,ic,jc,kc);
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    ! Global operations
    call MPI_ALLReduce(Beta_local,Beta_global,1,MPI_DOUBLE_PRECISION,&
                      &MPI_SUM,MPI_COMM_WORLD,ierr);

    Beta_global = sqrt(Beta_global) / Norm_cri;

    ! v1 = r0 / norm(r0)
    ! v1 = r0 / beta
    DO kc = 1, nz_local
     DO ic = 1, nx_local
       BSS_V(ic,kc,1) = ShockV(ic,kc) / Beta_global;
       BSS_A(ic,kc,1) = ShockAc(ic,kc)/ Beta_global;
     ENDDO
    ENDDO
    DO kc = 1, nz_local
     DO jc = 1, Ny-1
      DO ic = 1, nx_local
       DO iVar = 1, NumVar
         BSS_U(iVar,ic,jc,kc,1) = dUcons0(iVar,ic,jc,kc) / Beta_global;
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    ! Arnoldi Process
    HessenbergMat = 0.0_rk;
    ! ve = V * e1;
    VE = 0.0_rk;
    VE(1) = 1.0_rk;
    
    KrylovSubItr: DO JM = 1, MaxKrylovSubSpace
      ! yj = A * vj
      ! Define the increment eps_fre
      eps_fre = Fre_eps / Norm_cri;
      ! Q + eps_fre * yj
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        ShockH(ic,kc) = ShockH_old(ic,kc) + eps_fre * BSS_V(ic,kc,JM);
        ShockV(ic,kc) = ShockV_old(ic,kc) + eps_fre * BSS_A(ic,kc,JM); 
       ENDDO
      ENDDO
      DO kc = 1, nz_local
       DO jc = 1, Ny-1
        DO ic = 1, nx_local
         DO iVar = 1, NumVar
           Ucons0(iVar,ic,jc,kc) = UconsOld(iVar,ic,jc,kc) + eps_fre * BSS_U(iVar,ic,jc,kc,JM);
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      ! Calculate the R(Q + eps_fre * yi)
      call Calculate_Jaco
      call Update_Variables
      call Calculate_Flux
      call Calculate_ShockAC3D

      ! ( R(Q + eps_fre * yi) - R(Q) ) / eps_fre
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        !ShockV(ic,kc) = BSS_V(ic,kc,JM) / DT0 - (ShockV(ic,kc) - ShockV_old(ic,kc)) / eps_fre;
         ShockV(ic,kc) = BSS_V(ic,kc,JM) / DT0 - BSS_A(ic,kc,JM);
        ShockAc(ic,kc) = BSS_A(ic,kc,JM) / DT0 - (ShockAc(ic,kc) - ShockAc_old(ic,kc)) / eps_fre;
       ENDDO
      ENDDO
      DO kc = 1, nz_local
       DO jc = 1, Ny-1
        DO ic = 1, nx_local
         DO iVar = 1, NumVar
           dUcons0(iVar,ic,jc,kc) = BSS_U(iVar,ic,jc,kc,JM) / DT0 - (dUcons0(iVar,ic,jc,kc) - dUcons0_old(iVar,ic,jc,kc)) / eps_fre;
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      !select case(Select_Kind_Precon_Solver)
      !CASE(1) ! None PC
      !  continue
      !CASE(2) ! DPLR
      !  call DPLUR_Line_Solver
      !  dUcons0 = dUcons
      !   ShockV = dShockH
      !  ShockAc = dShockV
      !end select
      
      ! Gram Schimdt Orthogonalization
      DO IM = 1, JM
        ! h(i,j) = (wj,vi)
        Prod_local = 0.0_rk;
        DO kc = 1, nz_local
         DO ic = 1, nx_local
          Prod_local = Prod_local +  ShockV(ic,kc) * BSS_V(ic,kc,IM);
          Prod_local = Prod_local + ShockAc(ic,kc) * BSS_A(ic,kc,IM);
         ENDDO
        ENDDO
        DO kc = 1, nz_local
         DO jc = 1, Ny-1
          DO ic = 1, nx_local
           DO iVar = 1, NumVar
             Prod_local = Prod_local + dUcons0(iVar,ic,jc,kc) * BSS_U(iVar,ic,jc,kc,IM);
           ENDDO
          ENDDO
         ENDDO
        ENDDO

        call MPI_ALLReduce(Prod_local,Prod_global,1,MPI_DOUBLE_PRECISION,&
                          &MPI_SUM,MPI_COMM_WORLD,ierr);
        Prod_global = Prod_global / (Norm_cri * Norm_cri);

        HessenbergMat(IM,JM) = Prod_global;

        ! wj = wj - h(i,j) * vi
        DO kc = 1, nz_local
         DO ic = 1, nx_local
          ShockV(ic,kc) = ShockV(ic,kc) - Prod_global * BSS_V(ic,kc,IM);
          ShockAc(ic,kc) = ShockAc(ic,kc) - Prod_global * BSS_A(ic,kc,IM);
         ENDDO
        ENDDO
        DO kc = 1, nz_local
         DO jc = 1, Ny-1
          DO ic = 1, nx_local
           DO iVar = 1, NumVar
             dUcons0(iVar,ic,jc,kc) = dUcons0(iVar,ic,jc,kc) - Prod_global * BSS_U(iVar,ic,jc,kc,IM);
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDDO !MGSItr

        ! h(j+1,j) = norm(wj)
        Prod_local = 0.0_rk;
        DO kc = 1, nz_local
         DO ic = 1, nx_local
          Prod_local = Prod_local +  ShockV(ic,kc) *  ShockV(ic,kc);
          Prod_local = Prod_local + ShockAc(ic,kc) * ShockAc(ic,kc);
         ENDDO
        ENDDO
        DO kc = 1, nz_local
         DO jc = 1, Ny-1
          DO ic = 1, nx_local
           DO iVar = 1, NumVar
             Prod_local = Prod_local + dUcons0(iVar,ic,jc,kc) * dUcons0(iVar,ic,jc,kc);
           ENDDO
          ENDDO
         ENDDO
        ENDDO

        call MPI_ALLReduce(Prod_local,Prod_global,1,MPI_DOUBLE_PRECISION,&
                          &MPI_SUM,MPI_COMM_WORLD,ierr);

        HessenbergMat(JM+1,JM) = sqrt(Prod_global) / Norm_cri;

        if(JM < MaxKrylovSubSpace) then
          ! v(j+1) = wj / h(j+1,j)
          DO kc = 1, nz_local
           DO ic = 1, nx_local
            BSS_V(ic,kc,JM+1) =  ShockV(ic,kc) / HessenbergMat(JM+1,JM);
            BSS_A(ic,kc,JM+1) = ShockAc(ic,kc) / HessenbergMat(JM+1,JM);
           ENDDO
          ENDDO
          DO kc = 1, nz_local
           DO jc = 1, Ny-1
            DO ic = 1, nx_local
             DO iVar = 1, NumVar
               BSS_U(iVar,ic,jc,kc,JM+1) = dUcons0(iVar,ic,jc,kc) / HessenbergMat(JM+1,JM);
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        end if

        ! Apply Givens Rotation to obtaing results
        ! h(:,j) = prod(Givens) * h(:,j)
        DO IM = 1, JM-1
          HIJ1 = GIVENS(1,IM) * HessenbergMat(IM,JM) + GIVENS(2,IM) * HessenbergMat(IM+1,JM);
          HIJ2 =-GIVENS(2,IM) * HessenbergMat(IM,JM) + GIVENS(1,IM) * HessenbergMat(IM+1,JM);
          HessenbergMat(IM  ,JM) = HIJ1;
          HessenbergMat(IM+1,JM) = HIJ2;
        ENDDO
        ! h_sqrt = sqrt(h(j,j)^2 + h(j+1,j)^2)
        H_sqrt = sqrt(HessenbergMat(JM,JM)  *HessenbergMat(JM,JM) &
                   &+ HessenbergMat(JM+1,JM)*HessenbergMat(JM+1,JM));
        Givens(1,JM) = HessenbergMat(JM,JM) / H_sqrt;   ! cos
        Givens(2,JM) = HessenbergMat(JM+1,JM) / H_sqrt; ! sin

        ! h(:,j) = Givens_j * h(:,j)
        HessenbergMat(JM  ,JM) = H_sqrt;
        HessenbergMat(JM+1,JM) = 0.0_rk;

        ! Check the convergence
        VE(JM+1) = -Givens(2,JM) * VE(JM);
        VE(JM  ) =  Givens(1,JM) * VE(JM);

        IF(abs(VE(JM+1)) < InnerGmresTol) then
          M_DIM = JM;
          EXIT KrylovSubItr
        ENDIF

    ENDDO KrylovSubItr

    ! Least Square Solution
    ! ve = V * e1 * beta
    DO IM = 1, M_DIM + 1
      Ve(IM) = Beta_global * Ve(IM);
    ENDDO
    ! Solving the H*y = ve
    DO IM = M_DIM,1,-1
      DO JM = M_DIM, IM+1, -1
        Ve(IM) = Ve(IM) - HessenbergMat(IM,JM) * Ve(JM);
      ENDDO
      Ve(IM) = Ve(IM) / HessenbergMat(IM,IM);
    ENDDO
    ! Update the solution
    DO kc = 1, nz_local
     DO ic = 1, nx_local
       ShockV(ic,kc) = 0.0_rk;
       ShockAc(ic,kc) = 0.0_rk;
       DO Lsub = 1, M_DIM
            ShockV(ic,kc) = ShockV(ic,kc) + BSS_V(ic,kc,Lsub) * Ve(Lsub);
          ShockAc(ic,kc) = ShockAc(ic,kc) + BSS_A(ic,kc,Lsub) * Ve(Lsub);
       ENDDO
       dShockH(ic,kc) = ShockV(ic,kc);
       dShockV(ic,kc) = ShockAc(ic,kc);
       ShockH(ic,kc) = ShockH_old(ic,kc) + dShockH(ic,kc);
       ShockV(ic,kc) = ShockV_old(ic,kc) + dShockV(ic,kc);
     ENDDO
    ENDDO
    DO kc = 1, nz_local
     DO jc = 1, Ny-1
      DO ic = 1, nx_local
       DO iVar = 1, NumVar
         dUcons0(iVar,ic,jc,kc) = 0.0_rk;
         DO Lsub = 1, M_DIM
           dUcons0(iVar,ic,jc,kc) = dUcons0(iVar,ic,jc,kc) + BSS_U(iVar,ic,jc,kc,Lsub) * Ve(Lsub);
         ENDDO
         dUcons(iVar,ic,jc,kc) = dUcons0(iVar,ic,jc,kc);
         Ucons0(iVar,ic,jc,kc) = UconsOld(iVar,ic,jc,kc) + dUcons(iVar,ic,jc,kc);
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    call Calculate_Jaco
    call Update_Variables

  end subroutine Gmres_SF_Solver

END Module Krylov_SubSpace