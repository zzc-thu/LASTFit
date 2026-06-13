!! 这里需要额外注意对应的PETSc内存向量的管理顺序，以及对应的自己程序的MPI进程对应的
!! 数据分布方式，实际情况需要根据实际情况进行调整。为了保证不出错，我们选用最直接的方法
!! 我们根据我们的MPI进程分布方式，直接手动指定对应的PETSc向量分区，手动指定每一个对应
!! 的MPI进程对应的局部数据分布方式和原始CFD求解器中的分布方式一致。
!
!MODULE Eigen_Solver_Module
!  use petsc
!  use slepc
!  use iso_c_binding
!  implicit none
!
!  ! 用户上下文类型（包含流场数据）
!  type :: UserContext
!    real(8), pointer :: Q(:)       ! 流场变量（假设一维数组）
!    integer :: n_local             ! 当前进程的局部数据量
!    integer :: N_global            ! 全局数据量
!  end type UserContext
!
!contains
!
!  ! 封装 dRdQx 为 PETSc 回调函数
!  subroutine compute_dRdQx_petsc(x_petsc, y_petsc, ctx, ierr) bind(C)
!    use iso_c_binding
!    use petscVec
!    type(c_ptr), value :: x_petsc, y_petsc, ctx
!    integer, intent(out) :: ierr
!    type(UserContext), pointer :: user_ctx
!    real(8), pointer :: x_local(:), y_local(:)
!
!    ! 获取用户上下文
!    ! c_f_pointer 将 C 指针转换为 Fortran 指针 F2003
!    call c_f_pointer(ctx, user_ctx)
!
!    ! 获取 PETSc 向量的本地数据指针
!    call VecGetArrayReadF90(x_petsc, x_local, ierr)
!    call VecGetArrayF90(y_petsc, y_local, ierr)
!
!    ! 调用 CFD 代码的 dRdQx
!    call dRdQx(x_local(1:user_ctx%n_local), y_local(1:user_ctx%n_local))
!
!    ! 恢复 PETSc 向量
!    call VecRestoreArrayReadF90(x_petsc, x_local, ierr)
!    call VecRestoreArrayF90(y_petsc, y_local, ierr)
!  end subroutine
!
!  ! 初始化特征值求解器
!  subroutine eigen_solver_init(ctx, n_local, N_global, eps, A, x, y, ierr)
!    type(UserContext), target, intent(inout) :: ctx
!    integer, intent(in) :: n_local, N_global
!    EPS, intent(out) :: eps
!    Mat, intent(out) :: A
!    Vec, intent(out) :: x, y
!    integer, intent(out) :: ierr
!
!    ! 设置用户上下文
!    ctx%n_local = n_local
!    ctx%N_global = N_global
!
!    ! 创建 PETSc 向量（与 CFD 数据分布一致）
!    call VecCreateMPI(PETSC_COMM_WORLD, n_local, N_global, x, ierr)
!    call VecDuplicate(x, y, ierr)
!
!    ! 创建壳矩阵 A，关联到 compute_dRdQx_petsc
!    call MatCreateShell(PETSC_COMM_WORLD, n_local, n_local, N_global, N_global, c_loc(ctx), A, ierr)
!    call MatShellSetOperation(A, MATOP_MULT, c_funloc(compute_dRdQx_petsc), ierr)
!
!    ! 创建特征值求解器
!    call EPSCreate(PETSC_COMM_WORLD, eps, ierr)
!    call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)
!    call EPSSetProblemType(eps, EPS_NHEP, ierr)  ! 非对称问题
!    call EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE, ierr)
!    call EPSSetTolerances(eps, 1e-6, 100, ierr)
!    call EPSSetFromOptions(eps, ierr)
!  end subroutine
!
!  ! 执行特征值求解
!  subroutine eigen_solver_solve(eps, ierr)
!    EPS, intent(inout) :: eps
!    integer, intent(out) :: ierr
!    call EPSSolve(eps, ierr)
!  end subroutine
!
!  ! 提取特征值
!  subroutine eigen_solver_get_results(eps, eigenvalues, ierr)
!    EPS, intent(in) :: eps
!    complex(8), allocatable, intent(out) :: eigenvalues(:)
!    integer, intent(out) :: ierr
!    integer :: nconv, i
!    real(8) :: kr, ki
!    Vec :: vr, vi
!
!    call EPSGetConverged(eps, nconv, ierr)
!    allocate(eigenvalues(nconv))
!    do i = 1, nconv
!      call EPSGetEigenpair(eps, i-1, kr, ki, vr, vi, ierr)
!      eigenvalues(i) = cmplx(kr, ki, kind=8)
!    end do
!  end subroutine
!
!  ! 清理资源
!  subroutine eigen_solver_finalize(eps, A, x, y, ierr)
!    EPS, intent(inout) :: eps
!    Mat, intent(inout) :: A
!    Vec, intent(inout) :: x, y
!    integer, intent(out) :: ierr
!    call EPSDestroy(eps, ierr)
!    call MatDestroy(A, ierr)
!    call VecDestroy(x, ierr)
!    call VecDestroy(y, ierr)
!  end subroutine
!
!END MODULE Eigen_Solver_Module