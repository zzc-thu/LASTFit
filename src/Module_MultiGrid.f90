Module MultiGrid
  use SF_Constant, only: ik, rk
  ! This module is used to calculate the MultiGrid method
  implicit none

  type :: TwoLevelGrid
    ! 细网格 （Level 1）
    real( kind = rk ), allocatable :: CV0(:,:,:,:)
    real( kind = rk ), allocatable :: ShockH0(:,:)
    real( kind = rk ), allocatable :: ShockV0(:,:)
    ! 细网格上的网格数目
    integer( kind = ik ) :: Nxi0, Neta0, Nzeta0
    ! 细网格上的残差
    real( kind = rk ), allocatable :: Res0_CV(:,:,:,:)
    real( kind = rk ), allocatable :: Res0_H(:,:)
    real( kind = rk ), allocatable :: Res0_V(:,:)

    ! 粗网格 （Level 2）
    real( kind = rk ), allocatable :: CV1(:,:,:,:)
    real( kind = rk ), allocatable :: ShockH1(:,:)
    real( kind = rk ), allocatable :: ShockV1(:,:)
    ! 粗网格上的网格数目
    integer( kind = ik ) :: Nxi1, Neta1, Nzeta1
    ! 粗网格上的残差
    real( kind = rk ), allocatable :: Res1_CV(:,:,:,:)
    real( kind = rk ), allocatable :: Res1_H(:,:)
    real( kind = rk ), allocatable :: Res1_V(:,:)

  end type TwoLevelGrid

! Here we solving the system of equations:
!   d | Cv |   | R  |
! --- | H0 | = | V0 |
!  dt | V0 |   | As |

contains

  subroutine Calculate_Fine_Residual
    ! This subroutine calculates the residual on the fine grid
    ! Input: CV0, ShockH0, ShockV0
    ! Output: Res0_CV, Res0_H, Res0_V
    
    ! 整个计算过程分为如下几个步骤：
    ! 1. 计算对应的激波高度下的网格相关系数，比如说Jacobian矩阵，转换关系等等。
    ! 2. 计算在此状态下的通量函数获得此网格下的残差向量
    
    
    
    
  end subroutine Calculate_Fine_Residual

  subroutine Calculate_Coarse_Residual
    
  
  end subroutine Calculate_Coarse_Residual

END Module MultiGrid