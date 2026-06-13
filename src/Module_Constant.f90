MODULE SF_Constant
  ! This is the constant used for SF solver
  implicit none

  ! define the precisions
  INTEGER, PARAMETER :: ik = 4; ! ik = 8 for int64
  INTEGER, PARAMETER :: rk = 8; 

  ! Constant 
  REAL( kind = rk ), PARAMETER:: pi = 4.0_rk * ATAN(1.0_rk) ! Define the constant PI
  REAL( kind = rk ), PARAMETER::  Converge_Tol = 0.5e-14    ! Define the Converge Tol
  REAL( kind = rk ), PARAMETER::     eps_small = 1.D-11     ! Small Constant 
  REAL( kind = rk ), PARAMETER::   LaxFriSmall = 5.D-2      ! Small Constant for Lax-flux-split
  REAL( kind = rk ), PARAMETER::   Csutherland = 110.4_rk   ! Constant used for sutherland law
  REAL( kind = rk ), PARAMETER::       Fre_eps = 1.D-6      ! Small Constant for Frechet derivatives
  REAL( kind = rk ), PARAMETER:: InnerGmresTol = 0.1_rk     ! Inner tol for gmres solver
  REAL( kind = rk ), PARAMETER:: InnerBicgTol  = 0.1_rk     ! Inner tol for bicgstab solver

  INTEGER( kind = ik ),PARAMETER:: NumVar = 5_ik            ! Variables Density,U,V,W,T

  ! Define some constant for Krylov subspace method
  INTEGER( kind = ik ), PARAMETER:: MaxKrylovSubSpace = 10_ik  ! Krylov subspace dimensions
  INTEGER( kind = ik ), PARAMETER:: MaxRestartedNumber = 2_ik  ! Max restart numbers

  ! Define the filename for continue
  CHARACTER( len = 14 ), PARAMETER:: FilesForContinue = "SF_Results"
  CHARACTER( len = 14 ), PARAMETER:: PertForContinue  = "Pert_Results"

  ! Define the gost mesh layer
  INTEGER( kind = ik ), PARAMETER:: overLAP = 3_ik; 

  ! Define the scheme
  INTEGER( kind = ik ), PARAMETER::scheme_upwind = -1_ik
  INTEGER( kind = ik ), PARAMETER::scheme_center =  0_ik
  INTEGER( kind = ik ), PARAMETER::scheme_dowind =  1_ik 

  ! Define the Implicit Solver
  INTEGER( kind = ik ), PARAMETER::     RK4_Methods = 1_ik
  INTEGER( kind = ik ), PARAMETER::   DPLUR_Methods = 2_ik
  INTEGER( kind = ik ), PARAMETER::   GMRES_Methods = 3_ik
  INTEGER( kind = ik ), PARAMETER::MultiGridMethods = 4_ik
  
  INTEGER( kind = ik ), PARAMETER:: dplur_sub_iter  = 3_ik

  ! Define the Precondition Solver
  INTEGER( kind = ik ), PARAMETER::      No_precond = 0_ik
  INTEGER( kind = ik ), PARAMETER::    DPLR_precond = 1_ik
  INTEGER( kind = ik ), PARAMETER::MultGrid_precond = 2_ik

  ! Define the spatial scheme
  INTEGER( kind = ik ), PARAMETER::     first_order = 1_ik  ! 1st + 2nd order
  INTEGER( kind = ik ), PARAMETER::     fifth_order = 5_ik  ! 5st + 6nd order
  INTEGER( kind = ik ), PARAMETER::     weno5_order = 6_ik  ! weno5 + 6nd order
  INTEGER( kind = ik ), PARAMETER::     hybird_sch  = 7_ik  ! fifth + weno5 + 6nd order

  ! Define Analysis Type 
  INTEGER( kind = ik ), PARAMETER::  Steady_NS_Analysis = 1_ik
  INTEGER( kind = ik ), PARAMETER::Unsteady_NS_Analysis = 2_ik
  INTEGER( kind = ik ), PARAMETER::        LNS_Analysis = 3_ik
  INTEGER( kind = ik ), PARAMETER::       HLNS_Analysis = 4_ik

  ! Define the default number of harmonic for HLNS analysis
  ! 2m + 1 modes are used for analysis.  DefaultHarmonic = 4 means 9 modes are used
  ! based on given -4w0, -3w0, -2w0, -w0, 0, w0, 2w0, 3w0, 4w0
  INTEGER( kind = ik ), PARAMETER:: DefaultHarmonic = 4_ik 


END MODULE SF_Constant