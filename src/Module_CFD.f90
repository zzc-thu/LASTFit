MODULE SF_CFD_Global
  use SF_Constant, only: ik,rk,NumVar
  use SF_Constant, only: LNS_Analysis,Unsteady_NS_Analysis
  implicit none
  
  LOGICAL             :: IfConverge=.FALSE. ! logical judge of the converge
  INTEGER( kind = ik ):: ModelType          ! Model Type: 1: 2.5D with homogeneous direction z, 
                                            !             2: 3D model with singular lines
  INTEGER( kind = ik ):: Nx,Ny,Nz           ! Define the total grid numbers
  INTEGER( kind = ik ):: nx_local           ! local nx grid points
  INTEGER( kind = ik ):: nz_local           ! local nz grid points
  INTEGER( kind = ik ):: TotalTimeSteps     ! Total Calculation steps
  INTEGER( kind = ik ):: Iteration          ! Iteration index
  INTEGER( kind = ik ):: Iteration_start    ! Iteration index
  INTEGER( kind = ik ):: LNS_start    ! Iteration index  
  INTEGER( kind = ik ):: StepsOutputScreen  ! steps to output screen
  INTEGER( kind = ik ):: StepsWriteData     ! steps to write data to file
  INTEGER( kind = ik ):: StepChangeCFLStart ! StartSteps to change CFL
  INTEGER( kind = ik ):: StepChangeCFL      ! Steps to change CFL from CFL_Start to CFL_end   
  
  REAL( kind = rk ):: CFL                   ! The CFL numbers used for Calculation
  REAL( kind = rk ):: CFL_Start, CFL_End    ! The CFL_Start value and CFL_end value
  REAL( kind = rk ):: DT0                   ! The Delta Time Step
  REAL( kind = rk ):: Time_old,Time_new,Time_RK0,Time_RK1,Time_RK2,Time_RK3! Time variables
  
  REAL( kind = rk ):: AoA
  REAL( kind = rk ):: AoA_Deg               ! Angle of attack, in [rad] and Deg 
  REAL( kind = rk ):: AoS
  REAL( kind = rk ):: AoS_Deg               ! Angle of Sweep,  in [rad] and Deg
  
  INTEGER( kind = rk ):: Geometry_Type      ! Geometry_Type  
  REAL( kind = rk ):: theta_p_deg           ! wall angle in Deg
  REAL( kind = rk ):: theta_p_rad           ! wall angle in Rad
  REAL( kind = rk ):: Cp_Ref                ! Ref Cp
  REAL( kind = rk ):: xim                   ! Cluster coefficients for Round head
  REAL( kind = rk ):: ximax = 1.0_rk        ! 
  REAL( kind = rk ):: sigma_y               ! The first compress coefficients, along the wall normal direction
  REAL( kind = rk ):: yim                   ! The second compress parameters, along the wall normal direction
  REAL( kind = rk ):: ymax = 1.0_rk         ! 
  REAL( kind = rk ):: Lz_start              ! The start position of the z direction
  REAL( kind = rk ):: Lz_end                ! The end position of the z direction
  REAL( kind = rk ):: R_head,LModel         ! Radius of the head, and the length of the model
  INTEGER(kind =ik):: Nx_flate              !
  INTEGER(kind =ik):: Nx_head               !
  INTEGER(kind =ik):: Nx_head_hf            !


  REAL( kind = rk ):: Re_Ref                ! Ref Reynolds Numbers
  REAL( kind = rk ):: Te_Ref                ! Ref Temperatures
  REAL( kind = rk ):: Pr_Ref                ! Ref Prandtl Number
  REAL( kind = rk ):: TWall                 ! Ref Wall Temperatures
  REAL( kind = rk ):: Gamma                 ! Ratio of specific heat
  REAL( kind = rk ):: Mach_Ref              ! Ref Mach Numbers
  REAL( kind = rk ):: Mach_Normal           ! Ref Normal Mach Numbers
  REAL( kind = rk ):: Csthlnd_Ref           ! Ref Constant for Sutherland law

  ! Freestream variables
  REAL( kind = rk )::Rho_INF                ! Steady Freestream Density
  REAL( kind = rk ):: U_INF                 ! Steady Freestream Velocity along x direction
  REAL( kind = rk ):: V_INF                 ! Steady Freestream Velocity along y direction   
  REAL( kind = rk ):: W_INF                 ! Steady Freestream Velocity along z direction
  REAL( kind = rk ):: T_INF                 ! Steady Freestream Temperature
  REAL( kind = rk )::T0_INF                 ! Steady Stag Temperatures
  REAL( kind = rk ):: P_INF                 ! Steady Freestream Pressure
  REAL( kind = rk ):: C_INF                 ! Steady Freestream wave speed
  REAL( kind = rk )::  CV_INF(NumVar)       ! Freestream Conservative variables
  REAL( kind = rk )::FINV_INF(NumVar)       ! Freestream inviscid Flux along X directions 
  REAL( kind = rk )::GINV_INF(NumVar)       ! Freestream inviscid Flux along Y directions
  REAL( kind = rk )::HINV_INF(NumVar)       ! Freestream inviscid Flux along Z directions  
  
  ! Define the global residual 
  REAL( kind = rk ):: AVE_RES_Rho           ! average residual for all process
  REAL( kind = rk ):: AVE_RES_RhoU          ! average residual for RhoU           
  REAL( kind = rk ):: AVE_RES_RhoV          ! average residual for RhoV
  REAL( kind = rk ):: AVE_RES_RhoW          ! average residual for RhoW
  REAL( kind = rk ):: AVE_RES_RhoE          ! average residual for RhoE

  ! Define the output file No.
  INTEGER( kind = ik ):: OutputFileNo        ! Output file No. start from 0

  ! All the logical integer index is putted here
  INTEGER( kind = ik ):: IF_Continue_Calculate      ! Select new calculation or continue
  INTEGER( kind = ik ):: Select_Kind_Solver         ! Select flow solver
  INTEGER( kind = ik ):: Select_Kind_Precon_Solver  ! Select Precondition solver 
  INTEGER( kind = ik ):: Output_Residual_Mode       ! Select brief / detail residual output
  
  ! Define the Scheme selection
  INTEGER( kind = ik ):: Select_Scheme              ! Select the spatial scheme
  
  ! Define the Krylov parameters
  INTEGER( kind = ik ):: Nsize                      ! The size of the matrix

  ! Here we define the allocate arrays that used for SF Solver
  
  ! Mesh and Geometry Variables

  ! Define the Analysis type:
  INTEGER( kind = ik ):: AnalysisType

  INTEGER( kind = ik ):: IF_FixedDeltaTime

  REAL( kind = rk ):: FixedDeltaTime
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::Sxi_Surface,sxi_surface_total
  ! Wall Surface deg
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSX,WallSY,WallSZ 
  ! Wall Surface X, Y, Z
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSX_total,WallSY_total,WallSZ_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSXdxi,WallSYdxi,WallSZdxi                          
  ! Wall Surface X, Y, Z derivatives with xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSXdzeta,WallSYdzeta,WallSZdzeta                    
  ! Wall Surface X, Y, Z derivatives with zeta 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSXdxi_total,WallSYdxi_total,WallSZdxi_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallSXdzeta_total,WallSYdzeta_total,WallSZdzeta_total
  
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallNormalX,WallNormalY,WallNormalZ                    
  ! Wall normal vector (nwx,nwy,nwz)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallNormalX_total,WallNormalY_total,WallNormalZ_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallNormalXdxi,WallNormalYdxi,WallNormalZdxi           
  ! Wall normal vector (nwx,xi,nwy,xi,nwz,xi) derivatives with xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta     
  ! Wall normal vector (nwx,zeta,nwy,zeta,nwz,zeta) derivatives with zeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)  ::WallNormalXdxi_total,WallNormalYdxi_total,WallNormalZdxi_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)  ::WallNormalXdzeta_total,WallNormalYdzeta_total,WallNormalZdzeta_total
  
  ! 这里我们需要把对应的Heta定义为一个三维数组，因为对应的计算空间中的网格间距都是1，
  ! 但是对应的生成网格的时候，我们对应的eta并不是严格的等距离为1的网格
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Heta,HetaDeta,HetaDxi,HetaDzeta                                           
  ! distribution function and derivatives to parameters eta and xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::HetaDeta_total,HetaDxi_total,HetaDzeta_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::X_grid,Y_grid,Z_grid,Hdst,Hdst_steady          
  ! X,Y,Z Grid and the distance away from wall surface
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::X_grid_total,Y_grid_total,Z_grid_total
  
  ! Geom transformation coefficients: 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::dxdxi,dxdeta,dxdzeta,dxdtau            
  ! x(xi,eta,zeta,tau)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::dydxi,dydeta,dydzeta,dydtau            
  ! y(xi,eta,zeta,tau)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::dzdxi,dzdeta,dzdzeta,dzdtau            
  ! z(xi,eta,zeta,tau)

  ! Geom transformation coefficients:
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::dxidx,dxidy,dxidz,dxidt                
  ! xi(x,y,z,t)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::detadx,detady,detadz,detadt            
  ! eta(x,y,z,t)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::dzetadx,dzetady,dzetadz,dzetadt        
  ! zeta(x,y,z,t)     
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::nablaxi,nablaeta,nablazeta             
  ! sqrt(xix*xix+xiy*xiy+xiz*xiz)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::jaco,invJacodt                         
  ! Jacobian, and d/dt(1/J) 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::detadt_invJ,detadt_invJdeta 
  ! deta/dt/J, d/deta(deta/dt/J)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)  ::ax_tau,ay_tau,az_tau                   !

  ! Flow Variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::Ucons0,UconsOld,UconsNew,Ucons0_total                    
  ! Conservative variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::dUcons0,dUcons0_old                                     
  ! increament conservative variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::dUcons0_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::dUcons                                      
  ! increament  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rho,U,V,W,T,P,Mu                              
  ! Basic variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rho_total,U_total,V_total,W_total,T_total,P_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::BigU_Xi,BigU_Eta,BigU_Zeta                    
  ! BigU_Xi = u*xix + v*xiy + xit, BigV_Xi = u*etax + v*etay + etat
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::BigU_Xi_total,BigU_Eta_total,BigU_Zeta_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Cs_Xi,Cs_Eta,Cs_Zeta     
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rds_Xi,Rds_Eta,Rds_Zeta,Rds_Max               
  ! Spectral Radius along Xi,Eta and Max values
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::LocalDT                                       
  ! Local Delta Timesteps 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Cs                                            
  ! Speed of Sound
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhodx,Udx,Vdx,Wdx,Tdx                         
  ! derivatives with respect to x
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhody,Udy,Vdy,Wdy,Tdy                         
  ! derivatives with respect to y
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhodz,Udz,Vdz,Wdz,Tdz                         
  ! derivatives with respect to z
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhodxi,Udxi,Vdxi,Wdxi,Tdxi                    
  ! derivatives to xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhodeta,Udeta,Vdeta,Wdeta,Tdeta               
  ! derivatives to eta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Rhodzeta,Udzeta,Vdzeta,Wdzeta,Tdzeta          
  ! derivatives to zeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Udxi_total,Udeta_total,Udzeta_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Vdxi_total,Vdeta_total,Vdzeta_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Wdxi_total,Wdeta_total,Wdzeta_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::Tdxi_total,Tdeta_total,Tdzeta_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::invF,invG,invH,visF,visG,visH               
  ! inviscid flux F,G,H; viscous flux F,G,H
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::invFhat,invGhat,invHhat                     
  ! Fhat = F*xix/j + G*xiy/j + H*xiz/j;
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::visFhat,visGhat,visHhat   
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::visFhat_total,visGhat_total,visHhat_total

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::Fp,Fm,Gp,Gm,Hp,Hm                           
  ! Fhat = Fp + Fm; Ghat;Hhat;
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::Fpdxi,Fmdxi,Gpdeta,Gmdeta,Hpdzeta,Hmdzeta   
  ! Fpdxi,Fmdxi,Gpdeta,Gmdeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::invFhatdxi_new,invGhatdeta_new
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::visFcdxi,visGcdeta,visHcdzeta               
  ! viscous flux derivatives to xi and eta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::visFcdxi_total,visGcdeta_total,visHcdzeta_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::invFlux,visFlux                             
  ! invFlux = Fpdxi + Fmdxi + Gpdeta + Gmdeta + Hpdzeta + Hmdzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::invFlux_total,visFlux_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::Cf_Surface,Stan_Surface                         
  ! surface skin friction coefficients and Stan number
  
  ! Shock Dynamics
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockVN                                    
  ! Shock Velocity along shock normal directions
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockH,ShockV,ShockAc
  ! Shock Height, Velocity and Acceleration
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockH_total,ShockV_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockH_old,ShockV_old,ShockAc_old

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockH0,ShockV0
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::dShockH,dShockV                            
  ! delta Shock Height, Velocity
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockHdxi,ShockVdxi                        
  ! Shock Height, Velocity derivatives to xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockHdzeta,ShockVdzeta                    
  ! Shock Height, Velocity derivatives to zeta                      
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockHdxi_total,ShockVdxi_total
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockHdzeta_total,ShockVdzeta_total
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::FCV_SHK,GCV_SHK,HCV_SHK                    
  ! inviscid Flux behind the shock surface
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockNormalX,ShockNormalY,ShockNormalZ     
  ! Shock Normal Vector (X,Y,Z)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)::ShockXtau,ShockYtau,ShockZtau                
  ! Shock Velocity along X,Y,Z
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::LeftEigenVectorsEta                        
  ! Left EigenVectors for the largest EigenValues along Eta

  ! Flux Jacobian
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)::DIAG                                     
  ! Flux Jacobian DIAG values
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DFhatDU                              
  ! inviscid flux jacobian DFhat/DU, Fhat = (F *  xix + G *  xiy + H * xiz + U0 * xit)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DGhatDU                              
  ! inviscid flux jacobian DGhat/DU, Ghat = (F * etax + G * etay + H * etaz + U0 * etat)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DHhatDU                              
  ! inviscid flux jacobian DHhat/DU, Hhat = (F *zetax + G *zetay + H *zetaz + U0 * zetat)

  ! viscous flux Jacobian
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DvisFhatDp,DvisFhatDpxi,DvisFhatDpet,DvisFhatDpzt
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DvisGhatDp,DvisGhatDpxi,DvisGhatDpet,DvisGhatDpzt
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::DvisHhatDp,DvisHhatDpxi,DvisHhatDpet,DvisHhatDpzt
  
  ! The RK4 variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::dUconsRK1,dUconsRK2,dUconsRK3,dUconsRK4
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)::Ucons_RK1,Ucons_RK2,Ucons_RK3
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)    ::ShockAc_RK1,ShockAc_RK2,ShockAc_RK3,ShockAc_RK4
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)    ::ShockV_RK1,ShockV_RK2,ShockV_RK3

  ! The Full implicit variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::DRDH,DRDV,DRDHxi,DRDHzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      ::ShockAcdH,ShockAcdHdxi,ShockAcdHdzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      ::ShockAcdV,ShockAcdVdxi,ShockAcdVdzeta

  ! 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DUsDH,DUsDHdxi,DUsDHdzeta,DUsDV

  ! D(xix/J)/DH,D(xiy/J)/DH,D(xiz/J)/DH
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DxixJDH,DxiyJDH,DxizJDH
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DetaxJDH,DetayJDH,DetazJDH,DetatDH
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DzetaxJDH,DzetayJDH,DzetazJDH
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DJDH,DetatJDH,DJDHxi,DJDHzeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DetatJDV

  !d(etax/J)/dHxi,d(etay/J)/dHxi,d(etaz/J)/dHxi,d(etat/J)/dHxi
  !d(etax/J)/dHzeta,d(etay/J)/dHzeta,d(etaz/J)/dHzeta,d(etat/J)/dHzeta
  !d(xix/J)/dHxi,d(xiy/J)/dHxi,d(xiz/J)/dHxi
  !d(xix/J)/dHzeta,d(xiy/J)/dHzeta,d(xiz/J)/dHzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DetaxJDHxi,  DetayJDHxi,  DetazJDHxi,  DetatJDHxi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DetaxJDHzeta,DetayJDHzeta,DetazJDHzeta,DetatJDHzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DzetaxJDHxi, DzetayJDHxi, DzetazJDHxi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::DxixJDHzeta, DxiyJDHzeta, DxizJDHzeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::F_IM, G_IM, H_IM
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::DF_IMDxi, DF_IMDeta, DF_IMDzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::DG_IMDxi, DG_IMDeta, DG_IMDzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::DH_IMDxi, DH_IMDeta, DH_IMDzeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::DU_Deta

  ! The data parallel variables
  ! 这里存储的是三对角块矩阵的LU分解的结果，这里我们分解的形式是单位上三角矩阵和下三角矩阵
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:)::L_LUall,D_LUall,U_LUall
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    ::L_LUacUm,L_LUacUs
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      ::L_LUacH,L_LUacD
  ! IPIV是存储对应的lapack中的LU分解的整数指标矩阵
  INTEGER( kind = ik ),ALLOCATABLE,DIMENSION(:,:,:,:)::IPIV

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  ::ADU,CDU,DU
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  :: DPLUR_RHS
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: DH,DV
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: DH_xi,DH_zeta,DV_xi,DV_zeta
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: CO_AC_UT
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: AsUs_P1,AsUs_P3

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDH,BDUSDH,CDUSDH
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHxi,BDUSDHxi,CDUSDHxi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHzeta,BDUSDHzeta,CDUSDHzeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDH_xi,CDUSDH_xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDH_zeta,CDUSDH_zeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHxi_xi,CDUSDHxi_xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHxi_zeta,CDUSDHxi_zeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHzeta_xi,CDUSDHzeta_xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDHzeta_zeta,CDUSDHzeta_zeta
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDV,BDUSDV,CDUSDV
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDV_xi,CDUSDV_xi
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: ADUSDV_zeta,CDUSDV_zeta

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:)        :: RHS_Krylov
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:)        :: DeltaSolu
  
  ! Data used for gmres iteration
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: BSS_V,BSS_A
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:,:):: BSS_U
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: HessenbergMat,Givens
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:)        :: Ve

  ! Data used for LNS iterations
  ! 我们这里使用线性多步法来进行计算，这里我们需要对应的存储对应的多步法的数据
  ! 原则上我们需要保持对应的多步法求解精度和对应的RK4一至，因此我们选用了对应
  ! 的4阶精度的线性多步法进行计算。这里我们需要对应的存储对应的多步法的数据
  ! Qn0 = L1*Q^n + L2*Q^(n-1) + L3*Q^(n-2) + L4*Q^(n-3)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Rho0,U0,V0,W0,T0,P0 !steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Rho_pert,U_pert,V_pert,W_pert,T_pert,P_pert !perturbation state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Rho_inte,U_inte,V_inte,W_inte,T_inte,P_inte !perturbation state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  :: LNS_CVp0   ! Present conservative variables
  ! The following is the previous time perturbations
  ! Here p1,p2,p3,p4 means n-1,n-2,n-3,n-4 time
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  :: LNS_CVp1,LNS_CVp2,LNS_CVp3,LNS_CVp4

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: LNS_ShockHp0,LNS_ShockVp0
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: LNS_ShockHp1,LNS_ShockHp2,LNS_ShockHp3,LNS_ShockHp4
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: LNS_ShockVp1,LNS_ShockVp2,LNS_ShockVp3,LNS_ShockVp4

  ! Here we define the data used for adam-bashforth method
  ! The adam-bashforth method can be written as
  ! Qn+1 = Qn + dt * R(sum a_k Qn_k)
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  :: AB_Sum_CVp0,Rd
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: AB_Sum_SHp0
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: AB_Sum_SVp0
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockAcd
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockNormalxd,ShockNormalyd,ShockNormalzd,shockXtaud
  REAL( kind = rk )                                 :: LNS_dt

  ! parameters used for Perturbation
  INTEGER( kind = ik )                              :: Pert_Type                   ! Select mode of perturbation
  INTEGER( kind = rk )                              :: LNS_STEPS_MAX
  INTEGER( kind = ik )                              :: IF_Continue_LNS
  REAL( kind = rk )                                 :: k_infty
  REAL( kind = rk )                                 :: epsilon  
  
  ! parameters used for unsteady solver
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: Rho_free              ! Previous Freestream Density
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: U_free                ! Previous Freestream Velocity along x direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: V_free                ! Previous Freestream Velocity along y direction   
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: W_free                ! Previous Freestream Velocity along z direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: P_free                ! Previous Freestream Pressure
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: T_free                ! Previous Freestream Temperature
  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: Rho_pert_tau          ! Perturbation Density 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: U_pert_tau            ! Perturbation Velocity along x direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: V_pert_tau            ! Perturbation Velocity along y direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: W_pert_tau            ! Perturbation Velocity along z direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: P_pert_tau            ! Perturbation Pressure
   
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: Rho_pert_taud          ! Perturbation Density 
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: U_pert_taud            ! Perturbation Velocity along x direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: V_pert_taud            ! Perturbation Velocity along y direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: W_pert_taud            ! Perturbation Velocity along z direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)  :: P_pert_taud            ! Perturbation Pressure

  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Rho_interp            ! Interpolated Density
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: U_interp              ! Interpolated Velocity along x direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: V_interp              ! Interpolated Velocity along y direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: W_interp              ! Interpolated Velocity along z direction
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: P_interp              ! Interpolated Pressure
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: T_interp              ! Interpolated Temperature  
  
  REAl( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:,:)  :: Ucons0_steady       ! Steady state conservative variables
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: x_grid_steady         ! Steady state X-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: y_grid_steady         ! Steady state Y-coordinate  
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: z_grid_steady         ! Steady state Z-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockH_steady         ! Steady state Shock Height
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockV_steady         ! Steady state Shock Velocity
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockNormalX_steady   ! Steady state Shock Normal Vector X-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockNormalY_steady   ! Steady state Shock Normal Vector Y-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: ShockNormalZ_steady   ! Steady state Shock Normal Vector Z-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:)      :: shockXtau_steady      ! Steady state Shock X-tau
  !REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: i_t,i_t_tau,i_t_tau_tau
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: x_grid_interp         ! Interpolated X-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: y_grid_interp         ! Interpolated Y-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: z_grid_interp         ! Interpolated Z-coordinate
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Hdst_interp           ! Interpolated Height of mesh grid
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: Rho0_interp           ! Interpolated Density from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: U0_interp             ! Interpolated Velocity along x direction from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: V0_interp             ! Interpolated Velocity along y direction from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: W0_interp             ! Interpolated Velocity along z direction from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: P0_interp             ! Interpolated Pressure from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:,:,:)    :: T0_interp             ! Interpolated Temperature from steady state
  REAL( kind = rk ),ALLOCATABLE,DIMENSION(:)        :: etaSC               
    contains

SUBROUTINE Initilize_SFSolver
  use SF_Constant,only:overLap,MaxKrylovSubSpace
  implicit none
  integer( kind = ik )::nz_min,nz_max
  integer( kind = ik )::nx_min,nx_max
  integer( kind = ik )::nz_total_max
  
  !write(*,*)"Initilize_SFSolver"

  nx_min = 1 - overLap
  nx_max = nx_local + overLap

  nz_min = 1 - overLap
  nz_max = nz_local + overLap
  nz_total_max = Nz + overLap

  Nsize = Nx * (Ny - 1) * Nz * NumVar             ! 除开激波边界外的变量个数
  Nsize = Nsize + nx_local * nz_local             ! 激波高度的变量个数
  Nsize = Nsize + nx_local * nz_local             ! 激波速度的变量个数
  
  ALLOCATE(Sxi_Surface(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(Sxi_Surface_total(nx_min:nx_max,nz_min:nz_total_max))
  Sxi_Surface = 0.d0;
  Sxi_Surface_total = 0.d0;
  
  ALLOCATE(WallSX(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSY(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSZ(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSXdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSYdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSZdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSXdzeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSYdzeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallSZdzeta(nx_min:nx_max,nz_min:nz_max))

  ALLOCATE(WALLSX_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSY_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSZ_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSXdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSYdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSZdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSXdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSYdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WALLSZdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSX_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSY_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSZ_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSXdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSYdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSZdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSXdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSYdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  !ALLOCATE(WALLSZdzeta_total(nx_min:nx_max,nz_min:nz_total_max))  
  

  WALLSX_total = 0.d0;
  WALLSY_total = 0.d0;
  WALLSZ_total = 0.d0;
  WALLSXdxi_total = 0.d0;
  WALLSYdxi_total = 0.d0;
  WALLSZdxi_total = 0.d0;
  WALLSXdzeta_total = 0.d0;
  WALLSYdzeta_total = 0.d0;
  WALLSZdzeta_total = 0.d0;
  
  
  WallSX = 0.d0;
  WallSY = 0.d0;
  WallSZ = 0.d0;
  WallSXdxi = 0.d0;
  WallSYdxi = 0.d0;
  WallSZdxi = 0.d0;
  WallSXdzeta = 0.d0;
  WallSYdzeta = 0.d0;
  WallSZdzeta = 0.d0;

  
  ALLOCATE(WallNormalX(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalY(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalZ(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalXdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalYdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalZdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalXdzeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalYdzeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(WallNormalZdzeta(nx_min:nx_max,nz_min:nz_max))

  WallNormalX = 0.d0;
  WallNormalY = 0.d0;
  WallNormalZ = 0.d0;
  WallNormalXdxi = 0.d0;
  WallNormalYdxi = 0.d0;
  WallNormalZdxi = 0.d0;
  WallNormalXdzeta = 0.d0;
  WallNormalYdzeta = 0.d0;
  WallNormalZdzeta = 0.d0;

  ALLOCATE(WallNormalX_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalY_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalZ_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalXdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalYdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalZdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalXdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalYdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(WallNormalZdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  
  WallNormalX_total = 0.d0;
  WallNormalY_total = 0.d0;
  WallNormalZ_total = 0.d0;
  WallNormalXdxi_total = 0.d0;
  WallNormalYdxi_total = 0.d0;
  WallNormalZdxi_total = 0.d0;
  WallNormalXdzeta_total = 0.d0;  
  WallNormalYdzeta_total = 0.d0;
  WallNormalZdzeta_total = 0.d0;
  
  
  ALLOCATE(Heta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(HetaDxi(nx_min:nx_max,1:Ny,nz_min:nz_max))  
  ALLOCATE(HetaDeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(HetaDzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))

  ALLOCATE(HetaDxi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(HetaDeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(HetaDzeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  Heta = 0.d0;
  HetaDxi = 0.d0;
  HetaDeta = 0.d0;
  HetaDzeta = 0.d0;
  HetaDxi_total = 0.d0;
  HetaDeta_total = 0.d0;
  HetaDzeta_total = 0.d0;

  ALLOCATE(X_grid(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Y_grid(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Z_grid(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hdst(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hdst_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))

  X_grid = 0.d0;
  Y_grid = 0.d0;
  Z_grid = 0.d0;
  Hdst = 0.d0;
  Hdst_steady = 0.d0;
  
  ALLOCATE(X_grid_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Y_grid_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Z_grid_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  
  X_grid_total = 0.d0;
  Y_grid_total = 0.d0;
  Z_grid_total = 0.d0;
  
  ALLOCATE(dxdxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxdzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxdtau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dydxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dydeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dydzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dydtau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzdxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzdzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzdtau(nx_min:nx_max,1:Ny,nz_min:nz_max))

  dxdxi = 0.d0;
  dxdeta = 0.d0;
  dxdzeta = 0.d0;
  dxdtau = 0.d0;
  dydxi = 0.d0;
  dydeta = 0.d0;
  dydzeta = 0.d0;
  dydtau = 0.d0;
  dzdxi = 0.d0;
  dzdeta = 0.d0;
  dzdzeta = 0.d0;
  dzdtau = 0.d0;

  ALLOCATE(dxidx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxidy(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxidz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dxidt(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(detadx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(detady(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(detadz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(detadt(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzetadx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzetady(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzetadz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dzetadt(nx_min:nx_max,1:Ny,nz_min:nz_max))

  dxidx = 0.d0;
  dxidy = 0.d0;
  dxidz = 0.d0;
  dxidt = 0.d0;
  detadx = 0.d0;
  detady = 0.d0;
  detadz = 0.d0;
  detadt = 0.d0;
  dzetadx = 0.d0;
  dzetady = 0.d0;
  dzetadz = 0.d0;
  dzetadt = 0.d0;

  ALLOCATE(nablaxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(nablaeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(nablazeta(nx_min:nx_max,1:Ny,nz_min:nz_max))

  nablaxi = 0.d0;
  nablaeta = 0.d0;
  nablazeta = 0.d0;

  ALLOCATE(jaco(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invJacodt(nx_min:nx_max,1:Ny,nz_min:nz_max))

  jaco = 0.d0;
  invJacodt = 0.d0;
  
  ALLOCATE(detadt_invJ(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(detadt_invJdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  detadt_invJ = 0.d0;
  detadt_invJdeta = 0.d0; 

  ALLOCATE(ax_tau(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ay_tau(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(az_tau(nx_min:nx_max,nz_min:nz_max))

  ax_tau = 0.d0;
  ay_tau = 0.d0;
  az_tau = 0.d0;

  ALLOCATE(Ucons0(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(UconsOld(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(UconsNew(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUcons0(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUcons(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUcons0_old(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Ucons0_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(dUcons0_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  Ucons0 = 0.d0;
  UconsOld = 0.d0;
  UconsNew = 0.d0;
  dUcons0 = 0.d0;
  dUcons = 0.d0;
  dUcons0_old = 0.d0;
  Ucons0_total = 0.d0;
  dUcons0_total = 0.d0;

  ALLOCATE(Rho(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(U(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(V(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(W(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(T(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(P(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Mu(nx_min:nx_max,1:Ny,nz_min:nz_max))

  Rho = 0.d0;
  U = 0.d0;
  V = 0.d0;
  W = 0.d0;
  T = 0.d0;
  P = 0.d0;
  Mu = 0.d0;

  ALLOCATE(Rho_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(U_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(V_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(W_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(T_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(P_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  Rho_total = 0.d0;
  U_total = 0.d0;
  V_total = 0.d0;
  W_total = 0.d0;
  T_total = 0.d0;
  P_total = 0.d0;
  
  ALLOCATE(Rhodxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rhodeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rhodzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))

  Rhodxi = 0.d0;
  Udxi = 0.d0;
  Vdxi = 0.d0;
  Wdxi = 0.d0;
  Tdxi = 0.d0;
  Rhodeta = 0.d0;
  Udeta = 0.d0;
  Vdeta = 0.d0;
  Wdeta = 0.d0;
  Tdeta = 0.d0;
  Rhodzeta = 0.d0;
  Udzeta = 0.d0;
  Vdzeta = 0.d0;
  Wdzeta = 0.d0;
  Tdzeta = 0.d0;

  ALLOCAte(Udxi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Udeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Udzeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Vdxi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Vdeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Vdzeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Wdxi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Wdeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Wdzeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Tdxi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Tdeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(Tdzeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  Udxi_total = 0.d0;
  Udeta_total = 0.d0;
  Udzeta_total = 0.d0;
  Vdxi_total = 0.d0;
  Vdeta_total = 0.d0;
  Vdzeta_total = 0.d0;
  Wdxi_total = 0.d0;
  Wdeta_total = 0.d0;
  Wdzeta_total = 0.d0;
  Tdxi_total = 0.d0;
  Tdeta_total = 0.d0;
  Tdzeta_total = 0.d0;
  
  ALLOCATE(BigU_Xi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(BigU_Eta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(BigU_Zeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  BigU_Xi = 0.d0;
  BigU_Eta = 0.d0;
  BigU_Zeta = 0.d0;
  
  ALLOCATE(BigU_Xi_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(BigU_Eta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(BigU_Zeta_total(nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  
  BigU_Xi_total = 0.d0;
  BigU_Eta_total = 0.d0;
  BigU_Zeta_total = 0.d0;
  
  ALLOCATE(Cs_Xi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Cs_Eta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Cs_Zeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  Cs_Xi = 0.d0;
  Cs_eta = 0.d0;
  Cs_Zeta = 0.d0;

  ALLOCATE(Rds_Xi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rds_Eta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rds_Zeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rds_Max(nx_min:nx_max,1:Ny,nz_min:nz_max))

  Rds_Xi = 0.d0;
  Rds_Eta = 0.d0;
  Rds_Zeta = 0.d0;
  Rds_Max = 0.d0;

  ALLOCATE(LocalDT(nx_min:nx_max,1:Ny,nz_min:nz_max))

  LocalDT = 0.d0;

  ALLOCATE(Cs(nx_min:nx_max,1:Ny,nz_min:nz_max))

  Cs = 0.d0;

  ALLOCATE(Rhodx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdx(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rhody(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udy(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdy(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdy(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdy(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Rhodz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Udz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Vdz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Wdz(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Tdz(nx_min:nx_max,1:Ny,nz_min:nz_max))

  Rhodx = 0.d0;
  Udx = 0.d0;
  Vdx = 0.d0;
  Wdx = 0.d0;
  Tdx = 0.d0;
  Rhody = 0.d0;
  Udy = 0.d0;
  Vdy = 0.d0;
  Wdy = 0.d0;
  Tdy = 0.d0;
  Rhodz = 0.d0;
  Udz = 0.d0;
  Vdz = 0.d0;
  Wdz = 0.d0;
  Tdz = 0.d0;

  ALLOCATE(invF(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invG(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invH(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visF(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visG(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visH(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  invF = 0.d0;
  invG = 0.d0;
  invH = 0.d0;
  visF = 0.d0;
  visG = 0.d0;
  visH = 0.d0;

  ALLOCATE(invFhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invGhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invHhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visFhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visGhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visHhat(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  invFhat = 0.d0;
  invGhat = 0.d0;
  invHhat = 0.d0;
  visFhat = 0.d0;
  visGhat = 0.d0;
  visHhat = 0.d0;
  
  ALLOCATE(visFhat_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(visGhat_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(visHhat_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  visFhat_total = 0.d0;
  visGhat_total = 0.d0;
  visHhat_total = 0.d0;

  ALLOCATE(Fp(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Fm(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Gp(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Gm(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hp(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hm(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  Fp = 0.d0;
  Fm = 0.d0;
  Gp = 0.d0;
  Gm = 0.d0;
  Hp = 0.d0;
  Hm = 0.d0;

  ALLOCATE(Fpdxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Fmdxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Gpdeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Gmdeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hpdzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Hmdzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  Fpdxi = 0.d0;
  Fmdxi = 0.d0;
  Gpdeta= 0.d0;
  Gmdeta= 0.d0;
  Hpdzeta = 0.d0;
  Hmdzeta = 0.d0;
  
  ALLOCATE(invFhatdxi_new(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(invGhatdeta_new(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  invFhatdxi_new = 0.d0;
  invGhatdeta_new = 0.d0;
  
  ALLOCATE(visFcdxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visGcdeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visHcdzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  visFcdxi = 0.d0;
  visGcdeta= 0.d0;
  visHcdzeta = 0.d0;
  
  ALLOCATE(visFcdxi_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(visGcdeta_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(visHcdzeta_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  visFcdxi_total = 0.d0;
  visGcdeta_total = 0.d0;
  visHcdzeta_total = 0.d0;

  ALLOCATE(invFlux(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(visFlux(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  invFlux = 0.d0;
  visFlux = 0.d0;

  ALLOCATE(visFlux_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))
  ALLOCATE(invFlux_total(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_total_max))

  visFlux_total = 0.d0;
  invFlux_total = 0.d0;
  
  ALLOCATE(Cf_Surface(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(Stan_Surface(nx_min:nx_max,nz_min:nz_max))

  Cf_Surface = 0.d0;
  Stan_Surface = 0.d0;

  ALLOCATE(ShockVN(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockH(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAc(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockH0(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV0(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(dShockH(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(dShockV(nx_min:nx_max,nz_min:nz_max))

  ShockVN = 0.d0;
  ShockH  = 0.d0;
  ShockV  = 0.d0;
  ShockH0 = 0.d0;
  ShockV0 = 0.d0;
  dShockH  = 0.d0;
  dShockV  = 0.d0;  
  ShockAc = 0.d0;

  ALLOCATE(ShockH_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(ShockV_total(nx_min:nx_max,nz_min:nz_total_max))
  
  ShockH_total = 0.d0;
  ShockV_total = 0.d0;
  
  ALLOCATE(ShockH_old(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV_old(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAc_old(nx_min:nx_max,nz_min:nz_max))
  ShockH_old = 0.d0;
  ShockV_old = 0.d0;
  ShockAc_old = 0.d0;

  ALLOCATE(FCV_SHK(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(GCV_SHK(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(HCV_SHK(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  FCV_SHK = 0.d0;
  GCV_SHK = 0.d0;
  HCV_SHK = 0.d0;

  ALLOCATE(ShockHdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockVdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockHdzeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockVdzeta(nx_min:nx_max,nz_min:nz_max))

  ShockHdxi = 0.d0;
  ShockVdxi = 0.d0;
  ShockHdzeta = 0.d0;
  ShockVdzeta = 0.d0;

  ALLOCATE(ShockHdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(ShockVdxi_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(ShockHdzeta_total(nx_min:nx_max,nz_min:nz_total_max))
  ALLOCATE(ShockVdzeta_total(nx_min:nx_max,nz_min:nz_total_max))

  ShockHdxi_total = 0.d0;
  ShockVdxi_total = 0.d0;
  ShockHdzeta_total = 0.d0;
  ShockVdzeta_total = 0.d0;
  
  ALLOCATE(ShockNormalX(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockNormalY(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockNormalZ(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockXtau(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockYtau(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockZtau(nx_min:nx_max,nz_min:nz_max))

  ShockNormalX = 0.d0;
  ShockNormalY = 0.d0;
  ShockNormalZ = 0.d0;
  ShockXtau = 0.d0;
  ShockYtau = 0.d0;
  ShockZtau = 0.d0;

  ALLOCATE(LeftEigenVectorsEta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  LeftEigenVectorsEta = 0.d0;

  ALLOCATE(DIAG(nx_min:nx_max,1:Ny,nz_min:nz_max))

  DIAG = 0.d0;

  ALLOCATE(DFhatDU(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DGhatDU(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DHhatDU(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DFhatDU = 0.d0;
  DGhatDU = 0.d0;
  DHhatDU = 0.d0;

  ALLOCATE(DvisFhatDp(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisFhatDpxi(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisFhatDpet(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisFhatDpzt(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DvisFhatDp = 0.d0;
  DvisFhatDpxi = 0.d0;
  DvisFhatDpet = 0.d0;
  DvisFhatDpzt = 0.d0;

  ALLOCATE(DvisGhatDp(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisGhatDpxi(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisGhatDpet(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisGhatDpzt(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DvisGhatDp = 0.d0;
  DvisGhatDpxi = 0.d0;
  DvisGhatDpet = 0.d0;
  DvisGhatDpzt = 0.d0;

  ALLOCATE(DvisHhatDp(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisHhatDpxi(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisHhatDpet(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DvisHhatDpzt(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DvisHhatDp = 0.d0;
  DvisHhatDpxi = 0.d0;
  DvisHhatDpet = 0.d0;
  DvisHhatDpzt = 0.d0;
  
  ALLOCATE(dUconsRK1(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUconsRK2(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUconsRK3(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(dUconsRK4(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(ShockAc_RK1(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAc_RK2(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAc_RK3(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAc_RK4(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV_RK1(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV_RK2(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV_RK3(nx_min:nx_max,nz_min:nz_max))

  dUconsRK1 = 0.d0;
  dUconsRK2 = 0.d0;
  dUconsRK3 = 0.d0;
  dUconsRK4 = 0.d0;
  ShockAc_RK1 = 0.d0;
  ShockAc_RK2 = 0.d0;
  ShockAc_RK3 = 0.d0;
  ShockAc_RK4 = 0.d0;
  ShockV_RK1 = 0.d0;
  ShockV_RK2 = 0.d0;
  ShockV_RK3 = 0.d0;

  ALLOCATE(Ucons_RK1(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Ucons_RK2(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Ucons_RK3(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  Ucons_RK1 = 0.d0;
  Ucons_RK2 = 0.d0;
  Ucons_RK3 = 0.d0;

  ! Implicit variables
  ALLOCATE(DRDH(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DRDV(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DRDHxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DRDHzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DRDH = 0.d0;
  DRDV = 0.d0;
  DRDHxi = 0.d0;
  DRDHzeta = 0.d0;

  ALLOCATE(ShockAcdH(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAcdHdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAcdHdzeta(nx_min:nx_max,nz_min:nz_max))

  ShockAcdH = 0.d0;
  ShockAcdHdxi = 0.d0;
  ShockAcdHdzeta = 0.d0;

  ALLOCATE(ShockAcdV(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAcdVdxi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockAcdVdzeta(nx_min:nx_max,nz_min:nz_max))
  
  ShockAcdV = 0.d0;
  ShockAcdVdxi = 0.d0;
  ShockAcdVdzeta = 0.d0;

  ALLOCATE(DUsDH(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DUsDHdxi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DUsDHdzeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DUsDV(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  
  !
  ALLOCATE(DxixJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DxiyJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DxizJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetaxJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetayJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetazJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DzetaxJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DzetayJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DzetazJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))

  DxixJDH = 0.d0;
  DxiyJDH = 0.d0;
  DxizJDH = 0.d0;

  DetaxJDH = 0.d0;
  DetayJDH = 0.d0;
  DetazJDH = 0.d0;

  DzetaxJDH = 0.d0;
  DzetayJDH = 0.d0;
  DzetazJDH = 0.d0;

  ALLOCATE(DJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetatDH(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetatJDH(nx_min:nx_max,1:Ny,nz_min:nz_max))

  DJDH = 0.d0;
  DJDHxi = 0.d0;
  DJDHzeta = 0.d0;
  DetatDH = 0.d0;  
  DetatJDH = 0.d0;
  
  ALLOCATE(DetatJDV(nx_min:nx_max,1:Ny,nz_min:nz_max))
  DetatJDV = 0.d0;
  
  ALLOCATE(DetaxJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetayJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetazJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetatJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))

  ALLOCATE(DetaxJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetayJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetazJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DetatJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))

  ALLOCATE(DzetaxJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DzetayJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DzetazJDHxi(nx_min:nx_max,1:Ny,nz_min:nz_max))

  ALLOCATE(DxixJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DxiyJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DxizJDHzeta(nx_min:nx_max,1:Ny,nz_min:nz_max))
  
     DetaxJDHxi = 0.0_8;
     DetayJDHxi = 0.0_8;
     DetazJDHxi = 0.0_8;
     DetatJDHxi = 0.0_8;

   DetaxJDHzeta = 0.0_8;
   DetayJDHzeta = 0.0_8;
   DetazJDHzeta = 0.0_8;
   DetatJDHzeta = 0.0_8;

    DzetaxJDHxi = 0.0_8;
    DzetayJDHxi = 0.0_8;
    DzetazJDHxi = 0.0_8;

    DxixJDHzeta = 0.0_8;
    DxiyJDHzeta = 0.0_8;
    DxizJDHzeta = 0.0_8;
    
  ALLOCATE(F_IM(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(G_IM(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(H_IM(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  F_IM = 0.d0;
  G_IM = 0.d0;
  H_IM = 0.d0;

  ALLOCATE(DF_IMDxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DF_IMDeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DF_IMDzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DF_IMDxi = 0.d0;
  DF_IMDeta = 0.d0;
  DF_IMDzeta = 0.d0;

  ALLOCATE(DG_IMDxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DG_IMDeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DG_IMDzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DG_IMDxi = 0.d0;
  DG_IMDeta = 0.d0;
  DG_IMDzeta = 0.d0;

  ALLOCATE(DH_IMDxi(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DH_IMDeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DH_IMDzeta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DH_IMDxi = 0.d0;
  DH_IMDeta = 0.d0;
  DH_IMDzeta = 0.d0;

  ALLOCATE(DU_Deta(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DU_Deta = 0.d0;
 
  !这里存储的是三对角块矩阵的LU分解的结果，这里我们分解的形式是单位上三角矩阵和下三角矩阵
  ALLOCATE(L_LUall(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(D_LUall(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(U_LUall(1:NumVar,1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  L_LUall = 0.d0;
  D_LUall = 0.d0;
  U_LUall = 0.d0;

  ALLOCATE(L_LUacUm(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(L_LUacUs(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  L_LUacUm = 0.d0;
  L_LUacUs = 0.d0;

  ALLOCATE(L_LUacH(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(L_LUacD(nx_min:nx_max,nz_min:nz_max))

  L_LUacH = 0.d0;
  L_LUacD = 0.d0;
  
  ALLOCATE(IPIV(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ! IPIV是存储对应的lapack中的LU分解的整数指标矩阵
  IPIV = 0;

  ALLOCATE(ADU(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(CDU(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(DU(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  ADU = 0.d0;
  CDU = 0.d0;
  DU = 0.d0;
  
  ALLOCATE(DPLUR_RHS(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

  DPLUR_RHS = 0.d0;

  ALLOCATE(DH(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DV(nx_min:nx_max,nz_min:nz_max))

  DH = 0.d0;
  DV = 0.d0;
  
  ALLOCATE(DH_xi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DH_zeta(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DV_xi(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(DV_zeta(nx_min:nx_max,nz_min:nz_max))

  DH_xi = 0.d0;
  DH_zeta = 0.d0;
  DV_xi = 0.d0;
  DV_zeta = 0.d0;

  ALLOCATE(CO_AC_UT(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  CO_AC_UT = 0.d0;

  ALLOCATE(AsUs_P1(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(AsUs_P3(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  AsUs_P1 = 0.d0;
  AsUs_P3 = 0.d0;

  ALLOCATE(ADUSDH(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(BDUSDH(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDH(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDH = 0.d0;
  BDUSDH = 0.d0;
  CDUSDH = 0.d0;

  ALLOCATE(ADUSDHxi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(BDUSDHxi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHxi(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDHxi = 0.d0;
  BDUSDHxi = 0.d0;
  CDUSDHxi = 0.d0;

  ALLOCATE(ADUSDHzeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(BDUSDHzeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHzeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDHzeta = 0.d0;
  BDUSDHzeta = 0.d0;
  CDUSDHzeta = 0.d0;

  ALLOCATE(ADUSDH_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDH_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ADUSDH_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDH_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDH_xi = 0.d0;
  CDUSDH_xi = 0.d0;
  ADUSDH_zeta = 0.d0;
  CDUSDH_zeta = 0.d0;

  ALLOCATE(ADUSDHxi_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHxi_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ADUSDHxi_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHxi_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDHxi_xi = 0.d0;
  CDUSDHxi_xi = 0.d0;
  ADUSDHxi_zeta = 0.d0;
  CDUSDHxi_zeta = 0.d0;

  ALLOCATE(ADUSDHzeta_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHzeta_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ADUSDHzeta_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDHzeta_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDHzeta_xi = 0.d0;
  CDUSDHzeta_xi = 0.d0;
  ADUSDHzeta_zeta = 0.d0;
  CDUSDHzeta_zeta = 0.d0;

  ALLOCATE(ADUSDV(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ADUSDV_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ADUSDV_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  ADUSDV = 0.d0;
  ADUSDV_xi = 0.d0;
  ADUSDV_zeta = 0.d0;

  ALLOCATE(BDUSDV(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  BDUSDV = 0.d0;

  ALLOCATE(CDUSDV(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDV_xi(1:NumVar,nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(CDUSDV_zeta(1:NumVar,nx_min:nx_max,nz_min:nz_max))

  CDUSDV = 0.d0;
  CDUSDV_xi = 0.d0;
  CDUSDV_zeta = 0.d0;
  !write(*,*)"Initilize_SFSolver Finished"

  ALLOCATE(RHS_Krylov(1:Nsize))
  ALLOCATE(DeltaSolu(1:Nsize))

  RHS_Krylov = 0.d0;
   DeltaSolu = 0.d0;

  ! Gmres
  ALLOCATE(BSS_A(nx_min:nx_max,nz_min:nz_max,MaxKrylovSubSpace+1))
  ALLOCATE(BSS_V(nx_min:nx_max,nz_min:nz_max,MaxKrylovSubSpace+1))
  ALLOCATE(BSS_U(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max,MaxKrylovSubSpace+1))

  BSS_A = 0.d0;
  BSS_V = 0.d0;
  BSS_U = 0.d0;

  ALLOCATE(HessenbergMat(MaxKrylovSubSpace+1,MaxKrylovSubSpace))
  ALLOCATE(Givens(2,MaxKrylovSubSpace))
  ALLOCATE(Ve(MaxKrylovSubSpace+1))

  HessenbergMat = 0.d0;
  Givens = 0.d0;
  Ve = 0.d0;

  !

    ALLOCATE(Rho0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(U0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(V0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(W0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(P0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(T0(nx_min:nx_max,1:Ny,nz_min:nz_max))
    
    Rho0 = 0.d0;
    U0 = 0.d0;
    V0 = 0.d0;
    W0 = 0.d0;
    P0 = 0.d0;
    T0 = 0.d0;
    
    !ALLOCATE(Ucons0_steady(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    !ALLOCATE(x_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    !ALLOCATE(y_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    !ALLOCATE(z_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    !ALLOCATE(shockH_steady(nx_min:nx_max,nz_min:nz_max))
    !ALLOCATE(shockV_steady(nx_min:nx_max,nz_min:nz_max))
    
    if(AnalysisType == LNS_Analysis) then  
     ALLOCATE(Rho_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(U_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(V_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(W_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(P_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(T_pert(nx_min:nx_max,1:Ny,nz_min:nz_max))
     
     Rho_pert = 0.d0;
     U_pert = 0.d0;
     V_pert = 0.d0;
     W_pert = 0.d0;
     P_pert = 0.d0;
     T_pert = 0.d0;
    
     ALLOCATE(Rho_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(U_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(V_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(W_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(P_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(T_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
     
     Rho_free = 0.d0;
     U_free = 0.d0;
     V_free = 0.d0;
     W_free = 0.d0;
     P_free = 0.d0;
     T_free = 0.d0;
     
     ALLOCATE(Rho_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(U_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(V_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))    
     ALLOCATE(W_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(P_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
     
     Rho_pert_tau = 0.d0;
     U_pert_tau = 0.d0;
     V_pert_tau = 0.d0;
     W_pert_tau = 0.d0;
     P_pert_tau = 0.d0;
    
     ALLOCATE(Rho_pert_taud(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(U_pert_taud(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(V_pert_taud(nx_min:nx_max,1:Ny,nz_min:nz_max))    
     ALLOCATE(W_pert_taud(nx_min:nx_max,1:Ny,nz_min:nz_max))
     ALLOCATE(P_pert_taud(nx_min:nx_max,1:Ny,nz_min:nz_max))
     
     Rho_pert_taud = 0.d0;
     U_pert_taud = 0.d0;
     V_pert_taud = 0.d0;
     W_pert_taud = 0.d0;
     P_pert_taud = 0.d0;     
     
    ALLOCATE(LNS_CVp0(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(LNS_CVp1(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(LNS_CVp2(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(LNS_CVp3(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(LNS_CVp4(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(AB_Sum_CVp0(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))

    LNS_CVp0 = 0.d0;
    LNS_CVp1 = 0.d0;
    LNS_CVp2 = 0.d0;
    LNS_CVp3 = 0.d0;
    LNS_CVp4 = 0.d0;
    AB_Sum_CVp0 = 0.d0;

    ALLOCATE(LNS_ShockHp0(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockHp1(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockHp2(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockHp3(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockHp4(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(AB_Sum_SHp0(nx_min:nx_max,nz_min:nz_max))

    LNS_ShockHp0 = 0.d0;
    LNS_ShockHp1 = 0.d0;
    LNS_ShockHp2 = 0.d0;
    LNS_ShockHp3 = 0.d0;
    LNS_ShockHp4 = 0.d0;
    AB_Sum_SHp0 = 0.d0;
    
    ALLOCATE(LNS_ShockVp0(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockVp1(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockVp2(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockVp3(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(LNS_ShockVp4(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(AB_Sum_SVp0(nx_min:nx_max,nz_min:nz_max))

    LNS_ShockVp0 = 0.d0;
    LNS_ShockVp1 = 0.d0;
    LNS_ShockVp2 = 0.d0;
    LNS_ShockVp3 = 0.d0;
    LNS_ShockVp4 = 0.d0;
    AB_Sum_SVp0  = 0.d0;
    
    ALLOCATE(Rd(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(ShockAcd(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(ShockNormalxd(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(ShockNormalyd(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(ShockNormalzd(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(shockXtaud(nx_min:nx_max,nz_min:nz_max))
    
    Rd = 0.d0;
    ShockAcd = 0.d0;
    ShockNormalxd = 0.d0;
    ShockNormalyd = 0.d0;
    ShockNormalzd = 0.d0;
    ShockXtaud    = 0.d0;

    ALLOCATE(Ucons0_steady(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(x_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(y_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(z_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
    ALLOCATE(shockH_steady(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(shockV_steady(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(shockNormalX_steady(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(shockNormalY_steady(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(shockNormalZ_steady(nx_min:nx_max,nz_min:nz_max))
    ALLOCATE(ShockXtau_steady(nx_min:nx_max,nz_min:nz_max))
  
    Ucons0_steady = 0.d0;
    x_grid_steady = 0.d0;
    y_grid_steady = 0.d0;
    z_grid_steady = 0.d0;
    shockH_steady = 0.d0;
    shockV_steady = 0.d0;
    shockNormalX_steady = 0.d0;
    shockNormalY_steady = 0.d0;
    shockNormalZ_steady = 0.d0;
    ShockXtau_steady    = 0.d0;
    
    ALLOCATE(Rho_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(U_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(V_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(W_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(P_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(T_inte(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  
  Rho_inte = 0.d0;
    U_inte = 0.d0;
    V_inte = 0.d0;
    W_inte = 0.d0;
    P_inte = 0.d0;
    T_inte = 0.d0;

    ALLOCATE(x_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(y_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(z_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(Hdst_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(etaSC(1:2*Ny-1))
    
    x_grid_interp = 0.d0;
    y_grid_interp = 0.d0;
    z_grid_interp = 0.d0;
    Hdst_interp = 0.d0;
    etaSC = 0.d0;    
    
  endif

  if(AnalysisType == Unsteady_NS_Analysis) then
  ALLOCATE(Rho_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(U_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(V_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(W_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(P_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(T_free(nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  Rho_free = 0.d0;
  U_free = 0.d0;
  V_free = 0.d0;
  W_free = 0.d0;
  P_free = 0.d0;
  T_free = 0.d0;
  
    ALLOCATE(Rho_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(U_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(V_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(W_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(P_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    ALLOCATE(T_pert(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
    
    Rho_pert = 0.d0;
    U_pert = 0.d0;
    V_pert = 0.d0;
    W_pert = 0.d0;
    P_pert = 0.d0;
    T_pert = 0.d0;
  
  ALLOCATE(Rho_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(U_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(V_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))    
  ALLOCATE(W_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(P_pert_tau(nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  Rho_pert_tau = 0.d0;
  U_pert_tau = 0.d0;
  V_pert_tau = 0.d0;
  W_pert_tau = 0.d0;
  P_pert_tau = 0.d0;

  ALLOCATE(Rho_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(U_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(V_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(W_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(P_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(T_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  
  Rho_interp = 0.d0;
    U_interp = 0.d0;
    V_interp = 0.d0;
    W_interp = 0.d0;
    P_interp = 0.d0;
    T_interp = 0.d0;
  
  ALLOCATE(ShockH_steady(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(ShockV_steady(nx_min:nx_max,nz_min:nz_max))
  ALLOCATE(x_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(y_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(z_grid_steady(nx_min:nx_max,1:Ny,nz_min:nz_max))
  ALLOCATE(Ucons0_steady(1:NumVar,nx_min:nx_max,1:Ny,nz_min:nz_max))
  
  ShockH_steady = 0.d0;
  ShockV_steady = 0.d0;
  x_grid_steady = 0.d0;
  y_grid_steady = 0.d0;
  z_grid_steady = 0.d0;
  Ucons0_steady = 0.d0;
  
  ALLOCATE(x_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(y_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(z_grid_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(Hdst_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(etaSC(1:2*Ny-1))
  
  x_grid_interp = 0.d0;
  y_grid_interp = 0.d0;
  z_grid_interp = 0.d0;
  Hdst_interp = 0.d0;
  etaSC = 0.d0;
  
  ALLOCATE(Rho0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(U0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(V0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(W0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(P0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  ALLOCATE(T0_interp(nx_min:nx_max,1:2*Ny-1,nz_min:nz_max))
  
  Rho0_interp = 0.d0;
    U0_interp = 0.d0;
    V0_interp = 0.d0;
    W0_interp = 0.d0;
    P0_interp = 0.d0;
    T0_interp = 0.d0;
  
  endif
  
END SUBROUTINE Initilize_SFSolver

SUBROUTINE Finilize_SFSolver
  implicit none
  
  deallocate(Sxi_Surface)
  deallocate(Sxi_Surface_total)
  deallocate(WallSX,WallSY,WallSZ)
  deallocate(WallSX_total,WallSY_total,WallSZ_total)
  deallocate(WallSXdxi,WallSYdxi,WallSZdxi)
  deallocate(WallSXdzeta,WallSYdzeta,WallSZdzeta)
  deallocate(WallSXdxi_total,WallSYdxi_total,WallSZdxi_total)
  deallocate(WallSXdzeta_total,WallSYdzeta_total,WallSZdzeta_total)
  deallocate(WallNormalX,WallNormalY,WallNormalZ)
  deallocate(WallNormalX_total,WallNormalY_total,WallNormalZ_total)
  deallocate(WallNormalXdxi,WallNormalYdxi,WallNormalZdxi)
  deallocate(WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta)
  deallocate(WallNormalXdxi_total,WallNormalYdxi_total,WallNormalZdxi_total)
  deallocate(WallNormalXdzeta_total,WallNormalYdzeta_total,WallNormalZdzeta_total)
  deallocate(Heta,HetaDxi,HetaDeta,HetaDzeta)
  deallocate(X_grid,Y_grid,Z_grid,Hdst,Hdst_steady)
  deallocate(X_grid_total,Y_grid_total,Z_grid_total)
  deallocate(dxdxi,dxdeta,dxdzeta,dxdtau)
  deallocate(dydxi,dydeta,dydzeta,dydtau)
  deallocate(dzdxi,dzdeta,dzdzeta,dzdtau)
  deallocate(dxidx,dxidy,dxidz,dxidt)
  deallocate(detadx,detady,detadz,detadt)
  deallocate(dzetadx,dzetady,dzetadz,dzetadt)
  deallocate(nablaxi,nablaeta,nablazeta)
  deallocate(jaco,invJacodt)
  deallocate(detadt_invJ,detadt_invJdeta)
  deallocate(ax_tau,ay_tau,az_tau)
  deallocate(Ucons0,UconsOld,UconsNew,Ucons0_total)
  deallocate(dUcons0,dUcons,dUcons0_old)
  deallocate(dUcons0_total)
  deallocate(Rho,U,V,W,T,P,Mu)
  deallocate(Rho_total,U_total,V_total,W_total,T_total,P_total)
  deallocate(BigU_Xi,BigU_Eta,BigU_Zeta)
  deallocate(BigU_Xi_total,BigU_Eta_total,BigU_Zeta_total)
  deallocate(Cs_Xi,Cs_Eta,Cs_Zeta)
  deallocate(Rds_Xi,Rds_Eta,Rds_Zeta,Rds_Max)
  deallocate(LocalDT)
  deallocate(Cs)
  deallocate(Rhodx,Udx,Vdx,Wdx,Tdx)
  deallocate(Rhody,Udy,Vdy,Wdy,Tdy)
  deallocate(Rhodz,Udz,Vdz,Wdz,Tdz)
  deallocate(Rhodxi,Udxi,Vdxi,Wdxi,Tdxi)
  deallocate(Rhodeta,Udeta,Vdeta,Wdeta,Tdeta)
  deallocate(Rhodzeta,Udzeta,Vdzeta,Wdzeta,Tdzeta)
  deallocate(Udxi_total,Udeta_total,Udzeta_total)
  deallocate(Vdxi_total,Vdeta_total,Vdzeta_total)
  deallocate(Wdxi_total,Wdeta_total,Wdzeta_total)
  deallocate(Tdxi_total,Tdeta_total,Tdzeta_total)

  deallocate(invF,invG,invH,visF,visG,visH)
  deallocate(invFhat,invGhat,invHhat,visFhat,visGhat,visHhat)
  deallocate(visFhat_total,visGhat_total,visHhat_total)
  deallocate(Fp,Fm,Gp,Gm,Hp,Hm)
  deallocate(Fpdxi,Fmdxi,Gpdeta,Gmdeta,Hpdzeta,Hmdzeta)
  deallocate(invFhatdxi_new,invGhatdeta_new)
  deallocate(visFcdxi,visGcdeta,visHcdzeta)
  deallocate(visFcdxi_total,visGcdeta_total,visHcdzeta_total)
  deallocate(invFlux,visFlux)
  deallocate(invFlux_total,visFlux_total)
  deallocate(Cf_Surface,Stan_Surface)
  deallocate(ShockVN,ShockH,ShockV,ShockAc,ShockH0,ShockV0,dShockH,dShockV)
  deallocate(ShockH_total,ShockV_total)
  deallocate(ShockH_old,ShockV_old,ShockAc_old)
  deallocate(FCV_SHK,GCV_SHK,HCV_SHK)
  deallocate(ShockHdxi,ShockVdxi,ShockHdzeta,ShockVdzeta)
  deallocate(ShockNormalX,ShockNormalY,ShockNormalZ,ShockXtau,ShockYtau,ShockZtau)
  deallocate(LeftEigenVectorsEta)
  deallocate(DIAG,DFhatDU,DGhatDU,DHhatDU)
  deallocate(DvisFhatDp,DvisFhatDpxi,DvisFhatDpet,DvisFhatDpzt)
  deallocate(DvisGhatDp,DvisGhatDpxi,DvisGhatDpet,DvisGhatDpzt)
  deallocate(DvisHhatDp,DvisHhatDpxi,DvisHhatDpet,DvisHhatDpzt)
  deallocate(dUconsRK1,dUconsRK2,dUconsRK3,dUconsRK4)
  deallocate(Ucons_RK1,Ucons_RK2,Ucons_RK3)
  deallocate(ShockAc_RK1,ShockAc_RK2,ShockAc_RK3,ShockAc_RK4)
  deallocate(ShockV_RK1,ShockV_RK2,ShockV_RK3)
  deallocate(DRDH,DRDV,DRDHxi,DRDHzeta)
  deallocate(ShockAcdH,ShockAcdHdxi,ShockAcdHdzeta)
  deallocate(ShockAcdV,ShockAcdVdxi,ShockAcdVdzeta)
  deallocate(DUsDH,DUsDHdxi,DUsDHdzeta,DUsDV)
  deallocate(DxixJDH,DxiyJDH,DxizJDH)
  deallocate(DetaxJDH,DetayJDH,DetazJDH)
  deallocate(DzetaxJDH,DzetayJDH,DzetazJDH)
  deallocate(DJDH,DetatDH,DetatJDH,DetatJDV)
  deallocate(DJDHxi,DJDHzeta)
  deallocate(DetaxJDHxi,  DetayJDHxi,  DetazJDHxi,  DetatJDHxi)
  deallocate(DetaxJDHzeta,DetayJDHzeta,DetazJDHzeta,DetatJDHzeta)
  deallocate(DzetaxJDHxi, DzetayJDHxi, DzetazJDHxi)
  deallocate(DxixJDHzeta, DxiyJDHzeta, DxizJDHzeta)
  deallocate(F_IM, G_IM, H_IM)
  deallocate(DF_IMDxi, DF_IMDeta, DF_IMDzeta)
  deallocate(DG_IMDxi, DG_IMDeta, DG_IMDzeta)
  deallocate(DH_IMDxi, DH_IMDeta, DH_IMDzeta)
  deallocate(DU_Deta)
  deallocate(L_LUall,D_LUall,U_LUall)
  deallocate(L_LUacUm,L_LUacUs)
  deallocate(L_LUacH,L_LUacD)
  deallocate(IPIV)
  deallocate(ADU,CDU,DU)
  deallocate( DPLUR_RHS)
  deallocate( DH,DV)
  deallocate( DH_xi,DH_zeta,DV_xi,DV_zeta)
  deallocate( CO_AC_UT)
  deallocate( AsUs_P1,AsUs_P3)
  deallocate( ADUSDH,BDUSDH,CDUSDH)
  deallocate( ADUSDHxi,BDUSDHxi,CDUSDHxi)
  deallocate( ADUSDHzeta,BDUSDHzeta,CDUSDHzeta)
  deallocate( ADUSDH_xi,CDUSDH_xi)
  deallocate( ADUSDH_zeta,CDUSDH_zeta)
  deallocate( ADUSDHxi_xi,CDUSDHxi_xi)
  deallocate( ADUSDHxi_zeta,CDUSDHxi_zeta)
  deallocate( ADUSDHzeta_xi,CDUSDHzeta_xi)
  deallocate( ADUSDHzeta_zeta,CDUSDHzeta_zeta)
  deallocate( ADUSDV,ADUSDV_xi,ADUSDV_zeta)
  deallocate( BDUSDV)
  deallocate( CDUSDV,CDUSDV_xi,CDUSDV_zeta)
  deallocate( RHS_Krylov, DeltaSolu) 

    
    DEALLOCATE(Rho0)
    DEALLOCATE(U0)
    DEALLOCATE(V0)
    DEALLOCATE(W0)
    DEALLOCATE(P0)
    DEALLOCATE(T0)
    
  if(AnalysisType == LNS_Analysis) then     
        
    DEALLOCATE(Rho_pert)
    DEALLOCATE(U_pert)
    DEALLOCATE(V_pert)
    DEALLOCATE(W_pert)
    DEALLOCATE(P_pert)
    DEALLOCATE(T_pert)
    DEALLOCATE(Rho_pert_tau)
    DEALLOCATE(U_pert_tau)
    DEALLOCATE(V_pert_tau)
    DEALLOCATE(W_pert_tau)
    DEALLOCATE(P_pert_tau)
    DEALLOCATE(Rho_pert_taud)
    DEALLOCATE(U_pert_taud)
    DEALLOCATE(V_pert_taud)
    DEALLOCATE(W_pert_taud)
    DEALLOCATE(P_pert_taud)
    DEALLOCATE(Rho_free)
    DEALLOCATE(U_free)
    DEALLOCATE(V_free)
    DEALLOCATE(W_free)
    DEALLOCATE(P_free)
    DEALLOCATE(T_free)
    
    DEALLOCATE(LNS_CVp0)
    DEALLOCATE(LNS_CVp1)
    DEALLOCATE(LNS_CVp2)
    DEALLOCATE(LNS_CVp3)
    
    Deallocate(LNS_CVp0)
    DEALLOCATE(LNS_CVp1)
    DEALLOCATE(LNS_CVp2)
    DEALLOCATE(LNS_CVp3)
    DEALLOCATE(LNS_CVp4)
    DEALLOCATE(AB_Sum_CVp0)

    DEALLOCATE(LNS_ShockHp0)
    DEALLOCATE(LNS_ShockHp1)
    DEALLOCATE(LNS_ShockHp2)
    DEALLOCATE(LNS_ShockHp3)
    DEALLOCATE(LNS_ShockHp4)
    DEALLOCATE(AB_Sum_SHp0)
    
    DEALLOCATE(LNS_ShockVp0)
    DEALLOCATE(LNS_ShockVp1)
    DEALLOCATE(LNS_ShockVp2)
    DEALLOCATE(LNS_ShockVp3)
    DEALLOCATE(LNS_ShockVp4)
    DEALLOCATE(AB_Sum_SVp0)

    DEALLOCATE(Rd)
    DEALLOCATE(ShockAcd)
    DEALLOCATE(Ucons0_steady)
    DEALLOCATE(x_grid_steady)
    DEALLOCATE(y_grid_steady)
    DEALLOCATE(z_grid_steady)
    DEALLOCATE(ShockH_steady)
    DEALLOCATE(ShockV_steady)
    DEALLOCATE(ShockNormalX_steady)
    DEALLOCATE(ShockNormalY_steady)
    DEALLOCATE(ShockNormalZ_steady)
    DEALLOCATE(ShockXtau_steady)
    DEALLOCATE(ShockNormalxd)
    DEALLOCATE(ShockNormalyd)
    DEALLOCATE(ShockNormalzd)
    DEALLOCATE(ShockXtaud)
    
    DEALLOCATE(Rho_inte)
    DEALLOCATE(U_inte)    
    DEALLOCATE(V_inte)
    DEALLOCATE(W_inte)
    DEALLOCATE(P_inte)
    DEALLOCATE(T_inte)    
    
    DEALLOCATE(x_grid_interp)
    DEALLOCATE(y_grid_interp)
    DEALLOCATE(z_grid_interp)
    
    DEALLOCATE(Hdst_interp)
    DEALLOCATE(etaSC)
  endif

  if (AnalysisType == Unsteady_NS_Analysis) then
      
    DEALLOCATE(Rho_free)
    DEALLOCATE(U_free)
    DEALLOCATE(V_free)
    DEALLOCATE(W_free)
    DEALLOCATE(P_free)
    DEALLOCATE(T_free)
    
    DEALLOCATE(Rho_pert_tau)
    DEALLOCATE(U_pert_tau)
    DEALLOCATE(V_pert_tau)
    DEALLOCATE(W_pert_tau)
    DEALLOCATE(P_pert_tau)
    
    DEALLOCATE(Rho_pert)
    DEALLOCATE(U_pert)
    DEALLOCATE(V_pert)
    DEALLOCATE(W_pert)
    DEALLOCATE(P_pert)
    DEALLOCATE(T_pert)
    
    DEALLOCATE(ShockH_steady)
    DEALLOCATE(ShockV_steady)
    DEALLOCATE(x_grid_steady)
    DEALLOCATE(y_grid_steady)
    DEALLOCATE(z_grid_steady)
    DEALLOCATE(Ucons0_steady)
    
    DEALLOCATE(Rho0_interp)
    DEALLOCATE(U0_interp)    
    DEALLOCATE(V0_interp)
    DEALLOCATE(W0_interp)
    DEALLOCATE(P0_interp)
    DEALLOCATE(T0_interp)
    
    DEALLOCATE(Rho_interp)
    DEALLOCATE(U_interp)    
    DEALLOCATE(V_interp)
    DEALLOCATE(W_interp)
    DEALLOCATE(P_interp)
    DEALLOCATE(T_interp)    
    
    DEALLOCATE(x_grid_interp)
    DEALLOCATE(y_grid_interp)
    DEALLOCATE(z_grid_interp)
    
    DEALLOCATE(Hdst_interp)
    DEALLOCATE(Hdst_steady)
    DEALLOCATE(etaSC)
    
  endif
  
END SUBROUTINE Finilize_SFSolver

END MODULE SF_CFD_Global