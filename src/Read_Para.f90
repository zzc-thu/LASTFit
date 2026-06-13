subroutine Read_Para
  use SF_Constant,  only: ik,rk,NumVar,Csutherland, pi
  use SF_CFD_Global,only: ModelType,Nx,Ny,Nz,AoA_Deg,AoS_Deg,Re_Ref,Te_Ref,Pr_Ref,TWall,Gamma,Mach_Ref,&
                          xim,sigma_y,yim,Lz_start,Lz_end,R_head,LModel,Nx_head,theta_p_deg,theta_p_rad,&
                          Mach_Normal,Nx_head_hf,Nx_flate,T0_INF,StepChangeCFL,&
                          TotalTimeSteps,StepsOutputScreen,StepsWriteData,StepChangeCFLStart,&
                          CFL_Start, CFL_End,Output_Residual_Mode,Select_Scheme,AnalysisType,&
                          IF_FixedDeltaTime,FixedDeltaTime

  use SF_CFD_Global,only: IF_Continue_Calculate,Select_Kind_Solver,Select_Kind_Precon_Solver,OutputFileNo
  use MPI_Global,   only: MyId,npx0,npz0,MPI_INTEGER,MPI_COMM_WORLD,IERR,MPI_DOUBLE_PRECISION
  ! Arrays and matrix
  use SF_CFD_Global,only: Csthlnd_Ref,U_inf,V_inf,W_inf,P_inf,Rho_inf,T_inf,&
                          AoA,AoS,CV_INF,FINV_INF,GINV_INF,HINV_INF,C_inf
  use SF_CFD_Global,only: k_infty,epsilon,LNS_dt,LNS_STEPS_MAX,IF_Continue_LNS,Pert_Type,Geometry_Type
  ! funs and subs
  implicit none
  integer( kind = ik )::  Int_parameters(1:50)                             
  real(    kind = rk ):: Real_parameters(1:50)

  namelist /InitialSettings/ ModelType,Nx,Ny,Nz,npx0,npz0,AoA_Deg,AoS_Deg,Re_Ref,Te_Ref,Pr_Ref,TWall,Gamma,Mach_Ref,&
                             IF_Continue_Calculate,Select_Kind_Solver,Select_Kind_Precon_Solver,&
                             xim,sigma_y,yim,Lz_start,Lz_end,R_head,LModel,Nx_head,theta_p_deg,CFL_Start,CFL_End,&
                             TotalTimeSteps,StepsOutputScreen,StepsWriteData,StepChangeCFLStart,&
                             StepChangeCFL,Output_Residual_Mode,Select_Scheme,AnalysisType,&
                             k_infty,epsilon,LNS_dt,LNS_STEPS_MAX,IF_Continue_LNS,Pert_Type,Geometry_Type,&
                             IF_FixedDeltaTime,FixedDeltaTime
  

  ! Read Parameters Using Namelist, only the main process MyId = 0 read the information
  ! and then transfer it to all processes

  if( MyId == 0 ) then
    write(*,*)" Reading the input parameters "
     open(unit = 23,file="Config.cfg")
     read(unit = 23,NML=InitialSettings)
     close(unit = 23)
    write(*,*)" Reading finished "
  endif
  
    Int_parameters(1)           = Nx
    Int_parameters(2)           = Nx_head
    Int_parameters(3)           = Ny
    Int_parameters(4)           = Nz
    Int_parameters(5)           = npx0
    Int_parameters(6)           = npz0
    Int_parameters(7)           = IF_Continue_Calculate
    Int_parameters(8)           = Select_Kind_Solver
    Int_parameters(9)           = Select_Kind_Precon_Solver
    Int_parameters(10)          = TotalTimeSteps
    Int_parameters(11)          = StepsOutputScreen
    Int_parameters(12)          = StepsWriteData
    Int_parameters(13)          = StepChangeCFLStart
    Int_parameters(14)          = StepChangeCFL
    Int_parameters(15)          = Output_Residual_Mode
    Int_parameters(16)          = Select_Scheme
    Int_parameters(17)          = ModelType
    Int_parameters(18)          = AnalysisType
    Int_parameters(19)          = IF_FixedDeltaTime
    Int_parameters(20)          = Pert_Type
    Int_parameters(21)          = LNS_STEPS_MAX
    Int_parameters(22)          = IF_Continue_LNS
    Int_parameters(23)          = Geometry_Type
    
  CALL MPI_Bcast( Int_parameters,50_ik,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    
    Nx                          = Int_parameters(1)
    Nx_head                     = Int_parameters(2)
    Ny                          = Int_parameters(3)
    Nz                          = Int_parameters(4)
    npx0                        = Int_parameters(5)
    npz0                        = Int_parameters(6)
    IF_Continue_Calculate       = Int_parameters(7)
    Select_Kind_Solver          = Int_parameters(8)
    Select_Kind_Precon_Solver   = Int_parameters(9)
    TotalTimeSteps              = Int_parameters(10)
    StepsOutputScreen           = Int_parameters(11)
    StepsWriteData              = Int_parameters(12)
    StepChangeCFLStart          = Int_parameters(13)
    StepChangeCFL               = Int_parameters(14)
    Output_Residual_Mode        = Int_parameters(15)
    Select_Scheme               = Int_parameters(16)
    ModelType                   = Int_parameters(17)
    AnalysisType                = Int_parameters(18)
    IF_FixedDeltaTime           = Int_parameters(19)
    Pert_Type                   = Int_parameters(20)
    LNS_STEPS_MAX               = Int_parameters(21)
    IF_Continue_LNS             = Int_parameters(22)
    Geometry_Type               = Int_parameters(23)
 
  Real_parameters(1) =AoA_Deg                      
  Real_parameters(2) =AoS_Deg                      
  Real_parameters(3) =Re_Ref                       
  Real_parameters(4) =Te_Ref                       
  Real_parameters(5) =Pr_Ref                       
  Real_parameters(6) =TWall                        
  Real_parameters(7) =Gamma                        
  Real_parameters(8) =Mach_Ref                     
  Real_parameters(9) =theta_p_deg                  
  Real_parameters(10)=xim                          
  Real_parameters(11)=sigma_y                      
  Real_parameters(12)=yim                          
  Real_parameters(13)=Lz_start                     
  Real_parameters(14)=Lz_end                       
  Real_parameters(15)=R_head                      
  Real_parameters(16)=LModel                     
  Real_parameters(17)=CFL_Start                  
  Real_parameters(18)=CFL_End 
  Real_parameters(19)=FixedDeltaTime
  Real_parameters(20)=k_infty
  Real_parameters(21)=epsilon
  Real_parameters(22)=LNS_dt
  

  CALL MPI_Bcast( Real_parameters,50_ik,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  AoA_Deg       = Real_parameters(1) 
  AoS_Deg       = Real_parameters(2) 
  Re_Ref        = Real_parameters(3) 
  Te_Ref        = Real_parameters(4) 
  Pr_Ref        = Real_parameters(5) 
  TWall         = Real_parameters(6) 
  Gamma         = Real_parameters(7) 
  Mach_Ref      = Real_parameters(8) 
  theta_p_deg   = Real_parameters(9) 
  xim           = Real_parameters(10)
  sigma_y       = Real_parameters(11)
  yim           = Real_parameters(12)
  Lz_start      = Real_parameters(13)
  Lz_end        = Real_parameters(14)
  R_head        = Real_parameters(15)
  LModel        = Real_parameters(16)
  CFL_Start     = Real_parameters(17)
  CFL_End       = Real_parameters(18)
  FixedDeltaTime = Real_parameters(19)
  k_infty       = Real_parameters(20)
  epsilon       = Real_parameters(21)
  LNS_dt        = Real_parameters(22)
                      
  ! Define the Geometrical Information

  Nx_head_hf = Nx_head /2_ik;
 
  Nx_flate   = (Nx + 1_ik) - Nx_head;

  ! Transformation of Degree to Rad
  AoA = AoA_Deg * pi / 180.d0;
  AoS = AoS_Deg * pi / 180.d0;
  theta_p_rad = theta_p_deg * pi / 180.d0;

  ! The free stream variabels are used for nodimensionalization
  ! The primitive variables
  Rho_inf = 1.d0;
    U_inf = 1.d0 * COS( AoA ) * COS( AoS );
    V_inf = 1.d0 * SIN( AoA )
    W_inf = 1.d0 * SIN( AoS )
    T_inf = 1.d0;
    P_inf = 1.d0 / (Gamma * Mach_Ref * Mach_Ref)
    if (Pert_Type == 1) then
    C_inf = U_inf + sqrt( Gamma * P_inf / Rho_inf );
    else if (Pert_Type == 2) then
    C_inf = U_inf - sqrt( Gamma * P_inf / Rho_inf );
    else 
    C_inf = U_inf;
    endif
   T0_INF = T_inf * (1.d0 + 0.5d0 * (Gamma - 1.d0) * Mach_Ref * Mach_Ref);

  Mach_Normal = Mach_Ref * U_inf;
  
  
  ! The conservative variables
  call     UtoCV_CPG(Rho_inf,U_inf,V_inf,W_inf,P_inf,CV_INF(1:NumVar))
  ! The flux
  ! Here we already initialize the flux at the infinity
  call UtoFlux_CPG3D(Rho_inf,U_inf,V_inf,W_inf,P_inf,FINV_INF(1:NumVar),GINV_INF(1:NumVar),HINV_INF(1:NumVar))
  
  ! Surtherland Coefficient  
  Csthlnd_Ref = Csutherland / Te_Ref;

  ! Initialize the Output File Number;
  OutputFileNo = 0_ik;

end subroutine Read_Para