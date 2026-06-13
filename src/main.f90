program SFSolver
  !=========================================================================!
  !                                                                         !
  ! If your program takes more than 500 steps to reach a steady state,      !
  ! you are probably doing something wrong.                                 !
  !                                                                         !
  !                 Numerical Computation of Compressible and Viscous Flows !
  !                                                - Robert W. MacCormack   !
  !=========================================================================!
  use SF_Constant,       only: ik,RK4_Methods, DPLUR_Methods, GMRES_Methods, MultiGridMethods,&
                              & Steady_NS_Analysis,Unsteady_NS_Analysis,LNS_Analysis,HLNS_Analysis
  use SF_CFD_Global,     only: Iteration, TotalTimeSteps, Select_Kind_Precon_Solver,&
                              &Select_Kind_Solver,IfConverge,ModelType,AnalysisType,&
                              &IF_FixedDeltatime,FixedDeltatime,dT0,Rds_Max,StepsWriteData,&
                              &IF_Continue_Calculate,Iteration_start
  use MPI_Global,        only: MyId,MPI_COMM_WORLD,ierr
  ! matrix and arrays
  ! subroutine and function declaration
  use SF_CFD_Global,     only: Initilize_SFSolver, Finilize_SFSolver
  use MPI_GLOBAL,        only: Parallel_Initilize, Parallel_Finalize, &
                               Parallel_SplitXZ, Parallel_Allocate,&
                               MPI_Barrier
  use Runge_Kutta_4,     only: IncreamentRK4,IncreamentRK4_unsteady
  use DPLUR,             only: DPLUR_Solver
  use Krylov_SubSpace,   only: GMRES_iteration,Gmres_SF_Solver
  use Monitor,           only: Output_Brief_Information
  use OutputParaView,    only: Output_Results
  use SF_LNS,            only: LNS_Solver
  implicit none
  ! Local variables
  ! Parallel initialization
  call Parallel_Initilize
  ! Read parameters and Bcast to all processes
  call Read_Para
  ! Initialize the grid and the grid partition
  call Parallel_SplitXZ
  ! Allocate memory for the variables
  call Initilize_SFSolver
  ! Allocate memory for the parallel transfer
  call Parallel_Allocate
  if(ModelType == 1) then ! 2.5D cases
    ! Initialize the grid
    call Initial_Grid
    ! Jacobian calculation
    call Calculate_Jaco
    ! Initial Fields
    call Initial_Fields
  endif
  if(ModelType == 2) then ! slender body
    !! Initialize the grid
     call Singular_Initial_Grid
    ! Jacobian calculation
      call Singular_Calculate_Jaco
    ! Initial Fields
      call Singular_Initial_Fields
  endif
   if(ModelType == 3) then ! Elliptical cone 
      call Elliptical_Initial_Grid
    ! Jacobian calculation
      call Singular_Calculate_Jaco
    ! Initial Fields
      call Singular_Initial_Fields
   endif
   if(ModelType == 4) then ! Hifire5b   
      call Hifire5b_Initial_Grid
    ! Jacobian calculation
      call Singular_Calculate_Jaco
    ! Initial Fields
      call Singular_Initial_Fields
   endif   
  
  ! Here we start the main Time loop
  if(MyId == 0) then
    write(*,*)"=========================================================="
    write(*,*)"                                                          "
    write(*,*)" If your program takes more than 500 steps to reach a     "
    write(*,*)" steady state, you are probably doing something wrong.    "                      
    write(*,*)"                                                          "
    write(*,*)" Numerical Computation of Compressible and Viscous Flows  "
    write(*,*)"                                  - Robert W. MacCormack  "
    write(*,*)"                                                          "
    write(*,*)"=========================================================="
    write(*,*)"                    Start Time Loop                       "
  end if

  ! Select the kind analysis type for the solver
 SELECT CASE(AnalysisType)
 CASE(Steady_NS_Analysis) ! Steady Navier-Stokes Analysis
    ! In order to improve the performance, we can use the select case
    ! outside the time loop
  SELECT CASE(Select_Kind_Solver)
   CASE( RK4_Methods )
    ! Time loop
    Main_Time_LoopRK4: do Iteration = 1_ik, TotalTimeSteps
         ! Calculate the Time Steps, based on the CFL numbers
           IF(IF_FixedDeltatime == 1) THEN
             dT0 = FixedDeltatime
         ELSE
         ! Calculate the Time Steps, based on the CFL numbers
          call Calculate_TimeSteps(Iteration)
         END IF
         ! Here we have the time marching
          call IncreamentRK4
          !call test_complete_flux
          if(ModelType == 1) then
              call Calculate_Jaco
          else
              call Singular_Calculate_Jaco
          endif
          
          call Update_Variables
          call Output_Brief_Information(Iteration)
          call OutputResults(Iteration)
          !PAUSE
         ! Check the convergence
         If(IfConverge) then
          if(MyId == 0 ) then
           WRITE(*,*)
           WRITE(*,*)
           WRITE(*,"(' Iteration has Converged at step ',I9,'.')") Iteration
           WRITE(*,*)
           WRITE(*,*)
          endif
           CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
           EXIT Main_Time_LoopRK4
         endif
    end do Main_Time_LoopRK4
    
   CASE( DPLUR_Methods )
    ! Implicit methods need more variables
    !  call Calculate_Jaco_Implicit
    ! Time loop
    Main_Time_LoopDPLUR: do Iteration = 1_ik, TotalTimeSteps
         ! Calculate the Time Steps, based on the CFL numbers
          call Calculate_TimeSteps(Iteration)
         ! Here we have the time marching 
          call DPLUR_Solver

          call Output_Brief_Information(Iteration)
          call OutputResults(Iteration)

         ! Check the convergence
         If(IfConverge) then
          if(MyId == 0 ) then
           WRITE(*,*)
           WRITE(*,*)
           WRITE(*,"(' Iteration has Converged at step ',I9,'.')") Iteration
           WRITE(*,*)
           WRITE(*,*)
          endif
           CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
           EXIT Main_Time_LoopDPLUR
         endif
          
    end do Main_Time_LoopDPLUR

   CASE( GMRES_Methods )
    ! Time loop
    Main_Time_LoopGMRES: do Iteration = 1_ik, TotalTimeSteps
         ! Calculate the Time Steps, based on the CFL numbers
          call Calculate_TimeSteps(Iteration)
         ! Here we have the time marching
           call GMRES_iteration
         ! call Gmres_SF_Solver

          call Output_Brief_Information(Iteration)
          call OutputResults(Iteration)

         ! Check the convergence
         If(IfConverge) then
          if(MyId == 0 ) then
           WRITE(*,*)
           WRITE(*,*)
           WRITE(*,"(' Iteration has Converged at step ',I9,'.')") Iteration
           WRITE(*,*)
           WRITE(*,*)
          endif
           CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
           EXIT Main_Time_LoopGMRES
         endif
       
    end do Main_Time_LoopGMRES

   CASE( MultiGridMethods )
    ! Time loop
    Main_Time_LoopMG: do Iteration = 1_ik, TotalTimeSteps
         ! Calculate the Time Steps, based on the CFL numbers
          call Calculate_TimeSteps(Iteration)
         ! Here we have the time marching      
    end do Main_Time_LoopMG

   END SELECT

 CASE(Unsteady_NS_Analysis) ! Unsteady Navier-Stokes Analysis
     !Time loop
      !if(MyId == 0 ) write(*,*)"Reading Iteration Number"
      !open(113,file='RESU/IterationNumber.dat',form='unformatted',status='old')
      !read(113)Iteration_start
      !close(113)  

     Iteration_start = 1

    ! write(*,*)"Iteration_start=",Iteration_start
      Main_Time_LoopRK4_unsteady: do Iteration = Iteration_start, TotalTimeSteps     

         IF(IF_FixedDeltatime == 1) THEN
             dT0 = FixedDeltatime
         ELSE
         ! Calculate the Time Steps, based on the CFL numbers
          call Calculate_TimeSteps(Iteration)
         END IF
          !write(*,*)dT0
         ! Here we have the time marching
          call IncreamentRK4_unsteady
          !call test_complete_flux
          call Calculate_Jaco
          call Update_Variables
          if(mod(Iteration,StepsWriteData) == 0) then
          call Calculate_Pert(Iteration)
          endif
          call Calculate_CFL
          call Output_Brief_Information(Iteration)
          call OutputResults(Iteration)
          
         ! Check the convergence
         If(IfConverge) then
          if(MyId == 0 ) then
           WRITE(*,*)
           WRITE(*,*)
           WRITE(*,"(' Iteration has Converged at step ',I9,'.')") Iteration
           WRITE(*,*)
           WRITE(*,*)
          endif
           CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
           EXIT Main_Time_LoopRK4_unsteady
         endif 
    end do Main_Time_LoopRK4_unsteady
   
 CASE(LNS_Analysis) 
  ! Linearized Navier-Stokes Analysis
  ! In this case the linearized Navier-Stokes equations are solved
  ! also in the similar form as the Navier-Stokes equations, which
  ! means the conservation form of the equations are used.
  
  ! In linearized Navier-Stokes equations, we simulate the evolution
  ! of linearized perturbations of the flow field based on the solution
  ! of a fixed base flow. 

  call LNS_Solver

  call MPI_Barrier(MPI_COMM_WORLD,ierr);
  
 CASE(HLNS_Analysis)! Harmonic Linearized Navier-Stokes Analysis
 
 END SELECT


  ! Output Final Results
  if(MyId == 0) then
    write(*,*)"=========================================================="
    write(*,*)"                                                          "
    write(*,*)" If your program takes more than 500 steps to reach a     "
    write(*,*)" steady state, you are probably doing something wrong.    "                      
    write(*,*)"                                                          "
    write(*,*)" Numerical Computation of Compressible and Viscous Flows  "
    write(*,*)"                                  - Robert W. MacCormack  "
    write(*,*)"                                                          "
    write(*,*)"=========================================================="
    write(*,*)"                    End Time Loop                         "
    write(*,*)"                                                          "
    write(*,*)"              Output the final results                    "
  end if
  
!  call OutputResults(0_ik)

  ! Deallocate memory for the variables
  
  if( MyId == 0 ) then
    write(*,*)"=========================================================="
    write(*,*)"                                                          "
    write(*,*)"            Deallocate memory for the variables           "
    write(*,*)"                                                          "
    write(*,*)"=========================================================="
  end if

  call Finilize_SFSolver

  if( MyId == 0 ) then
    write(*,*)"=========================================================="
    write(*,*)"                                                          "
    write(*,*)"              Deallocate MPI Environment                  "
    write(*,*)"                                                          "
    write(*,*)"=========================================================="
  end if
  ! Finalize the mpi
  call Parallel_Finalize


end program SFSolver