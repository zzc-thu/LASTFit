#This is the makefile used for SF_3D
# GNU compiler
# Release mode
#FC = mpif90 
#FLFLAGS = -O2 -cpp -ffree-line-length-0 -fcray-pointer
# DEBUG MODE
#FC = mpif90 
#FLFLAGS = -O0 -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace -ffree-line-length-0 -fcray-pointer
# Intel compiler
## Release mode
#FC = mpiifort 
#FLFLAGS = -O3 -fPIC -xHost -cpp -heap-arrays -fp-model fast
## DEBUG MODE  -DDEBUG -check all bounds
FC = mpif90 
FLFLAGS = -g -O0 -cpp -ffree-line-length-0 -fcray-pointer

# LAPACK library settings

LIBS = -L/Users/xiyoucheng/petsc-3.24.1/arch-darwin-c-opt-complex/lib -lflapack -lfblas

#SYSNAME := $(shell uname -s)
#ifeq ($(SYSNAME), Darwin)
  #LIBS = -L./LAPACK/lapack-3.12.1/LIB_OSX -llapack -ltmglib -lrefblas
  # $(MKLROOT)/lib/libmkl_blas95_lp64.a
  # $(MKLROOT)/lib/libmkl_lapack95_lp64.a
#  LIBS = -L$(MKLROOT)/lib -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread \
#	-lmkl_core -liomp5 -lpthread -lstdc++ 
#endif

#ifeq ($(SYSNAME), Linux)
#  LIBS = -L./LAPACK/lapack-3.12.1/LIB_LINUX -llapack -ltmglib -lrefblas 
#endif

# MKL library for Intel compiler, you can used it for better performance maybe
#MKLLIBS = -L$(MKLROOT)/lib -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread \
#	-lmkl_core -liomp5 -lpthread -lstdc++ 
#if you are using AMD CPU, you can use the AOCL library for better performance


OBJS = util.o\
		Read_Para.o\
		Initial_Grid.o\
		Calculate_Jaco.o\
		Initial_Fields.o\
		Calculate_TimeSteps.o\
		Calculate_Flux.o\
		Update_Variables.o\
		OutputResults.o\
		Calculate_AcDeri.o\
		Calculate_UshkDeri_CV1.o\
		Calculate_UshkDeri_CV2.o\
		Calculate_UshkDeri_CV3.o\
		Calculate_UshkDeri_CV4.o\
		Calculate_UshkDeri_CV5.o\
		Complete_Flux_un_d.o\
		CompleteFlux_D_mpi.o\
		LNS_shockBC_d.o\
        Calculate_Pert.o\
		Singular_Initial_Grid.o\
        Singular_Initial_Fields.o\
		main.o

PROGRAM = SFSolver

all: create_dirs $(PROGRAM)

create_dirs:
	mkdir -p INIT
	mkdir -p CheckFiles
	mkdir -p JACO
	mkdir -p RESU
	mkdir -p LNSResults
	mkdir -p Pert

$(PROGRAM): Module_Constant.o Module_CFD.o Module_Parallel.o \
			Module_FD_5th.o Module_Dir.o Module_OutputParaView.o \
			Module_SFitting.o Module_RK4.o Module_DPLUR.o \
			Module_MultiGrid.o Module_Monitor.o Module_FD_1st.o \
			Module_Krylov.o Module_Griddistributions.o Module_LNS.o $(OBJS)
	$(FC) $(FLFLAGS) -o $(PROGRAM) $^ $(LIBS)

Module_Constant.o: Module_Constant.f90
	$(FC) $(FLFLAGS) -c Module_Constant.f90

Module_CFD.o: Module_CFD.f90 Module_Constant.o
	$(FC) $(FLFLAGS) -c Module_CFD.f90

Module_Parallel.o: Module_Parallel.f90 Module_Constant.o Module_CFD.o
	$(FC) $(FLFLAGS) -c Module_Parallel.f90	

Module_FD_5th.o: Module_FD_5th.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_FD_5th.f90

Module_Griddistributions.o: Module_Griddistributions.f90
	$(FC) $(FLFLAGS) -c Module_Griddistributions.f90

#Module_FD_Weno5th.o: Module_FD_Weno5th.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
#	$(FC) $(FLFLAGS) -c Module_FD_Weno5th.f90

Module_Dir.o: Module_Dir.f90
	$(FC) $(FLFLAGS) -c Module_Dir.f90

Module_Monitor.o: Module_Monitor.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_Monitor.f90

Module_OutputParaView.o: Module_OutputParaView.f90 Module_Constant.o Module_CFD.o Module_Parallel.o Module_Dir.o
	$(FC) $(FLFLAGS) -c Module_OutputParaView.f90

Module_SFitting.o: Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_SFitting.f90

Module_RK4.o: Module_RK4.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_RK4.f90

Module_DPLUR.o: Module_DPLUR.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_DPLUR.f90 $(LIBS)

Module_MultiGrid.o: Module_MultiGrid.f90 Module_Constant.o Module_CFD.o Module_Parallel.o	
	$(FC) $(FLFLAGS) -c Module_MultiGrid.f90

Module_FD_1st.o: Module_FD_1st.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_FD_1st.f90

Module_Krylov.o: Module_Krylov.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Module_Krylov.f90 $(LIBS)

Module_LNS.o: Module_LNS.f90 Module_Constant.o Module_CFD.o
	$(FC) $(FLFLAGS) -c Module_LNS.f90

util.o: util.f90 Module_Constant.o Module_CFD.o
	$(FC) $(FLFLAGS) -c util.f90

Read_Para.o: Read_Para.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Read_Para.f90

Initial_Grid.o: Initial_Grid.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Initial_Grid.f90

Singular_Initial_Grid.o: Singular_Initial_Grid.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Singular_Initial_Grid.f90

Calculate_Jaco.o: Calculate_Jaco.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Calculate_Jaco.f90

Initial_Fields.o: Initial_Fields.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Initial_Fields.f90

Singular_Initial_Fields.o: Singular_Initial_Fields.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Singular_Initial_Fields.f90

Calculate_TimeSteps.o: Calculate_TimeSteps.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Calculate_TimeSteps.f90

Calculate_Flux.o: Calculate_Flux.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Calculate_Flux.f90

Calculate_Pert.o: Calculate_Pert.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Calculate_Pert.f90

Update_Variables.o: Update_Variables.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c Update_Variables.f90

#Calculate_Pert.o: Calculate_Pert.f90 Module_Constant.o Module_CFD.o Module_Parallel.o
#    $(FC) $(FLFLAGS) -c Calculate_Pert.f90

OutputResults.o: OutputResults.f90 Module_Constant.o Module_CFD.o Module_Parallel.o Module_OutputParaView.o
	$(FC) $(FLFLAGS) -c OutputResults.f90

Calculate_AcDeri.o: Calculate_AcDeri.f90
	$(FC) $(FLFLAGS) -c Calculate_AcDeri.f90

Calculate_UshkDeri_CV1.o: Calculate_UshkDeri_CV1.f90
	$(FC) $(FLFLAGS) -c Calculate_UshkDeri_CV1.f90

Calculate_UshkDeri_CV2.o: Calculate_UshkDeri_CV2.f90
	$(FC) $(FLFLAGS) -c Calculate_UshkDeri_CV2.f90

Calculate_UshkDeri_CV3.o: Calculate_UshkDeri_CV3.f90
	$(FC) $(FLFLAGS) -c Calculate_UshkDeri_CV3.f90

Calculate_UshkDeri_CV4.o: Calculate_UshkDeri_CV4.f90
	$(FC) $(FLFLAGS) -c Calculate_UshkDeri_CV4.f90

Calculate_UshkDeri_CV5.o: Calculate_UshkDeri_CV5.f90
	$(FC) $(FLFLAGS) -c Calculate_UshkDeri_CV5.f90

Complete_Flux_un_d.o: Complete_Flux_un_d.f90 Module_Parallel.o
	$(FC) $(FLFLAGS) -c Complete_Flux_un_d.f90

CompleteFlux_D_mpi.o: CompleteFlux_D_mpi.f90 Module_Constant.o Module_Parallel.o
	$(FC) $(FLFLAGS) -c CompleteFlux_D_mpi.f90

LNS_shockBC_d.o: LNS_shockBC_d.f90 Module_Parallel.o
	$(FC) $(FLFLAGS) -c LNS_shockBC_d.f90
	
#TestCompleteFlux.o: TestCompleteFlux.f90 Module_Constant.o Module_CFD.o CompleteFlux_D_mpi.o
#	$(FC) $(FLFLAGS) -c TestCompleteFlux.f90

main.o: main.f90 Module_Constant.o Module_CFD.o Module_Parallel.o Module_RK4.o \
Module_DPLUR.o Module_Krylov.o Module_Monitor.o Module_LNS.o
	$(FC) $(FLFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(PROGRAM) INIT/* CheckFiles/* JACO/* 
