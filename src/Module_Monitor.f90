Module Monitor
! This module contain all the information that has 
! been used to monitoring the residuals of the 
! present shock fitting solver        
  ! basic variables
  use SF_Constant,  only: ik,rk
  use SF_CFD_Global,only: nx_local,Ny,nz_local,StepsOutputScreen,&
                          AVE_RES_Rho,AVE_RES_RhoU,AVE_RES_RhoV,AVE_RES_RhoW,AVE_RES_RhoE

  use MPI_GLOBAL,   only: MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,MyId,NumProcess,ierr,&
                          MPI_MAX,MPI_MIN,MPI_SUM
  ! arrays and matrix variables
  use SF_CFD_Global,only: dUcons,dShockH,dShockV
  ! subs and functions
  implicit none

  ! Define the Max and Min Residuals of the present variables
  ! Rho,RhoU,RhoV,RhoW,RhoE,ShockH,ShockV  
  real( kind = rk )::MaxRhoResidual,MaxRhoUResidual,MaxRhoVResidual,MaxRhoWResidual,MaxRhoEResidual
  real( kind = rk )::MinRhoResidual,MinRhoUResidual,MinRhoVResidual,MinRhoWResidual,MinRhoEResidual
  real( kind = rk )::AveRhoResidual,AveRhoUResidual,AveRhoVResidual,AveRhoWResidual,AveRhoEResidual
  real( kind = rk )::MaxShockHResidual,MaxShcokVResidual
  real( kind = rk )::MinShockHResidual,MinShockVResidual
  real( kind = rk )::AveShockHResidual,AveShockVResidual

    contains
   
   subroutine Output_Calculate_MinMaxAveResiduals
     use SF_CFD_Global,only: nx_local,Ny,nz_local,AVE_RES_Rho,AVE_RES_RhoE
     implicit none
     real( kind = rk )::ResidualRho, MaxRhoResidualLCL, MinRhoResidualLCL, AveRhoResidualLCL
     real( kind = rk )::ResidualRhoE,MaxRhoEResidualLCL,MinRhoEResidualLCL,AveRhoEResidualLCL
     real( kind = rk )::ResidualShockH,MaxShockHResidualLCL,MinShockHResidualLCL,AveShockHResidualLCL
     real( kind = rk )::ResidualShockV,MaxShockVResidualLCL,MinShockVResidualLCL,AveShockVResidualLCL
     integer(kind =ik)::ic,jc,kc,iVar
     real( kind = rk )::temp_NumPoints,temp_NumPointsSurf
     ! 边界处的变量更新由边界条件完全控制，因此，我们在这里，仅仅考虑内部计算区域中的变量变化
     ! 即计算对应的（1：nx_local，2:Ny-1，1:nz_local）的变量的残差
     ! Calculate the Max and Min Residuals of the present process     
     ! ==========================================================
     temp_NumPoints = REAL(nz_local*nx_local*(Ny-2_ik),rk)
     temp_NumPointsSurf = REAL(nz_local*nx_local,rk)

     ResidualRho = 0.0_rk;
     iVar = 1_ik;
     do kc = 1,nz_local
      do jc = 2,Ny-1
       do ic = 1,nx_local
         ResidualRho  = ResidualRho  + abs(dUcons(iVar,ic,jc,kc));  
       enddo
      enddo
     enddo
     MaxRhoResidualLCL = MAXVAL(dUcons(iVar,1:nx_local,2:Ny-1,1:nz_local));
     MinRhoResidualLCL = MINVAL(dUcons(iVar,1:nx_local,2:Ny-1,1:nz_local));
     AveRhoResidualLCL = (ResidualRho/temp_NumPoints)/REAL(NumProcess,rk);
     ! ==========================================================
     ! ==========================================================
     ResidualRhoE = 0.0_rk;
     iVar = 5_ik;
     do kc = 1,nz_local
      do jc = 2,Ny-1        
       do ic = 1,nx_local
         ResidualRhoE = ResidualRhoE + abs(dUcons(iVar,ic,jc,kc));  
       enddo
      enddo
     enddo
     MaxRhoEResidualLCL = MAXVAL(dUcons(iVar,1:nx_local,2:Ny-1,1:nz_local));
     MinRhoEResidualLCL = MINVAL(dUcons(iVar,1:nx_local,2:Ny-1,1:nz_local));
     AveRhoeResidualLCL = (ResidualRhoE/temp_NumPoints)/REAL(NumProcess,rk);
     ! ==========================================================
     ResidualShockH = 0.0_rk;
     do kc = 1, nz_local
      do ic = 1, nx_local
       ResidualShockH = ResidualShockH + dShockH(ic,kc)*dShockH(ic,kc);     
      enddo
     enddo
     MaxShockHResidualLCL = MAXVAL(dShockH(1:nx_local,1:nz_local));
     MinShockHResidualLCL = MINVAL(dShockH(1:nx_local,1:nz_local));
     AveShockHResidualLCL = SQRT(ResidualShockH/temp_NumPointsSurf)/REAL(NumProcess,rk);
     ! ==========================================================
     ResidualShockV = 0.0_rk;
     do kc = 1, nz_local
      do ic = 1, nx_local
       ResidualShockV = ResidualShockV + dShockV(ic,kc)*dShockV(ic,kc);
      enddo
     enddo
     MaxShockVResidualLCL = MAXVAL(dShockV(1:nx_local,1:nz_local));
     MinShockVResidualLCL = MINVAL(dShockV(1:nx_local,1:nz_local));
     AveShockVResidualLCL = SQRT(ResidualShockV/temp_NumPointsSurf)/REAL(NumProcess,rk);
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr);
     
     ! Get the global Max, Min and Ave Residuals
     call MPI_Allreduce(MaxRhoResidualLCL, MaxRhoResidual, 1_ik,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(MinRhoResidualLCL, MinRhoResidual, 1_ik,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(AveRhoResidualLCL, AveRhoResidual, 1_ik,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr);
     AVE_RES_Rho = AveRhoResidual;
     
     call MPI_Allreduce(MaxRhoEResidualLCL,MaxRhoEResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(MinRhoEResidualLCL,MinRhoEResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(AveRhoEResidualLCL,AveRhoEResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr);
     AVE_RES_RhoE = AveRhoEResidual;
     
     call MPI_Allreduce(MaxShockHResidualLCL,MaxShockHResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(MinShockHResidualLCL,MinShockHResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(AveShockHResidualLCL,AveShockHResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr);
     
     call MPI_Allreduce(MaxShockVResidualLCL,MaxShcokVResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(MinShockVResidualLCL,MinShockVResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr);
     call MPI_Allreduce(AveShockVResidualLCL,AveShockVResidual,1_ik,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr);
     
   end subroutine
    
   subroutine Output_Brief_Information(Iter)
     ! Output only the ave information
     use SF_Constant,       only: Converge_Tol
     use SF_CFD_Global,     only: CFL,DT0,IfConverge
     implicit none
     integer(kind = ik),intent(in):: Iter
     
     call Output_Calculate_MinMaxAveResiduals
     
     ! Check the Converge_Tol
     if((abs(MaxRhoResidual) < Converge_Tol) .AND. (abs(MinRhoResidual) < Converge_Tol)) then
      ! Flow field is converged
      if((abs(MaxShockHResidual) < Converge_Tol) .AND. (abs(MinShockHResidual) < Converge_Tol)) then
         IfConverge = .True.
      endif
     endif
     
     ! 
     if(MyId == 0 ) then
     ! Output the Residuals
      if(mod(Iter-1_ik,20_ik*StepsOutputScreen) == 0 ) then
        write(*,*)
        write(*,*)
        write(*,*)
        write(*,"(A,F12.5,A,ES12.5)")"    CFL  =  ", CFL ,",    TimeStep DT0 =  ", DT0     
        write(*,*)
        write(*,"(A)")"========================================================================================="
        write(*,"(A)")"   Its   |   Ave Rho Res   |   Ave RhoE Res   |   Ave ShockH Res   |   Ave ShockV Res   |"
      endif
      if(mod(Iter,StepsOutputScreen) == 0) then
        write(*,34274)Iter,AVE_RES_Rho,AVE_RES_RhoE,AveShockHResidual,AveShockVResidual
      endif
      if(mod(Iter,20_ik*StepsOutputScreen) == 0 ) then
        write(*,"(A)")"========================================================================================="
      endif  
    endif
     
34274 FORMAT(I8,2X,2(e14.4,5X),2X,e14.4,7X,e14.4)
      
   end subroutine
   
   subroutine Output_Detailed_Information
     ! Output the detailed information for all residuals
     implicit none
     
   end subroutine
    
END Module Monitor