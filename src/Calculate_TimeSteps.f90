Subroutine Calculate_TimeSteps(Iteration)
  use SF_Constant,                 only: ik,rk
  use SF_CFD_Global,               only: StepChangeCFLStart,CFL_End,StepChangeCFL,CFL,&
                                        CFL_Start,CFL_End,nz_local,Ny,nx_local,Gamma,Pr_Ref,&
                                        LocalDT,DT0

  use MPI_GLOBAL,                  only: ierr,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD
  ! matrix and some arrays
  use SF_CFD_Global,               only: Mu,Rds_Xi,Rds_Eta,Rds_Zeta,Rds_Max,&
                                        Cs,BigU_Xi,BigU_Eta,BigU_Zeta,Jaco,invJacodt,&
                                        nablaxi,nablaeta,nablazeta,Rho,DIAG
  use MPI_GLOBAL,                  only: Parallel_Exchange
  implicit none
  integer( kind = ik ), intent(in) :: Iteration
  ! Local variables
  integer( kind = ik ) :: ic,jc,kc
  real( kind = rk ) :: Rds_Vis 
  real( kind = rk ) :: Max_Values(2_ik)
  real( kind = rk ) :: LocalDT_PerThread
  
  
  ! Calculate the local CFL numbers, as the simulation always runs
  ! with small CFL numbers, we can use the local CFL numbers
    if( Iteration < StepChangeCFLStart ) then
        CFL = CFL_Start
    elseif( Iteration >= StepChangeCFLStart + StepChangeCFL - 1_ik ) then
        CFL = CFL_End
    else
        CFL = CFL_Start + (CFL_End - CFL_Start) * &
              REAL(Iteration - StepChangeCFLStart,kind=rk) / REAL(StepChangeCFL - 1_ik,kind=rk)
    endif
    
    ! Calculate the Spectral Radius, take the viscous effect into account
    Rds_Vis = 0.0_rk
    DO kc = 1, nz_local
        DO jc = 1, Ny
            DO ic = 1, nx_local
             Max_Values(1) = 4.d0 / 3.d0 * Mu(ic,jc,kc);
             Max_Values(2) = Gamma / Pr_Ref * Mu(ic,jc,kc);
             Rds_Vis = maxval(Max_Values) / Rho(ic,jc,kc);
             ! Obtaing the spectral radius along each directions
               Rds_Xi(ic,jc,kc)=abs(  BigU_Xi(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis*  nablaxi(ic,jc,kc))*  nablaxi(ic,jc,kc);
              Rds_Eta(ic,jc,kc)=abs( BigU_Eta(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis* nablaeta(ic,jc,kc))* nablaeta(ic,jc,kc);
             Rds_Zeta(ic,jc,kc)=abs(BigU_Zeta(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis*nablazeta(ic,jc,kc))*nablazeta(ic,jc,kc);
             ! Here we have the dxi = deta = dzeta = 1
              Rds_Max(ic,jc,kc)=Rds_Xi(ic,jc,kc)+Rds_Zeta(ic,jc,kc)+&
                                Rds_Eta(ic,jc,kc)+jaco(ic,jc,kc)*invJacodt(ic,jc,kc);
              LocalDT(ic,jc,kc)=CFL/Rds_Max(ic,jc,kc);
            END DO
        END DO
    END DO
    
    ! Obtain the minimum value of the local time steps
    LocalDT_PerThread = minval(LocalDT(1:nx_local,1:Ny,1:nz_local))

    ! Get the minimum value of the local time steps
    CALL MPI_ALLREDUCE(LocalDT_PerThread,DT0,1,&
                       MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)

   ! The following processes copy 
   ! For wall surfaces
   jc = 1_ik;
   DO kc = 1, nz_local
     DO ic = 1, nx_local
      DIAG(ic,jc,kc) = 1.0_rk / DT0 + Rds_Xi(ic,jc,kc) &
                                   &+ Rds_Zeta(ic,jc,kc) + jaco(ic,jc,kc)*invJacodt(ic,jc,kc);
     ENDDO
   ENDDO
   ! For inner points
   DO kc = 1, nz_local    
    DO jc = 2, Ny-1
     DO ic = 1, nx_local
      DIAG(ic,jc,kc) = 1.0_rk / DT0 + Rds_Max(ic,jc,kc)
     ENDDO
    ENDDO
   ENDDO
   ! For shock boundary condition
   ! 注意激波边界的变量由运动R-H关系式，精确确定，因此对应的时间步长理论上可以取到无穷大
   jc = Ny;
   DO kc = 1, nz_local
     DO ic = 1, nx_local
      DIAG(ic,jc,kc) = Rds_Xi(ic,jc,kc) &
                  &+ Rds_Zeta(ic,jc,kc) + jaco(ic,jc,kc)*invJacodt(ic,jc,kc);
     ENDDO
   ENDDO

   ! 由于Implicit方法理论上设计到全场的Jacobi矩阵及其运算
   ! 因此对应的谱半径应该在分块的边界处进行传输
   ! 但是由于我们的计算是在每个进程内部进行的，因此我们只需要在进程间的边界处进行传输
   ! 唯一需要注意的是
   call Parallel_Exchange(DIAG)
   call Parallel_Exchange(Rds_Xi)
   call Parallel_Exchange(Rds_Eta)
   call Parallel_Exchange(Rds_Zeta)
   call Parallel_Exchange(Rds_Max)

    end Subroutine Calculate_TimeSteps
    
    
  Subroutine Calculate_CFL
  use SF_Constant,                 only: ik,rk
  use SF_CFD_Global,               only: nz_local,Ny,nx_local,Gamma,Pr_Ref,&
                                         LocalDT,DT0,CFL

  use MPI_GLOBAL,                  only: ierr,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD
  ! matrix and some arrays
  use SF_CFD_Global,               only: Mu,Rds_Xi,Rds_Eta,Rds_Zeta,Rds_Max,&
                                        Cs,BigU_Xi,BigU_Eta,BigU_Zeta,Jaco,invJacodt,&
                                        nablaxi,nablaeta,nablazeta,Rho,DIAG
  use MPI_GLOBAL,                  only: Parallel_Exchange
  implicit none
  ! Local variables
  integer( kind = ik ) :: ic,jc,kc
  real( kind = rk ) :: Rds_Vis, Rds_max_global
  real( kind = rk ) :: Max_Values(2_ik)

    ! Calculate the Spectral Radius, take the viscous effect into account
    Rds_Vis = 0.0_rk
    DO kc = 1, nz_local
        DO jc = 1, Ny
            DO ic = 1, nx_local
             Max_Values(1) = 4.d0 / 3.d0 * Mu(ic,jc,kc);
             Max_Values(2) = Gamma / Pr_Ref * Mu(ic,jc,kc);
             Rds_Vis = maxval(Max_Values) / Rho(ic,jc,kc);
             ! Obtaing the spectral radius along each directions
               Rds_Xi(ic,jc,kc)=abs(  BigU_Xi(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis*  nablaxi(ic,jc,kc))*  nablaxi(ic,jc,kc);
              Rds_Eta(ic,jc,kc)=abs( BigU_Eta(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis* nablaeta(ic,jc,kc))* nablaeta(ic,jc,kc);
             Rds_Zeta(ic,jc,kc)=abs(BigU_Zeta(ic,jc,kc))+(Cs(ic,jc,kc)+2.d0*Rds_Vis*nablazeta(ic,jc,kc))*nablazeta(ic,jc,kc);
             ! Here we have the dxi = deta = dzeta = 1
              Rds_Max(ic,jc,kc)=Rds_Xi(ic,jc,kc)+Rds_Zeta(ic,jc,kc)+&
                                Rds_Eta(ic,jc,kc)+jaco(ic,jc,kc)*invJacodt(ic,jc,kc);
            END DO
        END DO
    END DO

    Rds_max_global = maxval(Rds_Max)
    
    CFL = dT0 * Rds_max_global

   ! 由于Implicit方法理论上设计到全场的Jacobi矩阵及其运算
   ! 因此对应的谱半径应该在分块的边界处进行传输
   ! 但是由于我们的计算是在每个进程内部进行的，因此我们只需要在进程间的边界处进行传输
   ! 唯一需要注意的是
   call Parallel_Exchange(Rds_Xi)
   call Parallel_Exchange(Rds_Eta)
   call Parallel_Exchange(Rds_Zeta)
   call Parallel_Exchange(Rds_Max)

end Subroutine Calculate_CFL