Module SF_LNS
   use SF_Constant,           only: LaxFriSmall,Csutherland,ik,rk,overLap,LaxFriSmall,NumVar
   use SF_Constant,           only: PertForContinue
   use SF_CFD_Global,         only: Ucons0,shockH,shockV,&  
                 & nx_local,Ny,nz_local,Gamma,Mach_Ref,Re_Ref,Csthlnd_Ref,Pr_Ref,&
                 & Finv_inf,Ginv_inf,Hinv_inf,CV_inf,WallSX,WallSY,WallSZ,&
                 & WallSXdxi,WallSYdxi,WallSZdxi,WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                 & Wallnormalx,Wallnormaly,Wallnormalz,&
                 & WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,&
                 & WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&    
                 & Heta,Hetadxi,Hetadeta,Hetadzeta,&
                 & dUcons0,shockAc,Nx,LNS_start
   use SF_CFD_Global,         only: LNS_CVp0,LNS_CVp1,LNS_CVp2,LNS_CVp3,LNS_CVp4,&
                                  & LNS_ShockHp0,LNS_ShockHp1,LNS_ShockHp2,LNS_ShockHp3,LNS_ShockHp4,&
                                  & LNS_ShockVp0,LNS_ShockVp1,LNS_ShockVp2,LNS_ShockVp3,LNS_ShockVp4,&
                                  & AB_Sum_CVp0,AB_Sum_SHp0,AB_Sum_SVp0
   use SF_CFD_Global,         only: Rho_inf,U_inf,V_inf,W_inf,P_inf,C_inf,T_inf,x_grid,LNS_dt,LNS_Steps_MAX,Rd,ShockAcd,&
                                  & Pert_Type,k_infty,epsilon,shockxtau,&
                                  & Rho_pert,U_pert,V_pert,W_pert,P_pert,T_pert,Rho0,U0,V0,W0,P0,T0,&
                                  & Rho_free,U_free,V_free,W_free,P_free,T_free,shockH_steady,shockV_steady,&
                                  & Rho_pert_tau,U_pert_tau,V_pert_tau,W_pert_tau,P_pert_tau,&
                                  & Rho_pert_taud,U_pert_taud,V_pert_taud,W_pert_taud,P_pert_taud,&
                                  & StepsWriteData,OutputFileNo,IF_Continue_LNS,StepsOutputScreen,dT0,&
                                  & x_grid,y_grid,z_grid,Ucons0_steady,x_grid_steady,y_grid_steady,z_grid_steady,&
                                  & x_grid_interp,y_grid_interp,z_grid_interp,Hdst_interp,etaSC,Hdst,&
                                  & Rho_inte,U_inte,V_inte,W_inte,P_inte,T_inte,&
                                  & shocknormalxd,shocknormalyd,shocknormalzd,&
                                  & shockNormalx_steady,shockNormaly_steady,shockNormalz_steady,&
                                  & shockNormalx,shockNormaly,shockNormalz,shockXtau_steady,shockXtau,shockXtaud
   
   use MPI_GLOBAL,            only: MyId,NumProcess,Parallel_Exchange,Parallel_Exchange_NumVar,&
                                  & Parallel_Exchange_Surface,ierr,MPI_COMM_WORLD,Parallel_Exchange_interp
   use OutputParaView,        only: output_LNSResults
   use SFitting,              only: ShockRelation3D_per,ShockRelation3D,ShockRelation3D_per0
   use grid_distribution,     only: GridDistributions
   implicit none

  contains

  subroutine LNS_Solver

    implicit none
    integer( kind = ik ):: LNS_Steps,ic,jc,kc
   
    ! read ShockH and ShockV and Ucons0 from the previous NS simulation
    call ReadSteadyState

    IF(IF_Continue_LNS == 1)then
        
       call READ_LNS_RESULTS
        
    else    
    
       LNS_Steps = 0;
       
       call Calculate_freestream(LNS_Steps)
       
       call LNS_BC_update(LNS_ShockHp0,LNS_ShockVp0,LNS_Steps) 
       
       LNS_CVp4     = LNS_CVp0
       LNS_ShockHp4 = LNS_ShockHp0
       LNS_ShockVp4 = LNS_ShockVp0
       
       ucons0 = ucons0_steady + lns_cvp4
       shockh = shockh_steady + lns_shockhp4
       shockv = shockv_steady + lns_shockvp4
       shockNormalx = shockNormalx_steady + shocknormalxd
       shockNormaly = shockNormaly_steady + shocknormalyd
       shockNormalz = shockNormalz_steady + shocknormalzd
       shockXtau    = shockXtau_steady    + shockxtaud
       
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
 
       LNS_Steps = 1;
       call Calculate_freestream(LNS_Steps)

       call LNS_Solver_RHS(LNS_CVp4,LNS_ShockHp4,LNS_ShockVp4,Rd,ShockAcd,shockAc,LNS_Steps,shocknormalxd,shocknormalyd,shocknormalzd,shockXtaud)
       ! Update the solution      
           
       LNS_CVp0     = LNS_CVp4     + Rd * LNS_dt
       LNS_ShockVp3 = LNS_ShockVp4 + ShockAcd * LNS_dt
       LNS_ShockHp3 = LNS_ShockHp4 + LNS_ShockVp4 * LNS_dt

       !write(*,*)shockacd(:,1)
       call LNS_BC_update(LNS_ShockHp3,LNS_ShockVp3,LNS_Steps)
   
       !write(*,*),'LNS_ShockVp3',LNS_ShockVp3(:,1)
        !output the linear solution
         LNS_CVp3 = LNS_CVp0
      
       Ucons0 = Ucons0_steady + LNS_CVp3
       shockH = shockH_steady + LNS_shockHp3
       shockV = shockV_steady + LNS_shockVp3
       shockNormalx = shockNormalx_steady + shocknormalxd
       shockNormaly = shockNormaly_steady + shocknormalyd
       shockNormalz = shockNormalz_steady + shocknormalzd      
       shockXtau    = shockXtau_steady    + shockxtaud
       
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
       
       
       LNS_Steps = 2;
       ! 2nd Adams-Bashforth method
       !====The second step for using n-1,n-2,to obtain n
       call Calculate_freestream(LNS_Steps)       
       
       AB_Sum_CVp0 = 1.5d0 * LNS_CVp3     - 0.5d0 * LNS_CVp4
       AB_Sum_SHp0 = 1.5d0 * LNS_ShockHp3 - 0.5d0 * LNS_ShockHp4
       AB_Sum_SVp0 = 1.5d0 * LNS_ShockVp3 - 0.5d0 * LNS_ShockVp4

       
       call LNS_Solver_RHS(AB_Sum_CVp0,AB_Sum_SHp0,AB_Sum_SVp0,Rd,ShockAcd,shockAc,LNS_Steps,shocknormalxd,shocknormalyd,shocknormalzd,shockXtaud)
       
        ! Update the solution 
     
       LNS_CVp0     = LNS_CVp3     + Rd * LNS_dt
       LNS_ShockVp2 = LNS_ShockVp3 + ShockAcd * LNS_dt
       LNS_ShockHp2 = LNS_ShockHp3 + LNS_ShockVp3 * LNS_dt
       
       call LNS_BC_update(LNS_ShockHp2,LNS_ShockVp2,LNS_Steps)
       !write(*,*),'LNS_ShockHp2',LNS_ShockHp2(:,1)
       !write(*,*),'LNS_ShockVp2',LNS_ShockVp2(:,1)
         LNS_CVp2 = LNS_CVp0
       
       Ucons0 = Ucons0_steady + LNS_CVp2  
       shockH = shockH_steady + LNS_shockHp2
       shockV = shockV_steady + LNS_shockVp2
       shockNormalx = shockNormalx_steady + shocknormalxd
       shockNormaly = shockNormaly_steady + shocknormalyd
       shockNormalz = shockNormalz_steady + shocknormalzd  
       shockXtau    = shockXtau_steady    + shockxtaud
       
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
     
      LNS_Steps = 3;
   !    3rd Adams-Bashforth method
       call Calculate_freestream(LNS_Steps)

       AB_Sum_CVp0 = 23.d0/12.d0 * LNS_CVp2     - 4.0d0/3.d0 * LNS_CVp3     + 5.d0/12.d0 * LNS_CVp4
       AB_Sum_SHp0 = 23.d0/12.d0 * LNS_ShockHp2 - 4.0d0/3.d0 * LNS_ShockHp3 + 5.d0/12.d0 * LNS_ShockHp4
       AB_Sum_SVp0 = 23.d0/12.d0 * LNS_ShockVp2 - 4.0d0/3.d0 * LNS_ShockVp3 + 5.d0/12.d0 * LNS_ShockVp4
       
       call LNS_Solver_RHS(AB_Sum_CVp0,AB_Sum_SHp0,AB_Sum_SVp0,Rd,ShockAcd,shockAc,LNS_Steps,shocknormalxd,shocknormalyd,shocknormalzd,shockXtaud)
       
        ! Update the solution 
       LNS_CVp0     = LNS_CVp2     + Rd * LNS_dt
       LNS_ShockVp1 = LNS_ShockVp2 + ShockAcd * LNS_dt
       LNS_ShockHp1 = LNS_ShockHp2 + LNS_ShockVp2  * LNS_dt
        
       call LNS_BC_Update(LNS_ShockHp1,LNS_ShockVp1,LNS_Steps) 

        LNS_CVp1 = LNS_CVp0    
       
       Ucons0 = Ucons0_steady + LNS_CVp1
       shockH = shockH_steady + LNS_shockHp1
       shockV = shockV_steady + LNS_shockVp1
       shockNormalx = shockNormalx_steady + shocknormalxd
       shockNormaly = shockNormaly_steady + shocknormalyd
       shockNormalz = shockNormalz_steady + shocknormalzd
       shockXtau    = shockXtau_steady    + shockxtaud
       
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
        !output the linear solution
        !call OutputLNS_Results(LNS_Steps)
    ENDIF     
    
    !write(*,*),'U_pert(1,:,1)',U_pert(1,:,1) 
    IF(IF_Continue_LNS == 1)then   
    if(MyId == 0 ) write(*,*)"Readind Iteration Number"
      open(113,file='LNSResults/IterationNumber.dat',form='unformatted',status='old')
      read(113)LNS_start
      close(113) 
     else
      LNS_start = 4;
     endif
     
     DO LNS_Steps = LNS_start, LNS_Steps_Max
        !====The fourth step for using n-1,n-2,n-3,n-4,to obtain n
       call Calculate_freestream(LNS_Steps)

       AB_Sum_CVp0 = 55.d0/24.d0 * LNS_CVp1     - 59.d0/24.d0 * LNS_CVp2     + 37.d0/24.d0 * LNS_CVp3     - 3.d0/8.d0 * LNS_CVp4
       AB_Sum_SHp0 = 55.d0/24.d0 * LNS_ShockHp1 - 59.d0/24.d0 * LNS_ShockHp2 + 37.d0/24.d0 * LNS_ShockHp3 - 3.d0/8.d0 * LNS_ShockHp4
       AB_Sum_SVp0 = 55.d0/24.d0 * LNS_ShockVp1 - 59.d0/24.d0 * LNS_ShockVp2 + 37.d0/24.d0 * LNS_ShockVp3 - 3.d0/8.d0 * LNS_ShockVp4

       
       call LNS_Solver_RHS(AB_Sum_CVp0,AB_Sum_SHp0,AB_Sum_SVp0,Rd,ShockAcd,shockAc,LNS_Steps,shocknormalxd,shocknormalyd,shocknormalzd,shockXtaud)

        ! Update the solution 
       
       LNS_CVp0     = LNS_CVp1     + Rd * LNS_dt
       LNS_ShockVp0 = LNS_ShockVp1 + ShockAcd * LNS_dt
       LNS_ShockHp0 = LNS_ShockHp1 + LNS_ShockVp1  * LNS_dt
       
       call LNS_BC_Update(LNS_ShockHp0,LNS_ShockVp0,LNS_Steps) 
         
       Ucons0 = Ucons0_steady + LNS_CVp0
       shockH = shockH_steady + LNS_shockHp0
       shockV = shockV_steady + LNS_shockVp0
       shockNormalx = shockNormalx_steady + shocknormalxd
       shockNormaly = shockNormaly_steady + shocknormalyd
       shockNormalz = shockNormalz_steady + shocknormalzd       
       shockXtau    = shockXtau_steady    + shockxtaud
       
  do kc = 1-overLAP,nz_local+overLAP
    do jc = 1,Ny
      do ic = 1-overLAP,nx_local+overLAP
        X_grid(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * ShockH(ic,kc) * Heta(ic,jc,kc);
      end do
    end do
  end do
       ! Update the fields for 1,2,3,4
        LNS_CVp4 = LNS_CVp3
        LNS_ShockHp4 = LNS_ShockHp3
        LNS_ShockVp4 = LNS_ShockVp3

        LNS_CVp3 = LNS_CVp2
        LNS_ShockHp3 = LNS_ShockHp2
        LNS_ShockVp3 = LNS_ShockVp2

        LNS_CVp2 = LNS_CVp1
        LNS_ShockHp2 = LNS_ShockHp1
        LNS_ShockVp2 = LNS_ShockVp1

        LNS_CVp1 = LNS_CVp0
        LNS_ShockHp1 = LNS_ShockHp0
        LNS_ShockVp1 = LNS_ShockVp0     
        
    
    do ic = 1-overLap, nx_local+overLap
        do jc = 1, Ny
            do kc = 1-overLap, nz_local+overLap
                call CVtoU_pert(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),LNS_CVp0(1:NumVar,ic,jc,kc),&
                             &   Rho_pert(ic,jc,kc), U_pert(ic,jc,kc), V_pert(ic,jc,kc), W_pert(ic,jc,kc), P_pert(ic,jc,kc))
                !T_pert(ic,jc,kc) = (Gamma * Mach_Ref * Mach_Ref *( P0(ic,jc,kc)+ P_pert(ic,jc,kc)) - Rho_pert(ic,jc,kc) * T0(ic,jc,kc)- Rho0(ic,jc,kc) * T0(ic,jc,kc)) / Rho0(ic,jc,kc)
                T_pert(ic,jc,kc) = (Gamma * Mach_Ref * Mach_Ref * P_pert(ic,jc,kc) - Rho_pert(ic,jc,kc) * T0(ic,jc,kc))/Rho0(ic,jc,kc)
            end do
        enddo
    enddo        
        !write(*,*),'U_pert(1,:,1)',U_pert(1,:,1) 
        if(mod(LNS_Steps,StepsOutputScreen) == 0) then
          if (MyId == 0) then
          write(*,*)"LNS_Steps =",LNS_Steps
          
          write(*,*)"Total time =",LNS_Steps*LNS_dt 
          endif
        endif
        
        !output the linear solution
        if(mod(LNS_Steps,StepsWriteData) == 0) then
            
        call OutputLNS_Results(LNS_Steps)
        endif
        !stop
     ENDDO
     
  end subroutine LNS_Solver

  ! 这里是一个对于计算线化SF_NS方程右端项的一个封装
  ! 主要是为了方便调用,这里直接给出对应的线化input的output
  !    call CompleteFlux_D_mpi(ucons0, ucons0d, shockh, shockhd, shockv, &
  !    & shockvd, nx_local, ny_local, nz_local, overlap, gamma, mach_ref, &
  !    & re_ref, csthlnd_ref, pr_ref, laxfrismall, finv_inf, ginv_inf, hinv_inf&
  !    & , cv_inf, wallsx, wallsy, wallsz, wallsxdxi, wallsydxi, wallszdxi, &
  !    & wallsxdzeta, wallsydzeta, wallszdzeta, wallnormalx, wallnormaly, &
  !    & wallnormalz, wallnormalxdxi, wallnormalydxi, wallnormalzdxi, &
  !    & wallnormalxdzeta, wallnormalydzeta, wallnormalzdzeta, heta, hetadxi, &
  !    & hetadeta, hetadzeta, r, rd, shockac_new, shockacd) 
  subroutine LNS_Solver_RHS(ucons0d_local,shockHd_local,shockvd_local,Rd_local,shockAcd_local,shockAc_local,Iter,shocknormalxd_local,shocknormalyd_local,shocknormalzd_local,shockXtaud_local)
    implicit none
    integer( kind = ik ), intent(in) :: Iter
    real( kind = rk ), dimension(1:NumVar,1-overlap:nx_local+overlap,1:Ny,1-overlap:nz_local+overlap),intent(in) :: ucons0d_local
    real( kind = rk ), dimension(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap), intent(in) :: shockHd_local,shockvd_local
    real( kind = rk ), dimension(1:NumVar,1-overlap:nx_local+overlap,1:Ny,1-overlap:nz_local+overlap),intent(out):: Rd_local
    real( kind = rk ), dimension(1-overlap:nx_local+overLap,1-overlap:nz_local+overLap), intent(out) :: shockAcd_local
    real( kind = rk ), dimension(1-overlap:nx_local+overLap,1-overlap:nz_local+overLap), intent(out) :: shocknormalxd_local,shocknormalyd_local,shocknormalzd_local,shockXtaud_local
    real (kind = rk ), dimension(1:NumVar,1-overlap:nx_local+overlap, ny, 1-overlap:nz_local+overlap) :: R_local
    real( kind = rk ), dimension(1-overlap:nx_local+overLap,1-overlap:nz_local+overLap), intent(out) :: shockAc_local

    R_local = 0.0_rk
    Rd_local = 0.0_rk
    shockAc_local = 0.0_rk
    shockAcd_local = 0.0_rk

     call COMPLETEFLUX_UN_D(Ucons0_steady,ucons0d_local,shockH_steady,shockHd_local,shockV_steady,shockVd_local,&
                        &Rho_free,Rho_free-Rho_inf,U_free,U_free-U_inf,V_free,V_free-V_inf,W_free,W_free-W_inf,P_free,P_free-P_inf,&
                        &Rho_pert_tau,Rho_pert_tau,U_pert_tau,U_pert_tau,V_pert_tau,V_pert_tau,W_pert_tau,W_pert_tau,&
                        &P_pert_tau,P_pert_tau,nx_local,Ny,nz_local,overLap,&
                        &Gamma,Mach_Ref,Re_Ref,Csthlnd_Ref,Pr_Ref,LaxFriSmall,&
                        &Finv_inf,Ginv_inf,Hinv_inf,CV_inf,WallSX,WallSY,WallSZ,&
                        &WallSXdxi,WallSYdxi,WallSZdxi,WallSXdzeta,WallSYdzeta,WallSZdzeta,&
                        &Wallnormalx,Wallnormaly,Wallnormalz,&
                        &WallNormalXdxi,WallNormalYdxi,WallNormalZdxi,WallNormalXdzeta,WallNormalYdzeta,WallNormalZdzeta,&    
                        &Heta,Hetadxi,Hetadeta,Hetadzeta,&
                        &R_local,Rd_local,shockAc_local,shockAcd_local,&
                        &shocknormalxd_local,shocknormalyd_local,shocknormalzd_local,shockXtaud_local)
    
    !call completeflux_d_mpi(ucons0_steady, ucons0d_local, shockh_steady, shockhd_local, shockv_steady, &
    !& shockvd_local, nx_local, ny, nz_local, overlap, gamma, mach_ref, &
    !& re_ref, csthlnd_ref, pr_ref, laxfrismall, finv_inf, ginv_inf, hinv_inf&
    !& , cv_inf, wallsx, wallsy, wallsz, wallsxdxi, wallsydxi, wallszdxi, &
    !& wallsxdzeta, wallsydzeta, wallszdzeta, wallnormalx, wallnormaly, &
    !& wallnormalz, wallnormalxdxi, wallnormalydxi, wallnormalzdxi, &
    !& wallnormalxdzeta, wallnormalydzeta, wallnormalzdzeta, heta, hetadxi, &
    !& hetadeta, hetadzeta, r_local, rd_local, shockac_local, shockacd_local) 
      
  end subroutine LNS_Solver_RHS

  
  subroutine LNS_shockBC(shockHd_local,shockvd_local,sRhod_local,sUd_local,sVd_local,sWd_local,sPd_local,sTd_local)
    implicit none
    
    real( kind = rk ), dimension(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap), intent(in) :: shockHd_local,shockvd_local
    real( kind = rk ), dimension(1-overlap:nx_local+overlap,1:Ny,1-overlap:nz_local+overlap),intent(out):: sRhod_local,sUd_local,sVd_local,sWd_local,sPd_local,sTd_local
    real (kind = rk ), dimension(1-overlap:nx_local+overlap,1:Ny,1-overlap:nz_local+overlap) :: sRho_local,sU_local,sV_local,sW_local,sP_local,sT_local
     
    call LINEARIZED_SHOCKBC_D(shockh_steady, shockhd_local, shockv_steady, shockvd_local, &
                            & Rho_free,Rho_free-Rho_inf,U_free,U_free-U_inf,V_free,V_free-V_inf,W_free,W_free-W_inf,&
                            & P_free,P_free-P_inf,nx_local, ny, nz_local, overlap, gamma, mach_ref, re_ref,&
                            & wallsx, wallsy, wallsz, k_infty, epsilon, pert_type, wallsxdxi, wallsydxi, &
                            & wallszdxi, wallsxdzeta, wallsydzeta, wallszdzeta, wallnormalx, &
                            & wallnormaly, wallnormalz, wallnormalxdxi, wallnormalydxi, &
                            & wallnormalzdxi, wallnormalxdzeta, wallnormalydzeta, wallnormalzdzeta, &
                            & heta, hetadxi, hetadzeta, hetadeta, srho_local, srhod_local, su_local, sud_local, sv_local, svd_local, sw_local&
                            & , swd_local, sp_local, spd_local, st_local, std_local)
    
 
    !write(*,*),'sRho_local(:,Ny,1)',sRho_local(:,Ny,1)
    !write(*,*),'sRhod_local(:,Ny,1)',sRhod_local(:,Ny,1)
    
  end subroutine LNS_shockBC
  
  subroutine ReadSteadyState
    ! This processes are very similar to previous one
    use SF_Constant,   only: FilesForContinue
    use MPI_GLOBAL,    only: MyId,ierr,MPI_COMM_WORLD
    use OutputParaView,only: output_results

    implicit none
    logical:: dir_exists
    integer( kind = ik ):: ic,jc,kc,iVar
    character(len=100)::FILENAME_FLOW
    character(len=100)::FILENAME_SHK    

    ! Here we read the initial fields from previous simulation
    ! in order to increase the efficient of the simulation
    ! all the processes read the sepearated files simultaneously
    if(MyId == 0) then
      write(*,'(A)')"We are reading the Steady fields from previous simulation"
      inquire(file="./RESU_STEADY/SF_Results00000.flowsfg", exist=dir_exists)
       if (.not. dir_exists) then
         write(*,'(A)')"The CONT file does not exist, you should have a steady state"
          stop
       else
         write(*,'(A)')"The CONT file exists, we will read the steady state from it"
       end if     
    endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
      
      ! Read the previous results
      write(filename_flow,'(A,A,I5.5,A)')'RESU_STEADY/',trim(FilesForContinue),MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow"
      open(112,file=filename_flow,form='unformatted',status='old')
       read(112)Ucons0_steady
       read(112)x_grid_steady
       read(112)y_grid_steady
       read(112)z_grid_steady
      close(112)
      
      ! Obtain the basic variables
    do ic = 1,nx_local
      do jc = 1,Ny
        do kc = 1,nz_local  
           call CVtoU_CPG(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),Ucons0_steady(1:NumVar,ic,jc,kc))
           T0(ic,jc,kc) = Gamma * Mach_ref * Mach_ref * P0(ic,jc,kc) / Rho0(ic,jc,kc) 
        end do
      end do
    end do
      
      ! Read the previous shock informations
      write(filename_shk,'(A,A,I5.5,A)')'RESU_STEADY/',trim(FilesForContinue),MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Reading ShockH and ShockV"
     open(113,file=filename_shk,form='unformatted',status='old')
      read(113)ShockH_steady
      read(113)ShockV_steady
     close(113)

     !Ucons0 = Ucons0_steady
     !shockH = ShockH_steady
     !shockV = ShockV_steady
      ! Update the Jacobian
      !call Calculate_Jaco
     
    if(MyId == 0) then
      write(*,'(A)')"We have read the steady fields from previous simulation"
      write(*,"(A)")"========================================================================================="
    end if    

    call Parallel_Exchange(Rho0)
    call Parallel_Exchange(U0)
    call Parallel_Exchange(V0)
    call Parallel_Exchange(W0)
    call Parallel_Exchange(P0)
    call Parallel_Exchange(T0)
    call Parallel_Exchange_NumVar(Ucons0_steady)
    !call Parallel_Exchange_NumVar(Ucons0)
    call Parallel_Exchange(X_grid_steady)
    call Parallel_Exchange(Y_grid_steady)
    call Parallel_Exchange(Z_grid_steady)
    ! Exchange the shock related variables
    call Parallel_Exchange_Surface(ShockH_steady)
    call Parallel_Exchange_Surface(ShockV_steady)
    !call Parallel_Exchange_Surface(shockH)
    !call Parallel_Exchange_Surface(shockV)
    
    !! Checking parts
    !!#ifdef DEBUG
    !if( MyId == 0 ) then
    !  write(*,'(A)')" Output the steady state results for checking "
    !end if
    !
    !do ic = 1,nx_local
    !  do jc = 1,Ny
    !      do kc = 1,nz_local
    !      	 call CVtoU_CPG(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
    !         T0(ic,jc,kc) = Gamma * Mach_ref * Mach_ref * P0(ic,jc,kc) / Rho0(ic,jc,kc) 
    !      end do
    !  end do
    !end do
    !OutputFileNo =
    !call output_results
    !
    call MPI_Barrier(MPI_COMM_WORLD,ierr);
    !stop
    !!#endif

  end subroutine ReadSteadyState

  subroutine Read_LNS_Results
   
    implicit none
    logical:: dir_exists
    character(len=200)::FILENAME_FLOW
    character(len=200)::FILENAME_SHK    
   
    ! Here we read the initial fields from previous simulation
    ! in order to increase the efficient of the simulation
    ! all the processes read the sepearated files simultaneously
    if(MyId == 0) then
      write(*,'(A)')"We are reading the LNS results from previous simulation"
      inquire(file="./LNSResults/Pert_Results_t0_00000.flowsfg", exist=dir_exists)
       if (.not. dir_exists) then
         write(*,'(A)')"The CONT file does not exist, you should have a LNS results"
          stop
       else
         write(*,'(A)')"The CONT file exists, we will read the LNS results from it"
       end if     
    endif
    
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr);
      
      ! Read the previous shock informations
      write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t0_',MyId,'.shksfg'
      if(MyId == 0 ) write(*,*)"Reading perturbaiton of ShockH and ShockV at Step = n" 
      open(1111,file=filename_shk,form='unformatted',status='old')
       read(1111)LNS_ShockHp1
       read(1111)LNS_ShockVp1
      close(1111)

      write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t1_',MyId,'.shksfg'
      if(MyId == 0 ) write(*,*)"Reading perturbaiton of ShockH and ShockV at Step = n-1" 
      open(1112,file=filename_shk,form='unformatted',status='old')
       read(1112)LNS_ShockHp2
       read(1112)LNS_ShockVp2
      close(1112) 
      
      write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t2_',MyId,'.shksfg'
      if(MyId == 0 ) write(*,*)"Reading perturbaiton of ShockH and ShockV at Step = n-2" 
      open(1113,file=filename_shk,form='unformatted',status='old')
       read(1113)LNS_ShockHp3
       read(1113)LNS_ShockVp3
      close(1113) 
      
      write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t3_',MyId,'.shksfg'
      if(MyId == 0 ) write(*,*)"Reading perturbaiton of ShockH and ShockV at Step = n-3" 
      open(1114,file=filename_shk,form='unformatted',status='old')
       read(1114)LNS_ShockHp4
       read(1114)LNS_ShockVp4
      close(1114) 
      
      !Read flow Perturbation
      write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t0_',MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow Perturbation at Step = n" 
      open(1115,file=filename_flow,form='unformatted',status='old')
       read(1115)LNS_CVp1
      close(1115)

      write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t1_',MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow Perturbation at Step = n-1" 
      open(1116,file=filename_flow,form='unformatted',status='old')
       read(1116)LNS_CVp2
      close(1116) 

      write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t2_',MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow Perturbation at Step = n-2" 
      open(1117,file=filename_flow,form='unformatted',status='old')
       read(1117)LNS_CVp3
      close(1117) 
      
      write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t3_',MyId,'.flowsfg'
      if(MyId == 0 ) write(*,*)"Reading flow Perturbation at Step = n-3" 
      open(1118,file=filename_flow,form='unformatted',status='old')
       read(1118)LNS_CVp4
      close(1118)  
      
      
    if(MyId == 0) then
      write(*,'(A)')"We have read the LNS Results from previous simulation"
      write(*,"(A)")"========================================================================================="
    end if    
  end subroutine Read_LNS_Results
  
  subroutine LNS_BC_update(LNS_shockHp_local,LNS_shockVp_local,Iter)
    implicit none
    ! Here we need to define the perturbations 
    integer( kind = ik ), intent(in) :: Iter
    real( kind = rk ), dimension(1-overlap:nx_local+overlap,1-overlap:nz_local+overlap), intent(in) :: LNS_shockHp_local,LNS_shockVp_local
    integer( kind = ik ):: ic,jc,kc      
    real( kind = rk ), dimension(:,:,:), allocatable :: sRho,sU,sV,sW,sP,sT

    allocate(sRho(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(sU(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))    
    allocate(sV(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(sW(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(sP(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(sT(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    
    sRho = 0.d0
    sU = 0.d0
    sV = 0.d0
    sW = 0.d0
    sP = 0.d0
    sT = 0.d0
    
    !transfer basic variables
    do ic = 1-overLap, nx_local+overLap
        do jc = 1, Ny
            do kc = 1-overLap, nz_local+overLap
                call CVtoU_pert(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),LNS_CVp0(1:NumVar,ic,jc,kc),&
                             &  Rho_pert(ic,jc,kc), U_pert(ic,jc,kc), V_pert(ic,jc,kc), W_pert(ic,jc,kc), P_pert(ic,jc,kc))
                !T_pert(ic,jc,kc) = (Gamma * Mach_Ref * Mach_Ref *( P0(ic,jc,kc)+ P_pert(ic,jc,kc)) - Rho_pert(ic,jc,kc) * T0(ic,jc,kc)- Rho0(ic,jc,kc) * T0(ic,jc,kc)) / Rho0(ic,jc,kc)
                T_pert(ic,jc,kc) = (Gamma * Mach_Ref * Mach_Ref * P_pert(ic,jc,kc) - Rho_pert(ic,jc,kc) * T0(ic,jc,kc))/Rho0(ic,jc,kc)
            end do
        enddo
    enddo
   
   call LNS_shockBC(LNS_shockHp_local,LNS_shockVp_local,sRho,sU,sV,sW,sP,sT)
   
   ! Here we update the perturbations along the shock surfaces
    jc = Ny
    DO kc = 1-overLap, nz_local+overLap
      DO ic = 1-overLap, nx_local+overLap
       Rho_pert(ic,jc,kc) = sRho(ic,jc,kc)
         U_pert(ic,jc,kc) = sU(ic,jc,kc)
         V_pert(ic,jc,kc) = sV(ic,jc,kc)
         W_pert(ic,jc,kc) = sW(ic,jc,kc)
         P_pert(ic,jc,kc) = sP(ic,jc,kc)
         T_pert(ic,jc,kc) = sT(ic,jc,kc)
             call UtoCV_pert(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),&
                      &  Rho_pert(ic,jc,kc), U_pert(ic,jc,kc), V_pert(ic,jc,kc), W_pert(ic,jc,kc), P_pert(ic,jc,kc),&
                      &  LNS_CVp0(1:NumVar,ic,jc,kc))
      ENDDO
    ENDDO
    
    jc = 1
    DO kc = 1-overLap, nz_local+overLap
      DO ic = 1-overLap, nx_local+overLap
      U_pert(ic,jc,kc) = 0.d0;
      V_pert(ic,jc,kc) = 0.d0;
      W_pert(ic,jc,kc) = 0.d0;
      T_pert(ic,jc,kc) = 0.d0;
      P_pert(ic,jc,kc) = (Rho_pert(ic,jc,kc)*T0(ic,jc,kc))/(gamma * Mach_ref * Mach_ref);
            call UtoCV_pert(Rho0(ic,jc,kc),U0(ic,jc,kc),V0(ic,jc,kc),W0(ic,jc,kc),P0(ic,jc,kc),&
                     &  Rho_pert(ic,jc,kc), U_pert(ic,jc,kc), V_pert(ic,jc,kc), W_pert(ic,jc,kc), P_pert(ic,jc,kc),&
                     &  LNS_CVp0(1:NumVar,ic,jc,kc))
         ENDDO
    ENDDO  
    
    
    call Parallel_Exchange(Rho_pert)
    call Parallel_Exchange(U_pert)
    call Parallel_Exchange(V_pert)
    call Parallel_Exchange(W_pert)
    call Parallel_Exchange(P_pert)
    call Parallel_Exchange(T_pert)
    call Parallel_Exchange_NumVar(LNS_CVp0) 
   
     deallocate(sRho)
     deallocate(sU)
     deallocate(sV)
     deallocate(sW)
     deallocate(sP)
     deallocate(sT)
     
  end subroutine LNS_BC_update
  
  subroutine OutputLNS_Results(Iter)
  
    implicit none
    integer( kind = ik ), intent(in) :: Iter
    integer( kind = ik ):: ic,jc,kc,loc_jc
    character(len=200)filename
    character(len=200)filename_shk,filename_flow
    
    !if(mod(Iter,StepsWriteData) == 0) then
        
    ! Here we output the results of the LNS results
    !if(MyId == 0) then
    !  write(*,"(A)")"========================================================================================="
    !endif
 
    ! Output the perturbation flow fields
     write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t0_',MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Saving perturbaiton of ShockH and ShockV at Step = n"
     open(111,file=filename_shk,form='unformatted',status='replace')
     write(111)LNS_ShockHp1
     write(111)LNS_ShockVp1
     close(111)
    
     write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t1_',MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Saving perturbaiton of ShockH and ShockV Step = n-1"
     open(112,file=filename_shk,form='unformatted',status='replace')
     write(112)LNS_ShockHp2
     write(112)LNS_ShockVp2
     close(112)
     
     write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t2_',MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Saving perturbaiton of ShockH and ShockV Step = n-2"
     open(113,file=filename_shk,form='unformatted',status='replace')
     write(113)LNS_ShockHp3
     write(113)LNS_ShockVp3
     close(113)
     
     write(filename_shk,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t3_',MyId,'.shksfg'
     if(MyId == 0 ) write(*,*)"Saving perturbaiton of ShockH and ShockV Step = n-3"
     open(114,file=filename_shk,form='unformatted',status='replace')
     write(114)LNS_ShockHp4
     write(114)LNS_ShockVp4
     close(114) 
     if(MyId == 0 ) write(*,*)"Saving perturbaiton of ShockH and ShockV finished"
     
     write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t0_',MyId,'.flowsfg'
     if(MyId == 0 ) write(*,*)"Saving flow Perturbation Step = n"
     open(115,file=filename_flow,form='unformatted',status='replace')
     write(115)LNS_CVp1
     close(115)
     
     write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t1_',MyId,'.flowsfg'
     if(MyId == 0 ) write(*,*)"Saving flow Perturbation Step = n-1"
     open(116,file=filename_flow,form='unformatted',status='replace')
     write(116)LNS_CVp2
     close(116)
  
     write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t2_',MyId,'.flowsfg'
     if(MyId == 0 ) write(*,*)"Saving flow Perturbation Step = n-2"
     open(117,file=filename_flow,form='unformatted',status='replace')
     write(117)LNS_CVp3
     close(117)
     
     write(filename_flow,'(A,A,A,I5.5,A)')'./LNSResults/',trim(PertForContinue),'_t3_',MyId,'.flowsfg'
     if(MyId == 0 ) write(*,*)"Saving flow Perturbation Step = n-3"
     open(118,file=filename_flow,form='unformatted',status='replace')
     write(118)LNS_CVp4
     close(118)
     
     !if(MyId == 0 ) write(*,*)"Saving Iteration Number"
     !open(113,file='./LNSResults/IterationNumber.dat',form='unformatted',status='replace')
     !write(113)Iter + 1
     !close(113)
     
     if(MyId == 0 ) write(*,*)"Saving flow field finished"
    ! Interp_LNS_Results
     
  call GridDistributions(Ny,(Ny+1)/2,(Ny+1)/2,1.05d0,etaSC)  
  
  do kc = 1-overLap, nz_local+overLap
    do jc = 1,Ny
     do ic = 1-overLap, nx_local+overLap
         Hdst(ic,jc,kc) = ShockH(ic,kc)* Heta(ic,jc,kc);
     enddo
    enddo
   enddo
  
  do kc = 1-overLap, nz_local+overLap
    do jc = 1,2*Ny-1
     do ic = 1-overLap, nx_local+overLap
          Hdst_interp(ic,jc,kc)   = shockH(ic,kc) * etaSC(jc)
          X_grid_interp(ic,jc,kc) = WallSx(ic,kc) + WallNormalX(ic,kc) * Hdst_interp(ic,jc,kc);
          Y_grid_interp(ic,jc,kc) = WallSy(ic,kc) + WallNormalY(ic,kc) * Hdst_interp(ic,jc,kc);
          Z_grid_interp(ic,jc,kc) = WallSz(ic,kc) + WallNormalZ(ic,kc) * Hdst_interp(ic,jc,kc);
       enddo
    enddo
  enddo
 
    do kc = 1, nz_local
      do ic = 1,nx_local
          do jc = 1,2*Ny-1
          	 if (Hdst_interp(ic,jc,kc) .GE. ShockH_steady(ic,kc))then
              loc_jc = jc
             exit 
             endif
          end do
          
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),Rho_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),Rho_inte(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),U_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),U_inte(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),V_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),V_inte(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),W_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),W_inte(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),P_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),P_inte(ic,1:loc_jc,kc),loc_jc)
    call CubicSplineInterp1D(Hdst(ic,1:Ny,kc),T_pert(ic,1:Ny,kc),Hdst_interp(ic,1:loc_jc,kc),T_inte(ic,1:loc_jc,kc),loc_jc)
          
          Rho_inte(ic,loc_jc:2*Ny-1,kc)= Rho_free(ic,Ny,kc)- Rho_inf
          U_inte(ic,loc_jc:2*Ny-1,kc)  = U_free(ic,Ny,kc)  - U_inf
          V_inte(ic,loc_jc:2*Ny-1,kc)  = V_free(ic,Ny,kc)  - V_inf
          W_inte(ic,loc_jc:2*Ny-1,kc)  = W_free(ic,Ny,kc)  - W_inf
          P_inte(ic,loc_jc:2*Ny-1,kc)  = P_free(ic,Ny,kc)  - P_inf
          T_inte(ic,loc_jc:2*Ny-1,kc)  = T_free(ic,Ny,kc)  - T_inf
      enddo 
    enddo
    
        if(MyId == 0) then
        write(filename_flow,'(A,A,"_",I5.5,A)')'LNSResults/',trim(PertForContinue),OutputFileNo,'.dat'
        open(55,file=filename_flow,form='unformatted',access = 'stream')
        endif
        call Parallel_OutputFunction(55,x_grid(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,y_grid(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,Rho_pert(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,U_pert(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,V_pert(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,W_pert(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,P_pert(:,:,4),Nx,Ny,nx_local,Ny)
        call Parallel_OutputFunction(55,T_pert(:,:,4),Nx,Ny,nx_local,Ny)    
        
        if(MyId == 0) then
        close(55)
        endif   
    ! Output the perturbation shock height/the perturbation shock velocity
     OutputFileNo = OutputFileNo + 1; 
     
    call Parallel_Exchange_interp(Rho_inte)
    call Parallel_Exchange_interp(U_inte)
    call Parallel_Exchange_interp(V_inte)
    call Parallel_Exchange_interp(W_inte)    
    call Parallel_Exchange_interp(P_inte)
    call Parallel_Exchange_interp(T_inte)
    call Parallel_Exchange_interp(X_grid_interp)
    call Parallel_Exchange_interp(Y_grid_interp)
    call Parallel_Exchange_interp(Z_grid_interp)
    call Parallel_Exchange_interp(Hdst_interp)
    call Parallel_Exchange(Hdst)
    
    call output_LNSResults
     
     !if(MyId == 0) then
     !write(*,"(A)")"=========================================================================================" 
     !endif


  end subroutine OutputLNS_Results

!  subroutine InterpToSteadyGrid(h_unsteady, phi_unsteady, h_steady, phi_interp)
!    !--- 输入输出参数声明 ---
!    use SF_Constant,  only: ik,rk,OverLap
!    use SF_CFD_Global,only: nx_local,nz_local,Ny
!    implicit none
!    real( kind = rk ), intent(in)  :: h_unsteady  (1-OverLap:nx_local+OverLap, 1:ny, 1-OverLap:nz_local+OverLap)   ! 非定常网格高度(3D)
!    real( kind = rk ), intent(in)  :: phi_unsteady(1-OverLap:nx_local+OverLap, 1:ny, 1-OverLap:nz_local+OverLap)  ! 非定常物理量(3D)
!    real( kind = rk ), intent(in)  :: h_steady    (1-OverLap:nx_local+OverLap, 1:2*ny-1, 1-OverLap:nz_local+OverLap)     ! 定常网格高度(3D)
!    real( kind = rk ), intent(out) :: phi_interp  (1-OverLap:nx_local+OverLap, 1:2*ny-1, 1-OverLap:nz_local+OverLap)  ! 输出物理量(3D)
!    
!    
!    !--- 局部变量 ---
!    integer( kind = ik ) :: ic, jc, kc
!    
!    !--- 三维网格遍历 ---
!        do kc = 1-OverLap, nz_local+OverLap
!            do ic = 1-OverLap, nx_local+OverLap
!                
!                !call LinearInterp1D(h_unsteady(ic,:,kc), phi_unsteady(ic,:,kc), h_steady(ic,:,kc), phi_interp(ic,:,kc))
!                call CubicSplineInterp1D(h_unsteady(ic,:,kc), phi_unsteady(ic,:,kc), h_steady(ic,:,kc), phi_interp(ic,:,kc))
!            end do
!        end do
!
!    end subroutine InterpToSteadyGrid
!    
!    subroutine LinearInterp1D(x, y, xi, yi) 
!    !--- 输入输出参数声明 --- 
!    use SF_Constant,  only: ik,rk 
!    use SF_CFD_Global,only: ny 
!    
!    implicit none 
!    real(kind = rk), intent(in) :: x(ny)    ! 原始网格点 
!    real(kind = rk), intent(in) :: y(ny)    ! 原始物理量 
!    real(kind = rk), intent(in) :: xi(ny)   ! 插值网格点 
!    real(kind = rk), intent(out) :: yi(ny)  ! 插值后的物理量 
!    
!  
!    !--- 局部变量 --- 
!    integer(kind = ik) :: i, j 
!    real(kind = rk) :: slope 
! 
!    !--- 遍历插值点 --- 
!    do i = 1, ny 
!        ! 寻找插值点所在的区间 
!        do j = 1, ny - 1 
!            if (xi(i) >= x(j) .and. xi(i) <= x(j+1)) then 
!                ! 计算线性插值的斜率 
!                slope = (y(j+1) - y(j)) / (x(j+1) - x(j)) 
!                ! 进行线性插值 
!                yi(i) = y(j) + slope * (xi(i) - x(j)) 
!                exit 
!            end if 
!        end do 
!    end do 
!    end subroutine LinearInterp1D 
!    
!    ! 导入必要的模块 
!subroutine CubicSplineInterp1D(x, y, xi, yi) 
!    !--- 输入输出参数声明 --- 
!    use SF_Constant,  only: ik,rk 
!    use SF_CFD_Global,only: ny 
!    implicit none 
!    real(kind = rk), intent(in) :: x(ny)    ! 原始网格点 
!    real(kind = rk), intent(in) :: y(ny)    ! 原始物理量 
!    real(kind = rk), intent(in) :: xi(2*ny-1)   ! 插值网格点 
!    real(kind = rk), intent(out) :: yi(2*ny-1)  ! 插值后的物理量 
! 
!    !--- 局部变量 --- 
!    integer(kind = ik) :: i, j 
!    real(kind = rk) :: h(ny-1), alpha(ny-1), l(ny), mu(ny), z(ny) 
!    real(kind = rk) :: a(ny), b(ny), c(ny), d(ny) 
!    real(kind = rk) :: dx, dy 
! 
!    ! 计算步长 
!    do i = 1, ny-1 
!        h(i) = x(i+1) - x(i) 
!    end do 
! 
!    ! 计算 alpha 
!    do i = 2, ny-1 
!        alpha(i) = 3.0_rk / h(i) * (y(i+1) - y(i)) - 3.0_rk / h(i-1) * (y(i) - y(i-1)) 
!    end do 
! 
!    ! 求解三对角线性方程组 
!    l(1) = 1.0_rk 
!    mu(1) = 0.0_rk 
!    z(1) = 0.0_rk 
!    do i = 2, ny-1 
!        l(i) = 2.0_rk * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1) 
!        mu(i) = h(i) / l(i) 
!        z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i) 
!    end do 
!    l(ny) = 1.0_rk 
!    z(ny) = 0.0_rk 
!    c(ny) = 0.0_rk 
!    do j = ny-1, 1, -1 
!        c(j) = z(j) - mu(j) * c(j+1) 
!        b(j) = (y(j+1) - y(j)) / h(j) - h(j) * (c(j+1) + 2.0_rk * c(j)) / 3.0_rk 
!        d(j) = (c(j+1) - c(j)) / (3.0_rk * h(j)) 
!        a(j) = y(j) 
!    end do 
! 
!    ! 遍历插值点 
!    do i = 1, 2*ny-1 
!        ! 寻找插值点所在的区间 
!        do j = 1, ny - 1 
!            if (xi(i) >= x(j) .and. xi(i) <= x(j+1)) then 
!                dx = xi(i) - x(j) 
!                ! 进行三次样条插值 
!                yi(i) = a(j) + b(j) * dx + c(j) * dx**2 + d(j) * dx**3 
!                exit 
!            end if 
!        end do 
!    end do 
!end subroutine CubicSplineInterp1D 

SUBROUTINE Calculate_freestream(Iter)
    IMPLICIT NONE
    integer( kind = ik ), intent(in) :: Iter
    integer( kind = ik ):: i,j,k
    real( kind = rk ), dimension(:,:,:),   allocatable :: Rho_per,U_per,V_per,W_per,P_per
    real( kind = rk ), dimension(:,:,:),   allocatable :: i_t,i_t_tau

    allocate(Rho_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(U_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(V_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(W_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(P_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(i_t(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    allocate(i_t_tau(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
    
    
   !CALL MPI_Barrier(MPI_COMM_WORLD,ierr); 
    
     j = ny
      do k = 1-overlap,nz_local+overlap
        do i = 1-overlap,nx_local+overlap

    if (pert_type .eq. 1)then !fast 
          i_t(i,j,k) = cos(k_infty*((x_grid(i,j,k)-1.d0) - C_inf * LNS_dt * Iter));
      i_t_tau(i,j,k) = -k_infty*(shockxtau(i,k)-c_inf)*sin(k_infty*((x_grid(i,j,k)-1.d0) - C_inf * LNS_dt * Iter));
       
          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);  
            u_per(i,j,k)    = epsilon * i_t(i,j,k);
            v_per(i,j,k)    = 0.d0;
            w_per(i,j,k)    = 0.d0;
            p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
           
        rho_pert_tau(i,j,k) = epsilon * mach_ref * i_t_tau(i,j,k);
          u_pert_tau(i,j,k) = epsilon * i_t_tau(i,j,k);
          v_pert_tau(i,j,k) = 0.d0;
          w_pert_tau(i,j,k) = 0.d0;
          p_pert_tau(i,j,k) = epsilon * mach_ref * gamma * i_t_tau(i,j,k)/ ( gamma * mach_ref * mach_ref );
            
           rho_free(i,j,k)  = rho_inf + rho_per(i,j,k);
             u_free(i,j,k)  = u_inf   + u_per(i,j,k);
             v_free(i,j,k)  = v_inf   + v_per(i,j,k);
             w_free(i,j,k)  = w_inf   + w_per(i,j,k);
             p_free(i,j,k)  = p_inf   + p_per(i,j,k);
             t_free(i,j,k)  = gamma * mach_ref * mach_ref * p_free(i,j,k) / rho_free(i,j,k);
          elseif (pert_type .eq. 2)then !slow
          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
            u_per(i,j,k)    = - epsilon * i_t(i,j,k);
            v_per(i,j,k)    = 0.d0;
            w_per(i,j,k)    = 0.d0;            
            p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
          elseif (pert_type .eq. 3)then
          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
            u_per(i,j,k)    = 0.d0;
            v_per(i,j,k)    = 0.d0;
            w_per(i,j,k)    = 0.d0;   
            p_per(i,j,k)    = 0.d0;
          !if (pert_type .eq. 4)then
          else
          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
            u_per(i,j,k)    = 0.d0;
            v_per(i,j,k)    = epsilon * i_t(i,j,k);
            w_per(i,j,k)    = 0.d0;            
            p_per(i,j,k)    = 0.d0;                
          endif    
        enddo
      enddo  

call Parallel_Exchange(Rho_free)     
call Parallel_Exchange(U_free)     
call Parallel_Exchange(V_free)     
call Parallel_Exchange(W_free)     
call Parallel_Exchange(P_free)     
call Parallel_Exchange(T_free)     
call Parallel_Exchange(rho_pert_tau)     
call Parallel_Exchange(u_pert_tau)     
call Parallel_Exchange(v_pert_tau)     
call Parallel_Exchange(w_pert_tau)     
call Parallel_Exchange(p_pert_tau)     

      
DEALLOCATE(i_t,i_t_tau)
DEALLOCATE(Rho_per,U_per,V_per,W_per,P_per)
      
END SUBROUTINE Calculate_freestream  
  
  
END Module SF_LNS