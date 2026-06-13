!SUBROUTINE Calculate_freestream(Iter)
!    use SF_Constant,   only: ik,rk,overLap
!    use SF_CFD_Global, only: Pert_type,k_infty,x_grid,C_inf,LNS_dt,shockXtau,epsilon,gamma,mach_ref,&
!                           & Rho_inf,U_inf,V_inf,W_inf,P_inf,T_inf,&
!                           & Rho_free,U_free,V_free,W_free,P_free,T_free,&
!                           & rho_pert_tau,u_pert_tau,v_pert_tau,w_pert_tau,p_pert_tau,&
!                           & nx_local,ny,nz_local,overlap
!    use MPI_GLOBAL,     only: MyId,ierr,MPI_COMM_WORLD
!    use MPI_GLOBAL,     only: Parallel_Exchange,MPI_Barrier
!    IMPLICIT NONE
!    integer( kind = ik ), intent(in) :: Iter
!    integer( kind = ik ):: i,j,k
!    real( kind = rk ), dimension(:,:,:),   allocatable :: Rho_per,U_per,V_per,W_per,P_per
!    real( kind = rk ), dimension(:,:,:),   allocatable :: i_t,i_t_tau
!
!    allocate(Rho_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(U_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(V_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(W_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(P_per(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(i_t(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    allocate(i_t_tau(1-overLap:nx_local+overLap,1:Ny,1-overLap:nz_local+overLap))
!    
!    
!   CALL MPI_Barrier(MPI_COMM_WORLD,ierr); 
!    
!     j = ny
!      do k = 1-overlap,nz_local+overlap
!        do i = 1-overlap,nx_local+overlap
!
!    if (pert_type .eq. 1)then !fast 
!          i_t(i,j,k) = cos(k_infty*((x_grid(i,j,k)-1.d0) - C_inf * LNS_dt * Iter));
!      i_t_tau(i,j,k) = -k_infty*(shockxtau(i,k)-c_inf)*sin(k_infty*((x_grid(i,j,k)-1.d0) - C_inf * LNS_dt * Iter));
!       
!          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);  
!            u_per(i,j,k)    = epsilon * i_t(i,j,k);
!            v_per(i,j,k)    = 0.d0;
!            w_per(i,j,k)    = 0.d0;
!            p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
!           
!        rho_pert_tau(i,j,k) = epsilon * mach_ref * i_t_tau(i,j,k);
!          u_pert_tau(i,j,k) = epsilon * i_t_tau(i,j,k);
!          v_pert_tau(i,j,k) = 0.d0;
!          w_pert_tau(i,j,k) = 0.d0;
!          p_pert_tau(i,j,k) = epsilon * mach_ref * gamma * i_t_tau(i,j,k)/ ( gamma * mach_ref * mach_ref );
!            
!           rho_free(i,j,k)  = rho_inf + rho_per(i,j,k);
!             u_free(i,j,k)  = u_inf   + u_per(i,j,k);
!             v_free(i,j,k)  = v_inf   + v_per(i,j,k);
!             w_free(i,j,k)  = w_inf   + w_per(i,j,k);
!             p_free(i,j,k)  = p_inf   + p_per(i,j,k);
!             t_free(i,j,k)  = gamma * mach_ref * mach_ref * p_free(i,j,k) / rho_free(i,j,k);
!          elseif (pert_type .eq. 2)then !slow
!          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
!            u_per(i,j,k)    = - epsilon * i_t(i,j,k);
!            v_per(i,j,k)    = 0.d0;
!            w_per(i,j,k)    = 0.d0;            
!            p_per(i,j,k)    = epsilon * mach_ref * gamma* i_t(i,j,k)/ ( gamma * mach_ref * mach_ref );
!          elseif (pert_type .eq. 3)then
!          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
!            u_per(i,j,k)    = 0.d0;
!            v_per(i,j,k)    = 0.d0;
!            w_per(i,j,k)    = 0.d0;   
!            p_per(i,j,k)    = 0.d0;
!          !if (pert_type .eq. 4)then
!          else
!          rho_per(i,j,k)    = epsilon * mach_ref * i_t(i,j,k);
!            u_per(i,j,k)    = 0.d0;
!            v_per(i,j,k)    = epsilon * i_t(i,j,k);
!            w_per(i,j,k)    = 0.d0;            
!            p_per(i,j,k)    = 0.d0;                
!          endif    
!        enddo
!      enddo  
!
!call Parallel_Exchange(Rho_free)     
!call Parallel_Exchange(U_free)     
!call Parallel_Exchange(V_free)     
!call Parallel_Exchange(W_free)     
!call Parallel_Exchange(P_free)     
!call Parallel_Exchange(T_free)     
!call Parallel_Exchange(rho_pert_tau)     
!call Parallel_Exchange(u_pert_tau)     
!call Parallel_Exchange(v_pert_tau)     
!call Parallel_Exchange(w_pert_tau)     
!call Parallel_Exchange(p_pert_tau)     
!
!      
!DEALLOCATE(i_t,i_t_tau)
!DEALLOCATE(Rho_per,U_per,V_per,W_per,P_per)
!      
!END SUBROUTINE Calculate_freestream