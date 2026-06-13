SUBROUTINE Update_Variables
  USE SF_Constant,    ONLY: ik, rk, NumVar,fifth_order,first_order,Unsteady_NS_Analysis
  USE SF_CFD_Global,  ONLY: nz_local,Ny,nx_local,Gamma,Mach_Ref,TWall,Select_Scheme,ModelType
  USE MPI_GLOBAL,     ONLY: npx,npx0
  
  ! Arrays and matrix
  USE SF_CFD_Global,  ONLY: U,V,W,P,Rho,T,Ucons0,&
                            ShockH0,ShockH,dShockH,&
                            ShockV0,ShockV,dShockV,dUcons,&
                            dxidx,dxidy,dxidz,dxidt,&
                            detadx,detady,detadz,detadt,&
                            dzetadx,dzetady,dzetadz,dzetadt,&
                            BigU_Xi,BigU_Eta,BigU_Zeta,&
                            AnalysisType,BigU_Xi_total,BigU_Eta_total,BigU_Zeta_total,&
                            Rho_total,U_total,V_total,W_total,P_total,T_total,shockH_total,shockV_total,Ucons0_total,Cs
  ! subroutines and functions
  USE SFitting,       ONLY: ShockRelation3D,ShockRelation3D_unsteady
  USE MPI_GLOBAL,     ONLY: Parallel_Exchange,Parallel_Exchange_NumVar,Parallel_Exchange_Surface
  implicit none
  ! Local variables
  integer( kind = ik ):: ic,jc,kc
  real( kind = rk ):: l_xi_abs,l_xi1,l_xi2,l_xi3,v_cdot_I,Mach_local
  ! The calculation of the variables is based on the conservative variables
  ! We need to update the primitive variables based on the conservative variables
  DO kc = 1,nz_local
    DO jc = 1,Ny
      DO ic = 1,nx_local
        ! Update the variables
        Call CVtoU_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
            T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc);
      ENDDO
    ENDDO
  ENDDO
  
  ! Boundary Conditions
   ! For Shock Boundary Conditions
    if(AnalysisType == Unsteady_NS_Analysis) then
        CALL ShockRelation3D_unsteady
    else
        CALL ShockRelation3D ! Here, Rho,U,V,W,P and T are updated based on the shock relation
        
    endif
    

    jc = Ny 
     DO kc = 1, nz_local
      DO ic = 1, nx_local
       CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
      ENDDO
     ENDDO

  select case(Select_Scheme)
  case(fifth_order)

   ! 这里的外插出口条件似乎是没有必要的，因为不管是粘性项还是对流项，在我们的边界处都是用
   ! 单侧差分来处理的，这意味着在ic = 0和ic = nx_local位置上的值都是外推的到的，所以不
   ! 需要再次外推（这是理论上的判断）。
   ! For outflow conditions
     !if( npx == npx0 - 1) then
     !  ic = nx_local
     !  do kc = 1, nz_local
     !    do jc = 2, Ny-1
     !        
     ! !l_xi_abs = 1.0_rk/sqrt(dxidx(ic,jc,kc)*dxidx(ic,jc,kc) + &
     ! !                      &dxidy(ic,jc,kc)*dxidy(ic,jc,kc) + &
     ! !                      &dxidz(ic,jc,kc)*dxidz(ic,jc,kc))
     ! !
     ! !l_xi1 = l_xi_abs*dxidx(ic,jc,kc)
     ! !l_xi2 = l_xi_abs*dxidy(ic,jc,kc)
     ! !l_xi3 = l_xi_abs*dxidz(ic,jc,kc)
     ! !
     ! !v_cdot_I = U(ic,jc,kc)*l_xi1+ V(ic,jc,kc)*l_xi2 + W(ic,jc,kc)*l_xi3
     ! ! 
     ! !Mach_local = v_cdot_I / Cs(ic,jc,kc)
     ! 
     !     !if (Mach_local >= 1.0_rk )then
     !            
     !          Rho(ic,jc,kc) = (4.d0 *Rho(ic-1,jc,kc) - 6.d0 *Rho(ic-2,jc,kc) &
     !                        &+ 4.d0 *Rho(ic-3,jc,kc) - 1.d0 *Rho(ic-4,jc,kc))
     !           U(ic,jc,kc) = (4.d0 *U(ic-1,jc,kc) - 6.d0 *U(ic-2,jc,kc) &
     !                        &+ 4.d0 *U(ic-3,jc,kc) - 1.d0 *U(ic-4,jc,kc))
     !           V(ic,jc,kc) = (4.d0 *V(ic-1,jc,kc) - 6.d0 *V(ic-2,jc,kc) &
     !                         &+ 4.d0 *V(ic-3,jc,kc) - 1.d0 *V(ic-4,jc,kc))
     !           W(ic,jc,kc) = (4.d0 *W(ic-1,jc,kc) - 6.d0 *W(ic-2,jc,kc) &
     !                         &+ 4.d0 *W(ic-3,jc,kc) - 1.d0 *W(ic-4,jc,kc))
     !           P(ic,jc,kc) = (4.d0 *P(ic-1,jc,kc) - 6.d0 *P(ic-2,jc,kc) &
     !                         &+ 4.d0 *P(ic-3,jc,kc) - 1.d0 *P(ic-4,jc,kc))
     !           T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc);
     !           CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),&
     !                           &W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
     !           
     !      !endif
     !    enddo
     !  enddo
     !endif
     !
     !!
     if (ModelType == 1)then
       if( npx == 0 ) then
       ic = 1
       do kc = 1, nz_local
         do jc = 2, Ny-1
           Rho(ic,jc,kc) = (4.d0 *Rho(ic+1,jc,kc) - 6.d0 *Rho(ic+2,jc,kc) &
                         &+ 4.d0 *Rho(ic+3,jc,kc) - 1.d0 *Rho(ic+4,jc,kc))
            U(ic,jc,kc) = ( 4.d0 *U(ic+1,jc,kc) - 6.d0 *U(ic+2,jc,kc) &
                         &+ 4.d0 *U(ic+3,jc,kc) - 1.d0 *U(ic+4,jc,kc))
            V(ic,jc,kc) = ( 4.d0 *V(ic+1,jc,kc) - 6.d0 *V(ic+2,jc,kc) &
                         &+ 4.d0 *V(ic+3,jc,kc) - 1.d0 *V(ic+4,jc,kc))
            W(ic,jc,kc) = ( 4.d0 *W(ic+1,jc,kc) - 6.d0 *W(ic+2,jc,kc) &
                         &+ 4.d0 *W(ic+3,jc,kc) - 1.d0 *W(ic+4,jc,kc))
            P(ic,jc,kc) = ( 4.d0 *P(ic+1,jc,kc) - 6.d0 *P(ic+2,jc,kc) &
                         &+ 4.d0 *P(ic+3,jc,kc) - 1.d0 *P(ic+4,jc,kc))
            T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc);
            CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),&
                            &W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
         enddo
       enddo       
     endif
     endif
     !
    
   !================================================================================================
   ! For Wall Surface Conditions
    jc = 1
     if(TWall < 0.d0) then ! Adiabatic Wall, dT/dN = 0;
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        T(ic,jc,kc) = 1.92D0 * T(ic,jc+1,kc) - 1.44D0 * T(ic,jc+2,kc) &
                   &+ 0.64D0 * T(ic,jc+3,kc) - 0.12D0 * T(ic,jc+4,kc);
       ENDDO
      ENDDO
     else   ! Isothermal Wall, T = TWall;
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        T(ic,jc,kc) = TWall;
       ENDDO
      ENDDO
     endif
     
     
     
     ! Density, Rho is calculated based on the continuity equation
     ! Other Variables for U, V, W;
     DO kc = 1, nz_local
      DO ic = 1, nx_local
       U(ic,jc,kc) = 0.d0;
       V(ic,jc,kc) = 0.d0;
       W(ic,jc,kc) = 0.d0;
       ! Here, the pressure is calculated based on dP/dN = 0
       P(ic,jc,kc) =  1.92D0 * P(ic,jc+1,kc) - 1.44D0 * P(ic,jc+2,kc) &
                   &+ 0.64D0 * P(ic,jc+3,kc) - 0.12D0 * P(ic,jc+4,kc);
       Rho(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / T(ic,jc,kc);
       ! Here, the pressure is calculated based on the continuity equation
       !P(ic,jc,kc) = Rho(ic,jc,kc) * T(ic,jc,kc) / (Gamma * Mach_Ref * Mach_Ref); 
       CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
      ENDDO
     ENDDO
    !================================================================================================    
  case(first_order)
   ! For outflow conditions
     if( npx == npx0 - 1) then
       ic = nx_local
       do kc = 1, nz_local
         do jc = 2, Ny-1
           Rho(ic,jc,kc)=Rho(ic-1,jc,kc)
            U(ic,jc,kc) =U(ic-1,jc,kc)
            V(ic,jc,kc) =V(ic-1,jc,kc)
            W(ic,jc,kc) =W(ic-1,jc,kc)
            P(ic,jc,kc) =P(ic-1,jc,kc)
            T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc);
            CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),&
                            &W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
         enddo
       enddo
     endif

     if (ModelType == 1)then
     if( npx == 0 ) then
       ic = 1
       do kc = 1, nz_local
         do jc = 2, Ny-1
           Rho(ic,jc,kc)=Rho(ic+1,jc,kc) 
            U(ic,jc,kc) =U(ic+1,jc,kc) 
            V(ic,jc,kc) =V(ic+1,jc,kc) 
            W(ic,jc,kc) =W(ic+1,jc,kc) 
            P(ic,jc,kc) =P(ic+1,jc,kc)
            T(ic,jc,kc) = Gamma * Mach_Ref * Mach_Ref * P(ic,jc,kc) / Rho(ic,jc,kc);
            CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),&
                            &W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
         enddo
       enddo       
     endif
     endif
   !================================================================================================
   ! For Wall Surface Conditions
    jc = 1
     if(TWall < 0.d0) then ! Adiabatic Wall
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        T(ic,jc,kc) = T(ic,jc+1,kc);
       ENDDO
      ENDDO
     else   ! Isothermal Wall
      DO kc = 1, nz_local
       DO ic = 1, nx_local
        T(ic,jc,kc) = TWall;
       ENDDO
      ENDDO
     endif
     ! Other Variables for U, V, W;
     DO kc = 1, nz_local
      DO ic = 1, nx_local
       U(ic,jc,kc) = 0.d0;
       V(ic,jc,kc) = 0.d0;
       W(ic,jc,kc) = 0.d0;
       ! Here, the pressure is calculated based on dP/dN = 0
        P(ic,jc,kc) =  P(ic,jc+1,kc);
        Rho(ic,jc,kc) = P(ic,jc,kc) * Gamma * Mach_Ref * Mach_Ref / T(ic,jc,kc);
       ! Here, the pressure is calculated based on the continuity equation
       !P(ic,jc,kc) = Rho(ic,jc,kc) * T(ic,jc,kc) / (Gamma * Mach_Ref * Mach_Ref); 
       CALL UtoCV_CPG(Rho(ic,jc,kc),U(ic,jc,kc),V(ic,jc,kc),W(ic,jc,kc),P(ic,jc,kc),Ucons0(1:NumVar,ic,jc,kc))
      ENDDO
     ENDDO      
    
  end select
    !================================================================================================
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
        BigU_Xi(ic,jc,kc) = U(ic,jc,kc)*  dxidx(ic,jc,kc) &
                         &+ V(ic,jc,kc)*  dxidy(ic,jc,kc) &
                         &+ W(ic,jc,kc)*  dxidz(ic,jc,kc) + dxidt(ic,jc,kc)

       BigU_Eta(ic,jc,kc) = U(ic,jc,kc)* detadx(ic,jc,kc) &
                         &+ V(ic,jc,kc)* detady(ic,jc,kc) &
                         &+ W(ic,jc,kc)* detadz(ic,jc,kc) + detadt(ic,jc,kc)

      BigU_Zeta(ic,jc,kc) = U(ic,jc,kc)*dzetadx(ic,jc,kc) &
                         &+ V(ic,jc,kc)*dzetady(ic,jc,kc) &
                         &+ W(ic,jc,kc)*dzetadz(ic,jc,kc) + dzetadt(ic,jc,kc)    
        ENDDO
      ENDDO
    ENDDO
    
    
    ! Parallel Transfer
if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then 
    call Parallel_Exchange_singularity_negative(BigU_Xi,BigU_Xi_total)
    call Parallel_Exchange_singularity(BigU_Eta,BigU_Eta_total)
    call Parallel_Exchange_singularity(BigU_Zeta,BigU_Zeta_total)
endif
    
    call Parallel_Exchange(BigU_Xi)
    call Parallel_Exchange(BigU_Eta)
    call Parallel_Exchange(BigU_Zeta)

    ! Update flow related variables
if (ModelType == 2 .or. ModelType == 3 .or. ModelType == 4) then
    CALL Parallel_Exchange_singularity(Rho,Rho_total)
    CALL Parallel_Exchange_singularity(U,U_total)
    CALL Parallel_Exchange_singularity(V,V_total)
    CALL Parallel_Exchange_singularity(W,W_total)
    CALL Parallel_Exchange_singularity(P,P_total)
    CALL Parallel_Exchange_singularity(T,T_total)
    CALL Parallel_Exchange_singularity_NumVar(Ucons0,Ucons0_total)
    ! Exchange the shock related variables
    CALL Parallel_Exchange_singularity_Surface(ShockH,ShockH_total)
    CALL Parallel_Exchange_singularity_Surface(ShockV,ShockV_total)
endif

    call Parallel_Exchange(Rho)
    call Parallel_Exchange(U)
    call Parallel_Exchange(V)
    call Parallel_Exchange(W)
    call Parallel_Exchange(P)
    call Parallel_Exchange(T)
    call Parallel_Exchange_NumVar(Ucons0)
    call Parallel_Exchange_Surface(ShockH)
    call Parallel_Exchange_Surface(ShockV)      
    
END SUBROUTINE Update_Variables