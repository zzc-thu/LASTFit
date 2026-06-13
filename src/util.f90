!===============================
!Transform relationship
!===============================
! This subroutine is used to calculate the 
! primitive variables from Consevative variables
SUBROUTINE CVtoU_CPG(Rho,U,V,W,P,CV)
    ! Transform Consevative variables to primitive
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL( kind = rk ),INTENT(OUT)::Rho,U,V,W,P
    REAL( kind = rk ),INTENT(IN)::CV(1:5)
  
    Rho=CV(1)
    U=CV(2)/CV(1)
    V=CV(3)/CV(1)
    W=CV(4)/CV(1)
    P=(GAMMA-1.D0)*(CV(5)-0.5D0*Rho*(U*U+V*V+W*W))

    if(P < 0.d0) then
        write(*,*)'Pressure is less than zero'
        stop
    endif

END SUBROUTINE CVtoU_CPG

! This subroutine is used to transform the 
! primitive variables to Consevative variables
SUBROUTINE UtoCV_CPG(Rho,U,V,W,P,CV)
    ! Transform primitive variables to Consevative variables
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL( kind = rk ),INTENT(IN)::Rho,U,V,W,P
    REAL( kind = rk ),INTENT(OUT)::CV(1:5)

    CV(1)=Rho
    CV(2)=Rho*U
    CV(3)=Rho*V
    CV(4)=Rho*W
    CV(5)=P/(GAMMA-1.D0)+0.5D0*Rho*(U*U+V*V+W*W)

END SUBROUTINE UtoCV_CPG

! This subroutine is used to calculate the flux in 2D
SUBROUTINE UtoFlux_CPG2D(Rho,U,V,W,P,F,G)
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL( kind = rk ),INTENT(IN)::Rho,U,V,W,P
    REAL( kind = rk ),INTENT(OUT)::F(1:5),G(1:5)

    F(1) = Rho * U
    F(2) = Rho * U * U + P
    F(3) = Rho * U * V
    F(4) = Rho * U * W
    F(5) = (P/(GAMMA - 1.D0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * U

    G(1) = Rho * V
    G(2) = Rho * V * U
    G(3) = Rho * V * V + P
    G(4) = Rho * V * W
    G(5) = (P/(GAMMA - 1.D0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * V

END SUBROUTINE UtoFlux_CPG2D

! This subroutine is used to calculate the flux in 3D
SUBROUTINE UtoFlux_CPG3D(Rho,U,V,W,P,F,G,H)
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL( kind = rk ),INTENT(IN)::Rho,U,V,W,P
    REAL( kind = rk ),INTENT(OUT)::F(1:5),G(1:5),H(1:5)

    F(1) = Rho * U
    F(2) = Rho * U * U + P
    F(3) = Rho * U * V
    F(4) = Rho * U * W
    F(5) = (P/(GAMMA - 1.D0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * U

    G(1) = Rho * V
    G(2) = Rho * V * U
    G(3) = Rho * V * V + P
    G(4) = Rho * V * W
    G(5) = (P/(GAMMA - 1.D0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * V

    H(1) = Rho * W
    H(2) = Rho * W * U
    H(3) = Rho * W * V
    H(4) = Rho * W * W + P
    H(5) = (P/(GAMMA - 1.D0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * W

    END SUBROUTINE UtoFlux_CPG3D
    
    
    SUBROUTINE CVtoU_pert(Rho0,U0,V0,W0,P0,CV1,Rho_dist,U_dist,V_dist,W_dist,P_dist)
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    ! CV1 is the conservative variables of perturbations
    REAL( kind = rk ),INTENT(IN) ::Rho0,U0,V0,W0,P0,CV1(1:5)
    REAL( kind = rk ),INTENT(OUT)::Rho_dist,U_dist,V_dist,W_dist,P_dist

    Rho_dist=CV1(1)
    U_dist=(CV1(2)-CV1(1)*U0)/Rho0
    V_dist=(CV1(3)-CV1(1)*V0)/Rho0
    W_dist=(CV1(4)-CV1(1)*W0)/Rho0
    P_dist=(CV1(5)-0.5D0*Rho_dist*(U0*U0+V0*V0+W0*W0)-Rho0*(U0*U_dist+V0*V_dist+W0*W_dist))*(GAMMA-1.D0)
    
    END SUBROUTINE CVtoU_pert
    
    SUBROUTINE UtoCV_pert(Rho0,U0,V0,W0,P0,Rho_dist,U_dist,V_dist,W_dist,P_dist,CV1)
    ! Transform primitive variables to Consevative variables
    USE SF_Constant,  only:rk
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL( kind = rk ),INTENT(IN) ::Rho0,U0,V0,W0,P0,Rho_dist,U_dist,V_dist,W_dist,P_dist
    REAL( kind = rk ),INTENT(OUT)::CV1(1:5)
    
    ! CV is the conservative variables of perturbations
    CV1(1)=Rho_dist
    CV1(2)=Rho_dist*U0+U_dist*Rho0
    CV1(3)=Rho_dist*V0+V_dist*Rho0
    CV1(4)=Rho_dist*W0+W_dist*Rho0
    CV1(5)=P_dist/(GAMMA-1.D0)+0.5D0*Rho_dist*(U0*U0+V0*V0+W0*W0)+Rho0*(U0*U_dist+V0*V_dist+W0*W_dist)

    END SUBROUTINE UtoCV_pert
    
    
    SUBROUTINE UtoCV_tau(Rho,U,V,W,Rho_tau,U_tau,V_tau,W_tau,P_tau,CV)
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL(8),INTENT(IN)::Rho,U,V,W,Rho_tau,U_tau,V_tau,W_tau,P_tau
    REAL(8),INTENT(OUT)::CV(5)
    
     CV(1)=Rho_tau
     CV(2)=Rho * U_tau + U * Rho_tau
     CV(3)=Rho * V_tau + V * Rho_tau
     CV(4)=Rho * W_tau + W * Rho_tau
     CV(5)=P_tau/(GAMMA-1.D0)+0.5D0*Rho_tau*(U*U+V*V+W*W)+Rho*(U_tau*U+V_tau*V+W_tau*W)
     
    END SUBROUTINE UtoCV_tau
    
    SUBROUTINE UtoCV_taud(Rho,U,V,W,Rho_tau,U_tau,V_tau,W_tau,P_tau,Rho_taud,U_taud,V_taud,W_taud,P_taud,CV)
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL(8),INTENT(IN)::Rho,U,V,W,Rho_tau,U_tau,V_tau,W_tau,P_tau,Rho_taud,U_taud,V_taud,W_taud,P_taud
    REAL(8),INTENT(OUT)::CV(5)
    
     CV(1)=Rho_taud
     CV(2)=Rho_tau * U_tau + Rho * U_taud + U_tau * Rho_tau + U * Rho_taud
     CV(3)=Rho_tau * V_tau + Rho * V_taud + V_tau * Rho_tau + V * Rho_taud
     CV(4)=Rho_tau * W_tau + Rho * W_taud + W_tau * Rho_tau + W * Rho_taud
     CV(5)=P_taud/(GAMMA-1.D0)+ Rho_tau*(U*U_tau+V*V_tau+W*W_tau)+Rho*(U_tau*U_tau+V_tau*V_tau+W_tau*W_tau+&
         U*U_taud+V*V_taud+W*W_taud)+0.5d0*Rho_taud*(U*U+V*V+W*W)+Rho_tau*(U*U_tau+V*V_tau+W*W_tau)
     
    END SUBROUTINE UtoCV_taud
    
    SUBROUTINE UtoFlux_tau(Rho,U,V,W,P,Rho_tau,U_tau,V_tau,W_tau,P_tau,F,G,H)
    USE SF_CFD_Global,only:GAMMA
    IMPLICIT NONE
    REAL(8),INTENT(IN)::Rho,U,V,W,P,Rho_tau,U_tau,V_tau,W_tau,P_tau
    REAL(8),INTENT(OUT)::F(5),G(5),H(5)
  
     F(1) = Rho * U_tau + U * Rho_tau
     F(2) = Rho_tau * U * U + Rho * 2 * U * U_tau + P_tau
     F(3) = Rho_tau * U * V + Rho * U_tau * V + Rho * U * V_tau
     F(4) = Rho_tau * U * W + Rho * U_tau * W + Rho * U * W_tau
     F(5) = (P_tau/(GAMMA - 1.0) + 0.5d0 * Rho_tau * (U*U + V*V + W*W) + Rho *(U_tau*U+V_tau*V+W_tau*W) + P_tau) * U &
         &+ (P/(GAMMA - 1.0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * U_tau
     
     G(1) = Rho * V_tau + V * Rho_tau
     G(2) = Rho_tau * V * U + Rho * V_tau * U + Rho * V * U_tau
     G(3) = Rho_tau * V * V + Rho * 2 * V * V_tau + P_tau
     G(4) = Rho_tau * V * W + Rho * V_tau * W + Rho * V * W_tau
     G(5) = (P_tau/(GAMMA - 1.0) + 0.5d0 * Rho_tau * (U*U + V*V + W*W) + Rho *(U_tau*U+V_tau*V+W_tau*W) + P_tau) * V &
         &+ (P/(GAMMA - 1.0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * V_tau
  
     H(1) = Rho * W_tau + W * Rho_tau
     H(2) = Rho_tau * W * U + Rho * W_tau * U + Rho * W * U_tau
     H(3) = Rho_tau * W * V + Rho * W_tau * V + Rho * W * V_tau
     H(4) = Rho_tau * W * W + Rho * 2 * W * W_tau + P_tau
     H(5) = (P_tau/(GAMMA - 1.0) + 0.5d0 * Rho_tau * (U*U + V*V + W*W) + Rho *(U_tau*U+V_tau*V+W_tau*W) + P_tau) * W &
         &+ (P/(GAMMA - 1.0) + 0.5d0 * Rho * (U*U + V*V + W*W) + P) * W_tau
     
    END SUBROUTINE UtoFlux_tau
    
SUBROUTINE UtoFlux_taud(Rho,U,V,W,P,Rho_tau,U_tau,V_tau,W_tau,P_tau,Rho_taud,U_taud,V_taud,W_taud,P_taud,F,G,H)          
USE SF_CFD_Global,only:GAMMA          
IMPLICIT NONE          
REAL(8),INTENT(IN)::Rho,U,V,W,P,Rho_tau,U_tau,V_tau,W_tau,P_tau,Rho_taud,U_taud,V_taud,W_taud,P_taud          
REAL(8),INTENT(OUT)::F(5),G(5),H(5)               
REAL(8) :: A, dA_dx, B, dB_dx      
REAL(8) :: A_G, dA_G_dx, B_G, dB_G_dx     
REAL(8) :: A_H, dA_H_dx, B_H, dB_H_dx       
F(1) = Rho_taud*U_tau + Rho_tau*U_taud + Rho*U_taud + U*Rho_taud      
F(2) = Rho_taud*U**2 + 4.0d0*Rho_tau*U*U_tau + 2.0d0*Rho*U_tau**2 + 2.0d0*Rho*U*U_taud + P_taud      
F(3) = Rho_taud*U*V + 2.0d0*Rho_tau*U_tau*V + 2.0d0*Rho_tau*U*V_tau + Rho*U_taud*V + 2.0d0*Rho*U_tau*V_tau + Rho*U*V_taud      
F(4) = Rho_taud*U*W + 2.0d0*Rho_tau*U_tau*W + 2.0d0*Rho_tau*U*W_tau + Rho*U_taud*W + 2.0d0*Rho*U_tau*W_tau + Rho*U*W_taud        
A = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_tau*(U**2 + V**2 + W**2) + Rho*(U_tau*U + V_tau*V + W_tau*W) + P_tau)     
dA_dx = (P_taud/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + &
         & Rho_tau*(U_tau*U + V_tau*V + W_tau*W) + Rho*(U_taud*U + U_tau**2 + V_taud*V + V_tau**2 + W_taud*W + W_tau**2) + P_taud)     
B = (P/(GAMMA-1.0d0) + 0.5d0*Rho*(U**2 + V**2 + W**2) + P)     
dB_dx = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + P_tau)     
F(5) = dA_dx*U + A*U_tau + dB_dx*U_tau + B*U_taud        

G(1) = Rho_taud*V_tau + Rho_tau*V_taud + Rho*V_taud + V*Rho_taud      
G(2) = Rho_taud*V*U + 2.0d0*Rho_tau*V_tau*U + 2.0d0*Rho_tau*V*U_tau + Rho*V_taud*U + 2.0d0*Rho*V_tau*U_tau + Rho*V*U_taud     
G(3) = Rho_taud*V**2 + 4.0d0*Rho_tau*V*V_tau + 2.0d0*Rho*V_tau**2 + 2.0d0*Rho*V*V_taud + P_taud      
G(4) = Rho_taud*V*W + 2.0d0*Rho_tau*V_tau*W + 2.0d0*Rho_tau*V*W_tau + Rho*V_taud*W + 2.0d0*Rho*V_tau*W_tau + Rho*V*W_taud              
A_G = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_tau*(U**2 + V**2 + W**2) + Rho*(U_tau*U + V_tau*V + W_tau*W) + P_tau)     
dA_G_dx = (P_taud/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + &
        & Rho_tau*(U_tau*U + V_tau*V + W_tau*W) + Rho*(U_taud*U + U_tau**2 + V_taud*V + V_tau**2 + W_taud*W + W_tau**2) + P_taud)
B_G = (P/(GAMMA-1.0d0) + 0.5d0*Rho*(U**2 + V**2 + W**2) + P)     
dB_G_dx = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + P_tau)   
G(5) = dA_G_dx*V + A_G*V_tau + dB_G_dx*V_tau + B_G*V_taud 

H(1) = Rho_taud*W_tau + Rho_tau*W_taud + Rho*W_taud + W*Rho_taud     
H(2) = Rho_taud*W*U + 2.0d0*Rho_tau*W_tau*U + 2.0d0*Rho_tau*W*U_tau + Rho*W_taud*U + 2.0d0*Rho*W_tau*U_tau + Rho*W*U_taud      
H(3) = Rho_taud*W*V + 2.0d0*Rho_tau*W_tau*V + 2.0d0*Rho_tau*W*V_tau + Rho*W_taud*V + 2.0d0*Rho*W_tau*V_tau + Rho*W*V_taud      
H(4) = Rho_taud*W**2 + 4.0d0*Rho_tau*W*W_tau + 2.0d0*Rho*W_tau**2 + 2.0d0*Rho*W*W_taud + P_taud        
A_H = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_tau*(U**2 + V**2 + W**2) + Rho*(U_tau*U + V_tau*V + W_tau*W) + P_tau)     
dA_H_dx = (P_taud/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + &
        & Rho_tau*(U_tau*U + V_tau*V + W_tau*W) + Rho*(U_taud*U + U_tau**2 + V_taud*V + V_tau**2 + W_taud*W + W_tau**2) + P_taud)     
B_H = (P/(GAMMA-1.0d0) + 0.5d0*Rho*(U**2 + V**2 + W**2) + P)     
dB_H_dx = (P_tau/(GAMMA-1.0d0) + 0.5d0*Rho_taud*(U**2 + V**2 + W**2) + 0.5d0*Rho_tau*(2*U*U_tau + 2*V*V_tau + 2*W*W_tau) + P_tau)     
H(5) = dA_H_dx*W + A_H*W_tau + dB_H_dx*W_tau + B_H*W_taud     

END SUBROUTINE UtoFlux_taud
    
