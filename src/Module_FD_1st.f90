Module FD1st_order
  use SF_Constant,  only: ik, rk, NumVar,overLAP
  use SF_CFD_Global,only: nx_local,Ny,nz_local,ModelType
  use MPI_Global,   only: npx,npx0,npz,npz0

  implicit none
  
  contains

  ! This subroutine calculate the 1st order derivative xi of the variables  
  subroutine Cal_Deri_Dxi_1st_up_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic-1,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    if(ModelType == 1)then
    ! BC Schemes
    if( npx == 0 ) then
      ic = 1
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic+1,jc,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    endif
    endif
    
    
    if( npx == npx0 - 1 ) then
      ic = nx_local
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic-1,jc,kc)
          ENDDO
        ENDDO
      ENDDO      
    endif

  end subroutine Cal_Deri_Dxi_1st_up_numvar

  subroutine Cal_Deri_Dxi_2nd_ce_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    
    df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = 0.5_rk*(f(iVar,ic+1,jc,kc) - f(iVar,ic-1,jc,kc))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    if( ModelType == 1)then 
    ! BC Schemes
    if( npx == 0 ) then
      ic = 1
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic+1,jc,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    endif
    endif
    
    if( npx == npx0 - 1 ) then
      ic = nx_local
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic-1,jc,kc)
          ENDDO
        ENDDO
      ENDDO      
    endif

  end subroutine Cal_Deri_Dxi_2nd_ce_numvar
  
  subroutine Cal_Deri_Dxi_1st_do_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic+1,jc,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    if( ModelType == 1)then 
    ! BC Schemes
    if( npx == 0 ) then
      ic = 1
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic+1,jc,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    endif
    endif
    
    
    if( npx == npx0 - 1 ) then
      ic = nx_local
      do kc = 1, nz_local
        do jc = 1, Ny
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic-1,jc,kc)
          ENDDO
        ENDDO
      ENDDO      
    endif

  end subroutine Cal_Deri_Dxi_1st_do_numvar
  
  subroutine Cal_Deri_Dxi_2nd_ce(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc

    df = 0.0_rk
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
            df(ic,jc,kc) = 0.5_rk*(f(ic+1,jc,kc) - f(ic-1,jc,kc))
        ENDDO
      ENDDO
    ENDDO

    if( ModelType == 1)then 
    ! BC Schemes
    if( npx == 0 ) then
      ic = 1
      do kc = 1, nz_local
        do jc = 1, Ny
            df(ic,jc,kc) = f(ic+1,jc,kc) - f(ic,jc,kc)
        ENDDO
      ENDDO
    endif
    endif
    
    
    if( npx == npx0 - 1 ) then
      ic = nx_local
      do kc = 1, nz_local
        do jc = 1, Ny
            df(ic,jc,kc) = f(ic,jc,kc) - f(ic-1,jc,kc)
        ENDDO
      ENDDO      
    endif

  end subroutine Cal_Deri_Dxi_2nd_ce

  ! This subroutine calculate the 1st order derivative eta of the variables
  subroutine Cal_Deri_Deta_1st_up_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar

    df = 0.0_rk
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 2, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic,jc-1,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! BC Schemes
      jc = 1
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc+1,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO

      jc = Ny
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic,jc-1,kc)
          ENDDO
        ENDDO
      ENDDO      

  end subroutine Cal_Deri_Deta_1st_up_numvar

  subroutine Cal_Deri_Deta_2nd_ce_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 2, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = 0.5_rk * (f(iVar,ic,jc+1,kc) - f(iVar,ic,jc-1,kc))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! BC Schemes
      jc = 1
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc+1,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO

      jc = Ny
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic,jc-1,kc)
          ENDDO
        ENDDO
      ENDDO      

  end subroutine Cal_Deri_Deta_2nd_ce_numvar
  
  subroutine Cal_Deri_Deta_1st_do_numvar(df,f)
      implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
      df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny-1
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc+1,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! BC Schemes
      jc = 1
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc+1,kc) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO

      jc = Ny
      do kc = 1, nz_local
        do ic = 1, nx_local
          do iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic,jc-1,kc)
          ENDDO
        ENDDO
      ENDDO      

  end subroutine Cal_Deri_Deta_1st_do_numvar
  
  subroutine Cal_Deri_Deta_2nd_ce(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc
    df = 0.d0;         
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 2, Ny-1
        DO ic = 1, nx_local
            df(ic,jc,kc) = 0.5_rk * (f(ic,jc+1,kc) - f(ic,jc-1,kc))
        ENDDO
      ENDDO
    ENDDO

    ! BC Schemes
      jc = 1
      do kc = 1, nz_local
        do ic = 1, nx_local
            df(ic,jc,kc) = f(ic,jc+1,kc) - f(ic,jc,kc)
        ENDDO
      ENDDO

      jc = Ny
      do kc = 1, nz_local
        do ic = 1, nx_local
            df(ic,jc,kc) = f(ic,jc,kc) - f(ic,jc-1,kc)
        ENDDO
      ENDDO      

  end subroutine Cal_Deri_Deta_2nd_ce

  ! This subroutine calculate the 1st order derivative zeta of the variables
  subroutine Cal_Deri_Dzeta_1st_up_per_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.0_rk
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc) - f(iVar,ic,jc,kc-1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO    

  end subroutine Cal_Deri_Dzeta_1st_up_per_numvar

  subroutine Cal_Deri_Dzeta_2nd_ce_per_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.0_rk
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = 0.5_rk*(f(iVar,ic,jc,kc+1) - f(iVar,ic,jc,kc-1))
          ENDDO
        ENDDO
      ENDDO
    ENDDO    

  end subroutine Cal_Deri_Dzeta_2nd_ce_per_numvar
  
  subroutine Cal_Deri_Dzeta_1st_do_per_numvar(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1:NumVar,1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc,iVar
    df = 0.d0;
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
          DO iVar = 1, NumVar
            df(iVar,ic,jc,kc) = f(iVar,ic,jc,kc+1) - f(iVar,ic,jc,kc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO    

  end subroutine Cal_Deri_Dzeta_1st_do_per_numvar

  subroutine Cal_Deri_Dzeta_2nd_ce_per(df,f)
    implicit none
    real( kind = rk ), intent(in):: f(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)
    real( kind = rk ), intent(out)::df(1-overLAP:nx_local+overLAP,1:Ny,1-overLAP:nz_local+overLAP)

    integer( kind = ik ):: ic,jc,kc

    df = 0.0_rk
    ! Here, we direct calculate the 1st order derivative of the variables
    DO kc = 1, nz_local
      DO jc = 1, Ny
        DO ic = 1, nx_local
            df(ic,jc,kc) = 0.5_rk*(f(ic,jc,kc+1) - f(ic,jc,kc-1))
        ENDDO
      ENDDO
    ENDDO    

  end subroutine Cal_Deri_Dzeta_2nd_ce_per

END Module FD1st_order