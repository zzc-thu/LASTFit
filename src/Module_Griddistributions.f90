module grid_distribution
    implicit none 
    contains 
 
    subroutine GridDistributions(N1, N2, N3, Hscale, eta) 
        integer, intent(in) :: N1, N2, N3 
        real(kind = 8), intent(in) :: Hscale 
        real(kind = 8), allocatable, dimension(:), intent(out) :: eta 
        integer :: N2s, N2e, N3s, N3e 
        integer :: i 
        real(kind = 8), allocatable, dimension(:) :: refX 
        real(kind = 8) :: ymax, y2max, y1max, y1i 
        real(kind = 8), allocatable, dimension(:) :: xi1_uni, y1, y3, TotalY 
        real(kind = 8) :: a, b 
        real(kind = 8) :: cdYdX1, cdYdX2 
        real(kind = 8), allocatable, dimension(:) :: XHermite, YHermite, Yderve 
 
        N2s = N1 
        N2e = N2s + N2 - 1 
        N3s = N2e 
        N3e = N3s + N3 - 1 
 
        allocate(refX(1:N3e)) 
        do i = 1, N3e 
            refX(i) = real(i, kind = 8) 
        end do 
 
        ymax = 1.0_8 * Hscale 
        y2max = 0.8_8 * Hscale 
        y1max = 0.2_8 * Hscale 
        y1i = 0.01_8 * Hscale 
 
        allocate(xi1_uni(1:N1)) 
        do i = 1, N1 
            xi1_uni(i) = -1.0_8 + 2.0_8 * (real(i - 1, kind = 8) / real(N1 - 1, kind = 8)) 
        end do 
 
        a = y1i * y1max / (y1max - 2.0_8 * y1i) 
        b = 1.0_8 + 2.0_8 * a / y1max 
 
        allocate(y1(1:N1)) 
        do i = 1, N1 
            y1(i) = a * (1.0_8 + xi1_uni(i)) / (b - xi1_uni(i)) 
        end do 
 
        allocate(y3(1:N3)) 
        do i = 1, N3 
            y3(i) = y2max + (ymax - y2max) * (real(i - 1, kind = 8) / real(N3 - 1, kind = 8)) 
        end do 
 
        allocate(TotalY(1:N3e)) 
        do i = 1, N1 
            TotalY(i) = y1(i) 
        end do 
        do i = 1, N3 
            TotalY(N3s + i - 1) = y3(i) 
        end do 
 
        ! Calculate the Derivatives at points TotalY(N2s) and TotalY(N2e) 
        cdYdX1 = (TotalY(N1) - TotalY(N1 - 1)) / (refX(N1) - refX(N1 - 1)) 
        cdYdX2 = (TotalY(N3s + 1) - TotalY(N3s)) / (refX(N3s + 1) - refX(N3s)) 
 
        allocate(XHermite(2)) 
        allocate(YHermite(2)) 
        allocate(Yderve(2)) 
        XHermite(1) = refX(N1) 
        XHermite(2) = refX(N3s) 
        YHermite(1) = TotalY(N1) 
        YHermite(2) = TotalY(N3s) 
        Yderve(1) = cdYdX1 
        Yderve(2) = cdYdX2 
 
        call Hermitezi(XHermite, YHermite, Yderve, refX(N2s:N2e), TotalY(N2s:N2e)) 
 
        allocate(eta(N1 + N2 + N3 - 2)) 
        eta = TotalY 
 
        deallocate(refX, xi1_uni, y1, y3, TotalY, XHermite, YHermite, Yderve) 
 
    end subroutine GridDistributions 
 
    subroutine Hermitezi(XHermite, YHermite, Yderve, X, Y) 
        real(kind = 8), dimension(2), intent(in) :: XHermite, YHermite, Yderve 
        real(kind = 8), dimension(:), intent(in) :: X 
        real(kind = 8), dimension(:), intent(out) :: Y 
        integer :: i 
        real(kind = 8) :: t 
 
        do i = 1, size(X) 
            t = (X(i) - XHermite(1)) / (XHermite(2) - XHermite(1)) 
            Y(i) = (2.0_8 * t**3 - 3.0_8 * t**2 + 1.0_8) * YHermite(1) & 
                 + (t**3 - 2.0_8 * t**2 + t) * (XHermite(2) - XHermite(1)) * Yderve(1) & 
                 + (-2.0_8 * t**3 + 3.0_8 * t**2) * YHermite(2) & 
                 + (t**3 - t**2) * (XHermite(2) - XHermite(1)) * Yderve(2) 
        end do 
 
    end subroutine Hermitezi 
 
end module grid_distribution