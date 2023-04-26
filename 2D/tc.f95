module tc
        use consts
        implicit none

        contains
! Converts between physical and conserved variables
                function ConsVar(rho, vx, vy, p)
                        implicit none
                        double precision, dimension(4) :: ConsVar
                        double precision :: rho, vx, vy, p
                        ConsVar(1) = rho
                        ConsVar(2) = rho*vx
                        ConsVar(3) = rho*vy
                        ConsVar(4) = p / (gamma-1) + 0.5 * (vx**2 + vy**2)*rho
                end function ConsVar

                function PhysVar(rho, rhovx, rhovy, E)
                        implicit none
                        double precision, dimension(4) :: PhysVar
                        double precision :: rho, rhovx, rhovy, E
                end function PhysVar

                function KHIS(x, y)
                        use consts
                        implicit none
                        double precision, dimension(4) :: KHIS
                        double precision :: x, y, eps = 0.01
                        if (abs(y-0.5) < 0.25) then
                                KHIS = ConsVar(2.D0, 0.5D0, sin(4*pi*x)*eps, 1.D0)
                        else
                                KHIS = ConsVar(1.D0, -0.5D0, sin(4*pi*x)*eps, 1.D0)
                        end if
                end function KHIS

! KHI saves initial conditions for Kelvin Helmholt into array uar
                subroutine KHI(uar)
                        use consts
                        implicit none
                        double precision, dimension(nx, ny, 4), intent(inout) :: uar
                        integer :: k, l
                       
                        do k=1, nx
                                do l=1, ny
                                        uar(k, l, :) = KHIS(dx*k, dy*l)
                                end do
                        end do
                end subroutine KHI

                function khiu0()
                        use consts
                        implicit none
                        double precision, dimension(nx, ny, 4) :: khiu0
                        call KHI(khiu0)
                end function khiu0

! Supmoll testcase - single peak
                function smtest()
                        use consts
                        implicit none
                        double precision, dimension(nx, ny) :: smtest
                        integer :: k,l
                        do k=1, nx
                                do l=1, ny
                                        smtest(k, l) = 0.0D0
                                end do
                        end do
                        smtest(128, 128) = 1.0D0

                end function smtest
end module tc

