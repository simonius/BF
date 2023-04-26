! File contains the helper routines to implement
! The domain of the forward facing step problem

module ffsdomain
        use tc
        use consts
        implicit none        

!       distance from the beginning of the rectangle torwrds the step in x
        integer, parameter :: kstep = nx/5

!       distance from the beginning of the rectangle torwards the step in y
        integer, parameter :: lstep = ny/5
!       reflective symmtry fix

        integer, parameter :: reflsymdist = 4

!       inflow variables (physical)
        double precision, parameter, dimension(4) :: uin = (/1.4D0, 3.D0, 0.D0, 1.D0/) 

!       initial condition (physical) 
        double precision, parameter, dimension(4) :: uinit = (/1.4D0, 3.D0, 0.D0, 1.D0/)

        contains

!       Initalizes domain
        subroutine initffs(u)
                implicit none
                double precision, intent(inout), dimension(nx, ny, ncons) :: u
                integer :: k, l
                
                do k=1, nx
                        do l=1, ny
                                u(k, l, :) = ConsVar(uinit(1), uinit(2), uinit(3), uinit(4))
                        end do
                end do
        end subroutine initffs

!       distance to the lower boundary
        function downdistfunc(k, l)
                implicit none
                integer, intent(in) :: k, l
                integer, dimension(2) :: downdistfunc

                if (k < kstep) then 
!                       first part of the domain, without step
                        downdistfunc(1) = k
                        downdistfunc(2) = l
                else
                        downdistfunc(1) = k
                        downdistfunc(2) = l - lstep + reflsymdist
                end if
        end function

!       distance to the upper boundary
        function updistfunc(k, l)
                implicit none
                integer, intent(in) :: k, l
                integer, dimension(2) :: updistfunc
               
                if (l < lstep) then
                        updistfunc(1) = kstep - k + reflsymdist
                else
                        updistfunc(1) = nx - k
                end if
                updistfunc(2) = ny - l
        end function

!       mirror speeds in x
        function mirrorx(u)
                implicit none
                double precision, dimension(ncons), intent(in) :: u
                double precision, dimension(ncons) :: mirrorx
                mirrorx(1) = u(1)
                mirrorx(2) = -u(2)
                mirrorx(3) = u(3)
                mirrorx(4) = u(4)
        end function

!       mirror speeds in y
        function mirrory(u)
                implicit none
                double precision, dimension(ncons), intent(in) :: u
                double precision, dimension(ncons) :: mirrory
                mirrory(1) = u(1)
                mirrory(2) = u(2)
                mirrory(3) = -u(3)
                mirrory(4) = u(4)
        end function


!       resets in and outflow boundaries and 
!       mirrors for reflective boundaries
!       in the x direction
        subroutine sanitizeffsx(u, t, nupper, kupper, lupper)
                double precision, dimension(nupper, kupper, lupper), intent(inout) :: u
                double precision, intent(in) :: t
                integer, intent(in) :: nupper, kupper, lupper
                double precision, dimension(ncons) :: uinc
                integer :: k, l

                uinc = ConsVar(uin(1), uin(2), uin(3), uin(4))
!               Copy the inflow values to the inflow boundary
                do k=1, kupper
                        u(1, k, :) = uinc
                end do
!               Copy the last grid values onto the outflow boundary
                do k=1, kupper
                        do l=1,lmrorder
                                u(nx-lmrorder+l, k, :) = u(nx-lmrorder, k, :)
                        end do
                end do

!               Reset the interrior of the step
                do k=kstep+2, nx
                        do l=1, lstep-2
                                u(k, l, :) = ConsVar(uinit(1), uinit(2), uinit(3), uinit(4))
                        end do
                end do

!               Mirror the values for the vertical part of the step in x
                do k=0, lmrorder
                        do l=1, lstep + lmrorder
                                u(kstep+k, l, :) = mirrorx(u(kstep-k-1, l, :))
                        end do
                end do 
 
        end subroutine sanitizeffsx

!       resets in and outflow boundaries and 
!       mirrors for reflective boundaries
!       in the y direction
        subroutine sanitizeffsy(u, t, nupper, kupper, lupper)
                double precision, dimension(nupper, kupper, lupper), intent(inout) :: u
                double precision, intent(in) :: t
                integer, intent(in) :: nupper, kupper, lupper
                double precision, dimension(ncons) :: uinc
                integer :: k, l
!               Mirror the values for the upper boundary                 
                do k=1, nupper
                        do l=0, lmrorder
                                u(k, ny-lmrorder+l, :) = mirrory(u(k, ny-lmrorder-l-1, :))
                        end do
                end do

!               Mirror the values for the first part of the lower boundary
                do k=1, kstep
                        do l=0, lmrorder
                                u(k, lmrorder-l+1, :) = mirrory(u(k, lmrorder+l+2, :))
                        end do
                end do 

!               Mirror the values for the second part of the lower boundary
                do k=kstep, nx
                        do l=0, lmrorder
                                u(k, lstep+lmrorder-l+1, :) = mirrory(u(k, lstep+lmrorder+l+2, :))
                        end do
                end do
        end subroutine sanitizeffsy
end module ffsdomain
