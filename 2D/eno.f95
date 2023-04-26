! ENO 2 Reconstruction - needed for the reference solution calculation and the 
! Entropy inequality predictor

module eno
        use consts
        implicit none
        contains

! Single ENO2 reconstruction, takes 3 arguments and decides which reconstruction has lower total variation
                subroutine ENO2(u, ul, ur)
                        implicit none
                        double precision, dimension(3), intent(in) :: u
                        double precision, intent(out) :: ul, ur
                        double precision :: a1, a2
                        a1 = u(2) - u(1)
                        a2 = u(3) - u(2)
                        if (a1**2 < a2**2) then
                                ul = (u(1) + u(2)) / 2.0D0
                                ur = u(2) + 0.5D0*a1
                        else
                                ul = u(2) - 0.5D0*a2
                                ur = (u(2) + u(3)) / 2.0D0
                        end if
                end subroutine ENO2

! ENO2 reconstruction for severall conserved variables in one space dimension
                subroutine ENO2Mat(u, ul, ur, nl)
                        implicit none
                        integer, intent(in) :: nl
                        double precision, dimension(nl, ncons), intent(in) :: u
                        double precision, dimension(nl, ncons), intent(out) :: ul, ur
                        integer :: i, j
                        do i=2,nl-1
                                do j=1, ncons
                                        call ENO2(u(i-1:i+1, j), ul(i, j), ur(i, j))
                                end do
                        end do
                        do j=1, ncons
                                call ENO2((/u(nl, j),  u(1, j), u(2, j)/), ul(1, j), ur(1, j))
                                call ENO2((/u(nl-1, j), u(nl, j), u(1, j)/), ul(nl, j), ur(nl, j))
                        end do
                end subroutine ENO2Mat


end module eno
