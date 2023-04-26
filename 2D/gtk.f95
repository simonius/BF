! Module gtk contains all needed functions for the implementation
! of linear combination fluxes

module gtk
        use consts
        use fluxes
        use fluxmatsparam
        implicit none

        contains

! Implementation of the LeFloch/Mercier/Rhode/Klein Matrix flux
! for finite domains, x component
! u: conserved variables
! alp: convex parameter
! k, l: indices in x and y direction
        function LMRKf(u, alp, k, l)
                implicit none
                double precision, dimension(ncons) :: LMRKf
                double precision, dimension(nx, ny, ncons) :: u
                double precision, dimension(nx, ny), intent(in) :: alp
                double precision, dimension(ncons) :: ul, ur
                double precision :: lalp
                integer, intent(in) :: k, l
                integer :: i, j
                integer, dimension(2) :: tolb, toub
                integer :: matindex
!               Lets first calculate the distance to the nearest boundary
                tolb = downdist(k, l)
                toub = updist(k, l)
!               We are only interested in x distances
                if (tolb(1) < lmrorder) then
                        matindex = tolb(1)
                else if (toub(1) < lmrorder) then ! We assume that there is enough space in the grid for one complete stencil
                        matindex = 2*lmrorder+1-toub(1)
                else
                        matindex = lmrorder+1
                end if


!               initialize
                do i=1, ncons
                        LMRKf(i) = 0.0D0
                end do

!               We return prematurely with zeros when asked outside of the domain
                if (tolb(1) < 1 .OR. tolb(2) < 1 .OR. toub(1) < 1 .OR. toub(2) < 1) then
                        return
                end if

!               and add together
                do i=1, 2*lmrorder+1
                        do j = 1, 2*lmrorder+1
                                if ((j > i) .AND. (fluxmats(matindex, i, j) .NE. 0.0D0)) then
!               We ASSUME that we stay inside the domain and array because tolb and toub is correct
                                        ul = u(k+i-matindex, l, :)
                                        ur = u(k+j-matindex, l, :)
                                        lalp = max(alp(k, l), alp(k+1, l))
                                        LMRKf = LMRKF + 2*fluxmats(matindex,i, j)*(lalp*fdisp(ul, ur)&
                                                        + (1-lalp)*fnum(ul, ur))
                                end if
                        end do
                end do



        end function LMRKf


! Implementation of the LeFloch/Mercier/Rhode/Klein Matrix flux
! for finite domains, y component

        function LMRKg(u, alp, k, l)
                implicit none
                double precision, dimension(ncons) :: LMRKg
                double precision, dimension(nx, ny, ncons) :: u
                double precision, dimension(nx, ny), intent(in) :: alp
                double precision, dimension(ncons) :: ul, ur
                double precision :: lalp
                integer, intent(in) :: k, l
                integer :: i, j
                integer, dimension(2) :: tolb, toub
                integer :: matindex
!               Lets first calculate the distance to the nearest boundary
                tolb = downdist(k, l)
                toub = updist(k, l)
!               We are only interested in x distances
                if (tolb(2) < lmrorder) then
                        matindex = tolb(2)
                else if (toub(2) < lmrorder) then ! We assume that there is enough space in the grid for one complete stencil
                        matindex = 2*lmrorder+1-toub(2)
                else
                        matindex = lmrorder+1
                end if
!               initialize
                do i=1, ncons
                        LMRKg(i) = 0.0D0
                end do

!               We return prematurely with zeros when asked outside of the domain
                if (tolb(1) < 1 .OR. tolb(2) < 1 .OR. toub(1) < 1 .OR. toub(2) < 1) then
                        return
                end if


!               and add together
                do i=1, 2*lmrorder+1
                        do j = 1, 2*lmrorder+1
                                if ((j > i) .AND. (fluxmats(matindex, i, j) .NE. 0.0D0)) then
!               We ASSUME that we stay inside the domain and array because tolb and toub is correct
                                        ul = u(k, l+i-matindex, :)
                                        ur = u(k, l+j-matindex, :)
                                        lalp = max(alp(k, l), alp(k, l+1))
                                        LMRKg = LMRKg + 2*fluxmats(matindex,i, j)*(lalp*gdisp(ul, ur)&
                                                        + (1-lalp)*gnum(ul, ur))
                                end if
                        end do
                end do



        end function LMRKg


end module gtk
