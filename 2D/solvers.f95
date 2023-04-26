! File contains Time integration routines and the Top 
! Solver routines that are called by particular test cases

module solvers
        use consts
        use fluxes
        use eiep
        use supmoll
        use eno
        use fluxmatsparam
        use gtk

        implicit none
        contains

!       Calculates Lax-Friedrichs scheme into du
!       For test/reference purposes
        subroutine LLF(du, u)
                double precision, dimension(nx, ny, 4), intent(out) :: du
                double precision, dimension(nx, ny, 4), intent(in) :: u
                double precision, dimension(4) :: deltax, deltay
                integer :: k, l

!               $OMP PARALLEL DO PRIVATE(k, deltax, deltay) SHARED(du, u)
                do l = 1, ny
                        do k=1, nx
                                deltax = fEulerLF(u(warp(k-1, nx), l, :), u(k, l, :)) - fEulerLF(u(k,l, :), u(warp(k+1, nx), l, :))
                                deltay = gEulerLF(u(k, warp(l-1, ny), :), u(k, l, :)) - gEulerLF(u(k, l, :), u(k, warp(l+1, ny), :))
                                du(k,l, :) = (deltax + deltay)/dx
                        end do
                end do
!               $OMP END PARALLEL DO
        end subroutine LLF

!       Calculates the ENO2 reconstruction combined with the LLF flux
!       Into $du$ for reference purposes
        subroutine LE2LF(du,u)
                double precision, dimension(nx, ny, 4), intent(out) :: du
                double precision, dimension(nx, ny, 4), intent(in) :: u
                double precision, dimension(nx, ny, 4) :: ul, ur, ub, ut
                double precision, dimension(4) :: deltax, deltay
                integer :: k, l

!               $OMP PARALLEL DO SHARED(u, ur, ul, ub, ut)
                do l=1, ny
                        call ENO2Mat(u(:, l, :), ul(:, l, :), ur(:, l, :), nx)
                end do
!               $OMP END PARALLEL DO

!               $OMP PARALLEL DO SHARED(u, ur, ul, ub, ut)
                do k=1, nx
                        call ENO2Mat(u(k, :, :), ub(k, :, :), ut(k, :, :), ny)
                end do
!               $OMP END PARALLEL DO


!               $OMP PARALLEL DO PRIVATE(k, deltax, deltay) SHARED(du, u)
                do l = 1, ny
                        do k=1, nx
                                deltax = fEulerLF(ur(warp(k-1, nx), l, :), ul(k, l, :)) &
                                        - fEulerLF(ur(k,l, :), ul(warp(k+1, nx), l, :))
                                deltay = gEulerLF(ut(k, warp(l-1, ny), :), ub(k, l, :)) &
                                        - gEulerLF(ut(k, l, :), ub(k, warp(l+1, ny), :))
                                du(k,l, :) = deltax/dx + deltay/dy
                        end do
                end do
!               $OMP END PARALLEL DO
        end subroutine LE2LF

!       calculates GT scheme into $du$.
 subroutine LEGT(du,u)
                double precision, dimension(nx, ny, 4), intent(out) :: du
                double precision, dimension(nx, ny, 4), intent(inout) :: u
                double precision, dimension(nx, ny) :: dm, alp
                double precision, dimension(4) :: deltax, deltay
                integer :: k, l


!               Reset BCs 
                call sanitizeux(u, 0.D0, nx, ny, ncons)

!               Calculate output of the entropy inequality predictors
                call E2LFeiep(u, dm)
                call NLmoll(dm, alp)

!               Calculate the changes in x
                !$OMP PARALLEL DO PRIVATE(k, deltax, deltay) SHARED(du, u, alp)
                do l = 1, ny
                        do k=1, nx
                                deltax = LMRKf(u, alp, k-1, l) - LMRKf(u, alp, k, l)
                                du(k,l, :) = (deltax)/dx
                        end do
                end do
                !$OMP END PARALLEL DO


                call sanitizeuy(u, 0.D0, nx, ny, ncons)
                call E2LFeiep(u, dm)
                call NLmoll(dm, alp)

!               Add the derivatives in y
                !$OMP PARALLEL DO PRIVATE(k, deltax, deltay) SHARED(du, u, alp)
                do l = 1, ny
                        do k=1, nx
                                deltay = LMRKg(u, alp, k, l-1) - LMRKg(u, alp, k, l)
                                du(k,l, :) = du(k, l, :) +  (deltay)/dy
                        end do
                end do
                !$OMP END PARALLEL DO


 end subroutine LEGT



! Calculates Forward Euler Step
! Using the Lax-Friedrichs Scheme
        subroutine fELF(unew, u)
                implicit none
                double precision, dimension(nx, ny, ncons), intent(out) :: unew
                double precision, dimension(nx, ny, ncons), intent(in) :: u
                double precision, dimension(nx, ny, ncons) :: du
                call LLF(du, u)
                unew = u + dt*du
        end subroutine fELF

! Calculates SSPRK33 step
! Using the ENO2 Lax-Friedrichs scheme
        subroutine SSPRK33ELF(unew, u)
                implicit none
                double precision, dimension (nx, ny, 4), intent(out) :: unew
                double precision, dimension (nx, ny, 4), intent(out) :: u
                double precision, dimension (nx, ny, 4) :: du, u1, u2
                call LE2LF(du, u)
                u1 = u + dt*du
                call LE2LF(du, u1)
                u2 = 0.75D0*u + 0.25D0*u1 + 0.25D0*dt*du
                call LE2LF(du, u2)
                unew = u/3.0D0 + 2.0D0/3.0D0*u2 + 2.0D0/3.0D0*dt*du 
        end subroutine SSPRK33ELF

! Calculates SSPRK33 step
! Using the GT scheme
        subroutine SSPRK33GT(unew, u)
                implicit none
                double precision, dimension (nx, ny, 4), intent(out) :: unew
                double precision, dimension (nx, ny, 4), intent(out) :: u
                double precision, dimension (nx, ny, 4) :: du, u1, u2
                call LEGT(du, u)
                u1 = u + dt*du
                call LEGT(du, u1)
                u2 = 0.75D0*u + 0.25D0*u1 + 0.25D0*dt*du
                call LEGT(du, u2)
                unew = u/3.0D0 + 2.0D0/3.0D0*u2 + 2.0D0/3.0D0*dt*du
        end subroutine SSPRK33GT



! Does steps until tend with ENO2 LF
        subroutine solE2LF(ures, u0, tend)
                use consts
                implicit none
                double precision, dimension(nx, ny, 4), intent(out) :: ures
                double precision, dimension(nx, ny, 4), intent(inout) :: u0
                double precision :: tend
                integer :: k, nsteps 
                nsteps = ceiling(tend/(2*dt))
                do k = 1, nsteps
                        call SSPRK33ELF(ures, u0)
                        call SSPRK33ELF(u0, ures)
                        if (mod(k, 100) == 1) then
                                write (*, *) k, " / ", nsteps
                        end if
                end do
                ures(:, :, :) = u0(:, :, :)
        end subroutine solE2LF

! Does steps until tend with SSPRK33GT
! New: with boundary conditions
        subroutine solEGT(ures, u0, tend)
                use consts
                implicit none
                double precision, dimension(nx, ny, 4), intent(out) :: ures
                double precision, dimension(nx, ny, 4), intent(inout) :: u0
                double precision :: tend
                integer :: k, nsteps

                nsteps = ceiling(tend/(2*dt))

                do k = 1, nsteps
                        call SSPRK33GT(ures, u0)
                        call SSPRK33GT(u0, ures)
                        if (mod(k, 100) == 1) then
                                write (*, *) k, " / ", nsteps
                        end if
                end do
                ures(:, :, :) = u0(:, :, :)
        end subroutine solEGT

! function version for python        
        function fsole2elf(u0, tend)
                use consts
                implicit none
                double precision, dimension(nx, ny, 4) :: fsole2elf
                double precision, dimension(nx, ny, 4), intent(inout) :: u0
                double precision :: tend
                call solE2LF(fsole2elf, u0, tend)
        end function fsole2elf

! function version for python        
        function fsolEGT(u0, tend)
                use consts
                implicit none
                double precision, dimension(nx, ny, 4) :: fsolEGT
                double precision, dimension(nx, ny, 4), intent(inout) :: u0
                double precision :: tend
                call solEGT(fsolEGT, u0, tend)
        end function fsolEGT

! pythonic init
        function init()
                integer :: init
                call Fluxmatsinit()
                fnum => fEulerCons
                fdisp => fEulerLLF
                gnum => gEulerCons
                gdisp => gEulerLLF
        end function init


end module solvers
