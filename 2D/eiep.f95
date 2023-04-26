! File contains the entropy inequality predictors

module eiep
        use consts
        use fluxes
        use eno

        contains

! Implementation of the smoothstep function
        elemental function smoothstep(x)
                double precision, intent(in) :: x
                double precision :: smoothstep

                if (x .LE. 0.0D0) then
                        smoothstep = 0.0D0
                else if (x .LE. 1.0D0) then
                        smoothstep = 6*x**5 - 15*x**4 + 10*x**3
                else
                        smoothstep = 1.0D0
                end if
        end function smoothstep

! Dissipative two points entropy flux in x direction
        function fdis(ua, ub)
        double precision, intent(in), dimension(ncons) :: ua, ub
        double precision :: fdis
        fdis = PUEuler(0.5D0*(ua + ub + lambda*(fEuler(ua) - fEuler(ub))))&
                                - 0.5D0*(PUEuler(ua) + PUEuler(ub) + lambda*(PFEuler(ua) - PFEuler(ub)))
        end function fdis

! Dissipative entropy two point flux in y direction
        function gdis(ua, ub)
        double precision, intent(in), dimension(ncons) :: ua, ub
        double precision :: gdis
        gdis = PUEuler(0.5D0*(ua + ub +lambda*(gEuler(ua) - gEuler(ub))))&
                                - 0.5D0*(PUEuler(ua) + PUEuler(ub) + lambda*(PGEuler(ua) - PGEuler(ub)))
        end function gdis

! ENO 2 Lax-Friedrichs Entropy inequality predictor
        subroutine E2LFeiep(u, dm)
                double precision, dimension(nx, ny, ncons), intent(in) :: u
                double precision, dimension(nx, ny), intent(out) :: dm
                double precision, dimension(nx, ny) :: dmt, dmr, Uar
                integer :: k, l
                integer, dimension(2) :: umax, umin
                double precision, dimension(nx, ny, 4) :: ul, ur, ub, ut
                double precision :: sref

                do k=1, nx
                        do l = 1, ny
                                Uar(k, l) = PUEuler(u(k, l, :))
                        end do
                end do
                
                umin = minloc(Uar)
                umax = maxloc(Uar)


! Calculation of the entropy dissipation reference
                sref = min(fdis(u(umin(1), umin(2), :), u(umax(1), umax(2), :)), &
                                fdis(u(umax(1), umax(2), :), u(umin(1), umin(2), :)),&
                                gdis(u(umin(1), umin(2), :), u(umax(1), umax(2), :)),&
                                gdis(u(umax(1), umax(2), :), u(umin(1), umin(2), :))) 

! ENO recovery of the conserved variables
!               $OMP PARALLEL SHARED(u, ur, ul, ub, ut, dmt, dmr, dm)
!               $OMP DO
                do l=1, nx
                        call ENO2Mat(u(l, :, :), ub(l, :, :), ut(l, :, :), ny)
                end do
                do l=1, ny
                        call ENO2Mat(u(:, l, :), ul(:, l, :), ur(:, l, :), nx)
                end do
!               $OMP END DO

! Calculation of the entropy dissipation of the recovered point values
!               $OMP DO PRIVATE(l)
                do k=1, nx
                        do l=1, ny
                                dmt(k, l) = PUEuler(0.5D0*(ut(k, l, :) + ub(k, warp(l+1, ny), :) +&
                                         lambda*(gEuler(ut(k, l, :)) - gEuler(ub(k, warp(l+1, ny), :)))))&
                                         - 0.5D0*(PUEuler(ut(k, l, :)) + PUEuler(ub(k, warp(l+1, ny), :)) &
                                         + lambda*(PGEuler(ut(k, l, :)) - PGEuler(ub(k, warp(l+1, ny), :))))
                                dmr(k, l) = PUEuler(0.5D0*(ur(k, l, :) + ul(warp(k+1, nx), l, :) +&
                                         lambda*(fEuler(ur(k, l, :)) - fEuler(ul(warp(k+1, nx), l, :)))))&
                                         - 0.5D0*(PUEuler(ur(k, l, :)) + PUEuler(ul(warp(k+1, nx), l, :)) &
                                         + lambda*(PFEuler(ur(k, l, :)) - PFEuler(ul(warp(k+1, nx), l, :))))  
                        end do
                end do
!               $OMP END DO

! Remapping into a convex parameter using the smoothstep function
!               $OMP DO PRIVATE(l)
                do k=1, nx
                        do l=1, ny
                                dm(k, l) = (min(dmt(k, l), dmr(k, l), dmt(k, warp(l-1, ny)),&
                                                 dmr(warp(k-1, nx), l))/sref - redc)/scalc
                        end do
                end do
!               $OMP END DO
!               $OMP DO 
                do l=1, ny
                        dm(:, l) = smoothstep(dm(:, l))
                end do
!               $OMP END DO
!               $OMP END PARALLEL                

        end subroutine


end module eiep
