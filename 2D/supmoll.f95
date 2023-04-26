! Efficient sup-mollification version using exponentiation of the operator. 

module supmoll
        use consts
        implicit none
        contains

!       subroutine does nonlinear mollification using sucessive small
!       Stencil supmollification on periodic grids
        subroutine NLmoll(dm, alp)
                double precision, dimension(nx, ny), intent(inout) :: dm
                double precision, dimension(nx, ny), intent(out) :: alp
                double precision :: ck 
                integer :: k, l, m


!               first round, form plateaus
!               $OMP PARALLEL SHARED(alp, dm, ck)
!               $OMP DO PRIVATE(k)
                do l=1, ny
                        do k=1, nx
                                alp(k, l) = max(dm(warp(k-1, nx), l), dm(k, l), dm(warp(k+1, nx), l), &
                                        dm(k, warp(l-1, ny)), dm(k, warp(l+1, ny)))
                        end do
                end do
!               $OMP END DO
!               next rounds, roll of viscosity
                do m= 2, mwidth-1, 2
                        ck = dble((m-mwidth)*(m-1))/dble(((m - 1 - mwidth)*m))
                        !write (*,*) m, ck
                        !$OMP DO  PRIVATE(k)
                        do l=1, ny
                                do k=1, nx
                                        dm(k, l) = max(alp(k, l), ck*alp(warp(k-1, nx), l),&
                                                ck*alp(warp(k+1, nx), l), ck*alp(k, warp(l-1, ny)), ck*alp(k, warp(l+1, ny)))
                                end do
                        end do
                        !$OMP END DO
                        ck = dble((m+1-mwidth)*(m))/dble(((m - mwidth)*(m+1)))
                        !$OMP DO  PRIVATE(k)
                        do l=1, ny
                                do k=1, nx
                                        alp(k, l) = max(dm(k, l), ck*dm(warp(k-1, nx), l),&
                                                ck*dm(warp(k+1, nx), l), ck*dm(k, warp(l-1, ny)), ck*dm(k, warp(l+1, ny)))
                                end do
                        end do
                        !$OMP END DO
                end do
!               $OMP END PARALLEL
        end subroutine NLmoll

! Pythonic version
        function pwNLmoll(dm)
                double precision, dimension(nx, ny), intent(inout) :: dm
                double precision, dimension(nx, ny) :: pwNLmoll
                call NLmoll(dm, pwNLmoll)
        end function pwNLmoll

end module supmoll
