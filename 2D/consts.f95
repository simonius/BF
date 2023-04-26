! File contains the parameters of the simulation
! Should be included in all other files

module consts
        implicit none
!       Fortran is not case sensitive 
!       Grid sizes
        integer, parameter :: nbase = 80
        integer, parameter :: nx = nbase*3
        integer, parameter :: ny = nbase

! Periodicity
        logical, parameter :: periodic = .FALSE.

!       # of conserved variables
        integer, parameter :: ncons = 4

!       snapshots that are saved
        integer, parameter :: nsnaps = 50
        double precision, parameter :: xmin = 0.0
        double precision, parameter :: xmax = 3.0
        double precision, parameter :: ymin = 0.0
        double precision, parameter :: ymax = 1.0

!       Used CFL number and assumed maximum speed
        double precision, parameter :: cbound = 6.0
        double precision, parameter :: CFL = 0.3
        double precision, parameter :: lambda = CFL/cbound

!       Grid constants
        double precision, parameter :: dx = (xmax-xmin) / nx
        double precision, parameter :: dy = (ymax-ymin) / ny
        double precision, parameter :: dt = lambda * min(dx, dy)

!       Ideal gas law constant
        double precision, parameter :: gamma = 1.4

!       Pi constant for initial conditions
        double precision, parameter :: PI=4.D0*DATAN(1.D0)

!       Scaling constants for the Entropy inequality predictors        
        double precision, parameter :: redc = 1.0D-3
        double precision, parameter :: scalc = 1.0D-3

!       the first one should always be smaller than the second one       
        integer, parameter :: lmrorder = 3
        integer, parameter :: mwidth = 10

!       f2py lacks support for abstract interfaces
        INTERFACE
                FUNCTION numflux(ul, ur)
                        double precision, dimension(4), INTENT(IN) :: ul, ur
                        DOUBLE PRECISION, dimension(4) :: numflux
                END FUNCTION numflux
        END INTERFACE

! Interface definitions of distance function
        INTERFACE
!               returns distance to the boundary of the domain
                FUNCTION distfunc(k, l)
                        integer, intent(in) :: k, l
                        integer, dimension(2) :: distfunc
                END FUNCTION
        END INTERFACE

! Pointers onto functions that give distance to the boundary 
        procedure(distfunc), pointer :: updist
        procedure(distfunc), pointer :: downdist

! Interface definition of sanitizer functions
        Interface
!               Resets ghost cells for reflective BC and so on
!               We can not use the global size variable here because interfaces can 
!               connect to other namespaces where this global could be missing
                subroutine sanitizer(u, t, nupper, kupper, lupper)
                        double precision, dimension(nupper, kupper, lupper), intent(inout) :: u
                        double precision, intent(in) :: t
                        integer, intent(in) :: nupper, kupper, lupper
                end subroutine
        end interface

! Pointers to routines that sanitize the boundary nodes after every timestep
        procedure(sanitizer), pointer :: sanitizeux
        procedure(sanitizer), pointer :: sanitizeuy

! Pointers for dissipative and conservative numerical fluxes
        procedure(numflux), pointer :: fnum !=> fEulerCons
        procedure(numflux), pointer :: fdisp !=> fEulerLF
        procedure(numflux), pointer :: gnum !=> gEulerCons
        procedure(numflux), pointer :: gdisp !=> gEulerLF

        contains

! Function recalculates the grid points outside the domain for values outside of the domain. 
! Assuming periodic boundaries
        function warp(i, imax)
                integer :: i, imax
                integer :: warp
                warp = MODULO(i-1, imax)+1
                
        end function warp

end module consts
