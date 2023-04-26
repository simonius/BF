! File contains the entropy conservative and entropy dissipative two-point fluxes

module fluxes
        use consts
        implicit none
        contains

!       calculates pressure from conserved variables
        function p(u)
                use consts
                implicit none
                double precision, dimension(4) :: u
                double precision :: p
                p = (gamma-1)*(u(4)-0.5*(u(2)**2 + u(3)**2)/u(1))
        end function p

!       calculates speed of sound from conserved variables
        function a(u)
                use consts
                implicit none
                double precision, dimension(4) :: u
                double precision :: a
                a = sqrt(gamma*p(u)/ u(1))
        end function a

!       calculated x component of the Euler flux
        function fEuler(u)
                use consts
                implicit none
                double precision, dimension(4) :: u
                double precision, dimension(4) :: fEuler
                double precision :: pu
                pu = p(u)
                fEuler(1) = u(2)
                fEuler(2) = u(2)**2/u(1) + pu
                fEuler(3) = u(2)*u(3)/u(1)
                fEuler(4) = u(2)/u(1)*(u(4) + pu)
        end function fEuler

!       calculated y component of Euler flux
        function gEuler(u)
                use consts
                implicit none
                double precision, dimension(4) :: u
                double precision, dimension(4) :: gEuler
                double precision :: pu
                pu = p(u)
                gEuler(1) = u(3)
                gEuler(2) = u(2)*u(3)/u(1)
                gEuler(3) = u(3)**2/u(1) + pu
                gEuler(4) = u(3)/u(1)*(u(4) + pu)
        end function gEuler

!       Calcualtes Lax-Friedrichs flux in x
        function fEulerLF(ul, ur)
                double precision, dimension(4) :: fEulerLF
                double precision, dimension(4), intent(in) :: ul, ur
                fEulerLF = (fEuler(ul) + fEuler(ur) + (ul - ur)/(2*lambda))/2
        end function fEulerLF

!       Lax-Friedrichs flux in y
        function gEulerLF(ub, ut)
                double precision, dimension(4) :: gEulerLF
                double precision, dimension(4), intent(in) :: ub, ut
                gEulerLF = (gEuler(ub) + gEuler(ut) + (ub - ut)/(2*lambda))/2
        end function gEulerLF

!       calculates the maximum speed for the conection of two given states
        function cmax(ul, ur)
                double precision, dimension(4) :: ul, ur
                double precision :: cmax
                double precision, dimension(2) :: vl, vr
                vl = ul(2:3) / ul(1)
                vr = ur(2:3) / ur(1)
                cmax = max(norm2(vl), norm2(vr)) + max(a(ul), a(ur))

        end function cmax

!       Calcualtes local Lax-Friedrichs flux in x
        function fEulerLLF(ul, ur)
                double precision, dimension(4) :: fEulerLLF
                double precision, dimension(4), intent(in) :: ul, ur
        
                fEulerLLF = (fEuler(ul) + fEuler(ur) + cmax(ul, ur)*(ul - ur))/2
        end function fEulerLLF

!       local Lax-Friedrichs flux in y
        function gEulerLLF(ub, ut)
                double precision, dimension(4) :: gEulerLLF
                double precision, dimension(4), intent(in) :: ub, ut

                gEulerLLF = (gEuler(ub) + gEuler(ut) + cmax(ub, ut)*(ub - ut))/2
        end function gEulerLLF



!       Entropy conservative two-point flux for the Euler equations of Gas-Dynamics
!       out of ARBITRARILY HIGH-ORDER ACCURATE ENTROPY STABLEESSENTIALLY NONOSCILLATORY SCHEMES FOR SYSTEMS OFCONSERVATION LAWS
!       By Fjordhlm, Mishra and Tadmor

        function fEulerCons(ul, ur)
                implicit none
                double precision, dimension(4) :: fEulerCons
                double precision, dimension(4), intent(in) :: ul, ur
                double precision :: pl, pr
                double precision :: betal, betar
                pl = p(ul)
                pr = p(ur)
                betal = ul(1)/(2*pl)
                betar = ur(1)/(2*pr)
                fEulerCons(1) = logmean(ul(1), ur(1))*meanv(ul(2)/ul(1), ur(2)/ur(1))
                fEulerCons(2) = meanv(ul(1), ur(1)) / (2*meanv(betal, betar)) + meanv(ul(2)/ul(1), ur(2)/ur(1))*fEulerCons(1)
                fEulerCons(3) = meanv(ul(3)/ul(1), ur(3)/ur(1))*fEulerCons(1)
                fEulerCons(4) = (1/(2*(gamma-1)*logmean(betal, betar))&
                                 - meanv((ul(2)/ul(1))**2+(ul(3)/ul(1))**2, (ur(2)/ur(1))**2 + (ur(3)/ur(1))**2)/2) &
                                *fEulerCons(1) + dot_product(meanv(ul(2:3)/ul(1), ur(2:3)/ur(1)), fEulerCons(2:3))
        end function fEulerCons

!       Version in y direction
        function gEulerCons(ub, ut)
                implicit none
                double precision, dimension(4) :: gEulerCons
                double precision, dimension(4), intent(in) :: ub, ut
                double precision :: pb, pt
                double precision :: betab, betat
                pb = p(ub)
                pt = p(ut)
                betab = ub(1)/(2*pb)
                betat =  ut(1)/(2*pt)
                gEulerCons(1) = logmean(ub(1), ut(1))*meanv(ub(3)/ub(1), ut(3)/ut(1))
                gEulerCons(2) =  meanv(ub(2)/ub(1), ut(2)/ut(1))*gEulerCons(1)
                gEulerCons(3) = meanv(ub(1), ut(1)) / (2*meanv(betab, betat)) + meanv(ub(3)/ub(1), ut(3)/ut(1))*gEulerCons(1)
                gEulerCons(4) = (1/(2*(gamma-1)*logmean(betab, betat))&
                                - meanv((ub(2)/ub(1))**2+(ub(3)/ub(1))**2, (ut(2)/ut(1))**2 + (ut(3)/ut(1))**2)/2)&
                                *gEulerCons(1) + dot_product(meanv(ub(2:3)/ub(1), ut(2:3)/ut(1)), gEulerCons(2:3))
        end function gEulerCons

!       Implementation of mean value 
        elemental function meanv(a, b)
                implicit none
                double precision :: meanv
                double precision, intent(in) :: a, b
                meanv = (a + b) / 2.0D0
        end function meanv

!       Implementation of logarithmic mean
        elemental function logmean(ul, ur)
                implicit none
                double precision, intent(in) :: ul, ur
                double precision :: logmean
                double precision, parameter :: eps = 10.0D-2
                double precision :: zeta
                double precision :: f
                double precision :: u
                double precision :: g
                zeta = ul/ur
                f = (zeta - 1)/(zeta + 1)
                u = f**2

                if (u < eps) then
                        g = 1 + u/3 + u**2/ 5 + u**3/7
                else
                        g = log(zeta)/(2*f)
                end if
                logmean = (ul + ur)/(2*g)
        end function logmean

!       Calculates Physical Entropy for the Euler System
        function PUEuler(u)
                implicit none
                double precision, dimension(ncons), intent(in) :: u
                double precision :: PUEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PUEuler = -u(1)*S
        end

!       Calculates entropy flux in x direction
        function PFEuler(u)
                implicit none
                double precision, dimension(ncons), intent(in) :: u
                double precision :: PFEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PFEuler = -u(2)*S
        end

!       Entropy flux in y direction
        function PGEuler(u)
                implicit none
                double precision, dimension(ncons), intent(in) :: u
                double precision :: PGEuler
                double precision :: S
                S = log(p(u)*u(1)**(-gamma))
                PGEuler = -u(3)*S
        end

!       Lax-Friedrichs numerical entropy flux
        function PFEulerLF(ul, ur)
                implicit none
                double precision, dimension(ncons), intent(in) :: ul, ur
                double precision :: PFEulerLF
                double precision :: S
                PFEulerLF = 0.5D0*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur)-PUEuler(ul))/(2*lambda))
        end

!       Lax-Friedrichs numerical entropy flux in y
        function PGEulerLF(ub, ut)
                implicit none
                double precision, dimension(ncons), intent(in) :: ub, ut
                double precision :: PGEulerLF
                double precision :: S
                PGEulerLF = 0.5D0*(PGEuler(ub) + PGEuler(ut) - (PUEuler(ut) - PUEuler(ub))/(2*lambda))
        end

end module fluxes
