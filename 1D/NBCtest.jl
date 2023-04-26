# Implementation of the boundary aware flux matrix construction procedure
# Numerical tests in 1D are also there.


using OffsetArrays
OA = OffsetArray

using LinearAlgebra
using GenericLinearAlgebra
using DifferentialEquations
using Plots
using LaTeXStrings
pyplot()

include("eno.jl")
disp = 1.0
setprecision(10000)

# Tadmors symmetric entropy conservative two-point flux for Burgers equation
function hBurgersCons(ul, ur)
        return (ul^2 + ul*ur + ur^2)/6
end

f(u) = u^2/2
U(u) = u^2/2
dU(u) = u
F(u) = u^3/3

# Equation 4.10b, Tadmor 1987
function HBurgersCons(ul, ur)
        g = hBurgersCons(ul, ur)
        vl, vr = dU(ul), dU(ur)
        G =  0.5*(vl + vr)*g + 0.5*(F(ul) + F(ur)) - 0.5*(vl*f(ul)) - 0.5*vr*f(ur)
end

function hLinCons(ul, ur)
	return (ul + ur)/2
end

# Lax-Friedrichs fluxes for reference
function hBurgersLF(ul, ur, lambda)
        return 0.5*((ul*ul+ur*ur)/2.0 - (ur-ul)/lambda)
end

function hLinLF(ul, ur, lambda)
        return 0.5*((ul+ur) - (ur-ul)/lambda)
end

# Riemann solvers 
function Riemann1D(xdbt, ul, ur, a, ainv, f)
        s = (f(ul) - f(ur))/(ul-ur)
        if ul > ur
                if xdbt < s
                        return ul
                else
                        return ur
                end
        else
                if xdbt <= a(ul)
                        return ul
                elseif xdbt >= a(ur)
                        return ur
                else
                        return ainv(xdbt)
                end
        end
end

# Godunov fluxes for reference
function RiemannBurg(xdbt, ul, ur)
        f(u) = u^2/2
        a(u) = u
        ainv(y) = y
        return Riemann1D(xdbt, ul, ur, a, ainv, f)
end

function hBurgersGod(ul, ur)
        f(u) = u^2/2
        return f(RiemannBurg(0, ul, ur))
end


# Coefficients for the High Order combination developed by LeFloch, Mercier and Rhode
function GetAlphaVec(k)
        b = zeros(BigFloat, k)
        b[1] = 1.0
        A = zeros(BigFloat,k,k)
        for r=1:k
                A[1, r] = r
        end
        for s=2:k
                for r = 1:k
                        A[s, r] = r^(2*s-1)
                end
        end
        alpha = A \ b
        return alpha
end


alphastore = [GetAlphaVec(1), GetAlphaVec(2), GetAlphaVec(3), GetAlphaVec(4), GetAlphaVec(5), GetAlphaVec(6), GetAlphaVec(7)]

# Calculates linear combination Coefficient matrix from Lefloch coefficient vector
# 2p is the width of the flux stencil.
function MakeAMatrix(p)
        A = OA(zeros(BigFloat, 2*p, 2*p), -p+1:p, -p+1:p)
        for r = 1:p
                for q = 0:r-1
                        A[-q, r-q] = A[-q, r-q] + (alphastore[p])[r]
                end
        end
        return 0.5*(A + transpose(A))
end

# Calculates vector $v$ as the difference between A and B
# 2p+1 is the width of the difference stencil
# z is the amount of points left of the cell, whose right edge should be calculated.
function OFV(p, z)
        b = zeros(BigFloat, 2*p+1)
        A = zeros(BigFloat, 2*p+1, 2*p+1)
        b[2] = 2.0
	for j=1:min(6, 2*p)
                for n=1:2*p+1
                        A[j, n] = (n-z-1)^(j-1)
                end
        end

	#b[2*p-3] = +1.0
	A[2*p+1, z+1] = 1.0
	beta = (A'*A) \ A'b
        return  OA(beta, -z:2p-z)
end

# converts v back to a matrix
function OFM(p, z)
	v = OFV(p, z)
	C = OA(zeros(BigFloat, 2*p+1,2*p+1), -z:2*p-z, -z:2*p-z)
	for k = -z:2*p-z
		C[0, k] = C[0, k] + v[k]/2
		C[k, 0] = C[k, 0] + v[k]/2
	end
	return C
end

function OFMw(p, z)
	v = OFV(p, z)
	C = OA(zeros(4*p, 4*p), -z+1-p:3p-z, -z+1-p:3p-z)
	for k = -z+1:2*p-z
		C[0, k] = C[0, k] + v[k]/2
		C[k, 0] = C[k, 0] + v[k]/2
	end
	return C
end

# T shifting operator to go one step furhter to the right boundary
function Tr(A)
	return OA(A, -1, -1)
end

# Embedding A into a stencil one point bigger to the left
function Er(A)
	kmax, lmax = size(OffsetArrays.no_offset_view(A))
	x, y = axes(A)
	BOA = OA(zeros(BigFloat, kmax+1, lmax+1), minimum(x)-1:maximum(x), minimum(y)-1:maximum(x))
	BOA[x, y] = A
	return BOA
end

# Embedding A into a stencil one point bigger to the left
function El(A)
	kmax, lmax = size(OffsetArrays.no_offset_view(A))
        x, y = axes(A)
        BOA = OA(zeros(BigFloat, kmax+1, lmax+1), minimum(x):maximum(x)+1, minimum(y):maximum(x)+1)
        BOA[x, y] = A
        return BOA
end

function Tl(A)
	return OA(A, 1, 1)
end

# Generates the complete set of matrices needed for the domain. 
function MatSet(p, doublify=true)
	A = MakeAMatrix(p)
	Mats = OA(Array{Any}(undef, 2p+1), -p:p)
	Mats[0] = A
	z = p
	localA = Tr(El(A))
	for r = 0:p-1
		localz = z + r
		OFVv = OFM(p, localz)
		Mats[r+1] = localA + OFVv
		localA = Tr(Mats[r+1])
	end
	localA = Er(A)
	for r = 0:-1:-p+1
		localz = z + r
		OFVm = OFM(p, localz)
		Mats[r-1] = Tl(localA - OFVm)
		localA = Mats[r-1]
	end
	
	if doublify
		# only after all calc we recalc as Float 64
		for s = -p:p
			Mats[s] = Float64.(Mats[s])
		end
	end
	return Mats
end

# generates »wide« stencil mat set with at least two stencils of each kind.
# 4p point wide
function WMatSet(p)
        A = OA(zeros(4*p, 4*p), -2p+1:2p, -2p+1:2p)
        A[-p+1:p, -p+1:p] = MakeAMatrix(p)
        Mats = OA(Array{Any}(undef, 4p-3), -2p+2:2p-2)
        Mats[0] = A
        z = p
        localA = Tr(A)
	c = 1
        for r = 1:p-1
		for s=1:2
                	localz = z + r
                	OFVm = OFMw(p, localz)
                	Mats[c] = localA + OFVm
                	localA = Tr(Mats[c])
			c = c+1
		end
        end
        localA = A
	c = 0
        for r = 0:-1:-p+2
		for s=1:2
                	localz = z + r
                	OFVm = OFMw(p, localz)
                	Mats[c-1] = Tl(localA - OFVm)
                	localA = Mats[c-1]
			c = c-1
		end
        end
        return Mats

end

function CalcFlux(u, k, M, p, lam; fnum = hBurgersCons, fdisp = hBurgersLF)
	N = size(u)[1]
	if k + p >= N
		ti = k-N + p+1
	elseif k-p <= 1
		ti = k-p-1
	else
		ti = 0
	end
	A = M[ti]
	lax, max = axes(M[ti])
	f = 0
	alp = 0.00#1
	for l in lax, m in max
		if m > l
			f = f + A[l, m]*(alp*fdisp(u[k+l], u[k+m], lam) + (1-alp)*fnum(u[k+l], u[k+m]))
		end
	end
	return 2*f
end

function CalcEFlux(u, k, M, p, lam; fnum = HBurgersCons)
        N = size(u)[1]
        if k + p >= N
                ti = k-N + p+1
        elseif k-p <= 1
                ti = k-p-1
        else
                ti = 0
        end
        A = M[ti]
        lax, max = axes(M[ti])
        f = 0
        
        for l in lax, m in max
                if m > l
                        f = f + A[l, m]*(fnum(u[k+l], u[k+m]))
                end
        end
        return 2*f
end


function BurgersFD!(du, u, p, t)
        N = length(u)
        dx = p[1]
        k = ceil(Int, p[2])
        ul = p[3]
        ur = p[4]
        M = p[5]
	dt = p[6]
	unew = zeros(N+2)
	unew[2:N+1] = u[:]
        unew[1] = ul(t)
        unew[N+2] = ur(t)
        for i = 2:N+1
                du[i-1] = (CalcFlux(unew, i-1, M, k, dt/dx) - CalcFlux(unew, i, M, k, dt/dx))/dx
        end
end

function CheckEIQ(u, p, t)
        N = length(u)
        du, dE = zeros(N), zeros(N)
        BurgersFD!(du, u, p, t)
        
        dx = p[1]
        k = ceil(Int, p[2])
        ul = p[3]
        ur = p[4]
        M = p[5]
        dt = p[6]
        unew = zeros(N+2)
        unew[2:N+1] = u[:]
        unew[1] = ul(t)
        unew[N+2] = ur(t)
        for i = 2:N+1
                dE[i-1] = dU(u[i-1])*du[i-1] - (CalcEFlux(unew, i-1, M, k, dt/dx) - CalcEFlux(unew, i, M, k, dt/dx))/dx
        end
        return dE
end

function PlotEIQ(xar, sol)
        N = length(sol(0.0))
        K = length(sol.t)
        ar = zeros(K, N)
        for k=1:K
                ar[k, :] = CheckEIQ(sol(sol.t[k]), sol.prob.p, sol.t[k])
        end
        heatmap(xar, sol.t, ar, c=cgrad([:black, :white]), xlabel="x", ylabel="t")
end

function BurgersENO2!(du, u, p, t)
	N = length(u)
        dx = p[1]
        k = ceil(Int, p[2])
        ul = p[3]
        ur = p[4]
        M = p[5]
	dt = p[6]
        unew = zeros(N+2)
	ulv = zeros(N+2)
	urv = zeros(N+2)
        unew[2:N+1] = u[:]
        unew[1] = ul(t)
        unew[N+2] = ur(t)
	for i=2:N+1
		ulv[i], urv[i] = classicENO2(unew[i-1:i+1])
	end
	ulv[1] = unew[1]
	urv[1] = 0.5*(unew[1] + unew[2])
	ulv[N+2] = 0.5*(unew[N+2] + unew[N+1])
	urv[N+2] = unew[N+2]

        for i = 2:N+1
		du[i-1] = (hBurgersLF(urv[i-1], ulv[i], dt/dx) - hBurgersLF(urv[i], ulv[i+1], dt/dx))/dx
        end
end



function numericaltest(N, tmax, u0func, xmax, Solver;order=2, CFL=0.5, ul=tzt, ur=tzt, DO = 0.0)
        dx = xmax / N
        dt = CFL*dx
        println("dx = " * string(dx))
        println("dt = " * string(dt))
        u0 = u0func.(collect(1:N) .* dx)
        tspan = (0, tmax)
	prob = ODEProblem(Solver, u0, tspan, [dx, order, ul, ur, MatSet(order), dt, DO])
        sol = solve(prob, SSPRK33(), dt=dt)
        return sol
end

# Tests the boundary fluxes
function BCtest(N, tmax, order, CFL, u0f = x->1.0, ul=t->1.0 + 0.02*exp(-(t-5)^2), ur=t->1.0)
	dx = 20.0 / (N+1)
        dt = CFL*dx
        println("dx = " * string(dx))
        println("dt = " * string(dt))
	u0 = zeros(N)
	xar = zeros(N)
	for k=1:N
		xar[k] = -10 + dx*k
		u0[k] = u0f(xar[k])
	end
        tspan = (0, tmax)
        
        prob = ODEProblem(BurgersFD!, u0, tspan, [dx, order, ul, ur, MatSet(order), dt])
        sol = solve(prob, SSPRK33(), dt=dt)
        return xar, sol
end

# Used for reference purposes
function BCRef(N, tmax, order, CFL)
	dx = 20.0 / (N+1)
        dt = CFL*dx
        println("dx = " * string(dx))
        println("dt = " * string(dt))
        u0 = zeros(N)
	xar = zeros(N)
        for k=1:N
                u0[k] = 1.0 #-sin(pi*(-10 + dx*k)/20)
		xar[k] = -10 + dx*k

        end
	ul(t) = 1.0 + 0.02*exp(-(t-5)^2)
        ur(t) = 1.0
        tspan = (0, tmax)
        prob = ODEProblem(BurgersENO2!, u0, tspan, [dx, order, ul, ur, 0, dt])
        sol = solve(prob, SSPRK33(), dt=dt)
        return xar, sol
end

# Piecewise constant interpolation
function Constf(Grid, uar, x)
	n = length(Grid)
        lb = 1
        ub = n
        im = 1
        while ub - lb > 1
                im = floor(Int, 0.5*(ub + lb))
                if Grid[im] < x
                        lb = im
                else
                        ub = im
                end
        end
        if abs(x - Grid[lb]) < abs(x - Grid[ub])
		return uar[lb]
	else
		return uar[ub]
	end
end

# Linear interpolation 
function Interp(Grid, Array, x)
	n = length(Grid)
	lb = 1
	ub = n
	im = 1
	while ub - lb > 1
		im = floor(Int, 0.5*(ub + lb))
		if Grid[im] < x
			lb = im
		else
			ub = im
		end
	end
	alp = (x - Grid[lb])/(Grid[ub] - Grid[lb])
	return alp*Array[ub] + (1-alp)*Array[lb]
end


# Function carries out the numerical tests for a fixes order, i.e. produces some optical
# tests of the entropy conservative boundary fluxes and a numerical converge analysis.
function NTBC(order)
	u0(x) = sinpi(x)
	tz(x, t) = 0
	tzt(t) = 0

	xar, sol = BCtest(100, 10.0, order, 0.25, x->-sin(pi*(x)/20), t->0.9 + 0.1*cos(0.5*pi*t)
, t->-(0.9 + 0.1*cos(0.5*pi*t)))
	for t in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
		scatter(xar, sol(t), xlabel="x", ylabel="u", label=L"u^n_k")
		savefig("BCtest"*string(Int(t))*"p"*string(order)*".pdf")
	end
        PlotEIQ(xar, sol)
        savefig("BCtestEIQ"*string(order)*".pdf")

	Ncont = 8192*2
	maxexp = 8
	minexp = 6
	steps = 8
	Nmin = 2^(minexp)
	Nmax = 2^(maxexp)
	ttest = 20.0
	eps = 0.0001
	N(k) = floor(Int, 2^(minexp + (maxexp-minexp)*(k-1)/(steps-1)))
	refxar, refsol = BCRef(Ncont, ttest, 2, 0.8)
	norms = zeros(steps)
	for k=1:steps
		Nloc = N(k)
		xar, sol = BCtest(Nloc, ttest, order, 0.25*64/Nloc)
		projref(x) = Interp(refxar, refsol(ttest), x)
		fh = 1:Nloc
		norms[k] = norm(projref.(xar[fh]) .- sol(ttest)[fh], 1)/Nloc
	end
	scatter(N.(1:steps), norms, yscale=:log10, xscale=:log10, xlabel="N", ylabel=L"L_1 error", color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(maxexp-minexp)], label="Order 1", linestyle=:dot, color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(2*(maxexp-minexp))], label="Order 2", linestyle=:dash, color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(3*(maxexp-minexp))], label="Order 3", linestyle=:solid, color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(4*(maxexp-minexp))], label="Order 4", linestyle=:dot, color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(5*(maxexp-minexp))], label="Order 5", linestyle=:dash, color=:black)
	plot!([2^minexp, 2^maxexp], [norms[1], norms[1]/2^(6*(maxexp-minexp))], label="Order 6", linestyle=:solid, color=:black)
	savefig("BCconvanap"*string(Int(order))*".pdf")
end


#functions for eigenvalue analysis of discrete operators
function LinMat(N, p)
	Dm = zeros(N+2, N+2)
	M = MatSet(p)
	for k=1:N+2
		u = zeros(N+2)
		u[k] = 1
		for l=2:N+1
                	Dm[l, k] = (CalcFlux(u, l-1, M, p, 0.5, fnum=hLinCons) - CalcFlux(u, l, M, p, 0.5, fnum=hLinCons, fdisp=hLinLF))
       		end
	end
	return Dm
end

function DiagEW(p)
	A = LinMat(128, p)
	v, w = eigen(A)
	scatter(v, label="p = "*string(p))
	savefig("EWsp"*string(p)*".pdf")
end

function ProjDm(N, p)
	Dm = LinMat(N, p)
	w, v = eigen(Dm)

end


