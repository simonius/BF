@inline minmod(a::Float64, b::Float64) = min(abs(a), abs(b))*(sign(a) + sign(b))/2.0 :: Float64
MUSCLr(u) = u[2, :] + 0.5*minmod.(u[2, :] - u[1, :], u[3, :] - u[2, :])
MUSCLl(u) = u[2, :] - 0.5*minmod.(u[2, :] - u[1, :], u[3, :] - u[2, :])



# Performs classic ENO reconstruction of order 2, u is the vector [u_{i-1}, u_i, u_{i+1}]
# returns mu, du/dx
function reconENO2(u)
        a1 = u[2] - u[1]
        a2 = u[3] - u[2]
        if a1^2 < a2^2
                return u[2], a1
        else
                return u[2], a2
        end
end

function classicENO2(u)
        mean, a = reconENO2(u)
        return mean - a/2, mean + a/2
end

function ENO2Mat(u)
	N, K = size(u)
	ul = zeros(N, K)
	ur = zeros(N, K)
	for i=2:N-1
		for k = 1:K
			ul[i, k], ur[i, k] = classicENO2(u[i-1:i+1, k])
		end
	end
	return ul, ur
end


function reconRC2(u)
	a1 = u[2] - u[1]
	a2 = u[3] - u[2]
#	if rand() < 0.5
#		return u[2] - a1/2, u[2] + a1/2
#	else
#		return u[2] - a2/2, u[2] + a2/2
#	end
	return u[2] - (a1 + a2)/8, u[2] +(a1 + a2)/8
end

function RC2Mat(u)
	 N, K = size(u)
        ul = zeros(N, K)
        ur = zeros(N, K)
        for i=2:N-1
                for k = 1:K
                        ul[i, k], ur[i, k] = reconRC2(u[i-1:i+1, k])
                end
        end
        return ul, ur
end

function LinU(ul, ur, U)
#	lam = range(0.0, stop=1.0, length=2^4+1)
	S(x) = U(x*ul+(1-x)*ur)
#	val, err = romberg(lam, S.(lam))
	val, err = quadgk(S, 0.0, 1.0)
	return val
end


# Entropy aware reconstruction. If reconstruction I and II are out of domain we return the mean cell value
# If only reconstruction I is out of domain we return reconstruction II and vice versa
# if both are in the domain we return the one with lower total entropy
function EMR2(u, U)
	a1 = u[2, :] - u[1, :]
	a2 = u[3, :] - u[2, :]
	U1 = 0
	U2 = 0
	try
		U1 = LinU(u[2, :] - a1/2, u[2, :] + a1/2, U)
	catch e
		try 
			U2 = LinU(u[2, :] - a2/2, u[2, :] + a1/2, U)
		catch
			return u[2, :], u[2, :]
		end
		return u[2, :] - a2/2, u[2, :] + a2/2
	end
	try 
		U2 = LinU(u[2, :] - a2/2, u[2, :] + a1/2, U)
	catch e
		return u[2, :] - a1/2, u[2, :] + a2/2
	end
	if U1 < U2
		return u[2, :] - a1/2, u[2, :] + a1/2
	else
		return u[2, :] - a2/2, u[2, :] + a2/2	
	end
end

function EMR2Mat(u, U)
	   N, K = size(u)
        ul = zeros(N, K)
        ur = zeros(N, K)
        for i=2:N-1
        	ul[i, :], ur[i, :] = EMR2(u[i-1:i+1, :], U) 
        end
        return ul, ur
end

