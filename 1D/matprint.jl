include("NBCtest.jl")
using Latexify


for p=2:3
	M = MatSet(p, false)
	for q=-p:p
		print(latexify(rationalize.(M[q])))
      	end
end


