include("NBCtest.jl")

# File prints out a given set of matrices as fortran code
# module name is also the name of the matrix.
# Todo: severall matrices in one module
function FortPrintMat(Ar, modname)
        # write a module
        ax = axes(Ar)
        dim = length(ax)
        
        println("module ", modname, "param")
        println("       implicit none")
        # write an arraz with the correcgt size
        print(  "       double precision, dimension(")
        for k=1:dim
                print(minimum(ax[k]), ":", maximum(ax[k]))
                if k < dim
                        print(", ")
                end     
        end
        println(") ::  ", modname)
        #write the contents
        println("contains")
        println("subroutine ", modname, "init()")

        ind = zeros(Int, dim)
        for a in CartesianIndices(Ar)
                for k=1:dim
                        ind[k] = a[k]
                end
                println("       ", modname, Tuple(ind), " = ", Ar[a])
        end

        println("end subroutine ", modname, "init")
        println("end module ", modname, "param")
end

p = 3
M = MatSet(p, false)
Mcomb = zeros(2*p+1, 2*p+1, 2*p+1)
for k=1:2*p+1
        if k-p-1 != 0
                Mcomb[k, :, :] = OffsetArrays.no_offset_view(M[k-p-1])
        else
                Mcomb[k, 2:end, 2:end] = OffsetArrays.no_offset_view(M[k-p-1])
        end
end

FortPrintMat(Float64.(Mcomb), "Fluxmats")
