"""

all steps combined into one: could pause and restart if this was too much.


"""

using HDF5
using BenchmarkTools

inputfile = "SellwoodStable.jl"
include(inputfile)

#####
# Linear Response
#####
# # Construct W matrices (basis FT)
# CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,d4ψ,FHT,basis,params)
    
# # Compute G(u) functions
# CallAResponse.RunGfunc(ndFdJ,FHT,params)

# # Compute matrix decomposition coefficient 
# CallAResponse.MakeaMCoefficients(FHT,params)


#####
# Secular response
#####

coupling = SecularResponse.BalescuLenardCouplingCreate(basis,FHT,params)
# coupling = SecularResponse.LandauBasisCouplingCreate(basis)

Jmin, Jmax, nJ = 0.1, 0.1, 5
Lmin, Lmax, nL = 1.2, 1.2, 5

tabJL = SecularResponse.Grid2D(Jmin,Jmax,nJ,Lmin,Lmax,nL)

@time tabJLcomp, totfric, totdiff, totflux = SecularResponse.GetSecular(tabJL,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,SRparams)

println(totflux[:,1])

# filename = "Jldomain.hf5"
# h5open(filename,"w") do file
#     write(file,"tabJL",tabJL)
#     write(file,"totfric",totfric)
#     write(file,"totdiff",totdiff)
#     write(file,"totflux",totflux)
# end