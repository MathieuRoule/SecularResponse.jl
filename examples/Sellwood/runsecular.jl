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
@time LinearResponse.RunLinearResponse(ψ,dψ,d2ψ,d3ψ,d4ψ,ndFdJ,FHT,basis,Linearparams)

#####
# Secular response
#####

coupling = SecularResponse.BalescuLenardCouplingCreate(basis,FHT,Linearparams)
# coupling = SecularResponse.LandauBasisCouplingCreate(basis)

Jmin, Jmax, nJ = 0.01, 0.3, 4
Lmin, Lmax, nL = 0.3, 1.5, 4

const tabJL = SecularResponse.Grid2D(Jmin,Jmax,nJ,Lmin,Lmax,nL)


totfric, totdiff, totflux, totdFdt = SecularResponse.GetSecular(tabJL,ψ,dψ,d2ψ,d3ψ,d4ψ,DF,ndFdJ,coupling,Secularparams)