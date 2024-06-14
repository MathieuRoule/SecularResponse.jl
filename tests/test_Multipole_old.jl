"""

all steps combined into one: could pause and restart if this was too much.


"""

using HDF5
using BenchmarkTools

inputfile = "SellwoodStable_old.jl"
include(inputfile)

#####
# Coupling structures
#####
couplingLM = SecularResponse.LandauMultipoleCoupling(Kw=1000,ε=0.0)
couplingLB = SecularResponse.LandauBasisCoupling(basis)
couplingBL = SecularResponse.BalescuLenardCoupling(basis,FHT,Linearparams)

#####
# Single coupling coefficient
#####

k1, k2 = -1, 2
a, e = 1.0, 0.05
Ω1, Ω2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e,Orbitalparams)

k1p, k2p = 0, 2
ap, ep = 1.0, 0.1
Ω1p, Ω2p = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,ap,ep,Orbitalparams)

ωres = k1*Ω1 + k2*Ω2

SecularResponse.CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,0.0+0.0*im,ψ,dψ,d2ψ,couplingLM,Linearparams)
ccLM = SecularResponse.CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,couplingLM,Linearparams)

SecularResponse.CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,0.0+0.0*im,ψ,dψ,d2ψ,couplingLB,Linearparams)
ccLB = SecularResponse.CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,couplingLB,Linearparams)

println("Landau-Multipole : psi = ",ccLM)
println("Landau-Basis : psi = ",ccLB)
println("Relative error = ",abs((ccLB-ccLM)/ccLM))

# #####
# # Secular response
# #####
J, L = 0.01, 1.2

k1, k2 = -1, 2
k1p, k2p = 0, 2

# First call, compilation allocations
_, _, dFdtLM = SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingLM,Secularparams)
_, _, dFdtLB = SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingLB,Secularparams)
_, _, dFdtBL = SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingBL,Secularparams)

println("Landau-Multipole : dF/dt = ",dFdtLM)
println("Landau-Basis : dF/dt = ",dFdtLB)
println("Balescu-Lenard : dF/dt = ",dFdtBL)

# # Time the second call
# @btime SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingLM,Secularparams)
# @btime SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingLB,Secularparams)
# @btime SecularResponse.GetSecularResContrib(J,L,k1,k2,k1p,k2p,lharmonic,ψ,dψ,d2ψ,DF,ndFdJ,couplingBL,Secularparams)

