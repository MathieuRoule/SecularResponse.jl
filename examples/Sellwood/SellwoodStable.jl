

include("./SellwoodStableCAR.jl")

import SecularResponse

const SRparams = SecularResponse.ParametersCreate(CARparams,basis,n1max=1,VERBOSE=0)

const βc(αc::Float64) = OrbitalElements.βcirc(αc,dψ,d2ψ,OEparams)

const DF(E::Float64,L::Float64) = OrbitalElements.mestel_Zang_DF(E,L,R0,Rin,Rout,Rmax,V0,xi,C0,q0,σ0,nu,mu)