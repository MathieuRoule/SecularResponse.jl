

include("./SellwoodStableCAR.jl")

import SecularResponse

SRparams = SecularResponse.ParametersCreate(params,basis,n1max=1,VERBOSE=0)

const βc(αc::Float64) = OrbitalElements.βcirc(αc,dψ,d2ψ,Ω₀,rmin=rmin,rmax=rmax)

const DF(E::Float64,L::Float64) = OrbitalElements.mestel_Zang_DF(E,L;
                                                                 R0=R0,Rin=Rin,Rmax=Rmax,V0=V0,
                                                                 xi=xi,C=C0,
                                                                 q=q0,sigma=σ0,
                                                                 nu=nu,mu=mu)