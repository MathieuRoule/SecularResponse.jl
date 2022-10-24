"""

all steps combined into one: could pause and restart if this was too much.


"""

inputfile = "SellwoodStable.jl"
include(inputfile)

#####
# Linear Response
#####
# Construct W matrices (basis FT)
CallAResponse.RunWmat(ψ,dψ,d2ψ,d3ψ,FHT,basis,params)
    
# Compute G(u) functions
CallAResponse.RunGfunc(ndFdJ,FHT,params)

# Compute matrix decomposition coefficient 
CallAResponse.MakeaMCoefficients(FHT,params)


#####
# Secular response
#####

Amin, Amax, nA = 0.3, 2.0, 100
Emin, Emax, nE = 0.0, 0.7, 101

tabAE = SecularResponse.AEgrid(Amin,Amax,nA,Emin,Emax,nE)

tabJL, totfric, totdiff, totflux = SecularResponse.GetSecular(tabAE,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,SRparams)


using HDF5

filename = "Jldomain.hf5"
h5open(filename,"w") do file
    write(file,"tabJL",tabJL)
    write(file,"totfric",totfric)
    write(file,"totdiff",totdiff)
    write(file,"totflux",totflux)
end
