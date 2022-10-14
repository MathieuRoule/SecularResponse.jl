"""

all steps combined into one: could pause and restart if this was too much.


"""

inputfile = "SellwoodStable.jl"
include(inputfile)

a, e = 0.8, 0.8

J, L, totfric, totdiff, totflux = SecularResponse.GetSecular(a,e,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,SRparams)

println("J, L = $J, $L")
println("Friction = ",totfric)
println("Diffusion = ",totdiff)
println("Flux = ",totflux)