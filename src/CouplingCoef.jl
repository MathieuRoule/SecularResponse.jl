

# make an abstract coupling type
# @WARNING: Should be defined before any coupling definition
abstract type AbstractCoupling end

include("coupling/LandauMultipole.jl")
include("coupling/LandauBasis.jl")
include("coupling/BalescuLenard.jl")
include("coupling/BLPdS.jl")