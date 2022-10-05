

include("coupling/LandauMultipole.jl")
include("coupling/LandauBasis.jl")
include("coupling/BalescuLenard.jl")
include("coupling/BLPdS.jl")

# make a generic Basis data type
CouplingType = Union{BalescuLenardCoupling,LandauBasisCoupling,LandauMultipoleCoupling,BLPdSCoupling}