module SecularResponse

# bring in the external dependencies
import OrbitalElements as OE
import AstroBasis as AB
import FiniteHilbertTransform as FHT
import CallAResponse as CAR
using LinearAlgebra
using HDF5

include("utils/ResonancesPair.jl")
include("utils/Parameters.jl")
include("utils/Grids.jl")

include("CouplingCoef.jl")

include("Prediction.jl")

end # module
