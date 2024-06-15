module SecularResponse

# bring in the external dependencies
import AstroBasis as AB
import FiniteHilbertTransform as FHT
using HDF5
import LinearResponse as LR
using LinearAlgebra
using OrbitalElements
using ProgressMeter
using SpecialFunctions

include("utils/ResonancesPair.jl")
include("utils/Parameters.jl")
include("utils/Grids.jl")
include("utils/IO.jl")

include("CouplingCoef.jl")
include("utils/ResonanceLines.jl")

include("Prediction.jl")

end # module
