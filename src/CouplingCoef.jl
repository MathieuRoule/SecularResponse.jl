

include("coupling/LandauMultipole.jl")
include("coupling/LandauBasis.jl")
include("coupling/BalescuLenard.jl")
include("coupling/BLPdS.jl")

"""
    CouplingCoefficient(k,k',J,J',ω)

"""
function CouplingCoefficient(k1::Int64,k2::Int64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             ω::Float64,
                             basis::AB.Basis_type;
                             COUPLING::String="Landau-Basis",
                             VERBOSE::Int64=0)

    if COUPLING == "Landau-Multipole"
        return LandauMultipoleCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p;VERBOSE=VERBOSE)
    elseif COUPLING == "Landau-Basis"
        return LandauBasisCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis;VERBOSE=VERBOSE)
    elseif COUPLING == "Balescu-Lenard"
        return BalescuLenardCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis,ω;VERBOSE=VERBOSE)
    elseif COUPLING == "BL-PdS"
        return BLPdSCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis,ω;VERBOSE=VERBOSE)
    else
        error("Unknow coupling")
    end
end