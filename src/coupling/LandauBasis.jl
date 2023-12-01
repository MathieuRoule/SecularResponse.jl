


struct LandauBasisCoupling{BT<:AB.AbstractAstroBasis} <: AbstractCoupling

    name::String 

    basis::BT

    Kw::Int64
    
    UFT::Vector{Float64}
    UFTp::Vector{Float64}
end

function LandauBasisCoupling(basis::BT;
                             name::String="LandauBasis",
                             Kw::Int64=32) where {BT<:AB.AbstractAstroBasis}

    return LandauBasisCoupling(name,basis,Kw,
                               zeros(Float64,basis.nradial),zeros(Float64,basis.nradial))
end


function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::ComplexF64,
                    ψ::F0,dψ::F1,d2ψ::F2,
                    coupling::LandauBasisCoupling,
                    Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k,J) 
    LR.WBasisFT(a,e,Ω1,Ω2,k1,k2,ψ,dψ,d2ψ,coupling.basis,coupling.UFT,Linearparams)
end

"""
    LandauBasisCC(k,k',J,J',ω)

Landau coupling coefficient using the basis elements
"""
function CouplingCoefficient(ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ψ::F0,dψ::F1,d2ψ::F2,
                             coupling::LandauBasisCoupling,
                             Linearparams::LR.LinearParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}
    #####
    # @ASSUMING the (k,J) part has been prepared
    #####

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1p *= -1
        k2p *= -1
    end
    # Computing the basis FT (k',J')
    LR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,ψ,dψ,d2ψ,coupling.basis,coupling.UFTp,Linearparams)

    res = 0.0
    for j = 1:Linearparams.nradial
        res -= coupling.UFT[j] * coupling.UFTp[j]
    end
    return res
end