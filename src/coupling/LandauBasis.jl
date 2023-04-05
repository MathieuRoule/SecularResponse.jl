


struct LandauBasisCoupling{BT<:AB.AbstractAstroBasis} <: AbstractCoupling

    name::String 

    basis::BT
    
    UFT::Vector{Float64}
    UFTp::Vector{Float64}
end

function LandauBasisCouplingCreate(basis::BT;
                                   name::String="LandauBasis") where {BT<:AB.AbstractAstroBasis}

    return LandauBasisCoupling(name,basis,
                               zeros(Float64,basis.nmax),zeros(Float64,basis.nmax))
end


function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::ComplexF64,
                    ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                    coupling::LandauBasisCoupling,
                    Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k,J) 
    LR.WBasisFT(a,e,Ω1,Ω2,k1,k2,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,Linearparams)
end

"""
    LandauBasisCC(k,k',J,J',ω)

Landau coupling coefficient using the basis elements
"""
function CouplingCoefficient(a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1::Int64,k2::Int64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ω::ComplexF64,
                             ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                             coupling::LandauBasisCoupling,
                             Linearparams::LR.LinearParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}
    """
    @ASSUMING the (k,J) part has been prepared
    """

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k',J')
    LR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFTp,Linearparams)

    res = 0.0
    for j = 1:Linearparams.nradial
        res -= coupling.UFT[j] * coupling.UFTp[j]
    end
    return res
end