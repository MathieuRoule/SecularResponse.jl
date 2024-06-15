


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


function CCPrepare!(
    a::Float64,
    e::Float64,
    Ω1::Float64,
    Ω2::Float64,
    k1::Int64,
    k2::Int64,
    lharmonic::Int64,
    ω::ComplexF64,
    model::Potential,
    coupling::LandauBasisCoupling,
    Linearparams::LR.LinearParameters
)

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k,J) 
    function basisft_integrand(r::Float64)
        # collect the basis elements (in place!)
        AB.tabUl!(coupling.basis, Linearparams.lharmonic, r)
        return coupling.basis.tabUl
    end
    _, L = EL_from_ae(a, e, model, Linearparams.Orbitalparams)
    LR.angle_fouriertransform!(coupling.UFT, basisft_integrand, a, e, k1, k2, model, Linearparams, L=L, Ω1=Ω1, Ω2=Ω2)
end

"""
    LandauBasisCC(k,k',J,J',ω)

Landau coupling coefficient using the basis elements
"""
function CouplingCoefficient(
    ap::Float64,
    ep::Float64,
    Ω1p::Float64,
    Ω2p::Float64,
    k1p::Int64,
    k2p::Int64,
    lharmonic::Int64,
    model::Potential,
    coupling::LandauBasisCoupling,
    Linearparams::LR.LinearParameters
)
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
    function basisft_integrand(r::Float64)
        # collect the basis elements (in place!)
        AB.tabUl!(coupling.basis, Linearparams.lharmonic, r)
        return coupling.basis.tabUl
    end
    _, Lp = EL_from_ae(ap, ep, model, Linearparams.Orbitalparams)
    LR.angle_fouriertransform!(coupling.UFTp, basisft_integrand, ap, ep, k1p, k2p, model, Linearparams, L=Lp, Ω1=Ω1p, Ω2=Ω2p)

    res = 0.0
    for j = 1:Linearparams.nradial
        res -= coupling.UFT[j] * coupling.UFTp[j]
    end
    return res
end