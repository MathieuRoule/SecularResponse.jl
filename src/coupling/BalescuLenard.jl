

struct BalescuLenardCoupling{BT<:AB.AbstractAstroBasis,HT<:FHT.AbstractFHT} <: AbstractCoupling

    name::String 

    basis::BT

    UFT::Vector{Float64}
    UFTp::Vector{Float64}

    fht::HT

    aMcoef::Array{Float64,4}
    tabωminωmax::Matrix{Float64}
    M::Matrix{ComplexF64}
    IMat::Matrix{ComplexF64}
    UFTXi::Vector{ComplexF64}

    ξ::Float64 # Active fraction
end

function BalescuLenardCoupling(
    basis::BT,
    fht::HT,
    params::LR.LinearParameters;
    name::String="BalescuLenard",
    ξ::Float64=1.0
) where {BT<:AB.AbstractAstroBasis, HT<:FHT.AbstractFHT}
    nradial = params.nradial
    tabaMcoef, tabωminωmax = LR.StageAXi(params)
    return BalescuLenardCoupling(
        name,
        basis,
        zeros(Float64,nradial),
        zeros(Float64,nradial),
        fht,
        tabaMcoef,
        tabωminωmax,
        zeros(ComplexF64,nradial,nradial),
        Matrix{ComplexF64}(I, nradial, nradial),
        zeros(ComplexF64,nradial),
        ξ
    )
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
    coupling::BalescuLenardCoupling,
    Linearparams::LR.LinearParameters
)

    if lharmonic < 0 
        @assert (Linearparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
        # M^{-l}(ω) = (M^{l}(-ω*))*
        ω = - conj(ω)
    end
    # Computing the basis FT (k,J) 
    function basisft_integrand(r::Float64)
        # collect the basis elements (in place!)
        AB.tabUl!(coupling.basis, Linearparams.lharmonic, r)
        return coupling.basis.tabUl
    end
    _, L = EL_from_ae(a, e, model, Linearparams.Orbitalparams)
    LR.angle_fouriertransform!(coupling.UFT, basisft_integrand, a, e, k1, k2, model, Linearparams, L=L, Ω1=Ω1, Ω2=Ω2)

    # Computing the response matrix
    LR.tabM!(ω,coupling.M,coupling.aMcoef,coupling.tabωminωmax,coupling.fht,Linearparams)

    if lharmonic < 0 
        # M^{-l}(ω) = (M^{l}(-ω*))*
        for j = 1:Linearparams.nradial
            for i = 1:Linearparams.nradial
                coupling.M[i,j] = conj(coupling.M[i,j])
            end
        end
    end

    tabXi = inv(Symmetric(coupling.IMat - coupling.ξ * coupling.M))
    for j = 1:Linearparams.nradial
        xiloc = 0.0 + im*0.0
        for i = 1:Linearparams.nradial
            xiloc += coupling.UFT[i] * tabXi[i,j]
        end
        coupling.UFTXi[j] = xiloc
    end
end


function CouplingCoefficient(
    ap::Float64,
    ep::Float64,
    Ω1p::Float64,
    Ω2p::Float64,
    k1p::Int64,
    k2p::Int64,
    lharmonic::Int64,
    model::Potential,
    coupling::BalescuLenardCoupling,
    Linearparams::LR.LinearParameters
)::ComplexF64

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

    res = 0.0 + im*0.0
    for j = 1:Linearparams.nradial
        res -= coupling.UFTXi[j] * coupling.UFTp[j]
    end
    return res
end