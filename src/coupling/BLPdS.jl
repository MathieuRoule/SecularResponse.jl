

struct BLPdSCoupling{BT<:AB.AbstractAstroBasis} <: AbstractCoupling

    name::String 

    LandauMultipole::LandauMultipoleCoupling
    LandauBasis::LandauBasisCoupling{BT}
    BalescuLenard::BalescuLenardCoupling{BT}
end

function BLPdSCoupling(LandauMultipole::LandauMultipoleCoupling,
                       LandauBasis::LandauBasisCoupling{BT},
                       BalescuLenard::BalescuLenardCoupling{BT};
                       name::String="BL-PdS") where {BT<:AB.AbstractAstroBasis}

    return BLPdSCoupling(name,LandauMultipole,LandauBasis,BalescuLenard)
end


function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::ComplexF64,
                    ψ::F0,dψ::F1,d2ψ::F2,
                    coupling::BLPdSCoupling,
                    Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,coupling.LandauMultipole,Linearparams)
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,coupling.LandauBasis,Linearparams)
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,coupling.BalescuLenard,Linearparams)
end

function CouplingCoefficient(ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ψ::F0,dψ::F1,d2ψ::F2,
                             coupling::BLPdSCoupling,
                             Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    psiMulti = CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,coupling.LandauMultipole,Linearparams)
    psiBasis = CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,coupling.LandauBasis,Linearparams)
    psidBL   = CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,coupling.BalescuLenard,Linearparams)

    return psidBL + (psiMulti - psiBasis) # "Pain de Sucre"
end