

struct BLPdSCoupling <: AbstractCoupling

    name::String 

    LandauMultipole::LandauMultipoleCoupling
    LandauBasis::LandauBasisCoupling
    BalescuLenard::BalescuLenardCoupling
end

function BLPdSCouplingCreate(LandauMultipole::LandauMultipoleCoupling,LandauBasis::LandauBasisCoupling,BalescuLenard::BalescuLenardCoupling;name::String="BL-PdS")

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

function CouplingCoefficient(a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1::Int64,k2::Int64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ω::ComplexF64,
                             ψ::F0,dψ::F1,d2ψ::F2,
                             coupling::BLPdSCoupling,
                             Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    psiMulti = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,coupling.LandauMultipole,Linearparams)
    psiBasis = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,coupling.LandauBasis,Linearparams)
    psidBL   = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,coupling.BalescuLenard,Linearparams)

    return psidBL + (psiMulti - psiBasis) # "Pain de Sucre"
end