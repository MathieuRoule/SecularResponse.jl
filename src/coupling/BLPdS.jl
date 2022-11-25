

struct BLPdSCoupling

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
                    ω::Float64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                    coupling::BLPdSCoupling,
                    params::Parameters)

    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.LandauMultipole,params)
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.LandauBasis,params)
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.BalescuLenard,params)
end

function CouplingCoefficient(a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1::Int64,k2::Int64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ω::Float64,
                             ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                             coupling::BLPdSCoupling,
                             params::Parameters)

    psiMulti = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.LandauMultipole,params)
    psiBasis = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.LandauBasis,params)
    psidBL   = CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ω,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.BalescuLenard,params)

    return psidBL + (psiMulti - psiBasis) # "Pain de Sucre"
end