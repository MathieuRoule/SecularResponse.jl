

struct LandauMultipoleCoupling <: AbstractCoupling

    name::String 

    tab::Array{Float64}
end

function LandauMultipoleCoupling(K::Int64;name::String="LandauMultipole")
    return LandauMultipoleCoupling(name,zeros(Float64,K))
end

function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::ComplexF64,
                    ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                    coupling::LandauMultipoleCoupling,
                    Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

end

function CouplingCoefficient(a::Float64,e::Float64,
                                Ω1::Float64,Ω2::Float64,
                                ap::Float64,ep::Float64,
                                Ω1p::Float64,Ω2p::Float64,
                                k1::Int64,k2::Int64,
                                k1p::Int64,k2p::Int64,
                                lharmonic::Int64,
                                ω::ComplexF64,
                                ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,
                                coupling::LandauMultipoleCoupling,
                                Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function}

    # TO DO !!!
    return 0.
end