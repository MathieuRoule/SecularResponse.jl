


struct LandauBasisCoupling

    name::String 

    basis::AB.Basis_type
    
    UFT::Array{Float64}
    UFTp::Array{Float64}
end

function LandauBasisCouplingCreate(basis::AB.Basis_type;name::String="BalescuLenard")

    return LandauBasisCoupling(name,basis,
                               zeros(Float64,basis.nmax),zeros(Float64,basis.nmax))
end


function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::Float64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    coupling::LandauBasisCoupling,
                    params::Parameters)

    # Computing the basis FT (k,J) 
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,lharmonic,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,params.CARparams)
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
                             ω::Float64,
                             ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                             coupling::LandauBasisCoupling,
                             params::Parameters)
    """
    @ASSUMING the (k,J) part has been prepared
    """

    # Computing the basis FT (k',J')
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,params.CARparams)

    res = 0.
    for i = 1:basis.nmax
        for j = 1:basis.nmax
            res -= coupling.UFT[i] * coupling.UFTp[j]
        end
    end
    return res
end