


struct LandauBasisCoupling

    name::String 

    basis::AB.BasisType
    
    UFT::Array{Float64}
    UFTp::Array{Float64}
end

function LandauBasisCouplingCreate(basis::AB.BasisType;name::String="LandauBasis")

    return LandauBasisCoupling(name,basis,
                               zeros(Float64,basis.nmax),zeros(Float64,basis.nmax))
end


function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::Complex{Float64},
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                    coupling::LandauBasisCoupling,
                    params::Parameters)

    if lharmonic < 0 
        @assert (params.CARparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k,J) 
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.basis,coupling.UFT,params.CARparams)
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
                             ω::Complex{Float64},
                             ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                             coupling::LandauBasisCoupling,
                             params::Parameters)
    """
    @ASSUMING the (k,J) part has been prepared
    """

    if lharmonic < 0 
        @assert (params.CARparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
    end
    # Computing the basis FT (k',J')
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling.basis,coupling.UFTp,params.CARparams)

    res = 0.
    for i = 1:params.CARparams.nradial
        for j = 1:params.CARparams.nradial
            res -= coupling.UFT[i] * coupling.UFTp[j]
        end
    end
    return res
end