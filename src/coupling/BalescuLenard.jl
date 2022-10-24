

struct BalescuLenardCoupling

    name::String 

    basis::AB.BasisType
    nradial::Int64

    UFT::Array{Float64}
    UFTp::Array{Float64}

    fht::FHT.FHTtype

    aMcoef::Array{Float64,4}
    M::Matrix{Complex{Float64}}
    IMat::Matrix{Complex{Float64}}
    Xi::Matrix{Complex{Float64}}

end

function BalescuLenardCouplingCreate(basis::AB.BasisType,fht::FHT.FHTtype,params::CAR.ResponseParameters;name::String="BalescuLenard")

    nradial = basis.nmax

    tabaMcoef = CAR.StageaMcoef(params)

    return BalescuLenardCoupling(name,basis,nradial,
                                zeros(Float64,nradial),zeros(Float64,nradial),
                                fht,
                                tabaMcoef,
                                zeros(Complex{Float64},nradial,nradial),
                                zeros(Complex{Float64},nradial,nradial),
                                zeros(Complex{Float64},nradial,nradial))
end

function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::Complex{Float64},
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    coupling::BalescuLenardCoupling,
                    params::Parameters)

    if lharmonic < 0 
        @assert (params.CARparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1 *= -1
        k2 *= -1
        # M^{-l}(ω) = (M^{l}(-ω*))*
        ω = - conj(ω)
    end
    # Computing the basis FT (k,J)
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,params.CARparams)

    # Computing the response matrix
    CAR.tabM!(ω,coupling.M,coupling.aMcoef,coupling.fht,params.CARparams)

    if lharmonic < 0 
        # M^{-l}(ω) = (M^{l}(-ω*))*
        for j = 1:params.CARparams.nradial
            for i = 1:params.CARparams.nradial
                coupling.M[i,j] = conj(coupling.M[i,j])
            end
        end
    end

    tabXi = inv(Symmetric(coupling.IMat - coupling.M)) # @ TO IMPROVE -- in place

    for j = 1:params.CARparams.nradial
        for i = 1:params.CARparams.nradial
            coupling.Xi[i,j] = tabXi[i,j]
        end
    end
end


function CouplingCoefficient(a::Float64,e::Float64,
                             Ω1::Float64,Ω2::Float64,
                             ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1::Int64,k2::Int64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ω::Complex{Float64},
                             ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                             coupling::BalescuLenardCoupling,
                             params::Parameters)

    """
    @ASSUMING the (k,J) part has been prepared
    """
    if lharmonic < 0 
        @assert (params.CARparams.lharmonic == -lharmonic) "Unexpected lharmonic"
        # ψ^{n,-l}_k = ψ^{n,l}_{-k}
        k1p *= -1
        k2p *= -1
    end
    # Computing the basis FT (k',J')
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFTp,params.CARparams)

    res = 0.
    for j = 1:params.CARparams.nradial
        for i = 1:params.CARparams.nradial
            res -= coupling.UFT[i] * coupling.Xi[i,j] * coupling.UFTp[j]
        end
    end
    return res
end