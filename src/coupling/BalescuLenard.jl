

struct BalescuLenardCoupling

    name::String 

    basis::AB.Basis_type
    UFT::Array{Float64}
    UFTp::Array{Float64}

    aMcoef::Array{Float64,4}
    M::Matrix{Complex{Float64}}
    Xi::Matrix{Complex{Float64}}

end

function BalescuLenardCouplingCreate(basis::AB.Basis_type, params::CAR.Parameters;name::String="BalescuLenard")

    return BalescuLenardCoupling()
end

function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::Float64,
                    ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                    coupling::BalescuLenardCoupling,
                    params::Parameters)

    # Computing the basis FT (k,J)
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,lharmonic,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,params.CARparams)

    # Computing the response matrix
    tabM!(ω,coupling.M,coupling.aMcoef,coupling.tabResVec,coupling.tabnpnq,RM.FHT,dψ,d2ψ,basis.nmax,params.CARparams)
    tabXi = inv(Symmetric(RM.IMat - RM.tabM)) # @ TO IMPROVE -- in place

    for j = 1:params.nradial
        for i = 1:params.nradial
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
                             ω::Float64,
                             ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                             coupling::BalescuLenardCoupling,
                             params::Parameters)

    """
    @ASSUMING the (k,J) part has been prepared
    """
    # Computing the basis FT (k',J')
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,ψ,dψ,d2ψ,d3ψ,coupling.basis,coupling.UFT,params.CARparams)

    res = 0.
    for j = 1:params.nradial
        for i = 1:params.nradial
            res -= FTk[i] * coupling.Xi[i,j] * FTkp[j]
        end
    end
    return res
end