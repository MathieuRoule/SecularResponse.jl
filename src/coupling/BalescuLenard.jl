



"""
    BalescuLenardCC(k,k',J,J',ω)

Balescu-Lenard coupling coefficient
"""
function BalescuLenardCC(k1::Int64,k2::Int64,
                         k1p::Int64,k2p::Int64,
                         lharmonic::Int64,
                         a::Float64,e::Float64,
                         Ω1::Float64,Ω2::Float64,
                         ap::Float64,ep::Float64,
                         Ω1p::Float64,Ω2p::Float64,
                         basis::Float64,
                         ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,
                         RM::RespMat,
                         ω::Float64;
                         Kw::Int64=50
                         EDGE::Float64=0.03,
                         VERBOSE::Int64=0,
                         KuTruncation::Int64=10000)

    FTk = zeros(Float64,basis.nmax)
    FTkp = zeros(Float64,basis.nmax)
    
    # Computing the basis FT (k,J) and (k',J')
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,FTk;Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,FTkp;Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)

    # Computing the response matrix
    tabM!(ω,RM.tabM,RM.tabaMcoef,RM.tabResVec,RM.tabnpnq,RM.FHT,dψ,d2ψ,basis.nmax,RM.Ω₀,rmin,rmax;VERBOSE=VERBOSE,KuTruncation=KuTruncation)
    tabXi = inv(Symmetric(RM.IMat - RM.tabM)) # @ TO IMPROVE -- in place
    res = 0.
    for j = 1:basis.nmax
        for i = 1:basis.nmax
            res -= FTk[i] * tabXi[i,j] * FTkp[j]
        end
    end
    return res
end