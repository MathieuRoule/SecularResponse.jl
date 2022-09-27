



"""
    LandauBasisCC(k,k',J,J',ω)

Landau coupling coefficient using the basis elements
"""
function LandauBasisCC(k1::Int64,k2::Int64,
                       k1p::Int64,k2p::Int64,
                       lharmonic::Int64,
                       a::Float64,e::Float64,
                       Ω1::Float64,Ω2::Float64,
                       ap::Float64,ep::Float64,
                       Ω1p::Float64,Ω2p::Float64,
                       basis::Float64;
                       Kw::Int64=50
                       EDGE::Float64=0.03,
                       VERBOSE::Int64=0)

    FTk = zeros(Float64,basis.nmax)
    FTkp = zeros(Float64,basis.nmax)
    
    # Computing the basis FT (k,J) and (k',J')
    CAR.WBasisFT(a,e,Ω1,Ω2,k1,k2,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,FTk;Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)
    CAR.WBasisFT(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,basis,ψ,dψ,d2ψ,d3ψ,FTkp;Kw=Kw,EDGE=EDGE,VERBOSE=VERBOSE)

    res = 0.
    for i = 1:basis.nmax
        for j = 1:basis.nmax
            res -= FTk[i] * FTkp[j]
        end
    end
    return res
end