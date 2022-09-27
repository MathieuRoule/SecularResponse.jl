



"""
    BLPdSCC(k,k',J,J',ω)

Balescu-Lenard coupling coefficient using 'Pain de Sucre' method
"""
function BLPdSCC(k1::Int64,k2::Int64,
                 k1p::Int64,k2p::Int64,
                 a::Float64,e::Float64,
                 ap::Float64,ep::Float64,
                 ω::Float64,
                 basis::AB.Basis_type;
                 VERBOSE::Int64=0)

    psidBL = BalescuLenardCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis,ω;VERBOSE=VERBOSE)
    psiMulti = LandauMultipoleCC(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis;VERBOSE=VERBOSE)
    psiBasis = LandauBasisCC(kk1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,basis;VERBOSE=VERBOSE)
    return psidBL + psiMulti - psiBasis # "Pain de Sucre"
end