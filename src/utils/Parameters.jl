

struct SecularParameters

    Linearparams::LR.LinearParameters

    lharmonic::Int64
    n1max::Int64
    dimension::Int64

    MTOT::Float64
    N::Int64

    nbResPair::Int64 
    tabResPair::Matrix{Int64}

    # Integration Parameters
    Kv::Int64
    VMAPN::Int64

    # Numerical derivative parameters (flux divergence)
    dJ::Float64
    dL::Float64

    secdir::String

    VERBOSE::Int64
end

function SecularParameters(basis::AB.AbstractAstroBasis;
                           Linearparams::LR.LinearParameters=LR.LinearParameters(),
                           n1max::Int64=10,
                           MTOT::Float64=1.,N::Int64=10^4,
                           Kv::Int64=200,VMAPN::Int64=2,
                           dJ::Float64=1.e-3,dL::Float64=1.e-3,
                           secdir::String="",VERBOSE::Int64=0)

    # Resonance vectors
    nbResPair, tabResPair = MakeTabResPair(Linearparams.lharmonic,n1max,basis.dimension)

    # !!! CHECH COMPATIBILITY !!! TO ADD

    return SecularParameters(Linearparams,
                      Linearparams.lharmonic,n1max,basis.dimension,
                      MTOT,N,nbResPair,tabResPair,Kv,VMAPN,dJ,dL,secdir,VERBOSE)
end
