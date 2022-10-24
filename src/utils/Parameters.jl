

struct Parameters

    CARparams::CAR.ResponseParameters

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

    VERBOSE::Int64
end

function ParametersCreate(CARparams::CAR.ResponseParameters,
                          basis::AB.BasisType;
                          lharmonic::Int64=2,n1max::Int64=10,
                          MTOT::Float64=1.,N::Int64=10^4,
                          Kv::Int64=200,VMAPN::Int64=2,
                          VERBOSE::Int64=0)

    # Resonance vectors
    nbResPair, tabResPair = MakeTabResPair(lharmonic,n1max,basis.dimension)

    # !!! CHECH COMPATIBILITY !!! TO ADD

    return Parameters(CARparams,
                      lharmonic,n1max,basis.dimension,
                      MTOT,N,nbResPair,tabResPair,Kv,VMAPN,VERBOSE)
end
