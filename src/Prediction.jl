
"""
    GetSecularResContrib(J,k,k')

Resonances (k,kp) contribution to the secular evolution.
"""
function GetSecularResContrib(a::Float64,e::Float64,
                                k1::Int64,k2::Int64,
                                k1p::Int64,k2p::Int64,
                                lharmonic::Int64,
                                ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                                DF::Function,ndFdJ::Function,
                                modedir::String,
                                Kv::Int64,
                                basis::AstroBasis.Basis_type,
                                Ω₀::Float64,
                                modelname::String,dfname::String,
                                rmin::Float64,rmax::Float64,
                                αmin::Float64,αmax::Float64,
                                βc::Function;
                                COUPLING::String="Landau-Basis",
                                VERBOSE::Int64=0,
                                NINT::Int64=32,
                                EDGE::Float64=0.03)

    # load a value of tabWmat, plus (a,e) values
    filename   = CallAResponse.AxiFilename(modedir,modelname,dfname,lharmonic,k1,k2,rb,Ku,Kv)
    file = h5open(filename,"r")
    close(file)

    # Considered J : associated frequencies
    Ω1, Ω2 = OE.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=false,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)
    kdotΩ = k1*Ω1 + k2*Ω2

    # k' related values
    ωminp, ωmaxp = OE.Findωminωmax(k1p,k2p,dψ,d2ψ,αmin,αmax,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
    # Resonance u
    ures = OE.Getϖ(kdotΩ/Ω0,ωminp,ωmaxp)

    # If no possible resonance (k,k')
    ( -1. < ures < 1.) || (return 0.)
        
    # Integration interval (along the resonance line)
    vmin, vmax = OE.FindVminVmax(ures,k1p,k2p,dψ,d2ψ,ωmin,ωmax,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
    # Integration step
    δv = (vmax - vmin)/(Kv)

    # Initialising the contribution
    fric, diff, flux = 0.0, 0.0, 0.0

    # Functions value in J
    Eval, Lval = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e)
    valF = DF(Eval,Lval)
    valkdFdJ = ndFdJ(k1,k2,Eval,Lval,kdotΩ)

    for kvval in 1:Kv
        vval = vmin + δv*(kvval-0.5)

        ####
        # (ures,v') -> (a',e')
        ####
        # (u',v') -> (α',β')
        αp, βp = OE.αβFromUV(ures,vval,k1p,k2p,ωminp,ωmaxp)
        # (α',β') -> (Ω1',Ω2')
        Ω1p, Ω2p = αp*Ω₀,αp*βp*Ω₀
        kdotΩp = k1p*Ω1p + k2p*Ω2p
        # (Ω1p,Ω2p) -> (ap,ep)
        ap, ep = OE.AEFromΩ1Ω2Brute(Ω1,Ω2,ψ,dψ,d2ψ,d3ψ,NINT=NINT,EDGE=EDGE,VERBOSE=VERBOSE)

        # need (E,L): this has some relatively expensive switches
        Evalp, Lvalp = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,ap,ep)
        valFp = DF(Evalp,Lvalp)
        valkdFdJp = ndFdJ(k1p,k2p,Evalp,Lvalp,kdotΩp)

        # compute Jacobians
        # (α,β) -> (u,v).
        # owing to the remapping of Ω, this has an extra 2/(ωmax-ωmin)
        Jacαβ = OE.JacαβToUV(k1p,k2p,ωminp,ωmaxp,vval)

        #(E,L) -> (α,β): this is the most expensive function here,
        # so we have pre-tabulated it
        JacEL = OE.JacELToαβAE(a,e,ψ,dψ,d2ψ,Ω₀)

        #(J) -> (E,L)
        JacJ = (1/Ω1)

        # Coupling coefficient
        SQpsid = (abs(CouplingCoefficient(k1,k2,k1p,k2p,lharmonic,a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,kdotΩ,basis,COUPLING=COUPLING,VERBOSE=VERBOSE)))^(2)

        commonpart = JacJ * JacEL * Jacαβ * SQpsid

        fric += commonpart * valkdFdJp
        diff += commonpart * valFp
        flux += commonpart * (valkdFdJp * valF - valkdFdJ * valFp)

    end

    # remove dimensionality from Ω mapping
    dimensionl = δv/Ω₀
    fric *= dimensionl
    diff *= dimensionl
    flux *= dimensionl
    
    return fric, diff, flux
end



"""SecularPrefactors
"""
function SecularPrefactors(Mtot::Float64,N::Int64,dim::Int64)
    
    pref = (dim == 2) ? 4*(pi^3) : 0.0
    μ = Mtot / N
    # Friction, Diffusion, Flux
    # Additional 2 for diffusion as Flux = Fric * F - 1/2 Diff * dFdJ
    return pref*μ, 2*pref*μ, pref*μ
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the hamonic l (decoupled harmonic numbers).
"""
function GetSecular(a::Float64,e::Float64,
                     ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                     DF::Function,ndFdJ::Function,
                     modedir::String,
                     Kv::Int64,
                     basis::AstroBasis.Basis_type,
                     n1max::Int64,lharmonic::Int64,
                     Ω₀::Float64,
                     modelname::String,dfname::String,
                     Ku::Int64,
                     rmin::Float64,rmax::Float64,
                     αmin::Float64,αmax::Float64,
                     βc::Function;
                     COUPLING::String="Landau-Basis",
                     VERBOSE::Int64=0,
                     NINT::Int64=32,
                     EDGE::Float64=0.03)

    """
    @IMPROVE: Make it a function of Jr,L (i.e. make a (Jr,L) -> (a,e) mapping
    """
    # Resonance vectors
    nbResPair, tabResPair = CallAResponse.MakeTabResPair(lharmonic,n1max,basis.dimension)
    (VERBOSE >= 0) && println("SecularResponse.GetSecular: Considering $nbResPair resonances.")

    # Initialising the contribution
    totfric = zeros(Float64,2)
    totdiff = zeros(Float64,2,2)
    totflux = zeros(Float64,2)

    for ires = 1:nbResPair

        # Considered resonance pair
        k1, k2, k1p, k2p = tabResPair[ires]
        (VERBOSE > 0) && println("SecularResponse.GetSecular: Starting on (($k1,$k2), ($k1p,$k2p)).")

        # Computing the associated contribution
        fric, diff, flux = GetSecularResContrib(a,e,
                                                k1,k2,k1p,k2p,lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,DF,ndFdJ,
                                                modedir,Kv,basis,Ω₀,modelname,dfname,
                                                rmin,rmax,αmin,αmax,βc;
                                                COUPLING=COUPLING,VERBOSE=VERBOSE,NINT=NINT,EDGE=EDGE)

        # Adding the contributions
        totfric[1] += k1 * fric
        totfric[2] += k2 * fric
        totdiff[1,1] += k1 * k1 * diff
        totdiff[1,2] += k1 * k2 * diff
        totdiff[2,1] += k2 * k1 * diff
        totdiff[2,2] += k2 * k2 * diff
        totflux[1] += k1 * flux
        totflux[2] += k2 * flux
    end

    return totfric, totdiff, totflux
end



struct RespMat
    tabResVec::Matrix{Int64}
    tabnpnq::Matrix{Int64}
    tabM::Array{Complex{Float64},2}
    IMat::Array{Complex{Float64},2}
    tabaMcoef::Array{Float64,4}
    FHT::FiniteHilbertTransform.FHTtype
end

nbResVec, tabResVec = CAR.MakeTabResVec(lharmonic,n1max,ndim) # Number of resonance vectors. ATTENTION, it is for the harmonics lharmonic
tabaMcoef = CAR.StageaMcoef(tabResVec,tab_npnq,Ku,Kv,nradial,
                                          modedir=modedir,modelname=modelname,dfname=dfname,lharmonic=lharmonic,rb=rb)