
"""
    GetSecularResContrib(J,k,k')

Resonances (k,kp) contribution to the secular evolution.
"""
function GetSecularResContrib(a::Float64,e::Float64,
                              Ω1::Float64,Ω2::Float64,
                              k1::Int64,k2::Int64,
                              k1p::Int64,k2p::Int64,
                              lharmonic::Int64,
                              ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,βc::Function,
                              DF::Function,ndFdJ::Function,
                              coupling::CouplingType,
                              params::Parameters)

    kdotΩ = k1*Ω1 + k2*Ω2
    ωres = kdotΩ + 0.0*im

    Ω₀ = params.CARparams.Ω₀
    rmin, rmax = params.CARparams.rmin, params.CARparams.rmax
    αmin, αmax = params.CARparams.αmin, params.CARparams.αmax
    
    # k' related values
    ωminp, ωmaxp = OE.Findωminωmax(k1p,k2p,dψ,d2ψ,Ω₀=Ω₀,rmin=rmin,rmax=rmax)
    # Resonance u
    ures = real(OE.Getϖ(ωres,ωminp,ωmaxp))

    # If no possible resonance (k,k')
    ( -1. < ures < 1.) || (return 0., 0., 0.)

    # Functions value in J
    Eval, Lval = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e)
    valF = DF(Eval,Lval)
    valkdFdJ = ndFdJ(k1,k2,Eval,Lval,kdotΩ)

    # Prepare the stable part of the coupling coefficients (not changing with J')
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,coupling,params)

    # Initialising the contribution
    fric, diff, flux = 0.0, 0.0, 0.0

    # Integration interval (along the resonance line)
    vmin, vmax = OE.FindVminVmax(ures,k1p,k2p,dψ,d2ψ,ωminp,ωmaxp,αmin,αmax,βc,Ω₀=Ω₀,rmin=rmin,rmax=rmax)

    # Integration step
    δvp = 1.0/params.Kv

    for kvval in 1:params.Kv

        # get the current v value
        vp   = δvp*(kvval-0.5)
        vval = CAR.vprime(vp,vmin,vmax,n=params.VMAPN)

        # vp -> v
        Jacvp = CAR.dvprime(vp,vmin,vmax,n=params.VMAPN)

        ####
        # (ures,v') -> (a',e')
        ####
        # (u',v') -> (α',β')
        αp, βp = OE.αβFromUV(ures,vval,k1p,k2p,ωminp,ωmaxp)
        # (α',β') -> (Ω1',Ω2')
        Ω1p, Ω2p = αp*Ω₀,αp*βp*Ω₀
        kdotΩp = k1p*Ω1p + k2p*Ω2p
        # (Ω1p,Ω2p) -> (ap,ep)
        ap, ep = OE.AEFromΩ1Ω2Brute(Ω1p,Ω2p,ψ,dψ,d2ψ,d3ψ,NINT=params.CARparams.NINT,EDGE=params.CARparams.EDGE,VERBOSE=params.CARparams.VERBOSE)

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
        JacEL = OE.JacELToαβAE(ap,ep,ψ,dψ,d2ψ,params.CARparams.Ω₀)

        #(J) -> (E,L)
        JacJ = (1/Ω1p)

        # Coupling coefficient
        SQpsid = (abs(CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,coupling,params)))^(2)

        commonpart = Jacvp * JacJ * JacEL * Jacαβ * SQpsid

        fric += commonpart * valkdFdJp
        diff += commonpart * valFp
        flux += commonpart * (valkdFdJp * valF - valkdFdJ * valFp)

    end

    # remove dimensionality from Ω mapping
    dimensionl = δvp/Ω₀
    fric *= dimensionl
    diff *= dimensionl
    flux *= dimensionl
    
    return fric, diff, flux
end

"""SecularPrefactors
"""
function SecularPrefactors(Mtot::Float64,N::Int64,dim::Int64)
    
    pref = (dim == 2) ? 4*(pi^3) : 0.0
    prefμ = (Mtot / N) * pref

    # Friction, Diffusion, Flux
    # Additional 2 for diffusion as Flux = Fric * F - 1/2 Diff * dFdJ
    return prefμ, 2*prefμ, prefμ
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).

@ATTENTION: We take into account both contribution from lharmonic and -lharmonic
"""
function GetSecular(a::Float64,e::Float64,
                     ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,βc::Function,
                     DF::Function,ndFdJ::Function,
                     coupling::CouplingType,
                     params::Parameters)

    """
    @IMPROVE: Make it a function of Jr,L (i.e. make a (Jr,L) -> (a,e) mapping
    """

    (params.VERBOSE > 0) && println("SecularResponse.GetSecular: Considering $(params.nbResPair) resonances pairs.")

    # Initialising the contribution
    totfric = zeros(Float64,2)
    totdiff = zeros(Float64,2,2)
    totflux = zeros(Float64,2)

    # Considered a, e : associated frequencies/actions
    Ω1, Ω2, J = OE.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e;action=true,TOLECC=params.CARparams.ELTOLECC,NINT=params.CARparams.NINT,EDGE=params.CARparams.EDGE)
    L = OE.LFromAE(ψ,dψ,d2ψ,d3ψ,a,e,TOLECC=params.CARparams.ELTOLECC)

    for ires = 1:params.nbResPair

        # Considered resonance pair
        k1, k2   = params.tabResPair[1,ires], params.tabResPair[2,ires]
        k1p, k2p = params.tabResPair[3,ires], params.tabResPair[4,ires]
        (params.VERBOSE > 0) && println("SecularResponse.GetSecular: Starting on (($k1,$k2), ($k1p,$k2p)).")

        # Computing the associated contribution for l = lharmonic
        fric, diff, flux = GetSecularResContrib(a,e,Ω1,Ω2,
                                                k1,k2,k1p,k2p,params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,βc,DF,ndFdJ,
                                                coupling,params)

        # Adding the contributions
        totfric[1] += k1 * fric
        totfric[2] += k2 * fric
        totdiff[1,1] += k1 * k1 * diff
        totdiff[1,2] += k1 * k2 * diff
        totdiff[2,1] += k2 * k1 * diff
        totdiff[2,2] += k2 * k2 * diff
        totflux[1] += k1 * flux
        totflux[2] += k2 * flux

        # Computing the associated contribution for l = -lharmonic
        fric, diff, flux = GetSecularResContrib(a,e,Ω1,Ω2,
                                                -k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,βc,DF,ndFdJ,
                                                coupling,params)

        # Adding the contributions
        totfric[1] += (-k1) * fric
        totfric[2] += (-k2) * fric
        totdiff[1,1] += k1 * k1 * diff
        totdiff[1,2] += k1 * k2 * diff
        totdiff[2,1] += k2 * k1 * diff
        totdiff[2,2] += k2 * k2 * diff
        totflux[1] += (-k1) * flux
        totflux[2] += (-k2) * flux
    end

    # Adding overall prefactors
    preffric, prefdiff,prefflux = SecularPrefactors(params.MTOT,params.N,params.dimension)
    totfric[1] *= preffric
    totfric[2] *= preffric
    totdiff[1,1] *= prefdiff
    totdiff[1,2] *= prefdiff
    totdiff[2,1] *= prefdiff
    totdiff[2,2] *= prefdiff
    totflux[1] *= prefflux
    totflux[2] *= prefflux

    return J, L, totfric, totdiff, totflux
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).
"""
function GetSecular(tabAE::Matrix{Float64},
                     ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,βc::Function,
                     DF::Function,ndFdJ::Function,
                     coupling::CouplingType,
                     params::Parameters)

    npts = (size(tabAE))[2]

    # Corresponding actions/frequencies
    tabJL = zeros(Float64,2,npts)

    # Initialising the contribution
    totfric = zeros(Float64,2,npts)
    totdiff = zeros(Float64,2,2,npts)
    totflux = zeros(Float64,2,npts)

    couplings = [deepcopy(coupling) for k=1:Threads.nthreads()]

    Threads.@threads for k = 1:npts

        t = Threads.threadid()

        a, e = tabAE[1,k], tabAE[2,k]

        J, L, fric, diff, flux = GetSecular(a,e,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,couplings[t],params)

        # Orbit
        tabJL[1,k] = J
        tabJL[2,k] = L

        # Evolution
        totfric[1,k] = fric[1]
        totfric[2,k] = fric[2]
        totdiff[1,1,k] = diff[1,1]
        totdiff[1,2,k] = diff[1,2]
        totdiff[2,1,k] = diff[2,1]
        totdiff[2,2,k] = diff[2,2]
        totflux[1,k] = flux[1]
        totflux[2,k] = flux[2]
    end

    return tabJL, totfric, totdiff, totflux
end