
"""
    GetSecularResContrib(J,k,k')

Resonances (k,kp) contribution to the secular evolution.
"""
function GetSecularResContrib(a::Float64,e::Float64,
                              k1::Int64,k2::Int64,
                              k1p::Int64,k2p::Int64,
                              lharmonic::Int64,
                              ψ::F0,dψ::F1,d2ψ::F2,d3ψ::F3,d4ψ::F4,βc::F5,
                              DF::F6,ndFdJ::F7,
                              coupling::CouplingType,
                              params::SecularParameters)  where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function, F5 <: Function, F6 <: Function, F7 <: Function}

    Linearparams = params.Linearparams
    Orbitalparams = Linearparams.Orbitalparams
    Ω₀ = Orbitalparams.Ω₀

    # Frequency
    Ω1, Ω2 = OE.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,Orbitalparams)
    kdotΩ = k1*Ω1 + k2*Ω2
    ωres = kdotΩ + 0.0*im
    
    # k' related values
    ωminp, ωmaxp = OE.Findωminωmax(k1p,k2p,dψ,d2ψ,Orbitalparams)
    # Resonance u'
    upres = real(OE.Getϖ(ωres/Ω₀,ωminp,ωmaxp))

    # If no possible resonance (k,k')
    ( -1. < upres < 1.) || (return 0., 0., 0.)

    # Functions value in J
    Eval, Lval = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,Orbitalparams)
    valF = DF(Eval,Lval)
    valkdFdJ = ndFdJ(k1,k2,Eval,Lval,kdotΩ)

    # Prepare the stable part of the coupling coefficients (not changing with J')
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling,Linearparams)

    # Initialising the contribution
    fric, diff, flux = 0.0, 0.0, 0.0

    # Integration interval (along the resonance line)
    vpmin, vpmax = OE.FindVminVmax(upres,k1p,k2p,dψ,d2ψ,ωminp,ωmaxp,βc,Orbitalparams)

    # Integration step
    δv2 = 1.0/params.Kv

    for kvval in 1:params.Kv

        # get the current v value
        v2   = δv2*(kvval-0.5)
        vpval = CAR.vFromvp(v2,vpmin,vpmax,params.VMAPN)

        ####
        # (ures',v') → (a',e')
        ####
        # (u',v') → (α',β')
        αp, βp = OE.αβFromUV(upres,vpval,k1p,k2p,ωminp,ωmaxp)
        # (α',β') → (Ω1',Ω2')
        Ω1p, Ω2p = OE.FrequenciesFromαβ(αp,βp,Ω₀)
        kdotΩp = k1p*Ω1p + k2p*Ω2p
        # (Ω1p,Ω2p) → (ap,ep)
        ap, ep = OE.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1p,Ω2p,Orbitalparams)

        # need (E,L): this has some relatively expensive switches
        Evalp, Lvalp = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,ap,ep,Orbitalparams)
        valFp = DF(Evalp,Lvalp)
        valkdFdJp = ndFdJ(k1p,k2p,Evalp,Lvalp,kdotΩp)

        #####
        # compute Jacobians
        #####

        # (u,v2) → (u,v') : dv'/ dv2
        Jacv2 = CAR.DvDvp(v2,vpmin,vpmax,params.VMAPN)

        # (u',v') → (α',β').
        # Renormalized. (2/(ωmax-ωmin) * |∂(α',β')/∂(u',v')|)
        RenormalizedJacαβ = OE.RenormalizedJacUVToαβ(k1p,k2p,upres,vpval)

        # (α',β') → (E',L') : |∂(E',L')/ ∂(α',β')| 
        JacEL = OE.JacαβToELAE(ψ,dψ,d2ψ,d3ψ,d4ψ,ap,ep,Orbitalparams)

        # (E',L') → (Jr',L') : |∂(Jr,L)/ ∂(E,L)|
        JacJ = (1/Ω1p)

        # Coupling coefficient
        SQpsid = (abs(CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling,Linearparams)))^(2)

        commonpart = Jacv2 * JacJ * JacEL * RenormalizedJacαβ * SQpsid

        fric += commonpart * valkdFdJp
        diff += commonpart * valFp
        flux += commonpart * (valkdFdJp * valF - valkdFdJ * valFp)

    end

    # remove dimensionality from Ω mapping
    dimensionl = δv2/Ω₀
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

    # Friction, Diffusion, Flux, dFdt
    # Additional 2 for diffusion as Flux = Fric * F - 1/2 Diff * dFdJ
    return prefμ, 2*prefμ, prefμ, -prefμ
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Friction, Diffusion, Flux) at given actions (Jr,L)
on the harmonic l (decoupled harmonic numbers).

@ATTENTION: We take into account both contribution from lharmonic and -lharmonic
"""
function GetSecular(J::Float64,L::Float64,
                     ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,βc::Function,
                     DF::Function,ndFdJ::Function,
                     coupling::CouplingType,
                     params::SecularParameters)

    # Initialising the contribution
    totfric = zeros(Float64,2)
    totdiff = zeros(Float64,2,2)
    totflux = zeros(Float64,2)

    totfluxdJ = zeros(Float64,2)
    totfluxdL = zeros(Float64,2)
    dFdt = 0.0

    Orbitalparams = params.Linearparams.Orbitalparams

    # (J,L) → (a,e)
    a, e = OE.ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,J,L,Orbitalparams)
    ((a <= 0.) || (e < 0.) || (e > 1.)) && error("Wrong domain ! For J = ",J," ; L = ",L," ; invertion gives a = ",a," ; e = ",e," ")
    adJ, edJ = OE.ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,J+params.dJ,L,Orbitalparams)
    ((adJ <= 0.) || (edJ < 0.) || (edJ > 1.)) && error("Wrong domain ! For J = ",J+params.dJ," ; L = ",L," ; invertion gives a = ",adJ," ; e = ",edJ," ")
    adL, edL = OE.ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,J,L+params.dL,Orbitalparams)
    ((adL <= 0.) || (edL < 0.) || (edL > 1.)) && error("Wrong domain ! For J = ",J," ; L = ",L+params.dL," ; invertion gives a = ",adL," ; e = ",edL," ")

    for ires = 1:params.nbResPair

        # Considered resonance pair
        k1, k2   = params.tabResPair[1,ires], params.tabResPair[2,ires]
        k1p, k2p = params.tabResPair[3,ires], params.tabResPair[4,ires]

        # Computing the associated contribution for l = lharmonic
        fric, diff, flux = GetSecularResContrib(a,e,k1,k2,k1p,k2p,params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)
        _, _, fluxdJ = GetSecularResContrib(adJ,edJ,k1,k2,k1p,k2p,params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)         
        _, _, fluxdL = GetSecularResContrib(adL,edL,k1,k2,k1p,k2p,params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)
        # Adding the contributions
        totfric[1] += k1 * fric
        totfric[2] += k2 * fric
        totdiff[1,1] += k1 * k1 * diff
        totdiff[1,2] += k1 * k2 * diff
        totdiff[2,1] += k2 * k1 * diff
        totdiff[2,2] += k2 * k2 * diff
        totflux[1] += k1 * flux
        totflux[2] += k2 * flux

        totfluxdJ[1] += k1 * fluxdJ
        totfluxdJ[2] += k2 * fluxdJ
        totfluxdL[1] += k1 * fluxdL
        totfluxdL[2] += k2 * fluxdL

        # Computing the associated contribution for l = -lharmonic
        fric, diff, flux = GetSecularResContrib(a,e,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)
        _, _, fluxdJ = GetSecularResContrib(adJ,edJ,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)
        _, _, fluxdL = GetSecularResContrib(adL,edL,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,coupling,params)

        # Adding the contributions
        totfric[1] += (-k1) * fric
        totfric[2] += (-k2) * fric
        totdiff[1,1] += k1 * k1 * diff
        totdiff[1,2] += k1 * k2 * diff
        totdiff[2,1] += k2 * k1 * diff
        totdiff[2,2] += k2 * k2 * diff
        totflux[1] += (-k1) * flux
        totflux[2] += (-k2) * flux

        totfluxdJ[1] += (-k1) * fluxdJ
        totfluxdJ[2] += (-k2) * fluxdJ
        totfluxdL[1] += (-k1) * fluxdL
        totfluxdL[2] += (-k2) * fluxdL
    end

    # Computing the flux divergence
    dFdt = (totfluxdJ[1]-totflux[1])/params.dJ + (totfluxdL[2]-totflux[2])/params.dL

    # Adding overall prefactors
    preffric, prefdiff, prefflux, prefdFdt = SecularPrefactors(params.MTOT,params.N,params.dimension)
    totfric *= preffric
    totdiff *= prefdiff
    totflux *= prefflux
    dFdt *= prefdFdt

    return a, e, totfric, totdiff, totflux, dFdt
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).
"""
function GetSecular(tabJL::Matrix{Float64},
                     ψ::Function,dψ::Function,d2ψ::Function,d3ψ::Function,d4ψ::Function,
                     DF::Function,ndFdJ::Function,
                     coupling::CouplingType,
                     params::SecularParameters)

    # check wmat directory before proceeding (save time if not.)
    CheckValidDirectory(params.secdir) || error("False directory")
    
    npts = (size(tabJL))[2]

    # Corresponding actions/frequencies
    tabAE = zeros(Float64,2,npts)

    # Initialising the contribution
    totfric = zeros(Float64,2,npts)
    totdiff = zeros(Float64,2,2,npts)
    totflux = zeros(Float64,2,npts)
    totdFdt = zeros(Float64,npts)

    couplings = [deepcopy(coupling) for k=1:Threads.nthreads()]

    # β circular curve
    βc(αc::Float64) = OE.βcirc(αc,dψ,d2ψ,params.Linearparams.Orbitalparams)

    p = Progress(npts,desc="Secular computations: ")

    Threads.@threads for k = 1:npts

        J, L = tabJL[1,k], tabJL[2,k]

        t = Threads.threadid()
        a, e, fric, diff, flux, dFdt = GetSecular(J,L,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,couplings[t],params)

        # Orbit
        tabAE[1,k] = a
        tabAE[2,k] = e

        # Evolution
        totfric[1,k] = fric[1]
        totfric[2,k] = fric[2]
        totdiff[1,1,k] = diff[1,1]
        totdiff[2,1,k] = diff[2,1]
        totdiff[1,2,k] = diff[1,2]
        totdiff[2,2,k] = diff[2,2]
        totflux[1,k] = flux[1]
        totflux[2,k] = flux[2]

        totdFdt[k] = dFdt

        next!(p)
    end

    outputfilename = SecularFilename(params,coupling.name)
    h5open(outputfilename, "w") do file
        # Basis results
        write(file,"tabJL",tabJL)
        write(file,"tabAE",tabAE)
        write(file,"totfric",totfric)
        write(file,"totdiff",totdiff)
        write(file,"totflux",totflux)
        write(file,"dFdt",totdFdt)
        # Parameters
        WriteParameters(file,params)
    end

    return tabAE, totfric, totdiff, totflux
end