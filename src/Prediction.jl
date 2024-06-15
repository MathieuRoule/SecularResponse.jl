
"""
    GetSecularResContrib(J,k,k')

Resonances (k,kp) contribution to the secular evolution on  harmonic l.
"""
function GetSecularResContrib(
    J::Float64,
    L::Float64,
    k1::Int64,
    k2::Int64,
    k1p::Int64,
    k2p::Int64,
    lharmonic::Int64,
    model::Potential,
    DF::F0,ndFdJ::F1,
    coupling::AbstractCoupling,
    params::SecularParameters
)  where {F0 <: Function, F1 <: Function}

    Linearparams = params.Linearparams
    Orbitalparams = Linearparams.Orbitalparams
    Ω₀ = frequency_scale(model)

    # Semi-major axis and eccentricity
    a, e = ae_from_actions(J,L,model,Orbitalparams)
    ((a <= 0.) || (e < 0.) || (e > 1.)) && error("Wrong domain ! For J = ",J," ; L = ",L," ; invertion gives a = ",a," ; e = ",e," ")

    # Frequency
    Ω1, Ω2 = frequencies_from_ae(a,e,model,Orbitalparams)
    kdotΩ = k1*Ω1 + k2*Ω2
    ωres = kdotΩ + 0.0*im
    
    # k' related values
    resp = Resonance(k1p,k2p,model,Orbitalparams)
    # Resonance u'
    upres = real(OrbitalElements.rescaled_ϖ(ωres/Ω₀,resp))

    # If no possible resonance (k,k')
    ( -1. < upres < 1.) || (return 0., 0., 0.)

    # Functions value in J
    Eval, Lval = EL_from_ae(a,e,model,Orbitalparams)
    valF = DF(Eval,Lval)
    valkdFdJ = ndFdJ(k1,k2,Eval,Lval,kdotΩ)

    # Prepare the stable part of the coupling coefficients (not changing with J')
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ωres,model,coupling,Linearparams)

    # Initialising the contribution
    fric, diff, flux = 0.0, 0.0, 0.0

    # Integration interval (along the resonance line)
    vpmin, vpmax = v_boundaries(upres,resp,model,Orbitalparams)

    # Integration step
    δv2 = 1.0/params.Kv

    for kvval in 1:params.Kv

        # get the current v value
        v2   = δv2*(kvval-0.5)
        vpval, dvpdv2 = v_from_vp_derivative(v2,vpmin,vpmax,n=params.VMAPN)

        ####
        # (ures',v') → (a',e')
        ####
        # (u',v') → (α',β')
        αp, βp = αβ_from_uv(upres,vpval,resp)
        # (α',β') → (Ω1',Ω2')
        Ω1p, Ω2p = frequencies_from_αβ(αp,βp,model)
        kdotΩp = k1p*Ω1p + k2p*Ω2p
        # (α',β') → (a',e')
        ap, ep = ae_from_αβ(αp,βp,model,Orbitalparams)

        # need (E,L): this has some relatively expensive switches
        Evalp, Lvalp = EL_from_ae(ap,ep,model,Orbitalparams)
        valFp = DF(Evalp,Lvalp)
        valkdFdJp = ndFdJ(k1p,k2p,Evalp,Lvalp,kdotΩp)

        #####
        # compute Jacobians
        #####

        # (u,v2) → (u,v') : dv'/ dv2
        Jacv2 = dvpdv2

        # (u',v') → (α',β').
        # Renormalized. (2/(ωmax-ωmin) * |∂(α',β')/∂(u',v')|)
        ωminp, ωmaxp = frequency_extrema(resp)
        RenormalizedJacαβ = (2/(ωmaxp-ωminp)) * uv_to_αβ_jacobian(upres,vpval,resp)

        # (α',β') → (E',L') : |∂(E',L')/ ∂(α',β')| = |∂(E',L')/ ∂(a',e')| / |∂(α',β')/∂(a',e')|
        # compute the Jacobian of the (α,β) ↦ (E,L) 
        JacEL = ae_to_EL_jacobian(ap,ep,model,Orbitalparams)/ae_to_αβ_jacobian(ap,ep,model,Orbitalparams)

        # (E',L') → (Jr',L') : |∂(Jr,L)/ ∂(E,L)|
        JacJ = (1/Ω1p)

        # Coupling coefficient
        # The (k,J) and ωres parts are already prepared and stored in coupling
        SQpsid = (abs(CouplingCoefficient(ap,ep,Ω1p,Ω2p,k1p,k2p,lharmonic,model,coupling,Linearparams)))^(2)

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
function GetSecular(
    J::Float64,
    L::Float64,
    model::Potential,
    DF::F0,
    ndFdJ::F1,
    coupling::AbstractCoupling,
    params::SecularParameters
)  where {F0 <: Function, F1 <: Function}

    # Initialising the contribution
    totfric = zeros(Float64,2)
    totdiff = zeros(Float64,2,2)
    totflux = zeros(Float64,2)

    totfluxdJ = zeros(Float64,2)
    totfluxdL = zeros(Float64,2)
    dFdt = 0.0

    # Derivative points
    JdJ, LdL = J + params.dJ, L + params.dL

    for ires = 1:params.nbResPair

        # Considered resonance pair
        k1, k2   = params.tabResPair[1,ires], params.tabResPair[2,ires]
        k1p, k2p = params.tabResPair[3,ires], params.tabResPair[4,ires]

        # Computing the associated contribution for l = lharmonic
        fric, diff, flux = GetSecularResContrib(J,L,k1,k2,k1p,k2p,params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)
        _, _, fluxdJ = GetSecularResContrib(JdJ,L,k1,k2,k1p,k2p,params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)         
        _, _, fluxdL = GetSecularResContrib(J,LdL,k1,k2,k1p,k2p,params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)
                                                
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
        fric, diff, flux = GetSecularResContrib(J,L,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)
        _, _, fluxdJ = GetSecularResContrib(JdJ,L,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)
        _, _, fluxdL = GetSecularResContrib(J,LdL,-k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                model,DF,ndFdJ,coupling,params)

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
    dFdt = (totfluxdJ[1]-totflux[1])/params.dJ + (totfluxdL[2]-totflux[2])/params.dJ

    # Adding overall prefactors
    preffric, prefdiff, prefflux, prefdFdt = SecularPrefactors(params.MTOT,params.N,params.dimension)
    totfric *= preffric
    totdiff *= prefdiff
    totflux *= prefflux
    dFdt *= prefdFdt

    return totfric, totdiff, totflux, dFdt
end

"""
    GetSecular(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).
"""
function GetSecular(
    tabJL::Matrix{Float64},
    model::Potential,
    DF::F0,
    ndFdJ::F1,
    coupling::AbstractCoupling,
    params::SecularParameters;
    store::Bool=false
)  where {F0 <: Function, F1 <: Function}

    # check wmat directory before proceeding (save time if not.)
    CheckValidDirectory(params.secdir) || error("False directory")
    
    npts = (size(tabJL))[2]

    # Initialising the contribution
    totfric = zeros(Float64,2,npts)
    totdiff = zeros(Float64,2,2,npts)
    totflux = zeros(Float64,2,npts)
    totdFdt = zeros(Float64,npts)

    couplings = [deepcopy(coupling) for k=1:Threads.nthreads()]

    p = Progress(npts,desc="Secular computations: ")

    Threads.@threads for k = 1:npts

        J, L = tabJL[1,k], tabJL[2,k]

        t = Threads.threadid()
        fric, diff, flux, dFdt = GetSecular(J,L,model,DF,ndFdJ,couplings[t],params)

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

    if store
        outputfilename = SecularFilename(params,coupling.name)
        h5open(outputfilename, "w") do file
            # Basis results
            write(file,"tabJL",tabJL)
            write(file,"totfric",totfric)
            write(file,"totdiff",totdiff)
            write(file,"totflux",totflux)
            write(file,"dFdt",totdFdt)
            # Parameters
            WriteParameters(file,params)
        end
    end

    return totfric, totdiff, totflux, totdFdt
end