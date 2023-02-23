
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
                              params::Parameters)  where {F0 <: Function, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function, F5 <: Function, F6 <: Function, F7 <: Function}

    CARparams = params.CARparams
    OEparams = CARparams.OEparams
    Ω₀ = OEparams.Ω₀

    # Frequency
    Ω1, Ω2 = OE.ComputeFrequenciesAE(ψ,dψ,d2ψ,d3ψ,d4ψ,a,e,OEparams)
    kdotΩ = k1*Ω1 + k2*Ω2
    ωres = kdotΩ + 0.0*im
    
    # k' related values
    ωminp, ωmaxp = OE.Findωminωmax(k1p,k2p,dψ,d2ψ,OEparams)
    # Resonance u
    ures = real(OE.Getϖ(ωres/Ω₀,ωminp,ωmaxp))

    # If no possible resonance (k,k')
    ( -1. < ures < 1.) || (return 0., 0., 0.)

    # Functions value in J
    Eval, Lval = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,a,e,OEparams)
    valF = DF(Eval,Lval)
    valkdFdJ = ndFdJ(k1,k2,Eval,Lval,kdotΩ)

    # Prepare the stable part of the coupling coefficients (not changing with J')
    CCPrepare!(a,e,Ω1,Ω2,k1,k2,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling,CARparams)

    # Initialising the contribution
    fric, diff, flux = 0.0, 0.0, 0.0

    # Integration interval (along the resonance line)
    vmin, vmax = OE.FindVminVmax(ures,k1p,k2p,dψ,d2ψ,ωminp,ωmaxp,βc,OEparams)

    # Integration step
    δv2 = 1.0/params.Kv

    for kvval in 1:params.Kv

        # get the current v value
        v2   = δv2*(kvval-0.5)
        vval = CAR.vFromvp(v2,vmin,vmax,params.VMAPN)

        # dv'/ dv2
        Jacvp = CAR.DvDvp(v2,vmin,vmax,params.VMAPN)

        ####
        # (ures,v') -> (a',e')
        ####
        # (u',v') -> (α',β')
        αp, βp = OE.αβFromUV(ures,vval,k1p,k2p,ωminp,ωmaxp)
        # (α',β') -> (Ω1',Ω2')
        Ω1p, Ω2p = OE.FrequenciesFromαβ(αp,βp,Ω₀)
        kdotΩp = k1p*Ω1p + k2p*Ω2p
        # (Ω1p,Ω2p) -> (ap,ep)
        ap, ep = OE.ComputeAEFromFrequencies(ψ,dψ,d2ψ,d3ψ,d4ψ,Ω1p,Ω2p,OEparams)

        # need (E,L): this has some relatively expensive switches
        Evalp, Lvalp = OE.ELFromAE(ψ,dψ,d2ψ,d3ψ,ap,ep,OEparams)
        valFp = DF(Evalp,Lvalp)
        valkdFdJp = ndFdJ(k1p,k2p,Evalp,Lvalp,kdotΩp)


        # compute Jacobians
        # (2/(ωmax-ωmin)) * ∂(α,β)/ ∂(u,v)
        Jacαβ = OE.JacαβToUV(k1p,k2p,vval)

        # ∂(E,L)/ ∂(α,β) 
        JacEL = OE.JacELToαβAE(ψ,dψ,d2ψ,d3ψ,d4ψ,ap,ep,OEparams)

        # ∂(Jr,L)/ ∂(E,L)
        JacJ = (1/Ω1p)

        # Coupling coefficient
        SQpsid = (abs(CouplingCoefficient(a,e,Ω1,Ω2,ap,ep,Ω1p,Ω2p,k1,k2,k1p,k2p,lharmonic,ωres,ψ,dψ,d2ψ,d3ψ,d4ψ,coupling,CARparams)))^(2)

        commonpart = Jacvp * JacJ * JacEL * Jacαβ * SQpsid

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

    # Friction, Diffusion, Flux
    # Additional 2 for diffusion as Flux = Fric * F - 1/2 Diff * dFdJ
    return prefμ, 2*prefμ, prefμ
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
                     params::Parameters)

    # Initialising the contribution
    totfric = zeros(Float64,2)
    totdiff = zeros(Float64,2,2)
    totflux = zeros(Float64,2)

    OEparams = params.CARparams.OEparams

    # Considered a, e : associated frequencies/actions
    a, e = OE.ComputeAEFromActions(ψ,dψ,d2ψ,d3ψ,J,L,OEparams)
    if (a <= 0.) || (e < 0.) || (e > 1.)
        error("Wrong domain ! For J = ",J," ; L = ",L," ; invertion gives a = ",a," ; e = ",e," ")
    end

    for ires = 1:params.nbResPair

        # Considered resonance pair
        k1, k2   = params.tabResPair[1,ires], params.tabResPair[2,ires]
        k1p, k2p = params.tabResPair[3,ires], params.tabResPair[4,ires]

        # Computing the associated contribution for l = lharmonic
        fric, diff, flux = GetSecularResContrib(a,e,
                                                k1,k2,k1p,k2p,params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,
                                                coupling,params)

        # println("k1, k2, k1p, k2p : ",k1," ",k2," ",k1p," ",k2p)
        # println("Flux : ",flux)
        
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
        fric, diff, flux = GetSecularResContrib(a,e,
                                                -k1,-k2,-k1p,-k2p,-params.lharmonic,
                                                ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,
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

    return a, e, totfric, totdiff, totflux
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
                     params::Parameters)

    # check wmat directory before proceeding (save time if not.)
    CheckValidDirectory(params.secdir) || error("False directory")
    
    npts = (size(tabJL))[2]

    # Corresponding actions/frequencies
    tabAE = zeros(Float64,2,npts)

    # Initialising the contribution
    totfric = zeros(Float64,2,npts)
    totdiff = zeros(Float64,2,2,npts)
    totflux = zeros(Float64,2,npts)

    couplings = [deepcopy(coupling) for k=1:Threads.nthreads()]

    # β circular curve
    βc(αc::Float64) = OE.βcirc(αc,dψ,d2ψ,params.CARparams.OEparams)

    p = Progress(npts,desc="Secular computations: ")

    Threads.@threads for k = 1:npts

        J, L = tabJL[1,k], tabJL[2,k]

        t = Threads.threadid()
        a, e, fric, diff, flux = GetSecular(J,L,ψ,dψ,d2ψ,d3ψ,d4ψ,βc,DF,ndFdJ,couplings[t],params)

        # Orbit
        tabAE[1,k] = a
        tabAE[2,k] = e

        # Evolution
        totfric[1,k] = fric[1]
        totfric[2,k] = fric[2]
        totdiff[1,1,k] = diff[1,1]
        totdiff[1,2,k] = diff[1,2]
        totdiff[2,1,k] = diff[2,1]
        totdiff[2,2,k] = diff[2,2]
        totflux[1,k] = flux[1]
        totflux[2,k] = flux[2]

        next!(p)
    end

    outputfilename = SecularFilename(params)
    h5open(outputfilename, "w") do file
        # Basis results
        write(file,"tabJL",tabJL)
        write(file,"tabAE",tabAE)
        write(file,"totfric",totfric)
        write(file,"totdiff",totdiff)
        write(file,"totflux",totflux)
        # Parameters
        WriteParameters(file,params)
    end

    return tabAE, totfric, totdiff, totflux
end