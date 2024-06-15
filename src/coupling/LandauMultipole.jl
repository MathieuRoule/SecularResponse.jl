

struct LandauMultipoleCoupling <: AbstractCoupling

    G::Float64  # Gravitational coupling strength
    ε::Float64  # Gravitational coupling softening length

    Kw::Int64

    tabr::Vector{Float64} # Positions of the first particle
    tabrp::Vector{Float64} # Positions of the second particle
    tabg::Vector{Float64} # Prefactors for the first particle
    tabgp::Vector{Float64} # Prefactors for the second particle
end

function LandauMultipoleCoupling(;G::Float64=1.0, ε::Float64=1.0e-3, Kw::Int64=32)
    return LandauMultipoleCoupling(
        G,
        ε,
        Kw,
        zeros(Float64,Kw),
        zeros(Float64,Kw),
        zeros(Float64,Kw),
        zeros(Float64,Kw)
    )
end

function CCPrepare!(
    a::Float64,
    e::Float64,
    Ω1::Float64,
    Ω2::Float64,
    k1::Int64,
    k2::Int64,
    lharmonic::Int64,
    ω::ComplexF64,
    model::Potential,
    coupling::LandauMultipoleCoupling,
    Linearparams::LR.LinearParameters
)
    tabr, tabg = coupling.tabr, coupling.tabg
    Orbitalparams = Linearparams.Orbitalparams
    tabrtabg!(a, e, Ω1, Ω2, k1, k2, model, tabr, tabg, Orbitalparams)
end

function CouplingCoefficient(
    ap::Float64,
    ep::Float64,
    Ω1p::Float64,
    Ω2p::Float64,
    k1p::Int64,
    k2p::Int64,
    lharmonic::Int64,
    model::Potential,
    coupling::LandauMultipoleCoupling,
    Linearparams::LR.LinearParameters
)
    #####
    # @ASSUMING the (k,J) part has been prepared
    #####
    tabrp, tabgp = coupling.tabrp, coupling.tabgp
    Orbitalparams = Linearparams.Orbitalparams
    tabrtabg!(ap, ep, Ω1p, Ω2p, k1p, k2p, model, tabrp, tabgp, Orbitalparams)

    pref = 4.0 / ((coupling.Kw*pi)^2)
    res = 0.0 
    for j = 1:coupling.Kw
        for i = 1:coupling.Kw
            Ulrirj = Ul(coupling.tabr[i],coupling.tabrp[j],lharmonic,coupling.G,coupling.ε)
            if isfinite(Ulrirj)
                res += coupling.tabg[i]*coupling.tabgp[j]*Ulrirj
            end
        end
    end
    res *= pref
    return res
end

"""
    Regularized hypergeometric function pF̃q([1/2,1/2,1], [1-l,1+l], 2a/(1+a))
"""
function regularized_pFq_special(l::Int64, a::Float64)

    @assert -1 < a <= 1 "a = $a must lie in (-1,1]"

    b = 2a/(1+a)
    eE, eK = ellipe(b), ellipk(b)

    if l == 0
        return 2*eK/pi
    elseif abs(l) == 1
        return 2*(eK - (1+a)*eE)/(pi*a)
    elseif abs(l) == 2
        return 2*((4-a^2)*eK - 4*(1+a)*eE)/(3*pi*(a^2))
    elseif abs(l) == 3
        return 2*((32-17a^2)*eK - (32-9a^2)*(1+a)*eE)/(15*pi*a^3)
    elseif abs(l) == 4
        return 2*((384-304a^2+25a^4)*eK - 16*(24-13a^2)*(1+a)*eE)/(105*pi*a^4)
    elseif abs(l) == 5
        return 2*((2048-2144a^2+411a^4)*eK - (2048-1632a^2+147a^4)*(1+a)*eE)/(315*pi*a^5)
    else
        error("Unknown regularized_pFq_special for l = ",l)
    end
end
"""
    Fourier tranform in (configuration) angle of the interaction potential
"""
function Ul(r::Float64, rp::Float64, lharmonic::Int64, G::Float64=1.,ε::Float64=1.e-5)

    if (r<=0.) || (rp<= 0.)
        return 0.
    end
    rm2 = r^2 + rp^2 + ε^2
    rm = sqrt(rm2)
    a = 2r*rp/(rm2)

    @assert abs(a) < 1 "a = $a must lie in ]-1,1[, does not work for r = $r and rp = $rp"

    return - (G / rm) * regularized_pFq_special(lharmonic,a) / sqrt(1+a)
end

"""
    rdθdw(u,a,e)

Integrand computation for FT of interaction potential
"""
function rdθdw(
    w::Float64,
    a::Float64,
    e::Float64,
    model::Potential,
    params::OrbitalParameters;
    kwargs...
)
    # Current location of the radius, r=r(w)
    rval = radius_from_anomaly(w,a,e,model,params)

    # Current value of the radial frequency integrand (almost dθ₁₂/dw)
    dθ1dw, dθ2dw = angles_gradient(w, a, e, model, params; kwargs...)

    # the velocity for integration (dθ1dw, dθ2dw)
    return rval, dθ1dw, dθ2dw
end

function tabrtabg!(
    a::Float64,
    e::Float64,
    Ω1::Float64,
    Ω2::Float64,
    k1::Int64,
    k2::Int64,
    model::Potential,
    tabr::Vector{Float64},
    tabg::Vector{Float64},
    params::OrbitalParameters
)

    # Number of middle points
    K = length(tabr)
    @assert length(tabg) == K "tabr and tabg have different length."

    # need angular momentum
    _, Lval = EL_from_ae(a,e,model,params)

    # Caution : Reverse integration (lower error at apocenter than pericenter)
    # -> Result to multiply by -1
    dw = -(2.0)/(K)

    #####
    # Warm-up for first point (Half a usual step)
    #####
    # Initialise the state vectors: w, θ1, (θ2-psi)
    # Reverse integration, starting at apocenter
    w, θ1, θ2 = 1.0, pi, 0.0
    # Step 1
    r, dθ1dw, dθ2dw = rdθdw(w, a, e, model, params; L=Lval, Ω1=Ω1, Ω2=Ω2)
    dθ1_1 = 0.5*dw*dθ1dw
    dθ2_1 = 0.5*dw*dθ2dw
    # Step 2
    w += 0.25*dw
    r, dθ1dw, dθ2dw = rdθdw(w, a, e, model, params; L=Lval, Ω1=Ω1, Ω2=Ω2)
    dθ1_2 = 0.5*dw*dθ1dw
    dθ2_2 = 0.5*dw*dθ2dw
    # Step 3
    dθ1_3 = dθ1_2
    dθ2_3 = dθ2_2
    # Step 4
    w += 0.25*dw
    r, dθ1dw, dθ2dw = rdθdw(w, a, e, model, params; L=Lval, Ω1=Ω1, Ω2=Ω2)
    dθ1_4 = 0.5*dw*dθ1dw
    dθ2_4 = 0.5*dw*dθ2dw
    #####
    # Better guess for the value of theta(u=1.0+dw/2)
    # with a first RK4 warm-up integration for a small step    
    θ1 += (dθ1_1 + 2.0*dθ1_2 + 2.0*dθ1_3 + dθ1_4)/(6.0)
    θ2 += (dθ2_1 + 2.0*dθ2_2 + 2.0*dθ2_3 + dθ2_4)/(6.0)
    # Store the values
    tabr[K] = r
    tabg[K] = dθ1dw * cos(k1*θ1 + k2*θ2)

    # Reverse integration
    for istep=(K-1):-1:1
        # Step 1
        # Same point as the previous 4th step
        dθ1_1 = dw*dθ1dw
        dθ2_1 = dw*dθ2dw
        # Step 2
        w += 0.5*dw
        r, dθ1dw, dθ2dw = rdθdw(w, a, e, model, params; L=Lval, Ω1=Ω1, Ω2=Ω2)
        dθ1_2 = dw*dθ1dw
        dθ2_2 = dw*dθ2dw
        # Step 3
        dθ1_3 = dθ1_2
        dθ2_3 = dθ2_2
        # Step 4
        w += 0.5*dw
        r, dθ1dw, dθ2dw = rdθdw(w, a, e, model, params; L=Lval, Ω1=Ω1, Ω2=Ω2)
        dθ1_4 = dw*dθ1dw
        dθ2_4 = dw*dθ2dw
        # Update the positions using RK4-like sum (Simpson's 1/3 rule)
        θ1 += (dθ1_1 + 2.0*dθ1_2 + 2.0*dθ1_3 + dθ1_4)/(6.0)
        θ2 += (dθ2_1 + 2.0*dθ2_2 + 2.0*dθ2_3 + dθ2_4)/(6.0)
        # Store the values
        tabr[istep] = r
        tabg[istep] = dθ1dw * cos(k1*θ1 + k2*θ2)
    end
end