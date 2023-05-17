

struct LandauMultipoleCoupling <: AbstractCoupling

    name::String 

    G::Float64

    Ku::Int64
    Kγ::Int64


    tabr ::Vector{Float64} # Positions of the first particle
    tabrp::Vector{Float64} # Positions of the second particle
    tabg ::Vector{Float64} # Prefactors for the first particle
    tabgp::Vector{Float64} # Prefactors for the second particle
end

function LandauMultipoleCoupling(;name::String="LandauMultipole",G::Float64=1.0,Ku::Int64=32,Kγ::Int64=32)
    return LandauMultipoleCoupling(name,G,Ku,Kγ,
                                   zeros(Float64,Ku),zeros(Float64,Ku),
                                   zeros(Float64,Ku),zeros(Float64,Ku))
end

function CCPrepare!(a::Float64,e::Float64,
                    Ω1::Float64,Ω2::Float64,
                    k1::Int64,k2::Int64,
                    lharmonic::Int64,
                    ω::ComplexF64,
                    ψ::F0,dψ::F1,d2ψ::F2,
                    coupling::LandauMultipoleCoupling,
                    Linearparams::LR.LinearParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    tabrtabg!(a,e,Ω1,Ω2,k1,k2,ψ,dψ,d2ψ,coupling.tabr,coupling.tabg,Linearparams.Orbitalparams)
end

function CouplingCoefficient(ap::Float64,ep::Float64,
                             Ω1p::Float64,Ω2p::Float64,
                             k1p::Int64,k2p::Int64,
                             lharmonic::Int64,
                             ψ::F0,dψ::F1,d2ψ::F2,
                             coupling::LandauMultipoleCoupling,
                             Linearparams::LR.LinearParameters)::Float64 where {F0 <: Function, F1 <: Function, F2 <: Function}

    #####
    # @ASSUMING the (k,J) part has been prepared
    #####
    tabrtabg!(ap,ep,Ω1p,Ω2p,k1p,k2p,ψ,dψ,d2ψ,coupling.tabrp,coupling.tabgp,Linearparams.Orbitalparams)

    pref = 4.0 / ((coupling.Ku*pi)^2)
    res = 0.0 
    for j = 1:coupling.Ku
        for i = 1:coupling.Ku
            Ulrirj = InteractionPotentialFT(coupling.tabr[i],coupling.tabrp[j],lharmonic,coupling.Kγ,coupling.G)
            if isnan(Ulrirj) || isinf(Ulrirj)
                #println("Infinite or Nan Ul : ri = ",coupling.tabr[i]," ; rj = ",coupling.tabr[j])
                res += 0.
            else
                res += coupling.tabg[i]*coupling.tabgp[j]*Ulrirj
            end
        end
    end
    res *= pref
    return res
end

function InteractionPotential(r::Float64,rp::Float64,
                              γ::Float64,
                              G::Float64=1.)::Float64

    return - G  / sqrt(r^2 + rp^2 - 2*r*rp*cos(γ))
end

function InteractionPotentialFT(r::Float64,rp::Float64,
                                lharmonic::Int64,
                                K::Int64,
                                G::Float64=1.)::Float64

    # u : -1 → 1
    # γ : 0 → π
    # u = 2γ/π - 1, γ = π(u+1)/2
    # dγ = π/2 du
    # 1/π ∫_0^π  dγ  = 1/2 ∫_-1^1  du
    function FTintegrand(u::Float64)::Float64
        # push integration forward on two different quantities: Θ(u),Θ(u)/r^2(u)
        γ = 0.5*pi*(u+1.0)

        return InteractionPotential(r,rp,γ,G) * cos(lharmonic*γ)
    end

    return 0.5 * OE.UnitarySimpsonIntegration(FTintegrand,K)
end

"""
    rdθdw(u,a,e)

Integrand computation for FT of interaction potential
"""
function rdθdw(w::Float64,
               a::Float64,e::Float64,L::Float64,
               Ω1::Float64,Ω2::Float64,
               ψ::F0,dψ::F1,d2ψ::F2,
               params::OE.OrbitalParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}


    # Current location of the radius, r=r(w)
    rval = OE.ru(w,a,e)

    # Current value of the radial frequency integrand (almost dθ₁₂/dw)
    gval = OE.ΘAE(ψ,dψ,d2ψ,w,a,e,params)

    # the velocity for integration (dθ1dw, dθ2dw)
    return rval, Ω1*gval, (Ω2 - L/(rval^(2)))*gval
end

function tabrtabg!(a::Float64,e::Float64,
                   Ω1::Float64,Ω2::Float64,
                   k1::Int64,k2::Int64,
                   ψ::F0,dψ::F1,d2ψ::F2,
                   tabr::Vector{Float64},tabg::Vector{Float64},
                   params::OE.OrbitalParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    # Number of middle points
    K = length(tabr)
    @assert length(tabg) == K "tabr and tabg have different length."

    # need angular momentum
    Lval = OE.LFromAE(ψ,dψ,a,e,params)

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
    r, dθ1dw, dθ2dw = rdθdw(w,a,e,Lval,Ω1,Ω2,ψ,dψ,d2ψ,params)
    dθ1_1 = 0.5*dw*dθ1dw
    dθ2_1 = 0.5*dw*dθ2dw
    # Step 2
    w += 0.25*dw
    r, dθ1dw, dθ2dw = rdθdw(w,a,e,Lval,Ω1,Ω2,ψ,dψ,d2ψ,params)
    dθ1_2 = 0.5*dw*dθ1dw
    dθ2_2 = 0.5*dw*dθ2dw
    # Step 3
    dθ1_3 = dθ1_2
    dθ2_3 = dθ2_2
    # Step 4
    w += 0.25*dw
    r, dθ1dw, dθ2dw = rdθdw(w,a,e,Lval,Ω1,Ω2,ψ,dψ,d2ψ,params)
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
    for istep=(K-1):1
        # Step 1
        # Same point as the previous 4th step
        dθ1_1 = dw*dθ1dw
        dθ2_1 = dw*dθ2dw
        # Step 2
        w += 0.5*dw
        r, dθ1dw, dθ2dw = rdθdw(w,a,e,Lval,Ω1,Ω2,ψ,dψ,d2ψ,params)
        dθ1_2 = dw*dθ1dw
        dθ2_2 = dw*dθ2dw
        # Step 3
        dθ1_3 = dθ1_2
        dθ2_3 = dθ2_2
        # Step 4
        w += 0.5*dw
        r, dθ1dw, dθ2dw = rdθdw(w,a,e,Lval,Ω1,Ω2,ψ,dψ,d2ψ,params)
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