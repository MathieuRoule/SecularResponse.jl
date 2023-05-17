
"""
    GetResonanceLines(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).
"""
function GetResonanceLines(tabωstart::Vector{ComplexF64},
                            ψ::F0,dψ::F1,d2ψ::F2,
                            coupling::BalescuLenardCoupling,
                            params::SecularParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    # check wmat directory before proceeding (save time if not.)
    CheckValidDirectory(params.secdir) || error("False directory")
    
    Linearparams = params.Linearparams
    Orbitalparams = Linearparams.Orbitalparams
    
    # Find the modes from starting values
    nmodes = length(tabωstart)
    ωModes = zeros(ComplexF64,nmodes)

    for imode = 1:nmodes
        ωModes[imode] = LR.FindDeterminantZero(tabωstart[imode],coupling.IMat,coupling.M,coupling.fht,coupling.aMcoef,coupling.tabωminωmax,Linearparams)
    end

    # Resonances (k1p,k2p)
    nbResVec, tabResVec = LR.MakeTabResVec(params.lharmonic,params.n1max,params.dimension)

    # Resonances lines store
    Kres = nbResVec*nmodes
    tabJLreslines = zeros(Float64,2,params.Kv,3,Kres)
    tabresonances = zeros(Float64,2,Kres)

    for ires = 1:nbResVec
        # Considered resonance pair
        k1, k2   = tabResVec[1,ires], tabResVec[2,ires]

        for imode = 1:nmodes
            # Considered mode
            ωM = ωModes[imode]
            ω0, γ = real(ωM), imag(ωM)

            ind = (ires-1)*nmodes+imode
            OE.GetResLineJL!(ω0-γ,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,1,ind),Orbitalparams)
            OE.GetResLineJL!(ω0,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,2,ind),Orbitalparams)
            OE.GetResLineJL!(ω0+γ,k1,k2,ψ,dψ,d2ψ,view(tabJLreslines,:,:,3,ind),Orbitalparams)
            tabresonances[1,ind], tabresonances[2,ind] = k1, k2
        end
    end

    outputfilename = SecularFilename(params,coupling.name)
    h5open(outputfilename, "r+") do file
        # Resonance vectors
        write(file,"tabResVec",tabresonances)
        # Resonance lines
        write(file,"tabJLreslines",tabJLreslines)
    end

    return tabresonances, tabJLreslines
end