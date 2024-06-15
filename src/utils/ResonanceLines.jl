
"""
    GetResonanceLines(J,lharmonic)

Secular evolution (Flux, Friction, Diffusion) at given actions (J_r,L)
on the harmonic l (decoupled harmonic numbers).
"""
function GetResonanceLines(
    tabωstart::Vector{ComplexF64},
    model::Potential,
    coupling::BalescuLenardCoupling,
    params::SecularParameters;
    store::Bool=false
)

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

        res = Resonance(k1,k2,model,Orbitalparams)

        for imode = 1:nmodes
            # Considered mode
            ωM = ωModes[imode]
            ω0, γ = real(ωM), imag(ωM)

            ind = (ires-1)*nmodes+imode
            actions_resonance_line!(view(tabJLreslines,:,:,1,ind),ω0-γ,res,model,Orbitalparams)
            actions_resonance_line!(view(tabJLreslines,:,:,2,ind),ω0,res,model,Orbitalparams)
            actions_resonance_line!(view(tabJLreslines,:,:,3,ind),ω0+γ,res,model,Orbitalparams)
            tabresonances[1,ind], tabresonances[2,ind] = k1, k2
        end
    end

    if store
        outputfilename = SecularFilename(params,coupling.name)
        h5open(outputfilename, "r+") do file
            # Resonance vectors
            write(file,"tabResVec",tabresonances)
            # Resonance lines
            write(file,"tabJLreslines",tabJLreslines)
        end
    end

    return tabresonances, tabJLreslines
end