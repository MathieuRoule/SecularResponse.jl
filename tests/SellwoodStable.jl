
using OrbitalElements
import AstroBasis
import FiniteHilbertTransform
import LinearResponse
import SecularResponse

using HDF5


##############################
# Basis
##############################
const G  = 1.

# Clutton-Brock (1972) basis
const basisname = "CluttonBrock"
const rb = 5.
const lmax,nradial = 2,10 # Usually lmax corresponds to the considered harmonics lharmonic
const basis = AstroBasis.CB72Basis(lmax=lmax,nradial=nradial,G=G,rb=rb) 

##############################
# Model Potential
##############################
const modelname = "Mestel"

const R0, V0 = 20., 1.
const ε0 = 1.0e-5
const model = TaperedMestel(R0,V0,ε0)

##############################
# Outputs directories
##############################
const dirname   = "./../data/KineticTheory/SecularResponse/2024_06_14/"
const wmatdir   = dirname
const gfuncdir  = dirname
const axidir    = dirname
const modedir   = dirname

##############################
# Model DF
##############################
const qDF = 11.4
const σDF = OrbitalElements.σMestelDF(R0,V0,qDF)
const CDF = OrbitalElements.NormConstMestelDF(G,R0,V0,qDF)

const Rin, Rout, Rmax = 1., 11.5, 20.   # Tapering radii
const ξDF = 0.5                         # Self-gravity fraction
const μDF, νDF = 5, 4                   # Tapering exponants

const dfname = "Sellwood_q_"*string(qDF)*"_xi_"*string(ξDF)*"_mu_"*string(μDF)*"_nu_"*string(νDF)

const DF(E::Float64,L::Float64) = ξDF * OrbitalElements.ZangDF(E,L,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF, νDF)

const ndFdJ(n1::Int64,n2::Int64,
            E::Float64,L::Float64,
            ndotΩ::Float64)   = ξDF * OrbitalElements.ZangndDFdJ(n1,n2,E,L,ndotΩ,R0,Rin,Rout,Rmax,V0,CDF,qDF,σDF,μDF,νDF)

#####
# Parameters
#####
# Radii for frequency truncations
const rmin = 0.1

const Orbitalparams = OrbitalParameters(;rmin=rmin)


# LinearResponse parameters

const Ku = 50    # number of u integration sample points
const FHT = FiniteHilbertTransform.LegendreFHT(Ku)
const Kv = 51    # number of allocations is directly proportional to this
const Kw = 52    # number of allocations is insensitive to this (also time, largely?

const VMAPN = 2
const ADAPTIVEKW = true
const KuTruncation=10000

const OVERWRITE = false

const lharmonic = 2
const n1max = 1  # maximum number of radial resonances to consider

const VERBOSE = 1

const Linearparams = LinearResponse.LinearResponse.LinearParameters(basis;
    Orbitalparams=Orbitalparams,
    Ω₀=frequency_scale(model),
    Ku=Ku,Kv=Kv,Kw=Kw,
    VMAPN=VMAPN,ADAPTIVEKW=ADAPTIVEKW,KuTruncation=KuTruncation,
    modelname=modelname,dfname=dfname,
    wmatdir=wmatdir,gfuncdir=gfuncdir,axidir=axidir,modedir=modedir,
    OVERWRITE=OVERWRITE,
    lharmonic=lharmonic,n1max=n1max,
    VERBOSE=VERBOSE
)

const secdir = "./../data/KineticTheory/SecularResponse/2024_06_14/"

const Secularparams = SecularResponse.SecularParameters(basis,Linearparams,secdir=secdir,n1max=1,VERBOSE=0,Kv=99)
