from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import Gauges as ga
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 0.9, "Height of free surface above bottom"),
    # tank
    ("tank_dim", (5., 1.2,), "Dimensions of the tank"),
    ("tank_sponge", (2., 2.), "Length of absorption zones (front/back, left/right)"),
    # waves
    ("waves", False, "Generate waves (True/False)"),
    ("wave_period", 0.8, "Period of the waves"),
    ("wave_height", 0.029, "Height of the waves"),
    ("wave_wavelength", None, "wavelength"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    ("wave_type", 'Linear', "wave type"),
    ("Bcoeff", None, "Bcoeff"),
    ("Ycoeff", None, "Ycoeff"),
    # caisson
    ("caisson_dim", (0.3, 0.1), "Dimensions of the caisson"),
    ("caisson_coords", None, "Dimensions of the caisson"),
    ("caisson_width", 0.9, "Width of the caisson"),
    ("free_x", (0.0, 0.0, 0.0), "Translational DOFs"),
    ("free_r", (0.0, 0.0, 1.0), "Rotational DOFs"),
    ("VCG", None, "vertical position of the barycenter of the caisson"),
    ("inertia", 0.236, "Inertia of the caisson"),
    ("rotation_angle", np.pi/12., "Initial rotation angle (in radians)"),
    # gauges
    ("PG", False, "Pressure gauge"),
    ("PG_dx", None, "Pressure gauge dx"),
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("he", 0, "characteristic mesh element size"),
    ("refinement_level", 0 ,"Set maximum element diameter to he/2**refinement_level"),
    ("rans", 0, "rans"),
    ("T", 10.0 ,"Simulation time"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.9,"Target cfl"),
    ("nsave",  20,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

# waves
if opts.waves is True:
    period = opts.wave_period
    height = opts.wave_height
    wavelength = opts.wave_wavelength
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave_type = opts.wave_type
    if wave_type == 'Fenton':
        Bcoeff = opts.Bcoeff
        Ycoeff = opts.Ycoeff
    else:
        Bcoeff = None
        Ycoeff = None
    wave = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                                 g=np.array([0., -9.81, 0.]), waveDir=direction, waveType=wave_type,
                                 Ycoeff=Ycoeff, Bcoeff=Bcoeff, wavelength=wavelength)

# tank options
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# caisson options
dim = opts.caisson_dim
VCG = opts.VCG
if VCG is None:
    VCG = dim[1]/2.
free_x = opts.free_x
free_r = opts.free_r
rotation = opts.rotation_angle
if opts.caisson_coords is None:
    coords = [tank_dim[0]/2., waterLevel]
else:
    coords = opts.caisson_coords
barycenter = (coords[0], coords[1]-dim[1]/2.+VCG)
inertia = opts.inertia
width = opts.caisson_width
caisson_mass = 1  # mass is actually not used (cancelled in inertia -> see below)



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if opts.waves is True:
    tank.setGenerationZones(x_n=left, waves=wave)
else:
    tank.setAbsorptionZones(x_nleft)
tank.setAbsorptionZones(x_p=right)

caisson3D = st.Rectangle(domain, dim=dim, coords=coords)
caisson3D.setRigidBody()
caisson3D.setMass(caisson_mass)
caisson3D.setConstraints(free_x=free_x, free_r=free_r)
if rotation:
    caisson3D.rotate(rotation)  # initial position for free oscillation
caisson3D.It = inertia/caisson3D.mass/width
if opts.waves is True:
    caisson3D.setRecordValues(pos=True, rot=True, F=True, M=True,
                            filename='jung_'+str(period)+'_'+str(height)+'_'+opts.wave_type+'_rans'+str(opts.rans))
else:
    caisson3D.setRecordValues(pos=True, rot=True, F=True, M=True,
                            filename='jung_freeosc_he'+str(opts.he)+'_tank'+str(tank.dim[0]))

# ----- BOUNDARY CONDITIONS ----- #

for bc in caisson3D.BC_list:
    bc.setNoSlip()

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setNoSlip()
if opts.waves is True:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, vert_axis=1)
else:
    tank.BC['x-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()
for bc in tank.BC_list:
    bc.setFixedNodes()






from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]


#refinement_level = opts.refinement_level
#he = (caisson3D.dim[-1])/12.0*(0.5**refinement_level)
he = opts.he
ecH = 3
domain.MeshOptions.he = he #coarse grid
st.assembleDomain(domain)
#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
if nDTout > 0:
    dt_out =  (T-dt_init)/nDTout
else:
    dt_out =  0#(T-dt_init)/nDTout
runCFL = opts.cfl



# GAUGES


if opts.PG is True:
    if opts.PG_dx is None:
        if opts.wave_wavelength is not None:
            gauge_dx=opts.wave_wavelength/20.
        else:
            gauge_dx=opts.tank_sponge[0]/20.
    else:
        gauge_dx = opts.PG_dx
    probes=np.linspace(-tank_sponge[0], tank_dim[0]+tank_sponge[1], ((tank_dim[0]+tank_sponge[0]+tank_sponge[1])/gauge_dx)+1)
    PG=[]
    zProbes=waterLevel*0.5
    for i in probes:
        PG.append((i, zProbes, 0.),)

    point_output=ga.PointGauges(gauges=((('p',),PG),
                                    ),
                            activeTime = (0., T),
                            sampleRate=0.,
                            fileName='point_gauges.csv')

    domain.auxiliaryVariables['twp'] += [point_output]

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True


#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = opts.rans # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
        basis=C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,3)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
        basis=C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
        basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
        basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = max(1.0e-12,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank.dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))

tank.BC['y+'].p_dirichlet.uOfXT = twpflowPressure_init
