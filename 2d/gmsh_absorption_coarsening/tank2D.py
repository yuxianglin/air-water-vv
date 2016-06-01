from proteus import Context
from proteus import Gauges as ga
import Domain
import SpatialTools as st
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 0.5, "Height of free surface above bottom"),
    # tank
    # ("tank_dim", (5., 1.2,), "Dimensions of the tank"),
    # ("tank_sponge", (2., 2.), "Length of absorption zones (front/back, left/right)"),
    # waves
    # ("waves", False, "Generate waves (True/False)"),
    ("wave_period", 0.8, "Period of the waves"),
    ("wave_height", 0.029, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("refinement_level", 0 ,"Set maximum element diameter to he/2**refinement_level"),
    ("T", 10.0 ,"Simulation time"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.9,"Target cfl"),
    ("nsave",  10,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# general options
water_level = opts.water_level

# waves
period = opts.wave_period
height = opts.wave_height
mwl = depth = opts.water_level
direction = opts.wave_dir
wave = wt.MonochromaticWaves(period, height, mwl, depth,
                            np.array([0., -9.81, 0.]), direction)
wavelength = wave.wavelength

sim_time = period*30.

# tank options
tank_dim = (wavelength*6, water_level+water_level/3.)
tank_sponge = (wavelength, wavelength*2)
# tank_dim = opts.tank_dim
# tank_sponge = opts.tank_sponge


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank = st.Tank2D(domain, tank_dim)
tank.setSponge(left=tank_sponge[0], right=tank_sponge[1])
# if tank_sponge[0]: left = True
# if tank_sponge[1]: right = True
# if opts.waves is True:
tank.setGenerationZones(left=True, waves=wave)
# else:
#     tank.setAbsorptionZones(left=left)
tank.setAbsorptionZones(right=True)


# ----- BOUNDARY CONDITIONS ----- #


tank.BC.top.setOpenAir()
tank.BC.bottom.setNoSlip()
# if opts.waves is True:
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave, vert_axis=1)
# else:
#     tank.BC.left.setNoSlip()
tank.BC.right.setNoSlip()
tank.BC.sponge.setNonMaterial()

tank.facets = np.array([[[0, 4, 5, 3]], [[4, 6, 7, 5]], [[6, 7, 2, 1]]])
tank.facetFlags = np.array([1, 2, 3])


he = wavelength/50.
tank.MeshOptions.refineFacet([0, 1], he)
tank.MeshOptions.refineSegment([0, 6, 7, 8, 1, 5, 9], he)
tank.MeshOptions.setRefinementFunction('{0}*1.1^(abs(x-4*{1})/{0})'.format(he, wavelength))
domain.MeshOptions.LcMax = tank_dim[1]/10.
# Mesh.CharacteristicLengthMax

# refinement_level = opts.refinement_level
# he = opts.he
# domain.MeshOptions.he = he #coarse grid
st.assembleDomain(domain)
domain.writeGeo('test1')



gauge_dx=wavelength/20.
probes=np.linspace(0., tank_dim[0], tank_dim[0]/gauge_dx+1)
PG=[]
zProbes = water_level/2.
for i in probes:
    PG.append((i, zProbes, 0.),)
point_output=ga.PointGauges(gauges=((('p'),PG),),
                            activeTime=(0.,sim_time),
                            sampleRate=0,
                            fileName='point_gauges.csv')

domain.auxiliaryVariables += [point_output]








##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=False
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = sim_time
nDTout = int(opts.T*opts.nsave)
dt_out =  (T-dt_init)/nDTout
runCFL = opts.cfl

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
useRANS = 0 # 0 -- None
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

tank.BC.top.p_dirichlet = twpflowPressure_init
