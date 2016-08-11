from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import MeshRefinement as mr
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 1., "Height of free surface above bottom"),
    # tank
    ("tank_x", None, "Dimensions of the tank"),
    ("tank_y", None, "Dimensions of the tank"),
    ("abs_x", 0, "Dimensions of the tank"),
    ("gen_x", 0, "Dimensions of the tank"),
    # waves
    ('waves', True, 'waves'),
    ("wave_period", 1.94, "Period of the waves"),
    ("wave_height", 0.025, "Height of the waves"),
    ("wave_dir", (1., 0., 0.), "Direction of the waves (from left boundary)"),
    # caisson
    # numerical options
    #("gen_mesh", True ,"Generate new mesh"),
    ("refinement_level", 50,"Set maximum element diameter to he/2**refinement_level"),
    ("use_gmsh", True, "gmsh"),
    ("refine_fs", True, ""),
    ("refinement_level_fs", 1., ""),
    ("T", None,"Simulation time"),
    ("he", None,"Simulation time"),
    ("grading", 1., 'grading'),
    ("refinement_level", 200., "refinement level"),
    ("dt_init", 0.001 ,"Initial time step"),
    ("cfl", 0.9 ,"Target cfl"),
    ("nsave",  10,"Number of time steps to save per second"),
    ("parallel", True ,"Run in parallel")])

# ----- CONTEXT ------ #

# general options
waterLevel = opts.water_level

# waves
if opts.waves is True:
    period = opts.wave_period
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                                g=np.array([0., -9.81, 0.]), waveDir=direction)
    wavelength = wave.wavelength


# tank options
if opts.waves is True:
    gen_length = opts.gen_x or wave.wavelength
    abs_length = opts.abs_x or 2*wave.wavelength
    tank_length = opts.tank_x or 3*wavelength
    tank_height = opts.tank_y or 1.5*opts.water_level
else:
    gen_length = opts.gen_x 
    abs_length = opts.abs_x 
    tank_length = opts.tank_x
    tank_height = opts.tank_y

tank_dim = (tank_length, tank_height)


# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()


# ----- SHAPES ----- #

tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=gen_length, x_p=abs_length)
if opts.waves is True:
    tank.setGenerationZones(x_n=True, waves=wave)
tank.setAbsorptionZones(x_p=True)


tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
if opts.waves is True:
    tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, vert_axis=1)
else:
    tank.BC['x-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()




##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]


# refinement_level = opts.refinement_level
# he = opts.he or wave.wavelength/opts.refinement_level
# ecH = 3.
# domain.MeshOptions.he = he #coarse grid


from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral

st.assembleDomain(domain)

he = he_max = wavelength/20.
he2 = height/4. #opts.he or wavelength/opts.refinement_level
# he2 = height/4.
ecH = 3
def mesh_grading(start, he, grading):
    return '{0}*{2}^(1+log((-1/{2}*(abs({1})-{0})+abs({1}))/{0})/log({2}))'.format(he, start, grading)

grading = opts.grading
tank.MeshOptions = mr.MeshOptions(tank)
# tank.MeshOptions.refineFacet([0, 1], he)
# tank.MeshOptions.setRefinementFunction('{0}*1.01^(abs(x-{1})/{0})'.format(he, tank.dim[0]))
# tank.MeshOptions.setRefinementFunction(mesh_grading(start='x-{0}'.format(tank.dim[0]), he=he, grading=grading))
# tank.MeshOptions.refineFacet([0, 1], he)
# tank.MeshOptions.setRefinementFunction('{0}*1.05^(abs(x-{1})/{0})'.format(he/4., 0))
# tank.MeshOptions.refineBox(he/15., he, x_min=-tank_sponge[0], x_max=tank.dim[0], y_min=waterLevel-0.1, y_max=waterLevel+0.1, z_min=0., z_max=0., restrict=[('facet', [0, 1])])
tank.MeshOptions.refineBox(he2, he*1000, x_min=-gen_length, x_max=tank.dim[0]+abs_length, y_min=waterLevel-height, y_max=waterLevel+height, z_min=0., z_max=0., restrict=[('facet', [0,1,2])])
tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-({0}+{1})'.format(waterLevel, height), he=he2, grading=grading))
tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-({0}-{1})'.format(waterLevel, height), he=he2, grading=grading))
# tank.MeshOptions.refineBox(he/4., he*1000, x_min=-tank_sponge[0], x_max=tank.dim[0]+tank_sponge[1], y_min=waterLevel-0.025, y_max=waterLevel+0.025, z_min=0., z_max=0., restrict=[('facet', [0])])
# tank.MeshOptions.setRefinementFunction('{0}*1.015^(abs(y-({1}+0.025))/{0})'.format(he/4., waterLevel), restrict=[('facet', [1, 2])])
# tank.MeshOptions.setRefinementFunction('{0}*1.015^(abs(y-({1}-0.025))/{0})'.format(he/4., waterLevel), restrict=[('facet', [1, 2])])
# for i in range(10):
#     tank.MeshOptions.setRefinementFunction('{0}*1.05^(sqrt((x-{3})^2+(y-({1}+0.1/({2}+1)))^2)/{0})'.format(he/4., waterLevel, i, tank.dim[0]), restrict=[('facet', [2])])
#     tank.MeshOptions.setRefinementFunction('{0}*1.05^(sqrt((x-{3})^2+(y-({1}-0.1/({2}+1)))^2)/{0})'.format(he/4., waterLevel, i, tank.dim[0]), restrict=[('facet', [2])])
# for i in range(10):
#     tank.MeshOptions.setRefinementFunction('{0}*1.05^(sqrt((x-{3})^2+(y-({1}+0.1/({2}+1)))^2)/{0})'.format(he/4., waterLevel, i, 0), restrict=[('facet', [1])])
#     tank.MeshOptions.setRefinementFunction('{0}*1.05^(sqrt((x-{3})^2+(y-({1}-0.1/({2}+1)))^2)/{0})'.format(he/4., waterLevel, i, 0), restrict=[('facet', [1])])
    # tank.MeshOptions.refineAroundPoint(coords=[tank.dim[0],0.8+0.2/10*i], lc_min=he/15, lc_max=he/5, dist_min=he/15, dist_max=he/15*1.1**5)
domain.MeshOptions.LcMax = he_max
# Mesh.CharacteristicLengthMax
domain.MeshOptions.he = he #coarse grid

mr._assembleRefinmentOptions(domain)
mr.writeGeo(domain, 'test1')
domain.MeshOptions.use_gmsh = opts.use_gmsh

from proteus import Gauges as ga
gauge_dx=0.25#wavelength/20. #0.25
nb_points = (tank.dim[0]+abs_length+gen_length)/gauge_dx
probes=np.linspace(-gen_length, tank.dim[0]+abs_length, int(nb_points))
PG0=[]
PG1=[]
for i in range(len(probes)):
   PG0.append((probes[i], waterLevel/2., 0.),)
   PG1.append(((probes[i], 0, 0.), (probes[i], tank.dim[1], 0.),))

# PG0.append((tank_dim[0]-0.00001, 0., 0.),)
# PG1.append(((tank_dim[0]-0.00001, 0, 0.), (tank_dim[0]-0.00001, tank.dim[1], 0.),))

point_output=ga.PointGauges(gauges=((('p',),PG0),
                                ),
                         #activeTime = (0., T),
                         sampleRate=0.,
                         fileName='point_gauges.csv')

domain.auxiliaryVariables['twp'] += [point_output]



fields = ('vof',)

columnGauge = ga.LineIntegralGauges(gauges=((fields, PG1),),
                                fileName='column_gauge.csv')
domain.auxiliaryVariables['vof'] += [columnGauge]

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
T = opts.T or 100*period
nDTout = int(T*opts.nsave)
if opts.nsave != 0:
    dt_out =  (T-dt_init)/nDTout
else:
    dt_out = 0
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
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.9
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
