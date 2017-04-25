from proteus import Domain, Context, AuxiliaryVariables
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # still water
    ("water_level", 4., "Height of free surface above seabed "),
    ("water_depth", 4., " water depth"),
    # tank
    ("tank_dim", (20., 15.,), "Dimensions of the tank"),
    ("tank_sponge", (0., 10.), "Length of relaxation zones zones (left, right)"),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", False, "Places Gauges in tank (5 per wavelength)"),
    #obst
    ("obst_dim",(5.,2.),"Dimensions of the oil tank"),
    ("obst_center",10.,"obstacle center at bottom"),
    # waves
    ("waves", False, "Generate waves (True/False)"),
    ("wave_height", 2., "Height of the waves"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("g", (0.,-9.81,0.), "Direction of the waves (from left boundary)"),
    ("trans",-1. ,"peak offset for solitary wave"),
    # mesh refinement
    ("he", 0.2, "Set characteristic element size"),
    # numerical options
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("T", 10.001, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave",  20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("movingDomain",True,"Switch on moving domain"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# waves
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = np.array(opts.wave_dir)
    g = np.array(opts.g)
    trans  = opts.trans * direction/sqrt(sum(direction[:]*direction[:]))
    wave = wt.SolitaryWave(waveHeight = height,
                           mwl = mwl,
                           depth = depth,
                           g = g,
                           waveDir = direction,
                           trans = trans)
else:
    wave=None

# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge
movingDomain = opts.movingDomain
# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options
boundaries=['left','right','top','bottom','obst','sponge']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
L=tank_dim
xSponge_1 = - tank_sponge[0]
xSponge_2=L[0] + tank_sponge[1]
obst_diameter=opts.obst_dim[0]
obst_height  =opts.obst_dim[1]
obst_center  =opts.obst_center
obst_point1 = obst_center-obst_diameter/2
obst_point2 = obst_center+obst_diameter/2
vertices=[[0.0,0.0],#0
          [obst_point1,0.0],#1
          [obst_point1,obst_height],#2
          [obst_point2,obst_height],#3
          [obst_point2,0.0],#4
          [L[0],0.0],#5
          [L[0],L[1]],#6
          [0.0,L[1]]#7
          ]

vertices=[[0.0,0.0],#0
          [obst_point1,0.0],#1
          [obst_point1,obst_height],#2
          [obst_point2,obst_height],#3
          [obst_point2,0.0],#4
          [L[0],0.0],#5
          [xSponge_2,0.0],#6
          [xSponge_2,L[1]],#7
          [L[0],L[1]],#8
          [0.0,L[1]]#9
          ]

vertexFlags=[
    boundaryTags['bottom'],#1
    boundaryTags['bottom'],#2
    boundaryTags['obst'],#3
    boundaryTags['obst'],#4
    boundaryTags['bottom'],#5
    boundaryTags['bottom'],#6
    boundaryTags['bottom'],#7
    boundaryTags['top'],#8
    boundaryTags['top'],#9
    boundaryTags['top']#10

]
segments=[]

for vb in range(len(vertices)-1):
    segments.append([vb,vb+1])
segments.append([9,0])
segments.append([5,8])


segmentFlags=[boundaryTags['bottom'],
              boundaryTags['obst'],
              boundaryTags['obst'],
              boundaryTags['obst'],
              boundaryTags['bottom'],
              boundaryTags['bottom'],
              boundaryTags['right'],
              boundaryTags['top'],
              boundaryTags['top'],
              boundaryTags['left'],
              boundaryTags['sponge']]



xRelaxCenter_1 = xSponge_1/2
xRelaxCenter_2 = (xSponge_2+L[0])/2
xCenter=L[0]/2
#regions=[[xCenter, 0.5*L[1]],
#         [xRelaxCenter_2, 0.5*L[1]]]
#regionFlags=[1,2]
regions=[[xCenter, 0.5*L[1]],
         ]
regionFlags=[1]

BCTags={       'left' : boundaryTags['left'],
               'right' : boundaryTags['right'],
               'bottom' : boundaryTags['bottom'],
               'top' : boundaryTags['top']  }
#'obst': 7,
#'sponge' : 8}
boundaryOrientations={'left' : np.array([-1., 0., 0.]),
                      'right' : np.array([ 1., 0., 0.]),
                      'bottom' : np.array([ 0., -1.,0.]),
                      'top' : np.array([ 0., 1., 0.])}

# ----- SHAPES ----- #
tank = st.CustomShape(domain,vertices=vertices,vertexFlags=vertexFlags,segments=segments,segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=BCTags, boundaryOrientations=boundaryOrientations)

#tank.BC['sponge'] =tank.BC_class(shape = tank, name = 'sponge')
tank.BC['obst'] =tank.BC_class(shape = tank, name = 'obst')
#tank.BC_list.append(tank.BC['sponge'])
tank.BC_list.append(tank.BC['obst'])
#regionIndex={'left':0,'right':1,'tank':2}


left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if opts.waves is True:
    smoothing = 0
    tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing)

#if left:
#    if opts.waves is True:
#tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
#        tank.setGenerationZones(flags=regionIndex['left'],epsFact_solid=xSponge_1/2.,center=np.array([xRelaxCenter_1,L[1]/2.]),orientation=boundaryOrientations['left'], waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
#    else:
#		pass

#if right:
#    tank.setAbsorptionZones(flags=regionIndex['right'], epsFact_solid=xSponge_2/2., center= np.array([xRelaxCenter_2, L[1]/2.]), orientation=boundaryOrientations['right'])



# ----- BOUNDARY CONDITIONS ----- #
#tank.BC['top'].setNoSlip()
tank.BC['top'].setAtmosphere()
#tank.BC['top'].setNoSlip()
tank.BC['right'].setNoSlip()
tank.BC['left'].setNoSlip()
tank.BC['bottom'].setFreeSlip()
#tank.BC['sponge'].setNonMaterial()
tank.BC['obst'].setNoSlip()

#tank.BC['obst'].setFixedNodes()
#tank.BC['right'].setFixedNodes()
#tank.BC['bottom'].setFixedNodes()
#tank.BC['top'].setFixedNodes()
#tank.BC['left'].setFixedNodes()
#for bc in tank.BC_list:
#    bc.setFixedNodes()
#for i in range(4):
#	tank.BC_list[i].setTank()
scriptmotion=True
class rigidpadal(AuxiliaryVariables.AV_base):
    def __init__(self,dt_init=0.001, scriptmotion=None):
        self.dt_init=dt_init
        self.scriptmotion=scriptmotion
        self.last_position=np.array([0.0,0.0,0.0])
        self.h=np.array([0.0,0.0,0.0])
        self.T=-dt_init

    def attachModel(self, model,ar):
        self.model=model
        self.ar=ar
        #import pdb;pdb.set_trace()
        self.writer=Archiver.XdmfWriter()
        return self

    def calculate_init(self):
        self.calculate()

    def calculate(self):
        import copy
        try:
            dt=self.model.levelModelList[-1].dt_last
        except:
            dt=self.dt_init
        if scriptmotion:
            self.T += dt
            self.new_position=np.array([self.scriptmotion(self.T),self.scriptmotion(self.T),0.0])
            self.h=self.new_position-self.last_position
            print "self.h=", self.h
            self.last_position=copy.deepcopy(self.new_position)
A=0.3
w=np.pi
motion=lambda t: A*np.sin(w*t)

padal=rigidpadal(dt_init=opts.dt_init,scriptmotion=motion)

quad_order=3


# ----- GAUGES ----- #

if opts.gauge_output:
    if left or right:
        gauge_dx = tank_sponge[0]/10.
    else:
        gauge_dx = tank_dim[0]/10.
    probes=np.linspace(-tank_sponge[0], tank_dim[0]+tank_sponge[1], (tank_sponge[0]+tank_dim[0]+tank_sponge[1])/gauge_dx+1)
    PG=[]
    PG2=[]
    LIG = []
    zProbes=waterLevel*0.5
    for i in probes:
        PG.append((i, zProbes, 0.),)
        PG2.append((i, waterLevel, 0.),)
        if i == probes[0]:
            LIG.append(((i, 0.+0.0001, 0.),(i, tank_dim[1]-0.0001,0.)),)
        elif i != probes[0]:
            LIG.append(((i-0.0001, 0.+0.0001, 0.),(i-0.0001, tank_dim[1]-0.0001,0.)),)
    tank.attachPointGauges(
        'twp',
        gauges = ((('p',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_pressure.csv'
    )
    tank.attachPointGauges(
        'ls',
        gauges = ((('phi',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_levelset.csv'
    )

    tank.attachLineIntegralGauges(
        'vof',
        gauges=((('vof',), LIG),),
        activeTime = (0., opts.T),
        sampleRate = 0,
        fileName = 'lineGauge.csv'
    )

# ----- ASSEMBLE DOMAIN ----- #

domain.MeshOptions.genMesh = opts.gen_mesh
domain.MeshOptions.he = opts.he
st.assembleDomain(domain)

# ----- REFINEMENT OPTIONS ----- #







##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]
mwl = depth = opts.water_level



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
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

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
useRANS = opts.useRANS # 0 -- None
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
sc = 0.5 # default: 0.5. Test: 0.25
sc_beta = 1.5 # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 1.0 # default: 1.0. Test: 0.1
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
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

ns_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-12,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-12,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-12,opts.mesh_tol*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    if opts.waves:
        wl = waterLevel + wave.eta(x,t)
    else:
        wl = waterLevel
    phi_L = tank_dim[1] - wl
    phi = x[1] - wl
    return p_L - g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                              -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
