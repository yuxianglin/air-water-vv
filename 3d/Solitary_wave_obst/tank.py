from proteus import Domain, Context, AuxiliaryVariables
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # still water
    ("water_level", 0.25, "Height of free surface above seabed "),
    ("water_depth", 0.25, " water depth"),
    # tank
    ("tank_dim", (1., 0.9, 0.75), "Dimensions of the tank"),
    ("tank_sponge", (0.5,1.), "Length of relaxation zones zones (left, right)"),
    ("tank_BC", 'freeslip', "boundary type of the tank"),
    ("gauge_output", False, "Places Gauges in tank (5 per wavelength)"),
	#obst
    ("obst_dim",(0.3,0.6),"Dimensions of the oil tank"),
	("obst_center",(0.5,0.45),"obstacle center at bottom"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_height", 0.1, "Height of the waves"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("g", (0.,0.,-9.81), "Direction of the waves (from left boundary)"),
    ("trans",-0.2 ,"peak offset for solitary wave"),
    # mesh refinement
    ("he", 0.01, "Set characteristic element size"),
    # numerical options
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("T", 5.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave",  20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("movingDomain",False,"Switch on moving domain"),
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





# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge
movingDomain = opts.movingDomain
# ----- DOMAIN ----- #

#domain = Domain.PlanarStraightLineGraphDomain()
domain = Domain.PiecewiseLinearComplexDomain()
# caisson options


boundaries=['left','right','bottom','top','front','back','sponge','obst']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

L=tank_dim
xSponge_1 = - tank_sponge[0]
xSponge_2=L[0] + tank_sponge[1]
obst_diameter=opts.obst_dim[0]
obst_height  =opts.obst_dim[1]
obst_center  =opts.obst_center

vertices=[[xSponge_1,0.0,0.0],#0
          [0.0,0.0,0.0],#1
          [L[0],0.0,0.0],#2
          [xSponge_2,0.0,0.0],#3
          [xSponge_2,L[1],0.0],#4
          [L[0],L[1],0.0],#5
          [0.0,L[1],0.0],#6
          [xSponge_1,L[1],0.0]]#7


vertexFlags=[boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom']]


for v,vf in zip(vertices,vertexFlags):
    vertices.append([v[0],v[1],L[2]])
    vertexFlags.append(boundaryTags['top'])
#adding obst
from math import ceil,pi,sin,cos
radius = obst_diameter/2.0
vbStart=len(vertices)
#points_on_circle=4*int(ceil(0.5*pi*(radius)/he))
dx=opts.he
points_on_circle=int(pi*obst_diameter/dx)
for cb in range(points_on_circle):
	vertices.append([obst_center[0]+radius*sin(float(cb)/float(points_on_circle)*2.0*pi),
		             obst_center[1]+radius*cos(float(cb)/float(points_on_circle)*2.0*pi),0.0])
	vertexFlags.append(boundaryTags['bottom'])
vtStart=len(vertices)
for cb in range(points_on_circle):
	vertices.append([obst_center[0]+radius*sin(float(cb)/float(points_on_circle)*2.0*pi),
		             obst_center[1]+radius*cos(float(cb)/float(points_on_circle)*2.0*pi),obst_height])
	vertexFlags.append(boundaryTags['obst'])
#print vertices
#print vertexFlags

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,0],
          [1,6],
          [2,5]]

segmentFlags=[boundaryTags['front'],
             boundaryTags['front'],
             boundaryTags['front'],
             boundaryTags['right'],
             boundaryTags['back'],
             boundaryTags['back'],
             boundaryTags['back'],
             boundaryTags['left'],
             boundaryTags['sponge'],
             boundaryTags['sponge'] ]


facets=[]
facetFlags=[]
facetHoles=[]

for s,sF in zip(segments,segmentFlags):
    facets.append([[s[0],s[1],s[1]+8,s[0]+8]])
    facetFlags.append(sF)
    facetHoles.append([])
bf=[[0,1,6,7],[2,3,4,5]]
tf=[]
for i in range(len(bf)):
    facets.append([bf[i]])
    facetFlags.append(boundaryTags['bottom'])
    facetHoles.append([])
    tf=[ss + 8 for ss in bf[i]]
    facets.append([tf])
    facetFlags.append(boundaryTags['top'])
    facetHoles.append([])
facets.append([[9,10,13,14]])
facetFlags.append(boundaryTags['top'])
facetHoles.append([])

obst_bottom_facet=[vbStart+cb for cb in range(points_on_circle)]
obst_top_facet=[vtStart+cb for cb in range(points_on_circle)]
facets.append([[1,2,5,6],obst_bottom_facet])
facetFlags.append(boundaryTags['bottom'])
facetHoles.append([[obst_center[0],obst_center[1],0.0]])

facets.append([obst_top_facet])
facets.append([obst_bottom_facet])
facetFlags.append(boundaryTags['obst'])#tank top
facetFlags.append(boundaryTags['obst'])#tank bottom
facetHoles.append([])
facetHoles.append([])
#sides of cylinder obst
for fN in range(len(obst_bottom_facet)):
    facets.append([[obst_bottom_facet[fN-1],obst_bottom_facet[fN],obst_top_facet[fN],obst_top_facet[fN-1]]])
    facetFlags.append(boundaryTags['obst'])
    facetHoles.append([])
holes=[[obst_center[0],obst_center[1],obst_height/2]]




#for i in range(0,3):
#	facetFlags.append(boundaryTags['bottom'])
#	facetFlags.append(boundaryTags['top'])

#print facets
#print facetFlags
xRelaxCenter_1 = xSponge_1/2
xRelaxCenter_2 = (xSponge_2+L[0])/2
xCenter=L[0]/2
regions=[[xRelaxCenter_1, 0.5*L[1], 0.5*L[2]],
         [xRelaxCenter_2, 0.5*L[1], 0.5*L[2]],
         [xCenter,        0.5*L[1], 0.5*L[2]]]
regionFlags=[1,2,3]
regionIndex={'left':1,'right':2,'tank':3}

BCTags={'left':boundaryTags['left'],
	   'right':boundaryTags['right'],
	   'front':boundaryTags['front'],
	    'back':boundaryTags['back'], 
	  'bottom':boundaryTags['bottom'], 
	    'top' :boundaryTags['top']}
              #'obst': 7,
              #'sponge' : 8}
boundaryOrientations={'left' :np.array([-1., 0., 0.]),
		              'right':np.array([ 1., 0., 0.]),
					  'front':np.array([ 0.,-1., 0.]),
                      'back' :np.array([ 0., 1., 0.]),
					 'bottom':np.array([ 0., 0.,-1.]),
					  'top'  :np.array([ 0., 0., 1.])}



# ----- SHAPES ----- #
#tank = st.Tank2D(domain, tank_dim)
tank = st.CustomShape(domain,vertices=vertices,vertexFlags=vertexFlags,facets=facets,facetFlags=facetFlags,
                      facetHoles=facetHoles, holes=holes, regions=regions, regionFlags=regionFlags,
                      boundaryTags=BCTags, boundaryOrientations=boundaryOrientations)
#tank = st.CustomShape(domain,vertices=vertices,vertexFlags=vertexFlags,facets=facets,facetFlags=facetFlags,
#                      holes=holes, regions=regions, regionFlags=regionFlags,
#                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)
#tank.facetHoles=np.array(facetHoles)
#tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
tank.BC['sponge'] =tank.BC_class(shape = tank, name = 'sponge')
tank.BC['obst'] =tank.BC_class(shape = tank, name = 'obst')
tank.BC_list.append(tank.BC['sponge'])
tank.BC_list.append(tank.BC['obst'])



left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if opts.waves is True:
    smoothing = opts.he*0.
    tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing,vert_axis=2)
if left:
    if opts.waves is True:
        #import pdb;pdb.set_trace()
        #tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
        tank.setGenerationZones(flags=regionIndex['left'],epsFact_solid=xSponge_1/2.,center=np.array([xRelaxCenter_1,L[1]/2.,L[2]/2.]),orientation=boundaryOrientations['left'], waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
    else:
		pass
		#tank.setAbsorptionZones(x_n=left, dragAlpha=0.5/1.004e-6)
else:
    tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(flags=regionIndex['right'], epsFact_solid=xSponge_2/2., center= np.array([xRelaxCenter_2, L[1]/2.,L[2]/2.]), orientation=boundaryOrientations['right'])

# ----- BOUNDARY CONDITIONS ----- #

tank.BC['top'].setAtmosphere()
if opts.tank_BC == 'noslip':
    tank.BC['front'].setNoSlip()
    tank.BC['back'].setNoSlip()
if opts.tank_BC == 'freeslip':
    tank.BC['front'].setFreeSlip()
    tank.BC['back'].setFreeSlip()
tank.BC['right'].setNoSlip()
tank.BC['bottom'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()
tank.BC['obst'].setFreeSlip()
for bc in tank.BC_list:
    bc.setFixedNodes()
for i in range(6):
	tank.BC_list[i].setTank()



FuncType='Lagrange'
Order=1
from SolidSolver import ElastoDynCylinder
class ElasticBar(AuxiliaryVariables.AV_base):
#class ElasticBar():
    def __init__(self,density=7800,E=207e9,nu=0.3,D=1.22,H=1.92,loc=np.array([10.,10.])):
        self.dt_init=opts.dt_init
        self.he=opts.he
        self.E=E
        self.nu=nu
        self.loc=loc
        self.density=density
        self.D=D
        self.H=H
        self.Solver=ElastoDynCylinder(E=E, nu=nu, g=opts.g, density=density, D=D, H=H, loc=loc)
        self.Solver.GenerateMesh(opts.he)
        self.Solver.set_dt(self.dt)
        self.Solver.set_functional(FuncType,Order)


    def attachModel(self,model,ar):
        self.model=model
        self.ar=ar
        self.writer=Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        #m = self.model.levelModelList[-1]
        return self

    def calculate(self):
        nExtEleBoundary =self.model.levelModelList[-1].mesh.nExteriorElementBoundaries_global

        EleBoundaryMTInd=self.model.levelModelList[-1].mesh.elementBoundaryMaterialTypes
        ExtEleBoundaryArray=self.model.levelModelList[-1].mesh.exteriorElementBoundariesArray
        #obstInd=[count if EleBoundaryMTInd[ind]==boundaryTags['obst'] for count, ind in enumerate(ExtEleBoundaryArray)]
        obstInd=[count for count, ind in enumerate(ExtEleBoundaryArray) if EleBoundaryMTInd[ind]==boundaryTags['obst'] ]
         
        
        Force_x=self.model.levelModelList[-1].coefficients.forces_x[obstInd]
        Force_y=self.model.levelModelList[-1].coefficients.forces_y[obstInd]
        Force_z=self.model.levelModelList[-1].coefficients.forces_z[obstInd]
        Force_coords=self.model.levelModelList[-1].ebqe['x'][obstInd]
        self.Solver.set_force(Force_x,Force_y,Force_z,Force_coords)

        self.Solver.set_variationalform()
        
        #zero=Constant(('0.0','0.0','0.0'))
        zero=('0.0','0.0','0.0')
        bottom=lambda x: x[2] < DOLFIN_EPS
        BC={bottom:zero}

        self.Solver.set_DBC(BC)

        self.Solver.solve()









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
g = opts.g




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
nd = 3
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
sc = 0.25 # default: 0.5. Test: 0.25
sc_beta = 1. # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 0.1 # default: 1.0. Test: 0.1
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
    wl = waterLevel + wave.eta(x,t)
    phi_L = tank_dim[nd-1] - wl
    phi = x[nd - 1] - wl
    print
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
