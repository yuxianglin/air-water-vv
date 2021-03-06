from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *   
from proteus.Profiling import logEvent
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import Gauges
from proteus.Gauges import PointGauges,LineGauges,LineIntegralGauges

#wave generator
windVelocity = (0.0,0.0,0.0)
inflowHeightMean = 1.0
currentVelocity=0.5
inflowVelocityMean = (cos(math.pi/4.0)*currentVelocity,cos(math.pi/4.0)*currentVelocity,0.0)
period = 1.94
omega = 2.0*math.pi/period
#waveheight = 0.020  #=0.025 without the current
#amplitude = waveheight/ 2.0
wavelength = 5.912
k = 2.0*math.pi/wavelength
rampTime=2.0*period
meanFrameVelocity= 2.694 #calculated from FFT
outflowHeight=inflowHeightMean
netcurrentVelocity=wavelength/period-meanFrameVelocity

#print inflowVelocityMean

#  Discretization -- input options  

genMesh=True
movingDomain=False
applyRedistancing=True
useOldPETSc=False
useSuperlu=False
timeDiscretization='be'#'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
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
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)     	 
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3) 	    
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
    
# Domain and mesh
L = (float(6.0*wavelength), float(3.0*wavelength), 1.50)

he = wavelength/100

GenerationZoneLength = wavelength*0.5
AbsorptionZoneLength= wavelength*1.0
spongeLayer = True #False  
levee=spongeLayer
slopingSpongeLayer=spongeLayer

ySponge_1 = GenerationZoneLength
ySponge_2 = L[1]-GenerationZoneLength
xSponge_1 = L[0]-AbsorptionZoneLength

yRelaxCenter_1 = ySponge_1/2.0
yRelaxCenter_2 = (ySponge_2+L[1])/2.0
xRelaxCenter_1 = (xSponge_1+L[0])/2.0

epsFact_solid = GenerationZoneLength/2.0
epsFact_solid_2 = AbsorptionZoneLength/2.0

nLevels = 1
weak_bc_penalty_constant = 100.0
quasi2D=False
if quasi2D:#make tank one element wide
    L = (L[0],he,L[2])

#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

structured=False  

gauge_dx=float(wavelength/10)
PGL=[]
LGL=[]
for i in range(0,int(L[0]/gauge_dx)): #+1 only if gauge_dx is an exact 
  PGL.append([gauge_dx*i,(ySponge_1+ySponge_2)/2.0,0.5])
  LGL.append([(gauge_dx*i,(ySponge_1+ySponge_2)/2.0,0),(gauge_dx*i,(ySponge_1+ySponge_2)/2.0,L[2])])
 

gaugeLocations=tuple(map(tuple,PGL)) 
columnLines=tuple(map(tuple,LGL)) 


pointGauges = PointGauges(gauges=((('u','v'), gaugeLocations),
                                (('p',),    gaugeLocations)),
                  activeTime = (0, 1000.0),
                  sampleRate = 0,
                  fileName = 'combined_gauge_0_0.5_sample_all.txt')


fields = ('vof',)

columnGauge = LineIntegralGauges(gauges=((fields, columnLines),),
                                 fileName='column_gauge.csv')

#lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_y=0':((0.0,0.0,0.0),(L[0],0.0,0.0))},linePoints=24)

#lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)


if useHex:   
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True    
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['empty','left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
        domain = Domain.RectangularDomain(L)
    elif spongeLayer:
        vertices=[[0.0,0.0,0.0],#0
                  [xSponge_1,0.0,0.0],#1
                  [L[0],0.0,0.0],#2
                  [L[0],ySponge_1,0.0],#3
                  [L[0],ySponge_2,0.0],#4
                  [L[0],L[1],0.0],#5
                  [xSponge_1,L[1],0.0],#6
                  [0.0,L[1],0.0],#7
                  [0.0,ySponge_2,0.0],#8
                  [0.0,ySponge_1,0.0],#9
                  [xSponge_1,ySponge_1,0.0],#10
                  [xSponge_1,ySponge_2,0.0]]#11

        vertexFlags=[]
        for i in range(0,int(len(vertices))):
            vertexFlags.append(boundaryTags['bottom'])
            

        for v,vf in zip(vertices,vertexFlags):
            vertices.append([v[0],v[1],L[2]])
            vertexFlags.append(boundaryTags['top'])

        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,4],
                  [4,5],
                  [5,6],
                  [6,7],
                  [7,8],
                  [8,9],
                  [9,0],
                 [9,10],
                 [10,3],
                 [8,11],
                 [11,4],
                 [1,10],
                 [10,11],
                 [11,6]]
                 
        segmentFlags=[boundaryTags['front'],
                     boundaryTags['front'],                                  
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['right'],
                     boundaryTags['back'],
                     boundaryTags['back'],               
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['left'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'],
                     boundaryTags['empty'] ]
        

        facets=[]
        facetFlags=[]
        #print int(len(vertices))
        for s,sF in zip(segments,segmentFlags):
            facets.append([[s[0],s[1],s[1]+int(len(vertices)/2),s[0]+int(len(vertices)/2)]])
            facetFlags.append(sF)

        bf=[[0,1,10,9],[1,2,3,10],[10,3,4,11],[11,4,5,6],[8,11,6,7],[9,10,11,8]]
        tf=[]
 
        for i in range(0,int(len(bf))):
         facets.append([bf[i]])
         tf=[ss + int(len(vertices)/2) for ss in bf[i]]
         facets.append([tf])

        for i in range(0,int(len(bf))):
         facetFlags.append(boundaryTags['bottom'])
         facetFlags.append(boundaryTags['top'])

        #print facets
        #print facetFlags

        regions=[[0.01*L[0], yRelaxCenter_1 ,0.0],#1
                 [xRelaxCenter_1, yRelaxCenter_1 ,0.0],#2
                 [xRelaxCenter_1, 0.5*L[1], 0.0],#3
                 [xRelaxCenter_1, yRelaxCenter_2 ,0.0],#4
                 [0.01*L[0], yRelaxCenter_2 ,0.0],#5
                 [0.5*L[0],0.5*L[1], 0.0]]#6

        regionFlags=[1,2,3,4,5,6]

        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags
                                                     )
        #go ahead and add a boundary tags member 
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="KVApq1.4q12feena%21.16e" % ((he**3)/6.0,)


        logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

        porosityTypes      = numpy.array([1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0,
                                          1.0])

        dragAlphaTypes = numpy.array([0.0,
                                      0.5/1.004e-6,
                                      0.5/1.004e-6,
                                      0.0,
                                      0.5/1.004e-6,
                                      0.0,
                                      0.5/1.004e-6,
                                      0.0,
                                      0.5/1.004e-6,
                                      0.0])

        dragBetaTypes = numpy.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        
        epsFact_solidTypes = np.array([0.0,epsFact_solid,epsFact_solid_2,0.0,epsFact_solid_2,0.0,epsFact_solid_2,0.0,epsFact_solid,0.0])

# Time stepping
T=40.0*period
dt_fixed = T
dt_init = min(0.1*dt_fixed,0.1*he)
runCFL=0.90
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False#True
backgroundDiffusionFactor=0.0
if useMetrics:
    ns_shockCapturingFactor  = 0.25
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.35
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.35
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.75
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 1.5
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
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
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-10,0.00001*he**2)
vof_nl_atol_res = max(1.0e-10,0.00001*he**2)
ls_nl_atol_res = max(1.0e-10,0.0001*he**2)
rd_nl_atol_res = max(1.0e-10,0.005*he)
mcorr_nl_atol_res = max(1.0e-10,0.0001*he**2)
kappa_nl_atol_res = max(1.0e-10,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-10,0.001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,0.0,-9.8]

# Initial condition
waterLine_x =  2*L[0]
waterLine_z =  inflowHeightMean
waterLine_y =  2*L[1]

def signedDistance(x):
    phi_z = x[2]-waterLine_z 
    return phi_z

def theta(x,t):
    return k*x[0] - omega*t + math.pi/2.0

def z(x):
    return x[2] - inflowHeightMean

def ramp(t):
  t0=10 #ramptime
  if t<t0:
    return 1
  else:
    return 1 

h = inflowHeightMean # - transect[0][1] if lower left hand corner is not at z=0

Y =  [0.01227860,     0.00006355]  #Surface elevation Fourier coefficients for non-dimensionalised solution, calculated from FFT
      
def waveHeight(x,t):
   waterDepth = inflowHeightMean 
   for i in range(0,int(len(Y))):  
       waterDepth += Y[i]*cos((i+1)*theta(x,t))/k
   return waterDepth*ramp(t)
 
B = [0.01089188,     0.00014504,     0.00000206,     0.00000003]  #Velocities Fourier coefficients for non-dimensionalised solution, calculated from FFT

def waveVelocity_u(x,t):
   wu = wavelength/period-meanFrameVelocity
   for i in range(0,int(len(B))): 
     wu += sqrt(abs(g[2])/k)*(i+1)*B[i]*cosh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*cos((i+1)*theta(x,t))
    
   return wu*ramp(t)

def waveVelocity_v(x,t):
   wv=0
   for i in range(0,int(len(B))): 
     wv += sqrt(abs(g[2])/k)*(i+1)*B[i]*sinh((i+1)*k*(z(x)+h))/cosh((i+1)*k*h)*sin((i+1)*theta(x,t)) 

   return wv*ramp(t)
   

#solution variables

def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def twpflowVelocity_u(x,t):
    waterspeed = waveVelocity_u(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def twpflowVelocity_v(x,t):
    waterspeed = waveVelocity_v(x,t)
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
    return H*windVelocity[1]+(1.0-H)*waterspeed

def twpflowVelocity_w(x,t):
    return 0.0

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - inflowHeightMean)

def outflowPressure(x,t):
  if x[2]>inflowHeightMean:
    return (L[2]-x[2])*rho_1*abs(g[2])
  else:
    return (L[2]-inflowHeightMean)*rho_1*abs(g[2])+(inflowHeightMean-x[2])*rho_0*abs(g[2])

def u_current(x,t):
    waterspeed = inflowVelocityMean[0]
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-inflowHeightMean-epsFact_consrv_heaviside*he)
    u = H*windVelocity[0] + (1.0-H)*waterspeed
    return u

def v_current(x,t):
    waterspeed = inflowVelocityMean[1]
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-inflowHeightMean-epsFact_consrv_heaviside*he)
    v = H*windVelocity[0] + (1.0-H)*waterspeed
    return v

def zeroVel(x,t):
    return 0.0


def phi_solid(x,t):
  if x[0]<=xSponge_1: 
      if x[1]<=ySponge_1:
          return yRelaxCenter_1-x[1] #region 1 
      elif x[1]>=ySponge_2:   
          return yRelaxCenter_2-x[1] #region 5
  elif x[0]>=xSponge_1:
  #    if x[1]<=ySponge_2:
          return xRelaxCenter_1-x[0] #region 2 and 3 
  #    elif x[1]>ySponge_2:   
  #        return epsFact_solid_2 -sqrt((x[0]-xSponge_1)**2 + (x[1]-ySponge_2)**2) #region 4


from collections import  namedtuple

RelaxationZone = namedtuple("RelaxationZone","phi_solid sign u v w")

class RelaxationZoneWaveGenerator(AV_base):
    """ Prescribe a velocity penalty scaling in a material zone via a Darcy-Forchheimer penalty
    
    :param zones: A dictionary mapping integer material types to Zones, where a Zone is a named tuple
    specifying the x coordinate of the zone center and the velocity components
    """
    def __init__(self,zones):
        assert isinstance(zones,dict)
        self.zones = zones
    def calculate(self):
        for l,m in enumerate(self.model.levelModelList):
            for eN in range(m.coefficients.q_phi.shape[0]):
                mType = m.mesh.elementMaterialTypes[eN]
                if self.zones.has_key(mType):
                    for k in range(m.coefficients.q_phi.shape[1]):
                        t = m.timeIntegration.t
                        x = m.q['x'][eN,k]
                        m.coefficients.q_phi_solid[eN,k] = self.zones[mType].sign*self.zones[mType].phi_solid(x,t)
                        m.coefficients.q_velocity_solid[eN,k,0] = self.zones[mType].u(x,t)
                        m.coefficients.q_velocity_solid[eN,k,1] = self.zones[mType].v(x,t)
                        m.coefficients.q_velocity_solid[eN,k,2] = self.zones[mType].w(x,t)
        m.q['phi_solid'] = m.coefficients.q_phi_solid
        m.q['velocity_solid'] = m.coefficients.q_velocity_solid

rzWaveGenerator = RelaxationZoneWaveGenerator(zones={1:RelaxationZone(phi_solid,
                                                                      1.0,#GENERATION
                                                                      u_current,
                                                                      v_current,  
                                                                      zeroVel),
                                                    2:RelaxationZone(phi_solid,
                                                                      -1.0,#ABSORPTION
                                                                      u_current,
                                                                      v_current,  
                                                                      zeroVel),
                                                    3:RelaxationZone(phi_solid,
                                                                      -1.0,#ABSORPTION
                                                                      u_current,
                                                                      v_current,  
                                                                      zeroVel),
                                                    4:RelaxationZone(phi_solid,
                                                                      -1.0,#ABSORPTION
                                                                      u_current,
                                                                      v_current,  
                                                                      zeroVel),
                                                    5:RelaxationZone(phi_solid,
                                                                      -1.0,#ABSORPTION
                                                                      u_current,
                                                                      v_current,  
                                                                      zeroVel)})
