from math import *
import numpy as np

from proteus import Domain
from proteus.mprans import SpatialTools as st

from proteus import MeshTools, AuxiliaryVariables

from proteus import Gauges as ga
from proteus import WaveTools as wt

from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral



# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()



# ----- WAVE CONDITIONS ----- #
period=1.94
waterLevel = 1.
waveDir=np.array([1, 0., 0.])
mwl=waterLevel #coordinate of the initial mean level of water surface
waveHeight=0.025
omega = float(2.0*pi/period)
wavelength = 4.997 #from dispersion law
inflowHeightMean=waterLevel
inflowVelocityMean =np.array([0.,0.,0.])
windVelocity = np.array([0.,0.,0.])
k = float(2.0*pi/wavelength)

#---------Domain Dimension
nd = 2


# ----- SHAPES ----- #

tank_dim = [20.0, 1.5]

L_leftSpo  = 5.# wavelength


boundaryOrientations = {'bottom': [0., -1.,0.],
                        'right': [1., 0.,0.],
                        'top': [0., 1.,0.],
                        'left': [-1., 0.,0.],
                        'sponge': None,
                       }
boundaryTags = {'bottom': 1,
                'right': 2,
                'top': 3,
                'left': 4,
                'sponge': 5,
               }


vertices=[[0.0,0.0],#0
          [L_leftSpo,0.0],#1
          [tank_dim[0],0.0],#2
          [tank_dim[0],tank_dim[1]],#3
          [L_leftSpo,tank_dim[1]],#4
          [0.0,tank_dim[1]],#5
          ]

vertexFlags=np.array([1, 1, 1,
                      3, 3, 3,
                     ])


segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,0],
          [1,4],
         ]

segmentFlags=np.array([1, 1, 2,
                       3, 3, 4,
                       5,
                      ])


regions = [ [ 0.1*tank_dim[0] , 0.1*tank_dim[1] ],
            [ 0.95*tank_dim[0] , 0.95*tank_dim[1] ] ]

regionFlags=np.array([1, 2])


tank = st.CustomShape(domain, vertices=vertices, vertexFlags=vertexFlags,
                      segments=segments, segmentFlags=segmentFlags,
                      regions=regions, regionFlags=regionFlags,
                      boundaryTags=boundaryTags, boundaryOrientations=boundaryOrientations)


##########################################
# Numerical Options and other parameters #
##########################################

nd = 2

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5

# Surface tension
sigma_01 = 0.0

# Gravity
g =np.array([0.,-9.8,0.])
gAbs=sqrt(sum(g**2))


#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------

weak_bc_penalty_constant = 10.0/nu_0 #100
dt_fixed = period/21.0
dt_init = min(0.1*dt_fixed,0.001)
T = 50*period
nDTout= int(round(T/dt_fixed))
runCFL = 0.9




# ----- WAVE input ----- #

"""
waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=wavelength,
                                  waveType='Fenton',
                                  Ycoeff = Ycoeff,
                                  Bcoeff = Bcoeff,
                                  )
"""
waveinput = wt.MonochromaticWaves(period=period,
                                  waveHeight=waveHeight,
                                  mwl=mwl,
                                  depth=waterLevel,
                                  g=g,
                                  waveDir=waveDir,
                                  wavelength=wavelength,
                                  waveType='Linear',
                                 )



# ----- BOUNDARY CONDITIONS ----- #

tank.BC.top.setOpenAir()
tank.BC.left.setUnsteadyTwoPhaseVelocityInlet(wave=waveinput, vert_axis=1, windSpeed=windVelocity)
#tank.BC.left.setFreeSlip()
tank.BC.bottom.setFreeSlip()

tank.BC.right.setFreeSlip()

#tank.BC.sponge.setParallelFlag0()
tank.BC.sponge.setNonMaterial()


# -----  GENERATION ZONE & ABSORPTION ZONE  ----- #


tank.setGenerationZones(flags=1, epsFact_solid=L_leftSpo/2, orientation=[1., 0., 0.], center=(L_leftSpo/2, 0., 0.),
                        waves=waveinput, windSpeed=windVelocity,
                        )

#
# ----- Output Gauges ----- #

gauge_dx=1.0 #0.25
probes=np.linspace(0., tank_dim[0], (tank_dim[0]/gauge_dx)+1)#
#LG=[]
#PG=[]
PG0=[]
for i in probes:
    #LG.append(((i,0.,0.), (i,tank_dim[1],0.)),)
    #PG.append((i,0.5,0.),)
    PG0.append((i, 0., 0.),)


point_output=ga.PointGauges(gauges=((('p'),PG0),
                                 ),
                          activeTime = (0., T),
                          sampleRate=0.,
                          fileName='point_gauges.csv')


domain.auxiliaryVariables += [point_output,
                              ]


# ----- Mesh ----- #

he=wavelength/100
domain.MeshOptions.elementSize(he)
st.assembleDomain(domain)

domain.writePoly("mesh")
triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

quad_order = 3


#  Discretization -- input options

checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True #False
useOnlyVF = False # if TRUE  proteus uses only these modules --> twp_navier_stokes_p + twp_navier_stokes_n
                  #                                              vof_p + vof_n
movingDomain=False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988

genMesh=True

# By DEFAULT on the other files.py -->  fullNewtonFlag = True 
#                                       multilevelNonlinearSolver & levelNonlinearSolver == NonlinearSolvers.Newton

useOldPETSc=False # if TRUE  --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.PETSc
                  # if FALSE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.KSP_petsc4py

useSuperlu = False #if TRUE --> multilevelLinearSolver & levelLinearSolver == LinearSolvers.LU

spaceOrder = 1
useHex     = False # used for discretization, if 1.0 --> CubeGaussQuadrature
                   #                          ELSE   --> SimplexGaussQuadrature

useRBLES   = 0.0 # multiplied with subGridError
useMetrics = 1.0 # if 1.0 --> use of user's parameters as (ns_shockCapturingFactor, ns_lag_shockCapturing, ecc ...)
useVF = 0.0 # used in the smoothing functions as (1.0-useVF)*smoothedHeaviside(eps_rho,phi) + useVF*fmin(1.0,fmax(0.0,vf))

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
ns_forceStrongDirichlet = False #True
backgroundDiffusionFactor=0.01


""" Default configuration by Tristan
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
    redist_Newton = True #False
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
"""


if useMetrics:
    ns_shockCapturingFactor  = 0.25 # magnifies numerical viscosity in NS (smoothening velocity fields)
    ns_lag_shockCapturing = True # lagging numerical viscosity speedsup Newton but destabilzes the solution
    ns_lag_subgridError = True # less nonlinear but less stable
    ls_shockCapturingFactor  = 0.25 # numerical diffusion of level set (smoothening phi)
    ls_lag_shockCapturing = True # less nonlinear but less stable
    ls_sc_uref  = 1.0 # reference gradient in numerical solution (higher=more diffusion)
    ls_sc_beta  = 1.0 # 1 is fully nonlinear, 2 is linear
    vof_shockCapturingFactor = 0.25 # numerical diffusion of level set (smoothening volume of fraction)
    vof_lag_shockCapturing = True # less nonlinear but less stable
    vof_sc_uref = 1.0 
    vof_sc_beta = 1.0 
    rd_shockCapturingFactor  = 0.25
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0 # control width of water/air transition zone
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0 # affects smoothing diffusion in mass conservation
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True # False
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
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
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
rd_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he)
mcorr_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
kappa_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*domain.MeshOptions.he**2)



#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4

# Initial condition
waterLine_x = 2*tank_dim[0]
waterLine_z = inflowHeightMean



def waveHeight(x,t):
   # waterDepth=waveinput.waterDepth(x[0], x[1], x[2], t)
    waterDepth=waveinput.eta(x, t)+waveinput.mwl
    return waterDepth


def wavePhi(x,t):
    if nd==2:
        return x[1]- waveHeight(x,t)
    else:
        return x[2]- waveHeight(x,t)



def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))


def signedDistance(x):
    phi_x = x[0]-waterLine_x # (tank_dim[0]-L_rightSpo) # x0_por
    phi_z = x[nd-1]-waterLine_z

    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*domain.MeshOptions.he,phi)))


