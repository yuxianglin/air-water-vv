from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context
import numpy as np

ct = Context.get()
domain = ct.domain
nd = domain.nd
genMesh = ct.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = RANS2P.LevelModel
if ct.useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if ct.useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if ct.useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if ct.movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

if hasattr (domain, 'porosityTypes'):
    porosityTypes = domain.porosityTypes
    dragAlphaTypes = domain.dragAlphaTypes
    dragBetaTypes = domain.dragBetaTypes
    epsFact_solid = domain.epsFact_solid
else:
    porosityTypes = None
    dragAlphaTypes = None
    dragBetaTypes = None
    epsFact_solid = None

coefficients = RANS2P.Coefficients(epsFact = ct.epsFact_viscosity,
                                   sigma = 0.0,
                                   rho_0 = ct.rho_0,
                                   nu_0 = ct.nu_0,
                                   rho_1 = ct.rho_1,
                                   nu_1 = ct.nu_1,
                                   g = ct.g,
                                   nd = nd,
                                   ME_model = int(ct.movingDomain)+0,
                                   VF_model = int(ct.movingDomain)+1,
                                   LS_model = int(ct.movingDomain)+LS_model,
                                   Closure_0_model = Closure_0_model,
                                   Closure_1_model = Closure_1_model,
                                   epsFact_density = ct.epsFact_density,
                                   stokes = False,
                                   useVF = ct.useVF,
                                   useRBLES = ct.useRBLES,
                                   useMetrics = ct.useMetrics,
                                   eb_adjoint_sigma = 1.0,
                                   eb_penalty_constant = ct.weak_bc_penalty_constant,
                                   forceStrongDirichlet = ct.ns_forceStrongDirichlet,
                                   turbulenceClosureModel = ct.ns_closure,
                                   movingDomain = ct.movingDomain,
                                   porosityTypes = porosityTypes,
                                   dragAlphaTypes = dragAlphaTypes,
                                   dragBetaTypes = dragBetaTypes,
                                   epsFact_solid = epsFact_solid,
                                   barycenters = ct.domain.barycenters)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet,
                       1: lambda x, flag: domain.bc[flag].u_dirichlet,
                       2: lambda x, flag: domain.bc[flag].v_dirichlet}


advectiveFluxBoundaryConditions =  {0: lambda x, flag: domain.bc[flag].p_advective,
                                    1: lambda x, flag: domain.bc[flag].u_advective,
                                    2: lambda x, flag: domain.bc[flag].v_advective}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: domain.bc[flag].u_diffusive},
                                   2: {2: lambda x, flag: domain.bc[flag].v_diffusive}}


class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if ct.signedDistance(x) < 0:
            return -(ct.L[1] - self.waterLevel)*ct.rho_1*ct.g[1] - (self.waterLevel - x[1])*ct.rho_0*ct.g[1]
        else:
            return -(ct.L[1] - self.waterLevel)*ct.rho_1*ct.g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class initialVelocity_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
      if ct.signedDistance(x) < 0:
        return ct.netcurrentVelocity
      else: 
        return 0.0

initialConditions = {0:PerturbedSurface_p(ct.waterLine_z),
                     1:initialVelocity_u(),
                     2:AtRest()}


