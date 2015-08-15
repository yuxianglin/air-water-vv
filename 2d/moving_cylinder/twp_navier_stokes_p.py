from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context
ct = Context.get()
genMesh = ct.genMesh
movingDomain = ct.movingDomain
L = ct.L
T = ct.T
nd = ct.nd
domain = ct.domain

LevelModelType = RANS2P.LevelModel
if ct.useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if ct.useRANS >= 1:
    Closure_0_model = 5
    Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2
        Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1
        Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS2P.Coefficients(epsFact=ct.epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = ct.rho_0,
                                   nu_0 = ct.nu_0,
                                   rho_1 = ct.rho_1,
                                   nu_1 = ct.nu_1,
                                   g=ct.g,
                                   nd=ct.nd,
                                   ME_model=int(ct.movingDomain)+0,
                                   VF_model=int(ct.movingDomain)+1,
                                   LS_model=int(ct.movingDomain)+LS_model,
                                   Closure_0_model=Closure_0_model,
                                   Closure_1_model=Closure_1_model,
                                   epsFact_density=ct.epsFact_density,
                                   stokes=False,
                                   useVF=ct.useVF,
                                   useRBLES=ct.useRBLES,
                                   useMetrics=ct.useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=ct.weak_bc_penalty_constant,
                                   forceStrongDirichlet=ct.ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ct.ns_closure,
                                   movingDomain=ct.movingDomain,
                                   barycenters=ct.barycenters)

domain = ct.domain

def getDBC_p(x,flag):
    if flag == ct.boundaryTags['right']:
        return lambda x,t: 0.0

def getDBC_u(x,flag):
    if flag == ct.boundaryTags['left']:
        return lambda x,t: -ct.speed*ct.inflowProfile(x,t)
    if flag in (ct.boundaryTags['front'],
                ct.boundaryTags['back'],
                ct.boundaryTags['top'],
                ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return lambda x,t: -ct.speed
        elif ct.wallBC=="no_slip_obstacle":
            return lambda x,t: 0.0
        elif ct.wallBC=="slip":
            return None
        else:
            raise RuntimeError
    elif flag == ct.boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == ct.boundaryTags['left']:
        return lambda x,t: 0.0
    elif flag in (ct.boundaryTags['front'],
                  ct.boundaryTags['back'],
                  ct.boundaryTags['top'],
                  ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return lambda x,t: 0.0
        elif ct.wallBC=="no_slip_obstacle":
            return lambda x,t: 0.0
        elif ct.wallBC=="slip":
            return None
        else:
            raise RuntimeError
    elif flag == ct.boundaryTags['obstacle']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

def getAFBC_p(x,flag):
    if flag == ct.boundaryTags['left']:
        return lambda x,t: ct.speed*ct.inflowProfile(x,t)
    elif flag in (ct.boundaryTags['top'],
                  ct.boundaryTags['bottom'],
                  ct.boundaryTags['front'],
                  ct.boundaryTags['back']):
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag in (ct.boundaryTags['front'],
                ct.boundaryTags['back'],
                ct.boundaryTags['top'],
                ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return None
        elif ct.wallBC=="no_slip_obstacle":
            return None
        elif ct.wallBC=="slip":
            return lambda x,t: 0.0
        else:
            raise RuntimeError

def getAFBC_v(x,flag):
    if flag in (ct.boundaryTags['front'],
                ct.boundaryTags['back'],
                ct.boundaryTags['top'],
                ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return None
        elif ct.wallBC=="no_slip_obstacle":
            return None
        elif ct.wallBC=="slip":
            return lambda x,t: 0.0
        else:
            raise RuntimeError

def getDFBC_u(x,flag):
    if flag == ct.boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag in (ct.boundaryTags['front'],
                  ct.boundaryTags['back'],
                  ct.boundaryTags['top'],
                  ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return None
        elif ct.wallBC=="no_slip_obstacle":
            return None
        elif ct.wallBC=="slip":
            return lambda x,t: 0.0
        else:
            raise RuntimeError

def getDFBC_v(x,flag):
    if flag == ct.boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag in (ct.boundaryTags['front'],
                  ct.boundaryTags['back'],
                  ct.boundaryTags['top'],
                  ct.boundaryTags['bottom']):
        if ct.wallBC=="no_slip_observer":
            return None
        elif ct.wallBC=="no_slip_obstacle":
            return None
        elif ct.wallBC=="slip":
            return lambda x,t: 0.0
        else:
            raise RuntimeError

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

class P_IC:
    def uOfXT(self,x,t):
        return ct.twpflowPressure_init(x,t)

class U_IC:
    def uOfXT(self,x,t):
        return 0.0

class V_IC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:P_IC(),
                     1:U_IC(),
                     2:V_IC()}
