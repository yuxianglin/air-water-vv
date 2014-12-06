from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from crump_weir import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
if useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2
coefficients = VOF.Coefficients(LS_model=LS_model,V_model=0,RD_model=RD_model,ME_model=1,
                                checkMass=False,useMetrics=useMetrics,
                                epsFact=epsFact_vof,sc_uref=vof_sc_uref,sc_beta=vof_sc_beta)
 
def getDBC_vof(x,flag):
   if flag == boundaryTags['top']:
       return lambda x,t: 1.0
   elif flag == boundaryTags['left'] and x[1] > waterLine_z:
       return lambda x,t: 1.0
   elif flag == boundaryTags['left'] and x[1] <=  waterLine_z:
       return lambda x,t: 0.0
   elif flag == boundaryTags['right'] and x[1] > downstream_water_level:
       return lambda x,t: 1.0
   elif flag == boundaryTags['right'] and x[1] <= downstream_water_level:
       return lambda x,t: 0.0

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag == boundaryTags['left']  and x[1] <= waterLine_z:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']  and x[1] <= downstream_water_level:
        return lambda x,t: 0.0

#    if flag == boundaryTags['top']:# or x[1] >= L[1] - 1.0e-12:
#        return None
#    elif flag == boundaryTags['left']  and x[1] > waterLine_z:
#        return lambda x,t: 1.0
#    elif flag == boundaryTags['right']:
#        return None
#    else:
#        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_H:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))
	    
initialConditions  = {0:PerturbedSurface_H()}