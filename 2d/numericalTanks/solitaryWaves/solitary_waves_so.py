"""
Split operator module for two-phase flow
"""

import os
from proteus.default_so import *
from proteus import Context


# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise NameError, 'Split operator module must end with "_so.py"'

try:
    case = __import__(name)
    Context.setFromModule(case)
    ct = Context.get()
except ImportError:
    raise ImportError, str(name) + '.py not found'

from proteus.BoundaryConditions import BC_Base
BC_Base.getContext() #[temp] resituate later

# List of p/n files
pnList = []

# moving mesh
if ct.movingDomain:
    pnList += [("moveMesh_p", "moveMesh_n")]

# Navier-Stokes and VOF
pnList += [("twp_navier_stokes_p", "twp_navier_stokes_n"),
           ("vof_p", "vof_n")]

# Level set
if not ct.useOnlyVF:
    pnList += [("ls_p", "ls_n"),
               ("redist_p", "redist_n"),
               ("ls_consrv_p", "ls_consrv_n")]

# Turbulence
if ct.useRANS > 0:
    pnList += [("kappa_p", "kappa_n"),
               ("dissipation_p", "dissipation_n")]

#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList=[0.0,ct.dt_init]+[ct.dt_init+ i*ct.dt_out for i in range(1,ct.nDTout+1)]
