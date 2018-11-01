
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
# add extra attributes here
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()

from InvadopodiaSteppables import InitializeMMPAndECM
initializeMMPAndECM=InitializeMMPAndECM(sim,_frequency=30000000)
steppableRegistry.registerSteppable(initializeMMPAndECM)

from InvadopodiaSteppables import InvadopodiaSteppable
steppableInstance=InvadopodiaSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)

from InvadopodiaSteppables import SolubleMMPSecretionSteppable
solubleMMPSecretionSteppable=SolubleMMPSecretionSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(solubleMMPSecretionSteppable)

from InvadopodiaSteppables import ECMDegradation
ecmDegradation=ECMDegradation(sim,_frequency=1)
steppableRegistry.registerSteppable(ecmDegradation)

from InvadopodiaSteppables import LogDataSteppable
logDataSteppable=LogDataSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(logDataSteppable)

from InvadopodiaSteppables import ExtraFieldVisualizationSteppable
extraFieldVisualizationSteppable=ExtraFieldVisualizationSteppable(sim,_frequency=30000000)
steppableRegistry.registerSteppable(extraFieldVisualizationSteppable)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        