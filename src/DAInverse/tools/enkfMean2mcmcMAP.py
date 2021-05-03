#!/aoe/bin/python27
import numpy as np
import os
import pdb

from utilities import readInputData, extractListFromDict
paramDict = readInputData('mcmcPost.in')
finalStep = float(paramDict['finalStep'])
    
if __name__ == "__main__":
   resultsFolder = './forMCMC'
   os.system('mkdir ' + resultsFolder)
   omegaEnsemble = np.loadtxt('./debugData/XC_' + str(finalStep))
   nDim, ns = omegaEnsemble.shape
   omegaInit = np.mean(omegaEnsemble, axis=1)
   assert len(omegaInit) == nDim, "Error, mean of omegaInit is not correct, try mean of omegaInit.T "
   np.savetxt(resultsFolder+'/omegaInit.data', omegaInit)

   
