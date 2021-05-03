#! /usr/bin/env python
# Jianxun Wang (vtwjx@vt.edu)
# Mar. 31, 2017

# SYSTEM MODULE
import numpy as np
import os
import pdb
import sys
import shutil
# LOCAL MODULE
import foamFileOperation as foamOp
from utilities import replace, readInputData

paramDict = readInputData('forwardModelInput.in')
caseName = paramDict['caseName']
nCell = int(paramDict['Ncell'])     
nModes = int(paramDict['Nmode'])

varVec = ['Xi', 'Eta', 'K', 'VA', 'VB', 'VC']
for var in varVec:
    print 'read foam file to KLModes: ', var
    klModes = np.zeros([nCell, nModes])
    for i in range(nModes):
        modeFilePath = os.path.join(caseName+'_KL', '0', var+'Mode'+str(i))
        mode = foamOp.readScalarFromFile(modeFilePath)
        klModes[:, i] = mode
    if not os.path.exists('randomData_'+var):
        os.system('mkdir randomData_'+var)
    klPath = os.path.join('randomData_'+var, 'KLModes.dat')
    np.savetxt(klPath, klModes)
    
        
        
