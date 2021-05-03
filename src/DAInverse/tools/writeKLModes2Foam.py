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
    print 'writing modes to foam file for variable: ', var
    varPath = os.path.join('randomData_'+var, 'KLModes.dat')
    klModes = np.loadtxt(varPath)
    if not ((klModes.shape[1] == nModes) and (klModes.shape[0] == nCell)):
        print "make sure number of modes and number of cell is correct in the forwardModelInput.in"
        sys.exit()
    else:
        for i in range(nModes):
            templateFilePath = os.path.join(caseName, '0', 'scalarFieldTemplate')
            modeFilePath = os.path.join(caseName, '0', var+'Mode'+str(i))
            shutil.copy(templateFilePath, modeFilePath)
            replace(modeFilePath, '<scalarField>', var+'Mode'+str(i))
            foamOp.writeScalarToFile(klModes[:, i], modeFilePath)
        
