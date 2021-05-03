#!/usr/bin/env python

# description        :Main Driver for Machine Learning for Turbulence Modelilng
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.13, 2016
# revision           :Oct.13, 2016
########################################################################################################################

#system import
import shutil
import subprocess
import os
import pdb
import time
import pp
import copy
import numpy as np
#local import
import peHillPlotGivenCases as peplot
from utilities import readInputData, replace, extractListFromDict

# specification read from dict
paramDict = readInputData('paraInput.in') 
mesh = paramDict['mesh']
turbulence_model = paramDict['turbulence_model']
templateCaseName = paramDict['templateCaseName']
openFoamSolver = paramDict['openFoamSolver']
FinalTimeStep = paramDict['FinalTimeStep']
ReVec = extractListFromDict(paramDict, 'ReVec')
featureNameList = ['QCriterion', 'TurbulenceIntensity','ReTurbulence', 'pGradAlongStream', 'TurbulenceTime', 
                   'PressureStress', 'UGradMisalignment','ConvectionTKE', 'TauRatio', 'Curvature']

########################################################################################################################
print('using mesh as, ', mesh)
print('Please check the template case is set as, ', turbulence_model)

caseGroupDir = os.path.join(os.getcwd(), mesh, 'cases_'+turbulence_model)  

caseNameList = []
timeDirList = []
for Re in ReVec:
    caseNameList.append('Re'+Re)
    timeDirList.append(FinalTimeStep)

os.chdir(caseGroupDir)
np.random.seed(800)
colors = peplot.genColorVec(caseNameList) 
peplot.plotScalar('q14', caseNameList, timeDirList, 1e5, colors) 

