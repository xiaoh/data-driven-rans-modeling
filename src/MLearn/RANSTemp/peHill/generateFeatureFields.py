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
import genFeature as gF
from utilities import readInputData, replace, extractListFromDict

# User specification
paramDict = readInputData('paraInput.in') 
mesh = paramDict['mesh']
turbulence_model = paramDict['turbulence_model']
templateCaseName = paramDict['templateCaseName']
openFoamSolver = paramDict['openFoamSolver']
FinalTimeStep = paramDict['FinalTimeStep']
ReVec = extractListFromDict(paramDict, 'ReVec')
flowType = 'peHill'
dataBaseDir = '../../database/' + flowType + '/fullField/'
########################################################################################################################
print('using mesh as, ', mesh)
print('Please check the template case is set as, ', turbulence_model)
caseGroupDir = os.path.join(os.getcwd(), mesh, 'cases_'+turbulence_model)  
print "Start to generate markers"


sampleNameListOrg = ['U', 'Tau', 'QCriterion', 'TurbulenceIntensity','ReTurbulence', 'pGradAlongStream', 'TurbulenceTime', 
                  'PressureStress', 'UGradMisalignment','ConvectionTKE', 'TauRatio', 'Curvature']



for Re in ReVec:
    ReCase = os.path.join(caseGroupDir, 'Re'+Re)
    ReCaseData = os. path.join(ReCase, FinalTimeStep+'.000000')
    shutil.copy(ReCase+'/system/controlDict.marker', ReCase+'/system/controlDict')
    subprocess.check_call('foamLog log.'+openFoamSolver, shell=True, cwd=ReCase)
    subprocess.check_call('writeCellCentres', shell=True, cwd=ReCase)
    subprocess.check_call('markerFoam > log.marker', shell=True, cwd=ReCase)
    # create a template feature file
    shutil.copy(ReCaseData+'/k', ReCaseData+'/scalarFieldTemplate')
    replace(ReCaseData+'/scalarFieldTemplate', 'k;', '<Feature>;')
    # construct all expanded features
    user_selection = {'strainRate', 'rotationRate', 'gradP', 'gradK'}
    featureMatrixSelected, featureLabels = gF.expandFeatures(ReCaseData, user_selection)

    shutil.rmtree(os.path.join(ReCase,'postProcessing'), ignore_errors=True)
    if os.path.exists(os.path.join(ReCase,'log.sample')):
        os.remove(os.path.join(ReCase,'log.sample'))
    for sampleName in sampleNameListOrg:
        shutil.copy(os.path.join(ReCase,'system/sampleDict.template'), os.path.join(ReCase,'system/sampleDict'))
        replace(os.path.join(ReCase,'system/sampleDict'), '<field>', sampleName)
        subprocess.check_call('sample >> log.sample', shell=True, cwd=ReCase)

    for i in np.arange(featureMatrixSelected.shape[1]):
        shutil.copy(os.path.join(ReCase,'system/sampleDict.template'), os.path.join(ReCase,'system/sampleDict'))
        replace(os.path.join(ReCase,'system/sampleDict'), '<field>', 'q'+str(i))
        subprocess.check_call('sample >> log.sample', shell=True, cwd=ReCase)
    newZeroFolder = os.path.join(ReCase, '0.000000')
    if os.path.exists(newZeroFolder):
        shutil.rmtree(newZeroFolder)
    shutil.copytree(ReCaseData, newZeroFolder)
    
  
