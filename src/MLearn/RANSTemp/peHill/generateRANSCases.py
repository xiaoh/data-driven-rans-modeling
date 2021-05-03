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

# User specification
paramDict = readInputData('paraInput.in') 
mesh = paramDict['mesh']
turbulence_model = paramDict['turbulence_model']
templateCaseName = paramDict['templateCaseName']
openFoamSolver = paramDict['openFoamSolver']
FinalTimeStep = paramDict['FinalTimeStep']
ReVec = extractListFromDict(paramDict, 'ReVec')
########################################################################################################################
print('using mesh as, ', mesh)
print('Please check the template case is set as, ', turbulence_model)
################################################   Functions   #########################################################
def callFoamSample(caseName, solver):
    """
    call openFOAM solver, sample 
    """
    subprocess.check_call(solver+ ' -case ' + caseName + ' > ' + caseName + '/log.'+solver, shell=True)
    subprocess.check_call('sample -case ' + caseName + ' > ' + caseName + '/log.sample', shell=True)
##############################################   Main Script   #########################################################
job_server = pp.Server()
job_server.set_ncpus()

caseGroupDir = os.path.join(os.getcwd(), mesh, 'cases_'+turbulence_model)  
shutil.rmtree(caseGroupDir, ignore_errors=True)
subprocess.check_call('mkdir '+caseGroupDir, shell=True)
templateCaseDir = os.path.join(os.getcwd(), mesh, templateCaseName)

jobs = []
tic = time.time()
for Re in ReVec:
    ReCase = os.path.join(caseGroupDir, 'Re'+Re)
    shutil.copytree(templateCaseDir, ReCase)
    # replace transport properties
    os.remove(ReCase+'/constant/transportProperties')
    shutil.copy(ReCase+'/constant/transportProperties_pool/transportProperties-Re'+Re, 
                ReCase+'/constant/transportProperties')
    jobs.append(
    job_server.submit(callFoamSample, (ReCase, openFoamSolver), modules = ("subprocess",))
    )
# Barrier
print "Running OpenFOAM ......"
part_sum1 = [job() for job in jobs]
jobs = []
job_server.print_stats()
toc = time.time()
print "Time used for caculating CFD claculation = ", toc-tic
# generate 10 features.
print "CFD Calculation is Done. Start to generate markers"
caseNameList = []
for Re in ReVec:
    ReCase = os.path.join(caseGroupDir, 'Re'+Re)
    caseNameList.append('Re'+Re)
    shutil.copy(ReCase+'/system/controlDict.marker', ReCase+'/system/controlDict')    
    subprocess.check_call('markerFoam > log.marker', shell=True, cwd=ReCase)
    for samplefile in os.listdir(ReCase+'/system/sampleDict_pool/'):
        shutil.copy(ReCase+'/system/sampleDict_pool/'+samplefile, ReCase+'/system/sampleDict')         
        subprocess.check_call('sample >> log.sample', shell=True, cwd=ReCase)

          

