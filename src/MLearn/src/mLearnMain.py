#!/usr/bin/env python

# description        :Main Driver for Machine Learning for Turbulence Modelilng
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.13, 2016
# revision           :Apr.28     , 2017
########################################################################################################################
# system module
import sys
import os
import time
import numpy as np
import pdb
import ast
import subprocess
import shutil
# local module
import genFeature as gF
import genDeltaTauTruth as gD
import regressionSolver as rg
import outputFormat as outF
import reconstructTau as recon
import plotLearningPerformance as postP
from utilities import readInputData, replace, extractListFromDict

########################################################################################################################
############################################ parse control parameters ##################################################
########################################################################################################################
paramDict = readInputData('mainInput.in') 
# test flow info
testFlow = paramDict['testFlow']
testFlowRe = paramDict['testFlowRe']
testFlownMesh = paramDict['testFlownMesh']
# training flows info
trainFlows = extractListFromDict(paramDict, 'trainFlows')
trainFlowsRes = extractListFromDict(paramDict, 'trainFlowsRes')
trainFlowsnMeshs = extractListFromDict(paramDict, 'trainFlowsnMeshs')
# turublence model used (LaunderSharmaKE, ...)
turbulence_model = paramDict['turbulence_model']
# combination keyword for feature
combinationKeyword_feature = paramDict['combinationKeyword_feature']
reGenerateFeatureFile = ast.literal_eval(paramDict['reGenerateFeatureFile'])
reGenerateDataFile = ast.literal_eval(paramDict['reGenerateDataFile'])
try: 
    labelList = extractListFromDict(paramDict, 'labelList')
except:
    labelList = []
# specifying parameters for genering training data
deltaOrienParamerization = paramDict['deltaOrienParamerization']
dataSparsenessList = extractListFromDict(paramDict, 'dataSparsenessList')
dataSourceList = extractListFromDict(paramDict, 'dataSourceList')
if (len(dataSourceList) != len(dataSparsenessList)) or (len(dataSourceList)-1 != len(trainFlows)) or \
   (len(dataSourceList)-1 != len(trainFlowsRes)) or (len(dataSourceList)-1 != len(trainFlowsnMeshs)):
   print "Error: check trainFlows, trainFlowsRes, trainFlowsnMeshs, dataSparsenessList, dataSourceList, \
          the should have same len as number of training flows"
   sys.exit(1)
# specifying regression model
regressionModel = paramDict['regressionModel']
responseCombinationFlag = paramDict['responseCombinationFlag']
# source file dir
srcDir = os.path.normpath(os.path.join(os.path.abspath(__file__), '../..'))
RANSDir = os.path.join(srcDir, 'RANSTemp')

# set random seeds
np.random.seed(2000)

outF.buildStarBlock("Case Parameters")
print "test flow: " 
print "   flow type: ", testFlow
print "   flow Re: ", testFlowRe
print "   flow nMesh", testFlownMesh
print ''
print "training flows:"
print "   flow types: ", trainFlows
print "   flow Res: ", trainFlowsRes
print "   flow nMeshes", trainFlowsnMeshs
print ''
print "Turbulence model: ", turbulence_model
print ''
print "Feature space: ", combinationKeyword_feature
print ''
print "orientation parameterization: ", deltaOrienParamerization


########################################################################################################################
############################################ Construct Input Space #####################################################
########################################################################################################################
outF.buildStarBlock("Start constructing features")
# combination keywords: 
# Full -- 10 scalar features + 47 expended features
# FullSelected -- subset of 57 features, a keywords list should also be included
# Scalar -- only 10 scalar features
# ScalarSelected -- subset of 10 scalar features, a keywords list should also be included
# Expanded -- only 47 expanded features
# ExpandedSelected -- subset of 47 expanded features, a keywords list should also be included

featureFolder = './featureSet_'+combinationKeyword_feature+'_userSelectedN'+str(len(labelList))
flowTypes = trainFlows + [testFlow]
flowRes = trainFlowsRes + [testFlowRe]
flownMeshs = trainFlowsnMeshs + [testFlownMesh]

if reGenerateFeatureFile or not os.path.exists(featureFolder):
    print "We will regenerate the feature files"    
    if os.path.exists(featureFolder):
        print "remove the original folder"
        shutil.rmtree(featureFolder)
    subprocess.check_call('mkdir '+featureFolder, shell=True)

    for flowType, flowRe, flownMesh in zip(flowTypes, flowRes, flownMeshs):
        print "We are generating flow feature file for: ", flowType, flowRe, 'nMesh-', flownMesh
        if combinationKeyword_feature in ['FullSelected', 'ScalarSelected', 'ExpandedSelected']:
            conFeature = gF.constructFeatureMain(flowType, flowRe, flownMesh, turbulence_model, 
                                                 combinationKeyword_feature, labelList)
        elif combinationKeyword_feature in ['Full', 'Scalar', 'Expanded']:
            conFeature = gF.constructFeatureMain(flowType, flowRe, flownMesh, turbulence_model, combinationKeyword_feature)
        else:
            print "Please provide valid combination keyword for features:\n \
                   FullSelected, ScalarSelected, ExpandedSelected\
                   Full, Scalar, Expanded"
        conFeature.genFeature()
else:
    print "I will read the feature files generated before."

########################################################################################################################
################################### Copy train/test RANS Cases To Current Folder #######################################
########################################################################################################################
currentRANSDir = 'RANSCases'
shutil.rmtree(currentRANSDir, ignore_errors=True)
for flowType, flowRe, flownMesh in zip(flowTypes, flowRes, flownMeshs):
    ransCaseDir = os.path.join(RANSDir, flowType, 'cell'+flownMesh, 'cases_'+turbulence_model, 'Re'+flowRe) 
    currentRansCaseDir = os.path.join(currentRANSDir, flowType+'Re'+flowRe)
    if not os.path.exists(currentRansCaseDir):
        shutil.copytree(ransCaseDir, currentRansCaseDir)
        shutil.move(os.path.join(currentRansCaseDir, '0',), os.path.join(currentRansCaseDir, '0.org'))
    
########################################################################################################################
############################################ Obtain Training Data ######################################################
########################################################################################################################
outF.buildStarBlock("Start generate training Tau discrepancy data")
currentTrainDataDir = './trainData'
for flowType, flowRe, flownMesh, dataSparseness, dataSource in zip(flowTypes, flowRes, flownMeshs, dataSparsenessList, dataSourceList):    
    # copy RANS Tau
    tauCurrentDir = os.path.join(currentRANSDir, flowType+'Re'+flowRe, '0.000000', 'Tau')
    tauRANSDir = os.path.join(currentRANSDir, flowType+'Re'+flowRe, '0.000000', 'TauRANS')
    shutil.copy(tauCurrentDir, tauRANSDir)
    replace(tauRANSDir, 'Tau', 'TauRANS')
    
    currCase = os.path.join(currentRANSDir, flowType+'Re'+flowRe)
    shutil.copy(os.path.join(currCase,'system/sampleDict.template'),os.path.join(currCase,'system/sampleDict'))
    replace(os.path.join(currCase,'system/sampleDict'), '<field>', 'TauRANS')
    subprocess.check_call('sample >> log.sampleTruth', shell=True, cwd=currCase)
    
    # copy DNS Tau
    if dataSparseness == 'fullField':
        dataCaseDir = os.path.join(srcDir, 'database', flowType, dataSparseness)
        tauDataDir = os.path.join(dataCaseDir, 'Tau_'+flowType+'-Re'+flowRe+'-'+dataSource+'-cell'+flownMesh)
        shutil.copy(tauDataDir, os.path.join(currentRANSDir, flowType+'Re'+flowRe, '0.000000', 'TauDNS'))     
        # sample
        shutil.copy(os.path.join(currCase,'system/sampleDict.template'),os.path.join(currCase,'system/sampleDict'))
        replace(os.path.join(currCase,'system/sampleDict'), '<field>', 'TauDNS')
        subprocess.check_call('sample >> log.sampleTruth', shell=True, cwd=currCase)

    elif dataSparseness == 'coarseData':
        #TODO : To be implemented
        pass
    if reGenerateDataFile or not os.path.exists(currentTrainDataDir):       
        print "Start generating training data, with"
        genData = gD.getTrainingDataMain(flowType, flowRe, flownMesh, dataSparseness, dataSource, 
                                     currentRANSDir, turbulence_model, deltaOrienParamerization)    
        deltaTauData_headerNames = genData.getDiscrepancy()
    else:
        print "I will read the data files generated before."
        if deltaOrienParamerization == 'Euler':
            deltaTauData_headerNames = gD.deltaTau_headerNames_Euler
        elif deltaOrienParamerization == 'Quaternion':
            deltaTauData_headerNames = gD.deltaTau_headerNames_Quaternion
# get data for test flow, only for comparison

########################################################################################################################
################################################## ML Learning #########################################################
########################################################################################################################
outF.buildStarBlock('Start Machine Learning Process')
rgModel = eval('rg.' + regressionModel + 'RegressionSolver(trainFlows, trainFlowsRes, testFlow, testFlowRe, \
           deltaOrienParamerization, featureFolder, currentTrainDataDir, deltaTauData_headerNames)')
predictResponseList = rgModel.learnPrediction(responseCombinationFlag)

########################################################################################################################
############################################## Reconstruction ##########################################################
########################################################################################################################
outF.buildStarBlock('Reconstructing Predicted Reynolds stress')
#recTau = recon.reconstructPredictedTau(currentRANSDir, testFlow, testFlowRe, int(testFlownMesh), deltaOrienParamerization, predictResponseList)
recTau = recon.reconstructPredictedTau(currentRANSDir, testFlow, testFlowRe, int(testFlownMesh), deltaOrienParamerization)
recTau.reconstructCorrectedTau()
########################################################################################################################
################################################### Plotting ###########################################################
########################################################################################################################
if not os.path.exists('figures'):
    subprocess.check_call('mkdir figures', shell=True)

outF.buildStarBlock('Post Plotting Results')
postPlt = postP.plotLearningPerformance(currentRANSDir, trainFlows, trainFlowsRes, testFlow, testFlowRe, deltaOrienParamerization)
postPlt.plotPhyDeltaTau_compareAll()
postPlt.plotTau_compareAll()



