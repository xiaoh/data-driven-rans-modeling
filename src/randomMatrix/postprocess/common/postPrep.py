#!/aoe/bin/python27

# description        :This file used in RMT case folder

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Jan.16, 2016
# revision           :Jan.16, 2016
####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.stats as ss

# system, file operation
import pdb
import sys
import os
import ast

# local import
import postFuns as pF
from utilities import readInputData, extractListFromDict, replace

def sampleAllComponents(caseName, caseRMT, casePhy):
    """
    sample all components in 0 folder in 8 lines. 
    
    """
    sampleCaseRMT = caseRMT + '/' + caseName + '_sample'
    sampleCasePhy = casePhy + '/' + caseName + '_sample'
    
    sampleFile_RMT = sampleCaseRMT + '/system/sampleDict'
    sampleFile_Phy = sampleCasePhy + '/system/sampleDict' 
    
    componentPool = ['c1', 'c2', 'c3', 'Xi', 'Eta', 'k', 'VA', 'VB', 'VC', 'deltaXi', 
                     'deltaEta', 'deltaK', 'deltaVA', 'deltaVB', 'deltaVC']
    componentEntPool = ['c1_SField', 'c2_SField', 'c3_SField', 'Xi_SField', 'Eta_SField', 
                        'TKE_SField', 'VA_SField', 'VB_SField', 'VC_SField',
                        'deltaXi_SField', 'deltaEta_SField', 'deltaK_SField', 
                        'deltaVA_SField', 'deltaVB_SField', 'deltaVC_SField']
                        
    componentKLDPool = ['c1_KLDistField', 'c2_KLDistField', 'c3_KLDistField', 'Xi_KLDistField', 
                        'Eta_KLDistField', 'TKE_KLDistField', 'VA_KLDistField', 'VB_KLDistField', 'VC_KLDistField',
                        'deltaXi_KLDistField', 'deltaEta_KLDistField', 'deltaK_KLDistField', 
                        'deltaVA_KLDistField', 'deltaVB_KLDistField', 'deltaVC_KLDistField']    
    # clear postProcessing folder
    print "cleaning sampling folder"
    os.system('rm -rf ' + sampleCaseRMT + '/postProcessing')
    os.system('rm -rf ' + sampleCaseRMT + '/log.sample')
    os.system('rm -rf ' + sampleCasePhy + '/postProcessing')
    os.system('rm -rf ' + sampleCasePhy + '/log.sample')
    for componentName in componentPool:
        print "processing sampling for " + componentName + " ..."
        replace(sampleFile_RMT, "<scalarField>", componentName); 
        replace(sampleFile_Phy, "<scalarField>", componentName);
        os.system('sample -case ' + sampleCaseRMT + ' >> log.sample')
        os.system('sample -case ' + sampleCasePhy + ' >> log.sample')
        os.system('cp ' + sampleCaseRMT+'/system/sampleDict.noRun ' + sampleFile_RMT)
        os.system('cp ' + sampleCasePhy+'/system/sampleDict.noRun ' + sampleFile_Phy)
        
    for componentName in componentEntPool:
        print "processing sampling for " + componentName + " ..."
        replace(sampleFile_RMT, "<scalarField>", componentName); 
        replace(sampleFile_Phy, "<scalarField>", componentName);
        os.system('sample -case ' + sampleCaseRMT + ' >> log.sample')
        os.system('sample -case ' + sampleCasePhy + ' >> log.sample')        
        os.system('cp ' + sampleCaseRMT+'/system/sampleDict.noRun ' + sampleFile_RMT)
        os.system('cp ' + sampleCasePhy+'/system/sampleDict.noRun ' + sampleFile_Phy)
            
    for componentName in componentKLDPool:
        print "processing sampling for " + componentName + " ..."
        replace(sampleFile_RMT, "<scalarField>", componentName); 
        os.system('sample -case ' + sampleCaseRMT + ' >> log.sample')      
        os.system('cp ' + sampleCaseRMT+'/system/sampleDict.noRun ' + sampleFile_RMT)           
    
# Main function
if __name__ == '__main__':
    mainInputFile = './mainInput.in'
    plotInputFile = './plotInfo.in'
    paramDict = readInputData(plotInputFile)
    paramDict_main = readInputData(mainInputFile)
    case_compared = paramDict['case_compared'] 
    caseName = paramDict_main['caseName']
    resultDir_RMT = 'resultData/'
    resultDir_Phy = case_compared + '/debugData/init/'
    RComponentDir = 'RComponent_samples/'
    deltaRComDir = 'deltaRComponent_samples/'  
    caseRMT = '.'
    casePhy = case_compared
    # parse plot control      
    allEntFlag = ast.literal_eval(paramDict['scalarEntFlag'])
    EntCFlag = ast.literal_eval(paramDict['EntCFlag'])
    EntXiEtaFlag = ast.literal_eval(paramDict['EntXiEtaFlag'])
    EntTKEFlag = ast.literal_eval(paramDict['EntTKEFlag'])
    EntVFlag = ast.literal_eval(paramDict['EntVFlag'])
    EntdeltaXiEtaFlag = ast.literal_eval(paramDict['EntdeltaXiEtaFlag'])
    EntdeltaKFlag = ast.literal_eval(paramDict['EntdeltaKFlag'])    
    EntdeltaVFlag = ast.literal_eval(paramDict['EntdeltaVFlag'])
    
    allMomentFlag = ast.literal_eval(paramDict['allMomentFlag'])
    
    if allEntFlag or EntCFlag:
        namePool = ['c1', 'c2', 'c3']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'RComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName)    
    if allEntFlag or EntXiEtaFlag:
        namePool = ['Xi', 'Eta']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'RComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName)                
    if allEntFlag or EntTKEFlag:
        namePool = ['TKE']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'RComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName) 
    if allEntFlag or EntVFlag:
        namePool = ['VA', 'VB', 'VC']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'RComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName) 
    if allEntFlag or EntdeltaXiEtaFlag:
        namePool = ['deltaXi', 'deltaEta']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'deltaRComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName)                                                  
    if allEntFlag or EntdeltaKFlag:
        namePool = ['deltaK']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'deltaRComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName) 
                                                     
    if allEntFlag or EntdeltaVFlag:
        namePool = ['deltaVA', 'deltaVB', 'deltaVC']
        for componentName in namePool:
            print "calclulate Entropy field for " + componentName + " ..."
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'deltaRComponent_samples/'+ componentName +'_s') 
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')
            SField_RMT, SField_Phy, SField_Phy2RMT = pF.estimateEntropyScalarField(scalarFields_RMT, 
                                                     scalarFields_Phy, caseRMT, casePhy, caseName, componentName)                                                     

    if allMomentFlag:
        namePool = ['deltaK']
        for componentName in namePool:
            print "calclulate moment field for " + componentName + " (RMT Cases) ... "
            scalarFields_RMT = np.loadtxt(resultDir_RMT + 'deltaRComponent_samples/'+ componentName +'_s') 
            pF.estimateFourMomentScalarField(scalarFields_RMT, caseRMT, caseName, componentName)
            print "calclulate moment field for " + componentName + " (Phy Cases) ... "            
            scalarFields_Phy = np.loadtxt(resultDir_Phy + componentName + '_s')           
            pF.estimateFourMomentScalarField(scalarFields_Phy, casePhy, caseName, componentName)

    #sample all components
    sampleAllComponents(caseName, caseRMT, casePhy)
    
       
