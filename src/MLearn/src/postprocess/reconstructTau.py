#!/usr/bin/env python

# description        :Reconstruct Tau - Learned deltaTau (different reconstructions based on different paramterization)
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.21, 2016
# revision           :Oct.22, 2016

########################################################################################################################

## Import system modules
# sci computing
import numpy as np
# system, file operation
import pdb
import os
import subprocess
import shutil
# sklearn importing
from sklearn.ensemble.forest import RandomForestRegressor
#from sklearn.ensemble.forest import RandomForestRegressor
# plotting
import matplotlib.pyplot as plt  # for plotting
import matplotlib as mp
# Import local modules
import utilities as uT
import outputFormat as outF
import genDeltaTauTruth as gD
import foamFileOperation as foamOp
import ReynoldsStressRF as ReRF



class reconstructPredictedTau:
    """
    reconstruct Reynolds stress field from predicted deltaTau component and baseline Reynolds stress field

    Arg:
    """

    def __init__(self, currentRANSDir, testFlow, testFlowRe, nMesh, deltaOrienParamerization, *predictResponseList):
        self.nCell = nMesh
        self.predictionsDir = './predictions'
        self.timeDir = '0.000000'
        self.currentRANSDir = currentRANSDir
        self.testFlow = testFlow
        self.testFlowRe = testFlowRe
        self.deltaOrienParamerization = deltaOrienParamerization
        if self.deltaOrienParamerization == 'Euler':
            self.deltaTauData_headerNames = gD.deltaTau_headerNames_Euler
        elif self.deltaOrienParamerization == 'Quaternion':
            self.deltaTauData_headerNames = gD.deltaTau_headerNames_Quaternion[:-1] #TODO: remove index list? but how to get index switching list???
                
        if len(predictResponseList) == 0:
            print "the prediction delta tau components will be read from dir: ./predictions/"
            self.predictResponseList = self._readFromPredictionFiles(self.deltaTauData_headerNames)
        else:
            self.predictResponseList = predictResponseList[0]
        # write predicted delta tau components to OpenFOAM field files
        self._writePredictedDiscrepancyFieldToOpenFOAM()

    def reconstructCorrectedTau(self):
        """
        reconstruct the corrected Tau with predicted discrepancies
        """
        caseZeroDir = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe, self.timeDir)
        mapTau = ReRF.ReynoldsStressRF(caseZeroDir, 'Tau', self.nCell, 1, correctInitTau='False')
        
        for component in self.deltaTauData_headerNames[4:]:
            exec(component + " = self.predictResponseList[component][:, -1]")
            exec(component + " = np.array([eval(component)]).T")
            
        if self.deltaOrienParamerization == 'Euler': 
            deltaV = np.hstack([deltaVA, deltaVB, deltaVC])
            tau_reconstruct = mapTau.perturbTauInBary(deltaXb, deltaYb, deltaLog2K, deltaV)               
        elif self.deltaOrienParamerization == 'Quaternion':

            tau_reconstruct = mapTau.perturbTauInBary(Tau_base,deltaXb,deltaYb,deltaLog2K,vx,vy,vz,theta,indexList)
        
        tauPath = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe, 
                                         self.timeDir, 'Tau')
        foamOp.writeTurbStressToFile(tau_reconstruct, tauPath)                                        
        ReCase = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe)
        shutil.copy(os.path.join(ReCase,'system/sampleDict.template'),os.path.join(ReCase,'system/sampleDict'))
        uT.replace(os.path.join(ReCase,'system/sampleDict'), '<field>', 'Tau')
        subprocess.check_call('sample >> log.samplePredict', shell=True, cwd=ReCase)         
        
    def _readFromPredictionFiles(self, deltaTauData_headerNames):
        """
        read predict delta tau component from files in ./predictions/
        
        Note: deltaTauData_headerNames is ['Index', 'ccx', 'ccy', 'ccz', .....], thus loop is from [4:]
        """
        predictResponseList = {}
        for component in deltaTauData_headerNames[4:]:
            componentDir = os.path.join(self.predictionsDir, 'predict_'+component)
            predictResponseList[component] = np.loadtxt(componentDir)
        return predictResponseList
         
    def _writePredictedDiscrepancyFieldToOpenFOAM(self):
        """
        write predicted delta tau components to OpenFOAM field files for further visualization
        """
        print "writting all learning components as OpenFOAM field files in 0.000000 folder" 
        
        for component in self.deltaTauData_headerNames[4:]:
            scalarField = self.predictResponseList[component][:, -1] # last column is the field
            componentPath = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe, 
                                         self.timeDir, component+'_predict')
            templateFilePath = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe, 
                                         self.timeDir, 'scalarFieldTemplate')
            if os.path.exists(componentPath):
                os.remove(componentPath)
            shutil.copy(templateFilePath, componentPath)
            uT.replace(componentPath, '<Feature>', component)
            foamOp.writeScalarToFile(scalarField, componentPath)
            
            # sample the truth
            ReCase = os.path.join(self.currentRANSDir, self.testFlow+'Re'+self.testFlowRe)
            shutil.copy(os.path.join(ReCase,'system/sampleDict.template'),os.path.join(ReCase,'system/sampleDict'))
            uT.replace(os.path.join(ReCase,'system/sampleDict'), '<field>', component+'_predict')
            subprocess.check_call('sample >> log.samplePredict', shell=True, cwd=ReCase)        

    
