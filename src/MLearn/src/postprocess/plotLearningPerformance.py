#!/usr/bin/env python

# description        :Plot learning performance (comparing learning results with the truth and training data)
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.24, 2016
# revision           :Oct.24, 2016

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



class plotLearningPerformance:
    """
    Plot learning performance (comparing learning results with the truth and training data)

    Arg:
    """

    def __init__(self, currentRANSDir, trainFlows, trainFlowsRes, testFlow, testFlowRe, deltaOrienParamerization):
        self.timeDir = '0'
        self.currentRANSDir = currentRANSDir
        self.testFlow = testFlow
        self.testFlowRe = testFlowRe
        self.deltaOrienParamerization = deltaOrienParamerization
        if self.deltaOrienParamerization == 'Euler':
            self.deltaTauData_headerNames = gD.deltaTau_headerNames_Euler
        elif self.deltaOrienParamerization == 'Quaternion':
            self.deltaTauData_headerNames = gD.deltaTau_headerNames_Quaternion[:-1] 
                
        self.testCase = self.testFlow+'Re'+self.testFlowRe
        self.trainCases = []
        for trainFlow, trainFlowRe in zip(trainFlows, trainFlowsRes):
            self.trainCases.append(trainFlow+'Re'+trainFlowRe)
        self.colors = ['blue', 'black', 'cyan', 'gray', 'darkblue', 'darkgreen']        
        self.currentDir = os.getcwd()
                
    def plotPhyDeltaTau_compareAll(self):
        """
        Plot learning prediction compared with benchmark and training data
        """
        exec('import ' + self.testFlow + 'PlotGivenCases as postPlt')

        os.chdir(self.currentRANSDir)
        for component in self.deltaTauData_headerNames[4:]:
            postPlt.fullCompareScalar(component, self.currentRANSDir, self.testCase, self.trainCases, self.timeDir, 
            1, self.colors, False, False)  
        subprocess.check_call('mv *.pdf ../figures', shell=True)
        os.chdir(self.currentDir)

    def plotTau_compareAll(self):
        """
        Plot learning prediction compared with benchmark and training data
        """
        exec('import ' + self.testFlow + 'PlotGivenCases as postPlt')

        os.chdir(self.currentRANSDir)
        postPlt.fullCompareSymmTensor('Tau', self.currentRANSDir, self.testCase, self.trainCases, self.timeDir,
                                      20, self.colors, False, False)  
        subprocess.check_call('mv *.pdf ../figures', shell=True)
        os.chdir(self.currentDir)    
    
        
        
    
