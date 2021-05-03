#!/usr/bin/env python

# description        :ML learning class (contain different regression function)
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.17, 2016
# revision           :Oct.21, 2016

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
import outputFormat as outF




class regressionSolver:
    """
    Build Machine Learning (Regression) Model, Superclass for other child classes

    Arg:
    """

    def __init__(self, trainFlows, trainFlowsRes, testFlow, testFlowRe, deltaOrienParamerization,
                 featureFolder, trainDataFolder, deltaTauData_headerNames):
        self.trainFlows = trainFlows
        self.trainFlowsRes = trainFlowsRes
        self.testFlow = testFlow
        self.testFlowRe = testFlowRe
        self.deltaOrienParamerization = deltaOrienParamerization
        self.featureFolder = featureFolder
        self.trainDataFolder = trainDataFolder
        self.deltaTauData_headerNames = deltaTauData_headerNames
        print "training: ", zip(trainFlows, trainFlowsRes)
        print "test: ", testFlow, testFlowRe
        self.rawTestFeatureM, self.rawTrainFeatureMergedM, self.rawTrainDataDict = self._readCaseFiles()
        self.nCluster = int(np.max(self.rawTrainDataDict['clusterMarker']) + 1)
        self.nCellTest = self.rawTestFeatureM.shape[0]
        self.nCellTrains = self.rawTrainFeatureMergedM.shape[0]
    
    def extractMatrixGivenCluster(self, M, idxCluster):
        """
        to extract submatrix from matrix M, with given cluster index: idxCluster. The first column
        of M should be index of cluster
        """
        cellIdx_idxCluster = np.where(M[:, 0] == idxCluster)
        cellIdx_idxCluster = cellIdx_idxCluster[0]
        subM = M[cellIdx_idxCluster , :]    
        return subM, cellIdx_idxCluster
    
    def splitRawFeatureM(self, rawFeatureM):
        """
        split raw feature matrix (clusterIdx, meshIdex, ccx, ccy, ccz, q1, q2, ...) to pure input features matrix
        """
        featureInputM = rawFeatureM[:, 5:]
        meshCoord = rawFeatureM[:, 2:5]
        return featureInputM, meshCoord
        
    def _readCaseFiles(self): 
        """
        Interface of reading feature and deltaTau files for training and test flows
        
        Arg:
            rawTestFeatureM: meshIdex, ccx, ccy, ccz, q1, q2, ...
            rawTrainFeatureMergedM: meshIdx, ccx, ccy, ccz, q1, q2, ...
            rawTrainDataDict: dict: meshIdx, ccx, ccy, ccz, deltaTau, ..
            
        """
        # read training raw input, output
        rawTrainFeatureMList = []
        rawTrainDataMList = []
        for trainFlow, trainFlowRe in zip(self.trainFlows, self.trainFlowsRes):            
            
            trainFlowFeatureDir = os.path.join(self.featureFolder, 'features-'+trainFlow+trainFlowRe)
            rawTrainFeatureM = np.loadtxt(trainFlowFeatureDir)
                        
            trainFlowDataDir = os.path.join(self.trainDataFolder, 'deltaTau_'+trainFlow+trainFlowRe+
                                               '_'+self.deltaOrienParamerization)
            rawTrainDataM = np.loadtxt(trainFlowDataDir)                                          
            # match training feature and data by mesh index            
            cellIdx_data = rawTrainDataM[:, 0].astype(int)
            rawTrainFeatureM_selected = rawTrainFeatureM[cellIdx_data, :]
            # put training feature and data to list
            rawTrainFeatureMList.append(rawTrainFeatureM_selected)
            rawTrainDataMList.append(rawTrainDataM)
        # aggregate the raw train feature and data
        rawTrainFeatureMergedM = rawTrainFeatureMList[0]
        rawTrainDataMergedM = rawTrainDataMList[0]
        if len(rawTrainFeatureMList) > 1:
            for i in range(len(rawTrainFeatureMList)-1):
                rawTrainFeatureMergedM = np.vstack([rawTrainFeatureMergedM, rawTrainFeatureMList[i+1]])
                rawTrainDataMergedM = np.vstack([rawTrainDataMergedM, rawTrainDataMList[i+1]])
        # separate training data (tau discrepancy) matrix
        rawTrainDataDict = self._separateDataMatrix(rawTrainDataMergedM)
        # read test raw input, output
        testFlowFeatureDir = os.path.join(self.featureFolder, 'features-'+self.testFlow+self.testFlowRe)
        rawTestFeatureM = np.loadtxt(testFlowFeatureDir)                                  
        # get clustering marker
        rawTestFeatureM, rawTrainFeatureMergedM,  rawTrainDataDict = \
        self._featureClustering(rawTestFeatureM, rawTrainFeatureMergedM,  rawTrainDataDict)
        
        return rawTestFeatureM, rawTrainFeatureMergedM,  rawTrainDataDict
   

    def _featureClustering(self, rawTestFeatureM, rawTrainFeatureMergedM,  rawTrainDataDict, *clusterArg):
        """
        classify feature space to sub-feature space for learning
        """
        
        if clusterArg == ():
            # no clustering
            nRawTestFeature = rawTestFeatureM.shape[0]
            nRawTrainFeatureMergedSelected = rawTrainFeatureMergedM.shape[0]
            nRawTrainData = len(rawTrainDataDict['Index'])
            if nRawTrainFeatureMergedSelected != nRawTrainData:
                print "Error: number of train feature data is not equal to # of train data"
                sys.exit(1)
            clusterIdx_test = np.ones((nRawTestFeature, 1))*0
            clusterIdx_train = np.ones((nRawTrainData, 1))*0
            
            rawTestFeatureM_clusterMarked = np.hstack([clusterIdx_test, rawTestFeatureM])
            rawTrainFeatureMergedM_clusterMarked = np.hstack([clusterIdx_train, rawTrainFeatureMergedM])
            rawTrainDataDict['clusterMarker'] = clusterIdx_train
        
        if len(clusterArg) > 0:
            pass
            #to be implemented
        
        return rawTestFeatureM_clusterMarked, rawTrainFeatureMergedM_clusterMarked, rawTrainDataDict
            
    def _separateDataMatrix(self, rawTrainDataM):
        """
        separate raw data matrix
        """
        print "separate raw training data matrix as variables:"
        print self.deltaTauData_headerNames
        i = 0
        rawTrainDataDict = {}
        for var in self.deltaTauData_headerNames:
            rawTrainDataDict[var] = rawTrainDataM[:, i]
            i = i + 1
        return rawTrainDataDict
        
                    

class randomForestRegressionSolver(regressionSolver):
    """
    Build Random Forest Machine Learning (Regression) Model, subclass of  for regressionSolver

    Arg:
    """
    def __init__(self, trainFlows, trainFlowsRes, testFlow, testFlowRe, deltaOrienParamerization,
                 featureFolder, trainDataFolder, deltaTauData_headerNames):
        #inherit init value from regression class.
        regressionSolver.__init__(self, trainFlows, trainFlowsRes, testFlow, testFlowRe, deltaOrienParamerization,
                                  featureFolder, trainDataFolder, deltaTauData_headerNames)
        self.nTree = 100 #TODO make it to be given outside
    
    def learnPrediction(self, responseCombinationFlag='combineLearning'):
        """
        """
        
        # clear prediction folder
        predictionDir = './predictions'
        if os.path.exists(predictionDir):
            print "remove previous predictions results"
            shutil.rmtree(predictionDir)            
        subprocess.check_call('mkdir '+predictionDir, shell=True)      
        # initialize random Forest regressor
        regModel = RandomForestRegressor(n_estimators=self.nTree, criterion='mse', max_depth=None, min_samples_split=1.0, 
                                         min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_features="sqrt", 
                                         max_leaf_nodes=None, bootstrap=True, oob_score=True, n_jobs=10, 
                                         random_state=None, verbose=0, warm_start=False)
        testPredictionClusterList = []
        for iCluster in range(self.nCluster):
            outF.buildSeparateLine('.')
            print "regression for cluster-", iCluster
               
            # prepare training input for iCluster
            rawTrainFeatureMergedM_iCluster, cellIdx_idxCluster_train = regressionSolver.extractMatrixGivenCluster\
                                                                        (self, self.rawTrainFeatureMergedM, iCluster)
            trainFeatureMergedM_iCluster, meshCoord_train = regressionSolver.splitRawFeatureM(self, 
                                                           rawTrainFeatureMergedM_iCluster)
            
            # prepare test input for iCluster
            rawTestFeatureM_iCluster, cellIdx_idxCluster_test = regressionSolver.extractMatrixGivenCluster\
                                                                (self, self.rawTestFeatureM, iCluster)            
            testFeatureM_iCluster, meshCoord_test = regressionSolver.splitRawFeatureM(self, rawTestFeatureM_iCluster)            
            
            # prepare training output for iCluster
            testPredicts = {}
            testPredicts['cellIdx_idxCluster'] = cellIdx_idxCluster_test
            testPredicts['meshCoord_test'] = meshCoord_test
            learningResponseList = self.deltaTauData_headerNames[4:]
            
            if responseCombinationFlag == 'separateLearning':                
                if self.deltaTauData_headerNames[-1] == 'indexList':
                    learningResponseList = learningResponseList[:-1]
                for learningResponse in learningResponseList:
                    trainData_iCluster = self.rawTrainDataDict[learningResponse]
                    print "    Learning for ", learningResponse
                    regModel.fit(trainFeatureMergedM_iCluster, trainData_iCluster)
                    # get corresponding feature importance
                    featureImportance = regModel.feature_importances_
                    np.savetxt('./predictions/featureImportance_'+learningResponse+'_cluster'+str(iCluster),  featureImportance)
                    print "    oob score for ", learningResponse, ' is ', regModel.oob_score_, '\n'
                    # prediction
                    testPredicts[learningResponse] = regModel.predict(testFeatureM_iCluster)                
                
            
            elif responseCombinationFlag == 'combineLearning':
                learningResponseList = self.deltaTauData_headerNames[4:]
                if self.deltaTauData_headerNames[-1] == 'indexList':
                    learningResponseList = learningResponseList[:-1]                
                
                responseDataCombine = np.zeros([self.nCellTrains, len(learningResponseList)]) 
                i = 0
                for learningResponse in learningResponseList:
                    responseDataCombine[:, i] = self.rawTrainDataDict[learningResponse]
                    i = i + 1
                print "    Learning for combined response"
                regModel.fit(trainFeatureMergedM_iCluster, responseDataCombine) 
                featureImportance = regModel.feature_importances_
                np.savetxt('./predictions/featureImportance_combined_cluster'+str(iCluster),  featureImportance)
                print '    oob score for combined learning is ', regModel.oob_score_, '\n'
                # prediction
                testPredicts_combined = regModel.predict(testFeatureM_iCluster)
                i = 0
                for learningResponse in learningResponseList:
                    testPredicts[learningResponse] = testPredicts_combined[:, i]
                    i = i + 1
            
            testPredictionClusterList.append(testPredicts)
        
        # reconstruct predictions of each cluster
        predictResponseDict = {}
        for predictResponseName in learningResponseList:
            
            predictionReconstruct = np.array([]) 
            
            for iCluster in range(self.nCluster):            
                
                predict_iCluster = testPredictionClusterList[iCluster]            
                predictM = predict_iCluster[predictResponseName]
                if predictM.ndim == 1:
                    predictM = np.array([predictM]).T
                meshCoordM = predict_iCluster['meshCoord_test']
                cellClusterIdxM = np.array([predict_iCluster['cellIdx_idxCluster']]).T
                predictResult_clusterI = np.hstack([cellClusterIdxM, meshCoordM, predictM])
                if iCluster == 0:
                    predictionReconstruct = predictResult_clusterI
                else:
                    predictionReconstruct = np.vstack([predictionReconstruct, predictResult_clusterI])
            # sort cell idx
            predictionReconstruct = predictionReconstruct[predictionReconstruct[:, 0].argsort()]
            header_iresponse = self.deltaTauData_headerNames[0:4] + [predictResponseName]
            headStr = ''
            for headerName in header_iresponse:
                headStr = headStr + headerName + '      '
            headStr = headStr + '\n'
            
            np.savetxt('./predictions/predict_'+predictResponseName,  predictionReconstruct, header = headStr)
            predictResponseDict[predictResponseName] = predictionReconstruct

        return predictResponseDict
    
        
class neuralNetworkRegressionSolver(regressionSolver):
    """
    Build Neural Network Machine Learning (Regression) Model, subclass of  for regressionSolver

    Arg:
    """
    def __init__(self, regression_caseVec, prediction_case):
        pass

class linearRegressionSolver(regressionSolver):
    """
    Build Linear Regression Machine Learning (Regression) Model, subclass of  for regressionSolver

    Arg:
    """
    def __init__(self, regression_caseVec, prediction_case):
        pass
    
