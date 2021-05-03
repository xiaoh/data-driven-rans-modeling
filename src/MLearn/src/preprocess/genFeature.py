#!/usr/bin/env python
# description        :Get Mean Flow Feature Files
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.13, 2016
# revision           :Oct.19, 2016
########################################################################################################################

# SYSTEM MODULE
import numpy as np
import matplotlib as mp
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pdb
import os
import shutil
import subprocess
import sys
# LOCAL MODULE
import foamFileOperation as foamOp
import utilities as uT
#sys.path.append("../../integrity-bases/")
import integrityBasis as base

coordName = ['ccx', 'ccy', 'ccz']
featureNameList_org = [ 'QCriterion', 'TurbulenceIntensity', 'ReTurbulence', 'pGradAlongStream', 
            'TurbulenceTime', 'PressureStress', 'UGradMisalignment', 'ConvectionTKE', 'TauRatio', 'Curvature']
featureLabel_org = [
                    set(['strainRate', 'rotationRate']), set(['K']), set(['wallDist', 'K']), set(['gradP']), 
                    set(['K/epsilon']), set(['nStress/sStress']), set(['UGradMisalignment']), set(['Convec/Product']), 
                    set(['tauN/tauS']), set(['curvature'])
                   ]

class constructFeatureMain():
    """
    main function to construct input feature space
    Arg:
        caseType:               flow type
        Re:                     Reynold number
        turbulence_model:       turbulence model used
        constructionKeyword:    Full, FullSelected, Scalar, ScalarSelected, Expended, ExpendedSelected
    """


    def __init__(self, caseType, Re, nMesh, turbulence_model, combinationKeyword, *labelList): 
        """
        initialization
        """
        self.caseType = caseType
        self.Re = Re
        self.nMesh = nMesh
        self.turbulence_model = turbulence_model
        self.combinationKeyword = combinationKeyword
        self.labelList = labelList
        srcDir = os.path.normpath(os.path.join(os.path.abspath(__file__), '../../..'))
        self.baseDir = os.path.join(srcDir, 'RANSTemp')
        nExpandedFeature = 47
        self.featureNameList_expand = ['q'+str(pn) for pn in range(nExpandedFeature)]
        self.featureNameList_full = featureNameList_org + self.featureNameList_expand
    
    def genFeature(self):
        """
        generate input feature files
        """
        dataDir = os.path.join(self.baseDir, self.caseType, 'cell'+self.nMesh, 'cases_'+self.turbulence_model, 
                               'Re'+self.Re, '0')  
        # obtain mesh information
        cellIdx = np.array([np.arange(int(self.nMesh))]).T
        coordMatrix = foamOp.readTurbCoordinateFromFile(dataDir)
        
        featureMatrix_expand, featureLabel_expand, featureMatrix_org, featureLabel_org, featureMatrix_full, \
        featureLabel_full = getFullFeatureSet(dataDir, int(self.nMesh))

        if self.combinationKeyword in ['FullSelected', 'Full']:
            featureMatrix = featureMatrix_full
            featureLabel = featureLabel_full
            featureNameList = self.featureNameList_full
        elif self.combinationKeyword in ['ScalarSelected', 'Scalar']:
            featureMatrix = featureMatrix_org
            featureLabel = featureLabel_org
            featureNameList = featureNameList_org            
        elif self.combinationKeyword in ['ExpandedSelected', 'Expanded']:
            featureMatrix = featureMatrix_expand
            featureLabel = featureLabel_expand
            featureNameList = self.featureNameList_expand

        if self.combinationKeyword in ['FullSelected', 'ScalarSelected', 'ExpandedSelected']:
            # select the feature out
            featureMatrixSelected = np.empty((nCell,0), float)
            featureLabelSelected = []
            featureNameListSelected = []            
            for i in range(featureMatrix.shape[1]):
                feature_description = featureLabel[i]
                if feature_description.issubset(user_selection):
                    print feature_description, ' IS SUBSET OF ', self.labelList
                    featureMatrixSelected = np.hstack([featureMatrixSelected, featureMatrix[:,[i]]])
                    featureLabelSelected.append(feature_description)
                    featureNameListSelected.append(featureNameList[i])
                else:
                    print feature_description, ' ISNOT A SUBSET OF ',  user_selection            
            print "I included ", i+1, "features in the learning process"
        
        elif self.combinationKeyword in ['Full', 'Scalar', 'Expanded']:
            featureMatrixSelected = featureMatrix
            featureLabelSelected = featureLabel
            featureNameListSelected = featureNameList        
        else:
            print "Please provide valid combination keyword for features:\n \
                   FullSelected, ScalarSelected, ExpandedSelected\
                   Full, Scalar, Expanded"
        #pdb.set_trace()
        # combine feature files
        M = np.hstack([cellIdx, coordMatrix, featureMatrixSelected])
        headerNames = ['Index'] + coordName + featureNameListSelected
        #pdb.set_trace()
        
        # add space between header
        headStr = ''
        for markerName in headerNames:
            headStr = headStr + markerName + '      '
        headStr = headStr + '\n'
            
        # write to txt files
        resultFolder = './featureSet_' + self.combinationKeyword + '_userSelectedN' + str(len(self.labelList))                
        outputName = os.path.join(resultFolder, 'features-'+self.caseType+self.Re)
        if os.path.exists(outputName):
            os.remove(outputName)
        np.savetxt(outputName, M, header = headStr)
            

########################################################################################################################
########################################## Utility Functions ###########################################################
########################################################################################################################

def getFullFeatureSet(dataDir, nCell):
    """
    get full feature set (original 10 features and expended 47 features) with their labels
    
    Arg:
        dataDir: dir containing raw feature openFOAM files
    """
    
    # get original features and corresponding labels
    featureMatrix_org = np.zeros([nCell, len(featureNameList_org)]) 
    i = 0
    for featureName in featureNameList_org:
        #print '#', i+1, 'read field', featureName
        fileName = os.path.join(dataDir, featureName)
        feature = foamOp.readScalarFromFile(fileName)
        featureMatrix_org[:, i] = feature
        i = i + 1

    # get expanded features and corresponding labels
    user_selection = {'strainRate', 'rotationRate', 'gradP', 'gradK'} # this label should be the same as in integrityBasis.py
    featureMatrix_expand, featureLabel_expand = expandFeatures(dataDir, user_selection, False)
    
    featureMatrix_full = np.hstack((featureMatrix_org, featureMatrix_expand))
    featureLabel_full = featureLabel_org + featureLabel_expand
    
    return featureMatrix_expand, featureLabel_expand, featureMatrix_org, featureLabel_org, \
    featureMatrix_full, featureLabel_full

def expandFeatures(caseDir, user_selection, reflectCoord=False):
    """
    Expanding features with given case dir
    
    
    Output:
        featureMatrixSelected
        featureLabels
        
    """
    S_vec, W_vec, dp_vec, dk_vec = generateFourField(caseDir, reflectCoord)
    nCell = S_vec.shape[0]
    featureMatrix, featureLabel = base.field_to_bases(S_vec, W_vec, dp_vec, dk_vec)
    featureMatrixSelected = np.empty((nCell,0), float)
    featureLabelSelected = []
    for i in range(featureMatrix.shape[1]):
        templateFilePath = caseDir+'/scalarFieldTemplate'
        featurei = featureMatrix[:, i]
        qiFilePath = caseDir+'/q'+str(i)
        shutil.copy(templateFilePath, qiFilePath)
        uT.replace(qiFilePath, '<Feature>', 'q'+str(i))
        foamOp.writeScalarToFile(featurei, qiFilePath)
        feature_description = featureLabel[i]
        if feature_description.issubset(user_selection):
            #print feature_description, ' IS SUBSET OF ', user_selection
            featureMatrixSelected = np.hstack([featureMatrixSelected, featureMatrix[:,[i]]])
            featureLabelSelected.append(feature_description)
        else:
            #print feature_description, ' ISNOT A SUBSET OF ',  user_selection
            pass
    
    #print "Finally, I included ", i+1, "expanded features"
    return featureMatrixSelected, featureLabelSelected

def generateFourField(caseDir, reflectCoord):
    """
        generate S, W, dp, dk
        """
    nDim = 3
    S_vec_raw = foamOp.readTurbStressFromFile(os.path.join(caseDir, 'STilda'))
    W_vec_raw = foamOp.readTensorFromFile(os.path.join(caseDir, 'WTilda'))
    dp_vec_raw = foamOp.readVectorFromFile(os.path.join(caseDir, 'gradPTilda'))
    dk_vec_raw = foamOp.readVectorFromFile(os.path.join(caseDir, 'gradKTilda'))
    nCell = S_vec_raw.shape[0]
    S_vec = np.zeros([nCell, nDim, nDim])
    W_vec = np.zeros([nCell, nDim, nDim])
    dp_vec = np.zeros([nCell, nDim, 1])
    dk_vec = np.zeros([nCell, nDim, 1])

    print "read four raw fields from OpenFOAM"
    for i in range(nCell):
        S_vec[i, :, :] = symmflat_To_tensor(S_vec_raw[i, :])
        W_vec[i, :, :] = nonSymmflat_To_tensor(W_vec_raw[i, :])
        dp_vec[i, :, :] = dp_vec_raw[[i], :].T
        dk_vec[i, :, :] = dk_vec_raw[[i], :].T
        if reflectCoord:
            S_vec[i, :, :] = reflectTensor(S_vec[i, :, :])
            W_vec[i, :, :] = reflectTensor(W_vec[i, :, :])
            dp_vec[i, :, :] = reflectVector(dp_vec[i, :, :])
            dk_vec[i, :, :] = reflectVector(dk_vec[i, :, :])
    return S_vec, W_vec, dp_vec, dk_vec

def symmflat_To_tensor(tauflat):
    """
        flatted tau (tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz) to
        matrix form of tau (11, 12, 13; 21, 22, 23; 31, 32, 33)
        """
    tau = np.zeros((3,3))
    
    tau[0,0] = tauflat[0]
    tau[0,1] = tauflat[1]
    tau[0,2] = tauflat[2]
    tau[1,1] = tauflat[3]
    tau[1,2] = tauflat[4]
    tau[2,2] = tauflat[5]
    
    tau[1,0] = tau[0,1]
    tau[2,0] = tau[0,2]
    tau[2,1] = tau[1,2]

    return tau
    
def reflectTensor(Matrix):
    """
    
    """
    Matrix_reflect = Matrix
    Matrix_reflect[0, 1] = -Matrix[0, 1]
    Matrix_reflect[1, 0] = -Matrix[1, 0]
    Matrix_reflect[1, 2] = -Matrix[1, 2]
    Matrix_reflect[2, 1] = -Matrix[2, 1]
    return Matrix_reflect
    
def reflectVector(Vector):
    """
    """
    Vector_reflect = Vector
    Vector_reflect[1, 0] = -Vector[1, 0]
    return Vector_reflect


def nonSymmflat_To_tensor(tauflat):
    """
        flatted tau (tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz) to
        matrix form of tau (11, 12, 13; 21, 22, 23; 31, 32, 33)
        """
    tau = np.zeros((3,3))
    
    tau[0,0] = tauflat[0]
    tau[0,1] = tauflat[1]
    tau[0,2] = tauflat[2]
    tau[1,0] = tauflat[3]
    tau[1,1] = tauflat[4]
    tau[1,2] = tauflat[5]
    tau[2,0] = tauflat[6]
    tau[2,1] = tauflat[7]
    tau[2,2] = tauflat[8]
    return tau



########################################################################################################################
################################################## Not USING ###########################################################
########################################################################################################################

def getFeature(caseType, Re):

    nameList = ['ccx', 'ccy', 'ccz', 'regionFlag', 'QCriterion', 'TurbulenceIntensity', 
                'ReTurbulence', 'pGradAlongStream', 'TurbulenceTime', 'PressureStress', 'UGradMisalignment', 
                'ConvectionTKE', 'TauRatio', 'Curvature']
    dataFolder = baseFolder + caseType + '/Re' + Re + '/0/'
    fileName = dataFolder + 'ccx'
    sizeCheck = foamOp.readScalarFromFile(fileName)
    nCell = sizeCheck.shape[0]
    M = np.zeros([nCell, len(nameList)])

    i = 0
    for featureName in nameList:
        print '#', i+1, 'read field', featureName
        fileName = dataFolder + featureName
        feature = foamOp.readScalarFromFile(fileName)
        M[:, i] = feature
        i = i + 1

    if expandFeaturesFlag:
        featureMatrixSelected = expandFeatures(dataFolder, user_selection)
        if nonSymm == 'exclude':
            print "exclude all nonsymmetric features"
            featureMatrixSelected = featureMatrixSelected[:, symmIdx]
        elif nonSymm == 'include':
            print "include all nonsymmetric features"
            featureMatrixSelected = featureMatrixSelected[:, antiIdx]
        else:
            print "include all features"
            pass
        M = np.hstack([M, featureMatrixSelected])
        headerNames = np.arange(M.shape[1])
        headerNames = [str(pn) for pn in headerNames]
        added = ['ccx', 'ccy', 'ccz', 'regionFlag']
        headerNames = added + headerNames
    else:    
        headerNames = ['ccx', 'ccy', 'ccz', 'regionFlag',
                   'QC', 'TI', 'Re_d',
                   'PGrad', 'TurbTime', 'PStress',
                   'UMisalign', 'ConvTKE', 'TauRatio', 'Curvature'] 

    # add space between header
    headStr = ''
    for markerName in headerNames:
        headStr = headStr + markerName + '      '
    headStr = headStr + '\n'
        
    # write to txt files
    outputName = resultFolder + 'markers-' + caseType + Re
    if os.path.exists(outputName):
        os.system('rm '+outputName)
    np.savetxt(outputName, M, header = headStr)




if __name__ == "__main__":
    reflectCoord = False    
    nonSymm = 'all'
    allIdx = np.arange(47)
    mask = np.ones(len(allIdx), dtype=bool)
    antiIdx = [10, 13, 14, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29, 35, 36, 37, 38, 39, 40]
    #antiIdx = [10]
    mask[antiIdx] = False
    symmIdx = allIdx[mask].tolist()
        
    expandFeaturesFlag = True
    user_selection = {'velocity', 'pressure', 'TKE'}
    resultFolder = 'genResults/'
    baseFolder = '../../markerGen/templateCases/RANS/'
    caseType = 'duct'
    
    for Re in ['2200', '2600', '2900', '3500']:
    #for Re in ['2200', '2600', '2900']:
    #for Re in ['3500']:
        print "We are generating mean flow features for training: case = ", caseType, 'with Re = ', Re
        getFeature(caseType, Re)





