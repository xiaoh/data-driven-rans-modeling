#!/aoe/bin/python27

# description        :Contain all functions for postprocessing the result of maxEnt

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.03, 2015
# revision           :Dec.05, 2015
####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.stats as ss
# system, file operation
import pdb
import os
import time
# plotting
import matplotlib.pyplot as plt
# Import local modules
import foamFileOperation as foamOp
from utilities import readInputData, extractListFromDict, replace
from ReynoldsStressRF import ReynoldsStressRF
def readMatrixSample(folderPath):
    """
    read matrix(2D-array) files(txt) in a folder, and combine them to be a (matrix ensemble) 3D-array
    
    Arg:
        folderPath: path of the directory containing matrix txt files
    return:
        matrixEns: matrix ensemble [nSample by nCell by 6]
    """

    nSample = 0
    for f in os.listdir(folderPath):
        filePath = os.path.join(folderPath, f)
        if os.path.isfile(filePath):
            # estimate number of the files = nSample
            nSample = nSample + 1
            # estimate number of the cells = nCell
            if nSample == 1:
                M = np.loadtxt(filePath)
                nCell, nTensor = M.shape
    
    # open space for storing of matrix ensemble [nSample by nCell by 6]
    matrixEns = np.zeros([nSample, nCell, nTensor])
    isample = 0
    for f in os.listdir(folderPath):
        filePath = os.path.join(folderPath, f)
        if os.path.isfile(filePath):
            matrixEns[isample, :, :] = np.loadtxt(filePath)
            isample = isample + 1
    return matrixEns

def statEva_scalarField(scalarFields):
    """
    perform statistical analysis of a scalar field ensemble (e.g., mean field, variance etc.)

    :arg:
        scalarFields: scalar field samples [nSample by nCell]
    :return:

    """
    meanField = np.mean(scalarFields, axis=0)
    varField = np.var(scalarFields, axis=0)
    stdField = np.std(scalarFields, axis=0)
    return meanField, varField, stdField

def statEva_2scalarFields(scalarField1, scalarField2, bw_method=0.3):
    """

    :param scalarField1: n by 1
    :param scalarField2: n by 1
    :return:
    """

    twoScalarFields = np.vstack((scalarField1, scalarField2))
    density2D = ss.gaussian_kde(twoScalarFields, bw_method)

    return density2D

def coord2cellIdx(coord, casefolder):
    """

    :param XYZ:
    :return:
    """

    foamOp.writeLocToFile(coord, casefolder+'/constant/obsLocations')
    os.system('getHostCellIdx -case '+ casefolder)
    cellidx = np.loadtxt(casefolder+'/constant/indexHost.txt')
    #pdb.set_trace()
    return cellidx

def getVerandHorCoord(x0, y0, z0, case_meshPath):
    """

    :param x0:
    :param y0:
    :return:
    """
    sampleDictPath = case_meshPath + '/system/sampleDict'
    sampleDictPathBck = case_meshPath + '/system/sampleDict.org'
    replace(sampleDictPath, "<x0>", str(x0))
    replace(sampleDictPath, "<y0>", str(y0))
    os.system('sample -case ' + case_meshPath)
    os.system('cp '+ sampleDictPathBck + ' ' + sampleDictPath)
    ylineU = np.loadtxt(case_meshPath+'/postProcessing/sets/0/line_x0_U.xy')
    yline = ylineU[:, 0]
    xlineU = np.loadtxt(case_meshPath+'/postProcessing/sets/0/line_y0_U.xy')
    xline = xlineU[:, 0]

    x0V_v = np.ones(yline.shape) * x0
    z0V_v = np.ones(yline.shape) * z0

    y0V_h = np.ones(xline.shape) * y0
    z0V_h = np.ones(xline.shape) * z0

    vCoord = np.vstack((x0V_v, yline, z0V_v)).T
    hCoord = np.vstack((xline, y0V_h, z0V_h)).T

    return hCoord, vCoord

def genKSamples(baselineDir, nCell):
    tauFile = '/0/Tau'
    mapTau = ReynoldsStressRF(baselineDir, tauFile, nCell, 1, 'True')
    Tau_base = mapTau.tau
    deltaFolder = './resultData/deltaRComponent_samples/'
    deltaLogkFile = deltaFolder+'deltaK_s'
    deltaLogks = np.loadtxt(deltaLogkFile)
    [nSample, ncell] = deltaLogks.shape
    k_base = Tau_base[:, 0] + Tau_base[:, 3] + Tau_base[:, 5]
    Logk_base = np.log(k_base)
    Logk_bases = np.tile(Logk_base, (nSample, 1))
    k_s = np.exp(Logk_bases + deltaLogks)
    # normalization
    Ub = 0.028
    k_snorm = k_s / Ub / Ub 
    k_basenorm = k_base / Ub / Ub   
    np.savetxt('./resultData/RComponent_samples/k_s', k_snorm)
    np.savetxt('./resultData/RComponent_samples/k_bar', k_basenorm)

def estimateEntropyScalarField(scalarFields_RMT, scalarFields_Phy, caseFolderRMT, caseFolderPhy, caseName, scalarName, kde_bw=0.2):
    """
    estimate entropy of scalarFields_RMT, scalarFields_Phy and
    estimate KL divergence of scalarFields_Phy with scalarFields_RMT
    """

    nSample, nCell = scalarFields_RMT.shape
    SField_RMT = np.zeros([1, nCell]) 
    SField_Phy = np.zeros([1, nCell]) 
    SField_Phy2RMT = np.zeros([1, nCell]) 
    tic = time.time()
    for icell in range(nCell):
        n_bins = 500
        
        sample_1Loc_RMT = scalarFields_RMT[:, icell]
        scalarMax_RMT = np.max(sample_1Loc_RMT)
        scalarMin_RMT = np.min(sample_1Loc_RMT)
        bins_RMT = np.linspace(scalarMin_RMT, scalarMax_RMT, n_bins)        
        
        sample_1Loc_Phy = scalarFields_Phy[:, icell]
        scalarMax_Phy = np.max(sample_1Loc_Phy)
        scalarMin_Phy = np.min(sample_1Loc_Phy)
        bins_Phy = np.linspace(scalarMin_Phy, scalarMax_Phy, n_bins)
                
        pdf_RMT = ss.kde.gaussian_kde(sample_1Loc_RMT, bw_method=kde_bw)
        pdf_Phy = ss.kde.gaussian_kde(sample_1Loc_Phy, bw_method=kde_bw)
        
        SField_RMT[0, icell] = ss.entropy(pdf_RMT(bins_RMT))
        SField_Phy[0, icell] = ss.entropy(pdf_Phy(bins_Phy))
        SField_Phy2RMT[0, icell] = ss.entropy(pdf_Phy(bins_Phy), pdf_RMT(bins_RMT))
    toc = time.time()
    print "collapse for estimating S of one scalar field is ", toc - tic, " s"
    foamOp.writeScalarToFile(SField_RMT.reshape((nCell, 1)), caseFolderRMT+'/'+caseName+'_sample/0/'+scalarName+'_SField')
    foamOp.writeScalarToFile(SField_Phy.reshape((nCell, 1)), caseFolderPhy+'/'+caseName+'_sample/0/'+scalarName+'_SField')
    foamOp.writeScalarToFile(SField_Phy2RMT.reshape((nCell, 1)), caseFolderRMT+'/'+caseName+'_sample/0/'+scalarName+'_KLDistField')
    
    return SField_RMT, SField_Phy, SField_Phy2RMT

def estimateFourMomentScalarField(scalarFields, caseFolder, caseName, scalarName):
    """
    estimate 1st moment (mean) of scalar fields samples
    estimate 2nd moment (variance) of scalar fields samples
    estimate 3rd moment (Skewness) of scalar fields samples
    estimate 4th moment (Kurtosis) of scalar fields samples
    """
    tic = time.time()
    nSample, nCell = scalarFields.shape
    # calculate 1st moment (mean) of scalar fields samples
    meanField = np.mean(scalarFields, axis=0)
    # calculate 2nd moment (variance) of scalar fields samples
    varField = np.var(scalarFields, axis=0)
    stdField = np.std(scalarFields, axis=0)
    # calculate 3rd moment (Skewness) of scalar fields samples    
    # calculate 4th moment (Kurtosis) of scalar fields samples
    SkewnessField = np.zeros([1, nCell])
    KurtosisField = np.zeros([1, nCell])     
    for icell in range(nCell):        
        sample_1Loc = scalarFields[:, icell]        
        SkewnessField[0, icell] = ss.skew(sample_1Loc)
        # Note, here we calculate the excess kurtosis, normal distribution is zero        
        KurtosisField[0, icell] = ss.kurtosis(sample_1Loc, fisher=True, bias=True)
    toc = time.time()
    
    print "collapse for estimating Moment of one scalar field is ", toc - tic, " s"        
    foamOp.writeScalarToFile(meanField.reshape((nCell, 1)), caseFolder+'/'+caseName+'_sample/0/'+scalarName+'_meanField')
    foamOp.writeScalarToFile(varField.reshape((nCell, 1)), caseFolder+'/'+caseName+'_sample/0/'+scalarName+'_varField')
    foamOp.writeScalarToFile(SkewnessField.reshape((nCell, 1)), caseFolder+'/'+caseName+'_sample/0/'+scalarName+'_SkewnessField')        
    foamOp.writeScalarToFile(KurtosisField.reshape((nCell, 1)), caseFolder+'/'+caseName+'_sample/0/'+scalarName+'_KurtosisField')

   
        
if __name__ == '__main__':
    #RPath = './resultData/R_samples'
    #matrixEns = readMatrixSample(RPath)
    baselineDir = './pehill_base'
    nCell = 3000
    genKSamples(baselineDir, nCell)
    pdb.set_trace()

