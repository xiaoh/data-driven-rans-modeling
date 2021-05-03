#!/usr/bin/env python

# description        :Generate delta field based on MLE or the estimation of variance of matrix G.

# author             :Jinlong Wu (jinlong@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.01, 2015
# revision           :Dec.01, 2015

####################################################################################################

# comments from HX
# The spacing around * is not consistent. I think you should generally leave space aroudn *, +, -, = etc.
# Please follow the python style document.
# Check with JX; He also as the similar problem.

## Import system modules
# sci computing
import numpy as np
import math
import numpy.linalg as LA

# system, file operation
import pdb
import fnmatch

# plotting
import matplotlib.pyplot as plt 

# local 
import foamFileOperation as foamOp
from utilities import readInputData

class GenerateDelta:
    """
    Generate delta field

    :public member variable:
        self.deltaField:    : generated deltaField (nCell by 1)

    """

    def __init__(self, nCell, nDim = 3):
        """
        Initialization

        :arg:
            nCell           : number of the cell points in one field
            nDim            : dimension of the matrix, default value is 3

        """

        self.nCell = nCell
        self.nDim = nDim
   

    def inferDelta(self, TauDNS, TauBase, nDelta):
        """
        Generate delta field based on MLE.
        
        :arg:
            TauDNS          : observation data (shape: Nobs x 6) 
            nDelta          : the number of delta used to search the maximum likelihood

        :variables:
            deltaPool:      pool of delta 
            pGPool:         calculated pG according to the pool of delta 
            Nobs:           number of observation
            maxIndex:       index to find the maximum pG in pGPool

        :return:
            deltaPool[maxIndex]:  delta which offer maximum likelihood pG
        """

        deltaPool = np.linspace(0, math.sqrt(2)/2, nDelta)
        deltaPool = deltaPool[1:-1]
        print nDelta, "delta is choosed from [", deltaPool[0], ", ", deltaPool[-1], "]"
        pGPool = np.zeros(deltaPool.shape)
        Nobs = TauDNS.shape[0]
        for iterDelta in range(len(deltaPool)):
            for iterN in range(Nobs):
                G = Tau2G(TauDNS[iterN,:], TauBase[iterN,:])
                pGPool[iterDelta] = pGPool[iterDelta] + \
                                    self.evalPG(deltaPool[iterDelta], G)
        maxIndex = np.argmax(pGPool)
        print "Maximum likelihood in log scale log(pG) is: ", pGPool[maxIndex]

        ## Set initial delta field
        self.deltaField = np.zeros([self.nCell, 1])

        ## Infer the delta field
        print "The infered delta based on MLE is: ", deltaPool[maxIndex]
        self.deltaField = self.deltaField * deltaPool[maxIndex]


    #def estimateDelta(self, caseName, Nsample):
    #    """
    #    Generate delta field based on the estimation of variance of G

    #    :arg:
    #        caseName        : The case name to locate the Tau sample files 
    #    """

    #    ## Implement the estimation of variance based on samples of Tau field
    #    ## TODO
    #    ## 1. continue to test the implementation
    #    ## 2. improve the efficiency

    #    ## Possible problem
    #    ## It seems that by capping the Reynolds stress to bottom line (two-component limit),
    #    ## the baseline Tau is nearly semi positive definite which leads to small det(L) and
    #    ## the reconstructed samples of matrix G may have very large determinant, and leads to
    #    ## a large estimation of variance. 

    #    Ncell = self.nCell
    #    ## Load baseline Tau
    #    baseFileName = caseName + '_base/0/TauCap'
    #    TauBase = readTurbStressFromFile(baseFileName)

    #    ## Construct samples of matrix G based on Tau samples
    #    Gsamples = np.zeros([Nsample, Ncell, self.nDim, self.nDim])
    #    for iterS in range(Nsample):
    #        sampleFileName = caseName + '_base-tmp_'+str(iterS+1)+'.0/0/Tau'
    #        TauSample = readTurbStressFromFile(sampleFileName)
    #        print "Currently processing the", iterS, "sample"
    #        for iterC in range(Ncell):
    #            Gsamples[iterS, iterC, :, :] = Tau2G(TauSample[iterC, :], TauBase[iterC, :])
    #            print "Currently processing the", iterC, "cell"

    #    Gmean = np.mean(Gsamples, axis = 0)

    #    ## Estimate the variance of samples of matrix G
    #    varG = np.zeros([Ncell, 1]) 
    #    for iterC in range(Ncell):
    #        for iterS in range(Nsample):
    #            varG[iterC, 0] = varG[iterC, 0] + \
    #                             LA.norm(Gsamples[iterS, iterC, :, :] - Gmean[iterC,:,:])**2
    #        varG[iterC, 0] = varG[iterC, 0] / Nsample

    #    pdb.set_trace()
    #    self.deltaField = varG**0.5

    #    # I think the 5 lines above should be like the following:
    #
    #    #for iterC in range(Ncell):
    #    #    for iterS in range(Nsample):
    #    #        varG[iterC, 0] = varG[iterC, 0] + (LA.norm(Gsamples[iterS, iterC, :] - Gmean[iterC, :])) ** 2

    #    # divide the entire array by Nsample instead of cell-by-cell: to improve efficiency
    #    # varG = varG / Nsample  

    #    # self.deltaField = varG**0.5 # or call sqrt  instead of using **0.5
            

    def estimateDelta(self, caseName, Nsample):
        """
        Generate delta field based on the estimation of variance of R

        :arg:
            caseName        : The case name to locate the Tau sample files 
        """

        ## Implement the estimation of variance based on samples of Tau field

        Ncell = self.nCell
        ## Load baseline Tau
        baseFileName = caseName + '/0/TauCap'
        TauBase = foamOp.readTurbStressFromFile(baseFileName)

        ## Construct samples of matrix G based on Tau samples
        Rsamples = np.zeros([Nsample, Ncell, self.nDim, self.nDim])
        Rmean = np.zeros([Ncell, self.nDim, self.nDim])
        for iterS in range(Nsample):
            sampleFileName = caseName + '-tmp_'+str(iterS+1)+'.0/0/Tau'
            TauSample = foamOp.readTurbStressFromFile(sampleFileName)
            print "Currently processing the", iterS, "sample"
            for iterC in range(Ncell):
                Rsamples[iterS, iterC, :, :] = Tau2R(TauSample[iterC, :])
                Rmean[iterC, :, :] = Tau2R(TauBase[iterC, :])
                #print "Currently processing the", iterC, "cell"

        #Rmean = np.mean(Rsamples, axis = 0)

        ## Estimate the variance of samples of matrix R
        varG = np.zeros([Ncell, 1]) 
        for iterC in range(Ncell):
            for iterS in range(Nsample):
                varG[iterC, 0] = varG[iterC, 0] + \
                                 LA.norm(Rsamples[iterS, iterC, :, :] - Rmean[iterC,:,:])**2
            varG[iterC, 0] = varG[iterC, 0] / Nsample / LA.norm(Rmean[iterC,:,:])**2

        self.deltaField = varG**0.5

        ## Save the estimation of delta field.
        deltaFileName = caseName + '/0/delta'
        foamOp.writeScalarToFile(self.deltaField, deltaFileName)

 

    def evalPG(self, delta, G): 
        """
            Evaluate the PDF of matrix G
        """

        nDim = self.nDim
        pG = self._cG(nDim,delta) + math.log(LA.det(G) * (nDim + 1.0) * (1 - delta**2)) - \
             math.log(2 * delta**2) + \
             (-(nDim + 1.0) / (2 * delta**2) * np.trace(G))
        return pG


    def _cG(self, nDim, delta):
        """
            Calculate the coefficient cG
        """

        nDim = self.nDim
        denominator = 0
        for j in range(nDim):
            j = j + 1
            denominator = denominator + math.lgamma((nDim + 1.0) / \
                                        (2 * delta**2) + (1.0 - j) / 2.0) 

        # HX: double checked against my writing and literature; OK
        numerator = math.log((2 * math.pi)) * (-nDim * (nDim - 1.0) / 4.0) + \
                    math.log((nDim + 1.0) / (2 * delta**2)) * \
                    (nDim * (nDim + 1.0) * (2 * delta**2)**(-1))
        cG = numerator - denominator
        return cG


def Tau2G(Tau, TauBase):
    """
    Get matrix G from Tau and L
    """
    
    L = R2L(Tau2R(TauBase))
    LT = L.transpose()
    G = np.dot(np.dot(LA.inv(LT), Tau2R(Tau)), LA.inv(L))

    return G


def Tau2R(Tau):
    """
    Reshape Tau [6, 1] to matrix [3, 3]
    """

    R = np.zeros([3, 3])
    idxs = np.triu_indices_from(R)
    R[idxs] = Tau
    R[idxs[1], idxs[0]] = Tau

    return R


def R2L(R):
    """
    Perform Cholesky decomposition to matrix R
    """
    
    L = LA.cholesky(R)

    return L
        

if __name__ == '__main__':
    ## Test the searching of delta based on MLE
    #M = np.loadtxt("testFiles/cavity_centerlineV_Tau.xy")
    #TauDNS = M[:,3:]
    #genDelta = GenerateDelta(16)
    #genDelta.inferDelta(TauDNS, TauDNS, 200)

    ## Test the estimation of delta based on Tau samples
    ModelInput = './forwardModelInput.in'
    paramDict = readInputData(ModelInput)        
    # baseline case folder name
    caseName = paramDict['caseName']    
    nCell = int(paramDict['Ncell'])
    MainInput = './MainInput.in'
    paramDict = readInputData(MainInput)
    nSample = int(paramDict['Ns'])
    genDelta = GenerateDelta(nCell)
    genDelta.estimateDelta(caseName, nSample)

    ## Test the implementation of pG
    #M = np.loadtxt("testFiles/testTau.xy")
    #TauDNS = M[:,3:]
    #genDelta = GenerateDelta(16)
    #pG = genDelta.evalPG(0.3, Tau2R(TauDNS[0,:]))
    #print pG
