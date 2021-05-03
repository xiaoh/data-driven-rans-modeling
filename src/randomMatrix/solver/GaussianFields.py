#!/usr/bin/env python

# description        :Generate gaussian random fields with zero mean and unit variance.

# author             :Jinlong Wu (jinlong@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.01, 2015
# revision           :Dec.01, 2015

####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.sparse as sp
import scipy.stats as ss
from scipy.stats import norm
import numpy.linalg as LA
# system, file operation
import pdb
import time
# plotting
#import seaborn as sns  # for statistical plotting
import matplotlib.pyplot as plt  # for plotting

## Import local modules
import StochasticProcess as mSP
import KLExpansion as KL


class GaussianF:
    """
    Generate nSample random fields with zero mean and unit gaussian distribution.

    :arg:
        nCell           : number of the cell points in one field
        xState          : spatial coordinate (More details in StochasticProcess)
        Arg_covGen      : dictionary to specify covariance matrix
        Arg_calKLModes  : dictionary to control KL expansion

    """

    def __init__(self, nCell, xState, Arg_covGen, Arg_calKLModes):
        self.nCell = nCell
        ## Set zero mean for random field construction
        self.meanField = np.zeros([self.nCell, 1])
        ## Initial a instance of GaussianProcess class
        self.gp = mSP.GaussianProcess(xState) 
        cov_sparse, covWeighted_sparse = self.gp.covGen(Arg_covGen)
        ## Initialize the KLExpansion class
        self.kl = KL.klExpansion(cov_sparse, covWeighted_sparse)
        ## Get KLModes
        [self.eigVal, self.KLModes] = self.kl.calKLModes(Arg_calKLModes)

    def sample(self, nSample):
        """
        Generate nSample gaussian random fields

        :arg:
            nSample:         number of samples of GP field

        :return:
            GaussianFields:  an ensemble of Gaussian Random fields [nSample by nCell by 1]
        """
        ## span storage space for GausianField
        GaussianFields = np.zeros((nSample, self.nCell, 1))
        [N, nKL] = self.KLModes.shape  # parse the dimension
        for iter in range(nSample):
            omegaVec = ss.norm.rvs(size=nKL); omegaVec = np.array([omegaVec]).T
            assert len(omegaVec) == nKL, \
            "Lengths of KL coefficient omega (%d) and KL modes (%d) differs!"% (len(omegaVec), nKL)
            GaussianFields[iter,:,:] = self.meanField + np.dot(self.KLModes, omegaVec)
        return GaussianFields

    def plotSinglePointPDF(self, GaussianFields, pIndex):
        """
        Plot PDF for a single point of a Gaussian Fields

        :arg:
            GaussianFields:  an ensemble of Gaussian Random fields [nSample by nCell by 1]
            pIndex:          index of point you want to study
        :return:
            None
        """
        fig = plt.figure()
        count, bins, ignored = plt.hist(GaussianFields[:,pIndex,:], 50, normed=1)
        plt.plot(bins, norm.pdf(bins,0,1),'r-', lw=5, alpha=0.6, label='gaussian pdf')
        plt.title('KLmodes ' + str(self.kl.nKL+1) + '; Point ' + str(pIndex))
        plt.tight_layout()
        fig.savefig('KLmodes_' + str(self.kl.nKL+1) + '-Point_' + str(pIndex) + '.pdf')
        print "Cell point test, gaussian distribution with zero mean and unit variance is expected"
        plt.show()




if __name__ == '__main__':

    # directory where the test data stored
    testDir = '../../src/RandomField/verificationData/klExpansion/mesh10/' 

    # Initialization 
    nSample = 10000
    nCell = 10
    xState = np.loadtxt(testDir + 'cellCenter3D.dat')
    lenXField = 0.5*np.ones(nCell)
    lenYField = 0.5*np.ones(nCell)
    lenZField = 0.5*np.ones(nCell)
    weightField = np.loadtxt(testDir + 'cellArea3D.dat')
    truncateTol = -np.log(1e-10)

    Arg_covGen = {
                    'lenXField': lenXField,
                    'lenYField': lenYField,
                    'lenZField': lenZField,
                    'weightField':weightField,
                    'truncateTol': truncateTol
                 }

    nKL = 9
    Arg_calKLModes = {
                        'nKL': nKL,
                        'weightField':weightField
                    }

    ## Generate gaussian random fields
    GaussianRF = GaussianF(nCell, xState, Arg_covGen, Arg_calKLModes)

    ## Test the distribtuion for cell point
    GaussianFields = GaussianRF.sample(nSample)

    GaussianRF.plotSinglePointPDF(GaussianFields, 0)
