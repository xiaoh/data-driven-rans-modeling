#!/usr/bin/env python

# description        :Generate Gamma fields based on Gaussian fields and PCE expansion.

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.03, 2015
# revision           :Dec.03, 2015

####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.sparse as sp
import scipy.stats as ss
import numpy.polynomial.polynomial as np3
from scipy.stats import norm
import numpy.linalg as LA
# system, file operation
import scipy.io
import pdb
import time
import os
# plotting
import matplotlib.pyplot as plt  # for plotting
## Import local modules
import GaussianFields as GP


class GammaF:
    """
    Generate Gamma fields based on Gaussian fields and PCE expansion.

    :arg:
        nCell:           number of the cell points in one field
        p_order:         Order of polynomials for PCE

    """

    def __init__(self, nCell_kl, nCell_pce, p_order, gammaMean, gammaMeanName, hostCellIdx=False):
        #TODO: nCell should be separated to nCell_pce and nCell_kl
        self.nCell_kl = nCell_kl
        self.nCell_pce = nCell_pce
        self.p_order = p_order
        self.gammaMean = gammaMean
        self.gammaMeanName = gammaMeanName
        self.hostCellIdx = hostCellIdx
        self.u_k = np.zeros([nCell_kl, p_order+1])
        ## find this source code path
        self.thisFilePath = os.path.dirname(os.path.realpath(__file__))
        ## save the parameters as m files for PCE MATLAB solver
        #gammaParaPath = self.thisFilePath+'/pce_matlab/para_gamma-'+gammaMeanName+'.mat'
        gammaParaPath = os.getcwd() +'/para_gamma-'+gammaMeanName+'.mat'
        scipy.io.savemat(gammaParaPath, {'mean_gamma':gammaMean, 'nCell': nCell_pce,
                                    'p_order': float(p_order), 'mean_gammaName': gammaMeanName})
        ## get basis and coefficient for PCE
        #self.u_k, self.P, self.highest_order, self.num_polys = self._PCECoefficients()

    def getPCECoefficient(self):
        """
        To get the PCE coefficients, which must be called before calling self.sample
        :return:
        """
        ## get basis and coefficient for PCE
        u_k, self.P, self.highest_order, self.num_polys = self._PCECoefficients()
        if self.hostCellIdx is not False:
            print "Mapping the u_k in PCE mesh to u_k in KL mesh"
            self._mapuk(u_k)
        else:
            self.u_k = u_k


    def _mapuk(self, u_k):
        """
        :param hostCellIdx:
        :return:
        """
        cellidx_kl = self.hostCellIdx[:, 0]
        cellidx_pce = self.hostCellIdx[:, 1]
        self.u_k[cellidx_kl.tolist(), :] = u_k[cellidx_pce.tolist(), :]

    def sample(self, nSample, GaussianFields):
        """
        draw nSample realizations (samples) of the Gamma field based on PCE

        :arg:
            nSample:        Number of realizations (samples) of the Gamma field
            GaussianFields: nSample of Gaussian fields [nSample, nCell]
        :return:
            GammaFields:    nSample of Gamma fields
        """
        GammaFields = np.zeros((nSample, self.nCell_kl, 1))
        for idx in np.arange(self.nCell_kl):
            for k in np.arange(self.P):
                GammaFields[:, idx, :] = GammaFields[:, idx, :] + self.u_k[idx, k] * \
                                         np3.polyval(GaussianFields[:, idx, :],
                                                     self.num_polys[k, 0:self.highest_order[k]+1]);

        return GammaFields


    def plotSinglePointPDF(self, GammaFields, pIndex, gammak_t):
        """
        Plot PDF for a single point of a Gamma Fields

        :arg:
            GammaFields:     an ensemble of Gamma Random fields [nSample by nCell by 1]
            pIndex:          index of point you want to study
            gammak_t:        analytical mean of Gamma
        :return:
            None
        """
        #pdb.set_trace()
        plt.figure()
        #pdb.set_trace()
        count, bins, ignored = plt.hist(GammaFields[:,pIndex,:], 100, normed=True)
        #pdb.set_trace()
        plt.plot(bins, ss.gamma.pdf(bins, gammak_t),'r-', lw=5, alpha=0.6, label='analytical Gamma pdf')
        plt.title('Point ' + str(pIndex))
        plt.tight_layout()
        #fig.savefig('Point_' + str(pIndex) + '.pdf')
        plt.show()



    def _PCECoefficients(self):
        """
        calculate PCE coefficients and PCE basis
        :arg:

        :return:
        """

        pathcwd = os.getcwd()
        pathofPCE_field = self.thisFilePath + '/pce_matlab/'
        #cmdString = 'echo \"PCE_field([\'' + self.gammaMeanName + '\'])\" ' + ' | matlab -nojvm -nodisplay -nosplash'
        cmdString = 'echo \"PCE_field([\'' + pathcwd + '\'], [\'' + self.gammaMeanName + '\'])\" ' + \
                    ' | matlab -nojvm -nodisplay -nosplash'
        #pdb.set_trace()
        # go to the PCE matlab src directory
        os.chdir(pathofPCE_field)
        os.system(cmdString)
        #pdb.set_trace()
        #os.system('matlab -nojvm -nodisplay -nosplash < PCE_field.m')
        # go back to the case directory
        os.chdir(pathcwd)
        poly_contents = scipy.io.loadmat('poly-' + self.gammaMeanName + '.mat')
        # reduce dimension (python import everything from matlab as 2D array)
        u_k = poly_contents['u_k']
        P = poly_contents['P'][0,0]
        # Account for zero-based indexing in python;
        highest_order = np.arange(P+1)
        # Note: python arange from low to high order
        # p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3
        # Matlab is the opposite; Also note one-based indexing
        # p[1] * x**3 + p[2] * x**2 + p[3] * x + p[4]
        # Revert coefficient arrangements
        num_polys = poly_contents['num_polys'][:, -1: : -1]

        return u_k, P, highest_order, num_polys

def calculateCovariance(fields, idx1, idx2):
    mean1 = np.mean(fields[:, idx1, 0])
    mean2 = np.mean(fields[:, idx2, 0])
    std1 = np.std(fields[:, idx1, 0])
    std2 = np.std(fields[:, idx2, 0])
    covar12 = np.mean((fields[:, idx1, 0] - mean1) * (fields[:, idx2, 0] - mean2))
    correlation12 = covar12/std1/std2 
    return correlation12    
     
if __name__ == '__main__':

    # directory where the test data stored
    #testDir = '../../RandomField/verificationData/klExpansion/mesh10/'
    testDir = './testMesh5/'
    # Initialization
    nSample = 1000
    nCell_kl = 6
    nCell_pce = 6
    xState = np.loadtxt(testDir + 'cellCenter3D.dat')
    lenXField = 0.2*np.ones(nCell_kl)
    lenYField = 0.2*np.ones(nCell_kl)
    lenZField = 0.2*np.ones(nCell_kl)
    weightField = np.loadtxt(testDir + 'cellArea3D.dat')
    truncateTol = -np.log(1e-10)

    Arg_covGen = {
                    'lenXField': lenXField,
                    'lenYField': lenYField,
                    'lenZField': lenZField,
                    'weightField':weightField,
                    'truncateTol': truncateTol
                 }

    nKL = 3
    Arg_calKLModes = {
                        'nKL': nKL,
                        'weightField':weightField
                    }

    ## Generate gaussian random fields
    GaussianRF = GP.GaussianF(nCell_kl, xState, Arg_covGen, Arg_calKLModes)

    ## Test the distribtuion for cell point
    GaussianFields = GaussianRF.sample(nSample)

    ## Plotting Gaussian field
    #GaussianRF.plotSinglePointPDF(GaussianFields, 0)

    ## initialize Gamma field class
    p_order = 3
    delta_x = np.zeros([nCell_kl, 1]) + 0.5;
    #gammaMean = 2/(delta_x**2)
    gammaMean = 3*np.ones([nCell_kl, 1])
    GammaRF = GammaF(nCell_kl, nCell_pce, p_order, gammaMean, 'xx')
    GammaRF.getPCECoefficient()
    GammaFields = GammaRF.sample(nSample, GaussianFields)
    
    # Test the covariance of Gaussian and Gamma
    idx1 = 0; idx2 = 2
    corrGaussian = calculateCovariance(GaussianFields, idx1, idx2)
    corrGamma = calculateCovariance(GammaFields, idx1, idx2)
    print "covariance of Gaussian is ",  corrGaussian
    print "covariance of Gamma is ",  corrGamma
    pdb.set_trace()
    #pdb.set_trace()
    GammaRF.plotSinglePointPDF(GammaFields, 0, 3.0)
