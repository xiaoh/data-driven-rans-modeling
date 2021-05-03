#!/usr/bin/env python

# description        :Generate random matrix field based on Gaussian and Gamma fields.

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
from scipy.stats import norm
import numpy.linalg as LA
# system, file operation
import pdb
import time
try:
    import pp
except ImportError, e:
    print e, 'Warning: no pp module, I only can use one core to do the calculation'
    hasPP = False
    pass
# plotting
import matplotlib.pyplot as plt  # for plotting
## Import local modules
import GaussianFields as Gau
import GammaFields as Gam

class RMF:
    """
    Generate random matrix field based on Gaussian and Gamma fields.

    Arg:
        Nsample         : Number of required fields

    """

    def __init__(self, nCell_kl, nCell_pce, Arg_covGen, Arg_calKLModes, deltaField, xCoord, hostCellIdx):
        self.d = 3  # generally dimension of random matrix is 3 for Reynolds stress problem
        self.nCell_kl = nCell_kl
        self.nCell_pce = nCell_pce
        self.hostCellIdx = hostCellIdx.astype(int)
        # TODO: deltaField need to be generated based on a comprehensive design
        # Note: deltaField must be [nCell_pce by 1]
        self.deltaField = deltaField
        self.sigmadField = self._getSigmad()
        self.GaussianRF = Gau.GaussianF(nCell_kl, xCoord, Arg_covGen, Arg_calKLModes)
        # Not work for PCE, temporarily shut down the PP
        self.hasPP = False
        if self.hasPP:
            self.job_server = pp.Server()
            self.job_server.set_ncpus(2)
            self.jobs = []

        #pdb.set_trace()

    def sample(self, nSample, p_order):
        """
        draw nSample realizations (samples) of the random matrix field [R]

        :arg:
            nSample:    Number of realizations (samples) of the random matrix
        :return:
            GFields:
        """
        self.nSample = nSample
        L11Fields, L22Fields, L33Fields = self._LDiagonal(p_order)
        L12Fields, L13Fields, L23Fields = self._LOffDiagonal()
        LFields = self._constructL(L11Fields, L22Fields, L33Fields, L12Fields, L13Fields, L23Fields)
        GFields = self._constructG(LFields)
        return GFields

    def _LDiagonal(self, p_order):
        """
        generate 3 Gamma field ensembles (nSamples samples), each is an ensemble of L_{a, a} field.


        :return:
            L11Fields: [nSample by nCell by 1]
            L22Fields: [nSample by nCell by 1]
            L33Fields: [nSample by nCell by 1]
        """

        # ## generate 3 ensembles of independent GP fields
        # GaussianFields_1 = self.GaussianRF.sample(self.nSample)
        # GaussianFields_2 = self.GaussianRF.sample(self.nSample)
        # GaussianFields_3 = self.GaussianRF.sample(self.nSample)
        # ## generate Gamma mean field
        # kField_1 = self._getGammaMean(1)
        # kField_2 = self._getGammaMean(2)
        # kField_3 = self._getGammaMean(3)
        # #TODO: nCell_kl need to input into GammaF
        # ## initialize Gamma PCE class
        # GammaRF_1 = Gam.GammaF(self.nCell_pce, p_order, kField_1)
        # GammaRF_2 = Gam.GammaF(self.nCell_pce, p_order, kField_2)
        # GammaRF_3 = Gam.GammaF(self.nCell_pce, p_order, kField_3)
        # ## get 3 Gamma fields
        # u1Fields = GammaRF_1.sample(self.nSample, GaussianFields_1)
        # u2Fields = GammaRF_2.sample(self.nSample, GaussianFields_2)
        # u3Fields = GammaRF_3.sample(self.nSample, GaussianFields_3)
        # ## get diagonal entry fields
        # L11Fields = self.sigmadField * (2*u1Fields)**0.5
        # L22Fields = self.sigmadField * (2*u2Fields)**0.5
        # L33Fields = self.sigmadField * (2*u3Fields)**0.5
        gammaMeanName = ['xx', 'yy', 'zz']
        ## generate 3 ensembles of independent GP fields
        GaussianFields = np.zeros([3, self.nSample, self.nCell_kl, 1])
        ## generate Gamma mean field
        kFields = np.zeros([3, self.nCell_pce, 1])
        ## initialize Gamma PCE class
        GammaRFs = np.ndarray((3,), dtype=np.object)
        ## get 3 Gamma fields
        uFields = np.zeros([3, self.nSample, self.nCell_kl, 1])
        for i in np.arange(3):

            GaussianFields[i, :, :, :] = self.GaussianRF.sample(self.nSample)
            kFields[i, :, :] = self._getGammaMean(i+1)
            GammaRFs[i] = Gam.GammaF(self.nCell_kl, self.nCell_pce, p_order, kFields[i, :, :], gammaMeanName[i], self.hostCellIdx)
        for i in np.arange(3):
            if(self.hasPP):
                self.jobs.append(
                    self.job_server.submit(GammaRFs[i].getPCECoefficient, (), modules = ("scipy.io", "os", ))
                                )
            else:
                GammaRFs[i].getPCECoefficient()

        #Barrier to make sure all cases are finished before moving on
        if(self.hasPP):
            part_sum1 = [job() for job in self.jobs]
            self.jobs = []
        part_sum1 = 0
        # finished all ensemble propagation. print stats
        if(self.hasPP):
            self.job_server.print_stats()

        for i in np.arange(3):
            uFields[i, :, :, :] = GammaRFs[i].sample(self.nSample, GaussianFields[i, :, :, :])

        ## get diagonal entry fields
        L11Fields = self.sigmadField * (2*uFields[0, :, :, :])**0.5
        L22Fields = self.sigmadField * (2*uFields[1, :, :, :])**0.5
        L33Fields = self.sigmadField * (2*uFields[2, :, :, :])**0.5

        return L11Fields, L22Fields, L33Fields

    def _LOffDiagonal(self):
        """
        generate 3 Gaussian field ensembles (nSamples samples), each is an ensemble of L_{a, b} field

        :return:
            L12Fields: [nSample by nCell]
            L13Fields: [nSample by nCell]
            L23Fields: [nSample by nCell]
        """
        ## generate 3 ensembles of independent GP fields
        GaussianFields_1 = self.GaussianRF.sample(self.nSample)
        GaussianFields_2 = self.GaussianRF.sample(self.nSample)
        GaussianFields_3 = self.GaussianRF.sample(self.nSample)
        ## generate sigma_d field
        L12Fields = self.sigmadField * GaussianFields_1
        L13Fields = self.sigmadField * GaussianFields_2
        L23Fields = self.sigmadField * GaussianFields_3

        return L12Fields, L13Fields, L23Fields

    def _mapSigmad(self, sigmadField_pce):
        """
        :param hostCellIdx:
        :return:
        """
        sigmadField_kl = np.zeros([self.nCell_kl, 1])
        cellidx_kl = self.hostCellIdx[:, 0]
        cellidx_pce = self.hostCellIdx[:, 1]
        sigmadField_kl[cellidx_kl.tolist(), :] = sigmadField_pce[cellidx_pce.tolist(), :]

        return sigmadField_kl

    def _getGammaMean(self, alpha):
        """
        get Gamma Mean Field (k(x) = 0.5*(d+1)/(delta(x)**2) + 0.5(1-alpha))

        :arg:
            alpha:          only be 1, 2 or 3, represent index of component
        :return:
            kField:         Gamma mean field (depending on delta(x) field)
        """
        kField = 0.5*(self.d + 1)/(self.deltaField**2) + 0.5*(1 - alpha)
        return kField

    def _getSigmad(self):
        """
        get sigmad field
        :return:
        """
        sigmadField_pce = self.deltaField * (self.d + 1)**(-0.5)
        sigmadField = self._mapSigmad(sigmadField_pce)

        return sigmadField

    def _constructL(self, L11Fields, L22Fields, L33Fields, L12Fields, L13Fields, L23Fields):
        """
        construct ensemble of upper triangular matrix field LFields

        :return:
            LFields: [nSample by (nCell by 3 by 3)] upper triangular matrix field
        """
        LFields = np.zeros([self.nSample, self.nCell_kl, self.d, self.d])
        LFields[:, :, 0, 0] = L11Fields[:, :, 0]
        LFields[:, :, 1, 1] = L22Fields[:, :, 0]
        LFields[:, :, 2, 2] = L33Fields[:, :, 0]
        LFields[:, :, 0, 1] = L12Fields[:, :, 0]
        LFields[:, :, 0, 2] = L13Fields[:, :, 0]
        LFields[:, :, 1, 2] = L23Fields[:, :, 0]

        return LFields

    def _constructG(self, LFields):
        """

        :param LFields:
        :return:
        """
        LFieldsT = np.transpose(LFields, (0, 1, 3, 2))
        GFields = np.zeros([self.nSample, self.nCell_kl, self.d, self.d])
        #TODO: this is not correct, how to do element-wise dot product need to be considered
        for isamp in np.arange(self.nSample):
            for icell in np.arange(self.nCell_kl):
                GFields[isamp, icell, :, :] = np.dot(LFieldsT[isamp, icell, :, :], LFields[isamp, icell, :, :])

        return GFields

