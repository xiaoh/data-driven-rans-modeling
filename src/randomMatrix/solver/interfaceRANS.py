#!/usr/bin/env python

# description        :To generated a group of Reynolds stress samples based on RMF
#                     and propagate this Reynolds stress ensemble to velocity.

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
import os
import os.path as ospt
import shutil
import pdb
import time
# plotting
import matplotlib.pyplot as plt  # for plotting
## Import local modules
import foamFileOperation as foamOp
from ReynoldsStressRF import ReynoldsStressRF
import utilities as Tool

class interfaceRANS:
    """
    To generated a group of Reynolds stress samples based on RMF,
    and propagate this Reynolds stress ensemble to velocity.

    Arg:
        nSample         : Number of required fields

    """

    def __init__(self, nSample, nCell_cfd):
        self.nSample = nSample
        self.nCell_cfd = nCell_cfd


    def getRbar(self, caseDir, tauFile):
        """

        :param dir:
        :return:
        """
        mapTau = ReynoldsStressRF(caseDir, tauFile, self.nCell_cfd, 1, 'True')
        #pdb.set_trace()
        # This suppose to be realizable Tau
        Tau = mapTau.tau
        # This is unrealizable tau
        #Tau = foamOp.readTurbStressFromFile(caseDir+tauFile)
        Rbar = self.tau2rm(Tau, self.nCell_cfd)
        return Rbar, Tau

    def tau2rm(self, tau, nCell):
        """
        reshape tau [nSample by 6] to random matrix [nSample, 3, 3]
        :param tau:
        :return:
        """
        Rbar = np.zeros([nCell, 3, 3])
        # TODO: this for loop needed to be replaced for efficiency
        for idx in np.arange(nCell):
            Rbar_idx = Rbar[idx, :, :]
            idxs = np.triu_indices_from(Rbar_idx)
            ## fill upper triangle (half symm)
            Rbar_idx[idxs] = tau[idx, :]
            ## fill full (full symm)
            Rbar_idx[idxs[1], idxs[0]] = tau[idx, :]
            ## fill into array of Rbar
            Rbar[idx, :, :] = Rbar_idx

        return Rbar

    def rm2tau(self, rm, nCell):
        """
        reshape random matrix [nSample, 3, 3] to tau
        :param rm:
        :return:

        """
        tau = np.zeros([nCell, 6])
        for icell in np.arange(nCell):
            rm_i = rm[icell, :, :]
            idxs = np.triu_indices_from(rm_i)
            tau[icell, :] = rm_i[idxs]
        return tau

    def genFolders(self, nSample, caseName, caseNameRun, TauSamples, writeInterval):
        """
        generate folder of different tau
        :param rm:
        :return:

        """
        # remove previous ensemble case files
        os.system('rm -fr '+caseName+'-tmp_*')
        os.system('rm -fr '+caseName+'-run')
        ii = 0
        caseCount = np.linspace(1, nSample, nSample)
        for case in caseCount:

            print "\n#", case, "/", nSample, " Creating folder for Case = ", case

            tmpCaseName = caseName + "-tmp_" + str(case)

            if(ospt.isdir(tmpCaseName)): #see if tmpCaseName's'directory is existed
                shutil.rmtree(tmpCaseName)
            shutil.copytree(caseName, tmpCaseName) # copy


            # Replace Tau ensemble for cases ensemble
            tauTemp = TauSamples[ii, :, :]
            tauFile = './' + tmpCaseName + '/0/Tau'
            foamOp.writeTurbStressToFile(tauTemp, tauFile)
            rasFile = ospt.join(os.getcwd(), tmpCaseName, "system", "controlDict")
            Tool.replace(rasFile, "<writeInterval>", str(writeInterval))
            ii += 1

        #generate observation folder
        if(ospt.isdir(caseNameRun)):
            shutil.rmtree(caseNameRun)
        shutil.copytree(caseName, caseNameRun) # copy
        rasFile = ospt.join(os.getcwd(), caseNameRun, "system", "controlDict")
        Tool.replace(rasFile, "<writeInterval>", str(writeInterval))





