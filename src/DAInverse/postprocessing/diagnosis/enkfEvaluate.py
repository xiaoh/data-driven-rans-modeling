#!/aoe/bin/python27

# description        :Evaluate algorithm of EnKS (EnKF)

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Mar.15, 2016
# revision           :Mar.15, 2016
####################################################################################################

## Import system modules
# sci computing
import numpy as np
# system, file operation
import pdb
import sys
import ast
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
## Local modules
from utilities import readInputData

def analysisCov_P():
    """
    analyze error covariance of state:
    :return:
    """
    PHT = np.loadtxt('./debugData/PHT_' + str(EnKFStep))
    HPHT = np.loadtxt('./debugData/HPHT_' + str(EnKFStep))
    K = np.loadtxt('./debugData/kalmanGain_' +  str(EnKFStep))
    pdb.set_trace()

if __name__ == '__main__':
    paramDict = readInputData('plotInfo.in')
    EnKFStep = float(paramDict['EnKFStep'])
    analysisCov_P()