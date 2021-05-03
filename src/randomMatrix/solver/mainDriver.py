#!/usr/bin/env python

# description        :Main drive of model-form uncertainty propagation by using Random matrix and maximum entropy theory

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
# plotting
import matplotlib.pyplot as plt  # for plotting
## Import local modules
from MaxEntSolver import MaxEntSolver
import KLExpansion as KL

np.random.seed(1000)
# initial MaxEntSolver class
maxEntSol = MaxEntSolver('mainInput.in')
tic = time.time()
maxEntSol.writeKLModesToOF()
maxEntSol.sampleTau()
elapsed = time.time() - tic
print("elapsed time for sample N random Tau = ", str(elapsed))
tic = time.time()
maxEntSol.mappingTau2Phy()
elapsed = time.time() - tic
print("elapsed time for mapping RMT to physical = ", str(elapsed))
#maxEntSol.plotTau1LocSamplesOnBay(4)
tic = time.time()
if maxEntSol.propagateFlag:
    maxEntSol.propagateTau()
elapsed = time.time() - tic
print("elapsed time for propagation = ", str(elapsed))
