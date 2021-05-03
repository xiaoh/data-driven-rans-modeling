#!/aoe/bin/python27

# description        :Plot of tensor ensemble

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.03, 2015
# revision           :Dec.06, 2015
####################################################################################################

## Import system modules
# sci computing
import numpy as np
# system, file operation
import pdb
import os
# plotting
import matplotlib.pyplot as plt
# Import local modules
import postFuns as post


RPath = './resultData/R_samples'
GPath = './resultData/G_samples'
# Get Ensemble of R Fields and Ensemble of G Fields
RFields = post.readMatrixSample(RPath)
GFields = post.readMatrixSample(GPath)
# Get size of Fields
nSample, nCell, nTensor = GFields.shape
# Ensemble means of R Fields and G Fields, respectively
RFieldMean = np.mean(RFields, axis=0)
GFieldMean = np.mean(GFields, axis=0)
# further collpase GFieldMean to [1 by 6] -> [G11, G12, G13, G22, G23, G33]
GMean = np.mean(GFieldMean, axis=0)
# G truth
GTruth = [1, 0, 0, 1, 0, 1]

# plot the 6 components with truth
plt.figure(1)
plt.clf()
xlabels = ['G11', 'G12', 'G13', 'G22', 'G23', 'G33']
x = np.arange(6)
lm, = plt.plot(x, GMean, 'o', c='red', markersize=10)
lt, = plt.plot(x, GTruth, 'x', c='blue', markersize=10)
plt.xticks(x, xlabels)
plt.xlim(-0.2, 6.2)
plt.legend((lm, lt), ('sample mean', 'truth'), loc=7)
plt.xlabel('G components')
plt.ylabel('Value')
plt.title("Mean of G Samples, nSampe = "+str(nSample)+" nCell = "+str(nCell))
fname = "GMean_nSample_"+str(nSample)+"_nCell_"+str(nCell)+'.pdf'
plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)


# plot the 6 components with truth
plt.figure(2)
plt.clf()
xlabels = ['G11', 'G12', 'G13', 'G22', 'G23', 'G33']
x = np.arange(6)
lm, = plt.plot(x, abs(GMean-GTruth), 'o', c='red', markersize=10)
plt.xticks(x, xlabels)
plt.xlim(-0.2, 6.2)
plt.xlabel('G components')
plt.ylabel('L1 Error')
plt.title("Error(L1)-Mean of G Samples, nSampe = "+str(nSample)+" nCell = "+str(nCell))
fname = "ErrGMean_nSample_"+str(nSample)+"_nCell_"+str(nCell)+'.pdf'
plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)

