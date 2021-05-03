#!/aoe/bin/python27

# description        :Plot of Scalar ensemble

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


deltaFolder = './resultData/deltaRComponent_samples/'
fileNameCommon = 'delta'
deltaRComponentPool = ['Xi', 'Eta', 'K', 'VA', 'VB', 'VC']
startCellIdx = 1000
endCellIdx = 1020
# loop 6 component of Reynolds stress discrepancy
i = 0
for dR_c in deltaRComponentPool:
    print 'processing the deltaR component: ' + dR_c
    #pdb.set_trace()
    scalarFields_temp = np.loadtxt(deltaFolder+fileNameCommon+dR_c+'_s')
    #pdb.set_trace()
    nSample, nCell = scalarFields_temp.shape
    meanField, varField, stdField = post.statEva_scalarField(scalarFields_temp)
    x = np.arange(nCell)
    x = x[startCellIdx:endCellIdx]
    # plot the mean field with std
    plt.figure(i+1)
    plt.clf()
    transparency = 0.2
    for s in range(nSample):
        asample, = plt.plot(x, scalarFields_temp[s, startCellIdx:endCellIdx],'o',color='#1b9e76',alpha=transparency)
    am, = plt.plot(x, meanField[startCellIdx:endCellIdx], 'o', c='red', markersize=10)
    #pdb.set_trace()
    # plt.errorbar(x, meanField[startCellIdx:endCellIdx], yerr=stdField, fmt=None, capsize=8,
    #                           elinewidth=2, zorder=2, ecolor='blue', label='sigma')
    plt.legend([am, asample],["Sample Mean", "Samples"],prop={'size':10},numpoints=1)
    xlabels = np.array([str(pn) for pn in x+1])
    plt.xticks(x, xlabels)
    plt.title('component for R discrepancy : '+dR_c)
    plt.xlabel('mesh cell index')
    plt.ylabel('Value')
    fname = 'deltaR_'+dR_c+'.pdf'
    plt.savefig(fname)






