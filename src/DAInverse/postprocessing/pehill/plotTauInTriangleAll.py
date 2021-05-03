#!/usr/bin/python
import sys, os
import os.path as osp
import numpy as np
import matplotlib.pyplot as plt
import pdb
from matplotlib.patches import Rectangle
import fnmatch
from functions import *
from foamFileOperation import *
from ReynoldsStressRF import *
from utilities import readInputData

paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir'] + '.000000'
mfuDir = os.getcwd()

xPosList = [0, 1, 2, 3, 4, 5, 6, 7, 8]
samplePlotList = range(1,101)

legendList = ['baseline', 'samples', 'DNS']

if not os.path.exists('figures/'):
    os.makedirs('figures/')

for numHill in range(0,1):
    for xPos in xPosList:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_'+str(sample)+'.0*'):
                        plotFlag = True
                if plotFlag == True:
                    Mbase = np.loadtxt('pehill-observation/postProcessing/sets/0/line_x' + str(xPos) + '_Tau.xy');
                    Msample = np.loadtxt(file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_Tau.xy');
                    Mdns = np.loadtxt('/l/jinlong/pehill-DA/dataDNS/sets/X/ransTau' + str(xPos));
                    
                    # Get Reynolds Stress Data
                    tauBase = Mbase[:,1:7] 
                    tauSample = Msample[:,1:7] 
                    tauDNS = Mdns[:,1:7]
                    
                    # Initialize R0
                    R0 = ReynoldsStressRF('None', tauBase, tauBase.shape[0],1);
                    
                    p1, p2, p3 = R0.plotTauOnBaycentric(tauBase,tauDNS,tauSample, 'o')
    lg = plt.legend([p1,p3,p2], legendList, loc = 0)
    lg.draw_frame(False)
    plt.savefig("figures/perturbCompAll.pdf")
    plt.clf()
