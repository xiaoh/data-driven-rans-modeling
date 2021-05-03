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

RComponentFolder = './resultData/RComponent_samples/'
deltaRFolder = './resultData/deltaRComponent_samples/'
VsName = ['VA_s', 'VB_s', 'VC_s']
VbarName = ['VA_bar', 'VB_bar', 'VC_bar']
deltaVsName = ['deltaVA_s', 'deltaVB_s', 'deltaVC_s']

VAs = np.loadtxt(RComponentFolder+VsName[0])
VBs = np.loadtxt(RComponentFolder+VsName[1])
VCs = np.loadtxt(RComponentFolder+VsName[2])

VAbar = np.loadtxt(RComponentFolder+VbarName[0])
VBbar = np.loadtxt(RComponentFolder+VbarName[1])
VCbar = np.loadtxt(RComponentFolder+VbarName[2])

deltaVAs = np.loadtxt(deltaRFolder+deltaVsName[0])
deltaVBs = np.loadtxt(deltaRFolder+deltaVsName[1])
deltaVCs = np.loadtxt(deltaRFolder+deltaVsName[2])

pdb.set_trace()
