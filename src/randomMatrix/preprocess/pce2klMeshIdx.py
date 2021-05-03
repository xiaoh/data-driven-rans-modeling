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
# system, file operation
import pdb
import os
# plotting
import matplotlib.pyplot as plt  # for plotting
## Import local modules
import foamFileOperation as foamOp
from utilities import readInputData


print "We are getting the pce mesh idx of locations of cell centers of kl mesh"
print "Usage:\n Please make sure ccx, ccy and ccz are available in caseName_base-kl"
# read mainInput file
mainInputFile = './mainInput.in'
paramDict = readInputData(mainInputFile)
caseName = paramDict['caseName']
nCell_kl = int(paramDict['nCell_kl'])
nCell_pce = int(paramDict['nCell_pce'])
# folder direct
cfdDir = os.path.join(os.getcwd(), caseName)
klDir = cfdDir + '-kl'
pceDir = cfdDir + '-pce'

xCoord_kl = foamOp.readTurbCoordinateFromFile(klDir+'/0/')
obsLocFile_pce = pceDir + '/constant/obsLocations'
foamOp.writeLocToFile(xCoord_kl, obsLocFile_pce)
# determine the existence of 'indexHost.txt'
path_IdxH = pceDir + '/constant/indexHost.txt'
pdb.set_trace()
if os.path.isfile(path_IdxH):
    os.remove(path_IdxH)
os.system('getHostCellIdx -case ' + pceDir + ' &> log.pce2kl')


