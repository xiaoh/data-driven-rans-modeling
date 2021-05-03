#!/aoe/bin/python27
import pdb
import numpy as np
import os
import sys
import KLReducedModel as kl
import ReynoldsStressRF as map
from sigmaFieldOperations import computeSigmaField
from utilities import readInputData, extractListFromDict

paramDict = readInputData('forwardModelInput.in')
kernelType = paramDict['kernelType']
m = int(paramDict['NmodeMax'])
KLRatio = float(paramDict['KLRatioXi'])
baseCaseName = paramDict['caseName']
rbfKernel = paramDict['rbfKernel']
rbfLengthScale = float(paramDict['rbfLengthScale'])
lenVec = extractListFromDict(paramDict, 'lenVecXi')
lenXi = np.array([float(pn) for pn in lenVec])
hyperPara = np.hstack((6, lenXi))
#pdb.set_trace()

baseDir = "./" + baseCaseName
# generate 2D meshing and cell areas
scatteredSigmaFile = baseDir + "/constant/scatSigma.dat"
meshCoords = np.loadtxt('./klExpansionDataXi3D/cellCenter3D.dat')
sigmaCoords = np.absolute(computeSigmaField(scatteredSigmaFile, meshCoords, \
                                       rbfKernel, rbfLengthScale))

areaCoords = np.loadtxt('./klExpansionDataXi3D/cellArea3D.dat')


# Initial KL class and call the solver
KL3D = kl.KLReducedModel(meshCoords, kernelType, hyperPara, 
                         m, KLRatio, sigmaCoords, areaCoords)

KLModes = np.loadtxt('./klExpansionDataXi3D/KLmodes.dat')
(ncell, nModes) = KLModes.shape
nModes = nModes - 2

XiorgM = np.loadtxt('./debugData/init/deltaXiM')
nSamples = 1000
omegaPrime = np.zeros([nModes, nSamples])
for i in np.arange(nSamples):
    num = float(i+1) 
    sampleCaseName = baseDir + '-tmp_' + str(num) + '/0/'
    mapTau = map.ReynoldsStressRF(sampleCaseName, 'Tau', ncell, 1)
    Xiorg = np.array([XiorgM[i, :]])
    Xiorg = Xiorg.T
    XiFac = mapTau.getXiFactor()
    XiFac = np.array([XiFac])
    XiFac = XiFac.T
    Xinew = Xiorg * XiFac
    omegaPrimeTemp = KL3D.projectFieldReduced(0, Xinew, None, KLModes)
    omegaPrime[:, i] = omegaPrimeTemp[:, 0]
    print "process sample :", i+1
np.savetxt("./omegaXiPrime", omegaPrime)
