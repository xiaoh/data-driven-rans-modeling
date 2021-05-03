#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import fnmatch
from utilities import readInputData

paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir']+'.000000'
xPosList = [0, 1, 2, 3, 4, 5, 6, 7, 8]
Ub = 0.028

for xPos in xPosList:
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, 'pehill_base-*'):
            filenames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_U.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_U.xy']
            normalizedFileNames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_UNorm.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_UNorm.xy']

            # Normalize 
            for i in range(len(filenames)):
                filename = filenames[i]
                M = np.loadtxt(filename)
                M[:,1] = M[:,1] / Ub
                M[:,2] = M[:,2] / Ub
                M[:,3] = M[:,3] / Ub
                np.savetxt(normalizedFileNames[i],M)
    filenames = ['pehill_truth' + '/postProcessing/sets/0/line_x' + str(xPos) + '_U.xy']
    normalizedFileNames = ['pehill_truth' + '/postProcessing/sets/0/line_x' + str(xPos) + '_UNorm.xy']

    # Normalize 
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        M[:,1] = M[:,1] / Ub
        M[:,2] = M[:,2] / Ub
        M[:,3] = M[:,3] / Ub
        np.savetxt(normalizedFileNames[i],M)

for xPos in xPosList:
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, 'pehill_base-*'):
            filenames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_Tau.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_Tau.xy']
            normalizedFileNames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_TauNorm.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_TauNorm.xy']
            # Normalize 
            for i in range(len(filenames)):
                filename = filenames[i]
                M = np.loadtxt(filename)
                M[:,1] = M[:,1] / Ub / Ub
                M[:,2] = M[:,2] / Ub / Ub
                M[:,3] = M[:,3] / Ub / Ub
                M[:,4] = M[:,4] / Ub / Ub
                M[:,5] = M[:,5] / Ub / Ub
                M[:,6] = M[:,6] / Ub / Ub
                np.savetxt(normalizedFileNames[i],M)
    filenames = ['pehill_truth' + '/postProcessing/sets/0/line_x' + str(xPos) + '_Tau.xy']
    normalizedFileNames = ['pehill_truth' + '/postProcessing/sets/0/line_x' + str(xPos) + '_TauNorm.xy']
    # Normalize 
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        M[:,1] = M[:,1] / Ub / Ub
        M[:,2] = M[:,2] / Ub / Ub
        M[:,3] = M[:,3] / Ub / Ub
        M[:,4] = M[:,4] / Ub / Ub
        M[:,5] = M[:,5] / Ub / Ub
        M[:,6] = M[:,6] / Ub / Ub
        np.savetxt(normalizedFileNames[i],M)
