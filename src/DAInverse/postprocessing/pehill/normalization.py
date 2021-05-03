#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import fnmatch
from utilities import readInputData
import ast
import pdb

paramDict = readInputData('plotInfo.in')
UorgFlag = ast.literal_eval(paramDict['UorgFlag'])
doSampling = ast.literal_eval(paramDict['doSampling'])
timeDir = paramDict['timeDir'] + '.000000'
lastTimeDir = paramDict['lastTimeDir'] + '.000000'
xPosList = [0, 1, 2, 3, 4, 5, 6, 7, 8]
Ub = 0.028
pathroot = os.getcwd()

if doSampling:
    for folder in os.listdir('.'):
        if fnmatch.fnmatch(folder, 'pehill_base-*'):
            print "do sampling for " + folder
            os.chdir(folder)
            os.system('rm -rf 0/uniform')
            os.system('rm -rf 0/dz.dat')
            os.system('cp 0/U 0/Uorg')
            os.system('cp '+ lastTimeDir+'/U ' + lastTimeDir+'/Uorg')
            if UorgFlag:
                os.system('cp ../sampleDict-org ./system/sampleDict' )
                os.system('sample -time 0 > newsample.log')
                os.system('sample -time '+timeDir + '>> newsample.log')
            else:
                os.system('cp ../sampleDict-nonorg ./system/sampleDict')
                os.system('sample -time 0 >> newsample.log')
                os.system('sample -time '+timeDir + '>> newsample.log')        
            os.chdir(pathroot)

if UorgFlag:
    print "We are using propagated velocity Uorg"
else:
    print "We are using updated velocity U"
    
for xPos in xPosList:
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, 'pehill_base-*'):
            if fnmatch.fnmatch(file, 'pehill_base-observation'):
                filenames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_U.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_U.xy']            
            else:
                if UorgFlag:
                    filenames = [file + '/postProcessing/sets/0/line_x' + str(xPos) + '_Uorg.xy', \
                         file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_Uorg.xy']
                else:
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
