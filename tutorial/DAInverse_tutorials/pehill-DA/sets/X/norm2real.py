#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import fnmatch

xPosList = [0, 1, 2, 3, 4, 5, 6, 7, 8]
Ub = 0.028

for xPos in xPosList:
     filenames = ['x' + str(xPos) + '.0_U2.xy']
     normalizedFileNames = ['x' + str(xPos) + '.0_Ureal.xy']

     # Normalize 
     for i in range(len(filenames)):
         filename = filenames[i]
         M = np.loadtxt(filename)
         M[:,1] = M[:,1] * Ub
         M[:,2] = M[:,2] * Ub
         M[:,3] = M[:,3] * Ub
         np.savetxt(normalizedFileNames[i],M)
