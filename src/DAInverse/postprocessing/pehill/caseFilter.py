#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from pylab import *
from utilities import readInputData

paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir']
resCtrl = float(paramDict['resCtrl'])

def main():
    if not os.path.exists('divergedCases/'):
        os.mkdir('divergedCases/')
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, 'pehill-tmp*'):
            Mres = np.loadtxt(file+'/logs/Ux_0')
            res = Mres[-1,-1]
            if res >= resCtrl:
                os.system('mv '+file+' divergedCases/')
                print file, '\n'

if __name__ == "__main__":
   main()
