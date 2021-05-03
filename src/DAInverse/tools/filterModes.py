#!/aoe/bin/python27

import numpy as np
import math
import os
import pdb

dirLists = ['klExpansionDataXi3D', 'klExpansionDataEta3D', 'klExpansionDatak3D']
paraNames = ['Xi', 'Eta', 'k']

def main():
    index = 0
    for dir in dirLists:
        filenameMode = dir + '/KLmodes.dat'
        filenameEig = dir + '/eig.dat'
        os.system('cp ' + filenameMode + ' ' + dir + '/KLmodes.org')
        os.system('cp ' + filenameEig + ' ' + dir + '/eig.org')
        Mmode = np.loadtxt(filenameMode)
        Meig = np.loadtxt(filenameEig)
        MmodeFilter = Mmode[:,0:2]
        MeigFilter = np.array([])
        for col in range(Mmode.shape[1]-2):
            mode = Mmode[:,col+2]
            eig = Meig[col]
            if np.abs(np.sum(mode)) > 1e-5:
                MmodeFilter = np.vstack((MmodeFilter.T, mode)).T
                MeigFilter = np.append(MeigFilter, eig)
        np.savetxt(dir+'/KLmodes.dat', MmodeFilter)
        np.savetxt(dir+'/eig.dat', MeigFilter)
        print "The number of current modes for", paraNames[index], "is", len(MeigFilter), "\n"
        index = index + 1

if __name__ == "__main__":
    main()

