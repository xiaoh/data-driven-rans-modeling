
import numpy as np
import foamFileOperation as utils

outfname = './pehill_base/constant/sigmaZero.dat'
xfname = './pehill_base/0/ccx'
xmin = 9

x = utils.readScalarFromFile(xfname)
np.savetxt(outfname, x>xmin, fmt='%i')


