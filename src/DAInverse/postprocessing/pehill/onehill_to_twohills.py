
import numpy as np 
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

from foamFileOperation import writeScalarToFile

def onehill_to_twohills(mask_1h, h1, foamfile, foamname, fill=0.0):
    ''' Expand scalar field h1 to two hills by padding.
    '''
    h2 = np.ones(h1.shape[0]*2) * fill
    h2[mask_1h] = h1
    writeScalarToFile(h2, foamfile)
    change_foam_name(foamfile, foamname)

def change_foam_name(foamfile, name):
    fh, abspath = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(foamfile) as old_file:
            for line in old_file:
                if 'object' in line:
                    new_file.write('    object      ' + name + ';')
                else: 
                    new_file.write(line)
    remove(foamfile)
    move(abspath, foamfile)
