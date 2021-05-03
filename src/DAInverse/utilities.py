### <ownDA>
### Version 1.0
### DAmain.model.openFoamCavity.tool
################################################################################
### This Function has such fuctions:
### For replacing character strings in files


##import for handle files
import tempfile
import os
import os.path as ospt
import shutil
import numpy as np
import pdb
import sys
######################  Function description  ########################

# Parameter description:
# (1)filePath: the directory of the file which you want to change
# (2)pattern: the old character string which you want replace
# (3)subst: the new character string

def replace(filePath, pattern, subst):
    #Create temp file
    fh, abs_path = tempfile.mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(filePath)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    #Remove original file
    os.remove(filePath)
    #Move new file
    shutil.move(abs_path, filePath)


def readInputData(InputFile):
    """Function is to split user input parameters make it to a dict
    
    Arg: 
    InputFile: user specifed input file
    
    Return:
    paramDict: parameter dictionary
    """
    paramDict = {}
    try:
        f = open(InputFile, 'r')
    except IOError:
        print 'I cannot open the file you give me: ', InputFile
        sys.exit(1)
    else:
        finDA = open(InputFile, "r")
        for line in finDA:
            # ignore empty or comment lines
            if line.isspace() or line.startswith("#"):
                continue

            # allow colon ":" in input dictionary for convenience
            param, value = line.strip().replace(':',' ').split(None, 1)
            paramDict[param] = value
            # print "in the constructor, param, value, type(value): ", param, value, type(value)
            f.close()
    return paramDict

def extractListFromDict(paramDict, key):
    """ Extract a list from a dictionary item corresponding to key
    
    Arg:
    paramDict: parameter dictionary
    key: key in paramter dictionary

    """
    valueList_ = paramDict[key].split(',')
    valueList = [pn.strip() for pn in valueList_]

    return valueList


def parameterSampling(Ns, Npara, paraVorg, cmu, csigma):
    """ Function is to sample parameters
    
    Args:
    Ns: number of total ensemble members
    Npara: number of parameters
    paraVorg: vector of parameters   

    Returns:
    XCorg
    """

    XCorg = np.zeros([Ns, Npara])
    # To add random error to the ensembleC
    np.random.seed(1000);       # Fix the pseudo random number state
    for i in np.arange(Ns):
        XCorg[i,:] = paraVorg + np.random.normal(cmu, csigma, Npara)
    XCorgMean = np.mean(XCorg, axis=0)
    print "XCorgMean =", XCorgMean
    return XCorg




