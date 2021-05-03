#!/usr/bin/env python

# description        :Functions of formating output
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.17, 2016
# revision           :Oct.17, 2016

########################################################################################################################

## Import system modules
# sci computing
import numpy as np
# system, file operation
import pdb
import os
# sklearn importing
import sklearn as sk 
#from sklearn.ensemble.forest import RandomForestRegressor
# plotting
import matplotlib.pyplot as plt  # for plotting
import matplotlib as mp
# Import local modules



def buildStarBlock(outStr, nStar = 100):
    """
    This is output section title
    """
    
    nStr = len(outStr)  
    nDiff = nStar - nStr
    if nDiff > 10:
        nSpace = int(nDiff/2.0)
        startStr = ['*'] * nStar
        startStr = ''.join(startStr)
        blockStr = [' '] * nSpace
        blockStr = ''.join(blockStr)
        buildSeparateLine('-', nStar)
        print startStr
        print blockStr, outStr
        print startStr

    else:
        outStrSplit = outStr.split()
        nStrSplit = len(outStrSplit)
        outStrUp = ' '.join(outStrSplit[:int(nStrSplit/2)])
        outStrDown = ' '.join(outStrSplit[int(nStrSplit/2):])
        nDiff = nStar - len(outStrUp)
        nSpace = int(nDiff/2.0)
        startStr = ['*'] * nStar
        startStr = ''.join(startStr)
        blockStr = [' '] * nSpace
        blockStr = ''.join(blockStr)
        buildSeparateLine('-', nStar)
        print startStr
        print blockStr, outStrUp
        print blockStr, outStrDown
        print startStr
    #blockFull = ''.join([' '] * nStar)
    #print blockFull  

def buildSeparateLine(marker, length=100):
    """
    Build a separate line with given marker and length
    """
    markerLine = [marker] * length
    markerLine = ''.join(markerLine)
    print markerLine
    
def addSpace(outstr, length=20):
    """
    add space befor out str
    """
    blockStr = [' '] * length
    blockStr = ''.join(blockStr)
    print blockStr, outStr        


if __name__ == "__main__": 
    
    outStr = 'Building the solver'
    buildStarBlock(outStr)
    buildSeparateLine('^')
