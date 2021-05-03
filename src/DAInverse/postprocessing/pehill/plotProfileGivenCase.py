#!/usr/bin/env python

# description        :Comparing Tau and velocity profiles for given cases.

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :May.28, 2016
# revision           :May.28, 2016
########################################################################################################################
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from pylab import *
import foamFileOperation as foam
from utilities import readInputData
import ast


######################################### Global Constant ##############################################################
Ub = 0.028
nu = 1e-5
lineIdx = [3, 4, 5, 6, 7, 8, 9, 10]
xPosList = [1, 2, 3, 4, 5, 6, 7, 8] 
########################################################################################################################

def genColorVec(caseNameList):
    """
    Generate colors vector corresponding to different cases
    """
    colors = []
    for i in range(len(caseNameList)):
        colors.append('#%06X' % randint(0, 0xFFFFFF))
    return colors


def plotU(caseNameList, timeDirList, scaleV):
    """
    Compare U profiles for a case list
    """
    plt.clf()
    ax1=plt.subplot(111)
    
    i = 0
    pList = range(len(caseNameList))
    for caseName in caseNameList:
        print "plot U profiles for case: ", caseName, "at timeDir: ", timeDirList[i]   
        maxRelNut = _calculateMaxNut(caseName, timeDirList[i])
        print "Maximum of Nut / Nu = ", maxRelNut
        if maxRelNut > 1e3:
            print "WARNING: the ratio of nut/nu is very large!!!"
        uProfileDict = _extractProfilesDict(caseName, timeDirList[i], 'U')
        # plot profiles
        legendLabel = caseNameList[i] + '_' + timeDirList[i]
        idx = 0
        for xPos in xPosList:
            M = uProfileDict[str(lineIdx[idx])]
            y = M[:, 0]
            x = M[:, 1] * scaleV / Ub + float(xPos)
            pList[i], = plt.plot(x, y, color=colors[i], alpha=1.0, lw = 2.0, label=legendLabel)
            legendLabel = None
            if i == 0:
                xzero = np.zeros(x.shape) + float(xPos)
                plt.plot(xzero, y, '--', color='black', alpha=1.0, lw = 1.0, mfc = 'none')
            idx = idx + 1                   
        i = i + 1   
    plt.legend(bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)    
    plt.axis([-0.5, 11, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+'$ U_x/U_b+x/H$')         
    _plotDomain()
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':8})    
    figureName = './U_compare'+caseName+'.pdf'
    plt.savefig(figureName)


def plotTau(caseNameList, timeDirList, scaleV):
    """
    Compare Tau profiles for a case list
    """
    plt.clf()
    ax1=plt.subplot(111)
    
    i = 0
    pList = range(len(caseNameList))
    for caseName in caseNameList:
        print "plot Tauxx profiles for case: ", caseName, "at timeDir: ", timeDirList[i]   
        maxRelNut = _calculateMaxNut(caseName, timeDirList[i])
        print "Maximum of Nut / Nu = ", maxRelNut
        if maxRelNut > 1e3:
            print "WARNING: the ratio of nut/nu is very large!!!"
        uProfileDict = _extractProfilesDict(caseName, timeDirList[i], 'Tau')
        # plot profiles
        legendLabel = caseNameList[i] + '_' + timeDirList[i]
        idx = 0
        for xPos in xPosList:
            M = uProfileDict[str(lineIdx[idx])]
            y = M[:, 0]
            x = M[:, 1] * scaleV / Ub /Ub + float(xPos)
            pList[i], = plt.plot(x, y, color=colors[i], alpha=1.0, lw = 2.0, label=legendLabel)
            legendLabel = None
            if i == 0:
                xzero = np.zeros(x.shape) + float(xPos)
                plt.plot(xzero, y, '--', color='black', alpha=1.0, lw = 1.0, mfc = 'none')
            idx = idx + 1                   
        i = i + 1   
    plt.legend(bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)    
    plt.axis([-0.5, 11, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+r'$\tau_{xx}/U_b^2+x/H$')         
    _plotDomain()
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':8})    
    figureName = './Tauxx_compare'+caseName+'.pdf'
    plt.savefig(figureName)

    plt.clf()
    ax1=plt.subplot(111)
    
    i = 0
    pList = range(len(caseNameList))
    for caseName in caseNameList:
        print "plot Tauxy profiles for case: ", caseName, "at timeDir: ", timeDirList[i]   
        maxRelNut = _calculateMaxNut(caseName, timeDirList[i])
        print "Maximum of Nut / Nu = ", maxRelNut
        if maxRelNut > 1e3:
            print "WARNING: the ratio of nut/nu is very large!!!"
        uProfileDict = _extractProfilesDict(caseName, timeDirList[i], 'Tau')
        # plot profiles
        legendLabel = caseNameList[i] + '_' + timeDirList[i]
        idx = 0
        for xPos in xPosList:
            M = uProfileDict[str(lineIdx[idx])]
            y = M[:, 0]
            x = M[:, 2] * scaleV / Ub /Ub + float(xPos)
            pList[i], = plt.plot(x, y, color=colors[i], alpha=1.0, lw = 2.0, label=legendLabel)
            legendLabel = None
            if i == 0:
                xzero = np.zeros(x.shape) + float(xPos)
                plt.plot(xzero, y, '--', color='black', alpha=1.0, lw = 1.0, mfc = 'none')            
            idx = idx + 1                   
        i = i + 1   
    plt.legend(bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)    
    plt.axis([-0.5, 11, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+r'$\tau_{xy}/U_b^2+x/H$')         
    _plotDomain()
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':8})    
    figureName = './Tauxy_compare'+caseName+'.pdf'
    plt.savefig(figureName)


    
def _calculateMaxNut(caseName, timeDir):
    """
    calculate maximum of nu_T/nu to check the convergence
    """
    filename = caseName + '/' + timeDir + '.000000/nut'
    nutVel = foam.readScalarFromFile(filename)
    relNutVel = nutVel / nu
    maxRelNut = np.max(relNutVel)
    return maxRelNut
    
        
def _extractProfilesDict(caseName, timeDir, var='U'):
    """
    Extract profiles of U or Tau with given xPos vector
    """
    
    profileDict = {}
    for idx in lineIdx:
        filename = caseName + '/postProcessing/sets/' + timeDir + '.000000' + '/line' + str(idx) + '_' + var + '.xy'        
        M = np.loadtxt(filename, comments = '%')
        profileDict[str(idx)] = M    
    return profileDict

def _plotDomain():
    """
    Plot domain of periodic hill
    """
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    plt.plot(y,h,'g-')
    
if __name__ == "__main__":
    caseNameList = ['pehill_TwoHill2800_cell3000', 'pehill_TwoHill2800_enkf']
    timeDirList = ['10000', '10000']
    np.random.seed(1000)
    colors = genColorVec(caseNameList) 
    plotU(caseNameList, timeDirList, 2.0)
    #plotTau(caseNameList, timeDirList, 20)    
