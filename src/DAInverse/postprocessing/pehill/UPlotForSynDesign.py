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
baseLine = paramDict['baseLine']
perturbLine = paramDict['perturbLine']
benchLine = paramDict['benchLine']
timeDir = paramDict['timeDir']
plotUx = paramDict['plotUx']
plotUy = paramDict['plotUy']
plotTau = paramDict['plotTau']
plotDeltaTau = paramDict['plotDeltaTau']
caseName = paramDict['caseName']
designCaseLoc = paramDict['designCaseLoc']

if not os.path.exists('figures/'):
    os.mkdir('figures/')

xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]

def plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    plt.plot(y,h,'g-')

def plotSamples(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = designCaseLoc + '/postProcessing/sets/' + timeDir + '.000000' + \
                    '/line_x' + str(xPos) + '_' + para + '.xy'
        M = np.loadtxt(filename, comments = '%')
        y = M[:,yCol]
        x = M[:,xCol] * scaleV + float(xPos)
        p1, = plt.plot(x, y, 'k-', alpha = 1, lw = 2, mfc = 'none')
    return p1

def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
    M = np.loadtxt(filename, comments = '%')
    y = M[:,yCol]
    x = M[:,xCol] * scaleV + float(xPos)
    if dash == 'None':
        p, = plt.plot(x, y, ls, lw=2, color=cl)
    else:
        p, = plt.plot(x, y, ls, lw=2, color=cl, dashes=dash)
    return p

def plotBaseline(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = designCaseLoc + '/postProcessing/sets/0' + \
                  '/line_x' + str(xPos) + '_' + para + '.xy'
        p3 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'darkred', (12,3), para)
    return p3

def plotDNS(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        if para == 'U':
            filename = 'sets/X/x' + str(xPos) + '.0_U2.xy'
        else:
            filename = 'sets/X/dnsTau' + str(xPos)
        p4 = plotProfile(filename, xPos, 1, xCol, yCol, '-', 'black', (9,2,2,2), para)
    return p4

def plotProfiles(scaleV, xCol, yCol, para = 'U'):
    p1 = plotSamples(scaleV, xCol, yCol, para)
    p3 = plotBaseline(scaleV, xCol, yCol, para)
    p4 = plotDNS(scaleV, xCol, yCol, para)
    return p1, p3, p4

def main(iShow=False):

    #plot Ux
    if plotUx == 'True':
        scaleV = 1/0.028
        plt.figure();
        ax1=plt.subplot(111)
        plotDomain()
        p1, p3, p4 = plotProfiles(scaleV, 1, 0)
        plt.axis([-0.5, 14, 0, 3.05])
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad  U_x/U_b+x/H$')
        ax1.set_aspect(1.5)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 4)
        lg = plt.legend([p1,p3,p4],["Syn-truth","Baseline","DNS"],loc = 7)
        lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':12})
        plt.savefig("design.pdf")

    if iShow:
        plt.show()

if __name__ == "__main__":
   main(True)
