#!/aoe/bin/python27

# description        :Calculate the max nut/nu

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :May.12, 2016
# revision           :May.12, 2016
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

nu = 1e-5

def obtainMaxNut():
    """
    """
    maxRelNut = []
    for sampleFile in os.listdir('.'): 
        if fnmatch.fnmatch(sampleFile, 'pehill_base-tmp*'):
            filename = sampleFile + '/' + timeDir + '.000000/nut'
            nutVel = foam.readScalarFromFile(filename)
            relNutVel = nutVel / nu
            maxRelNut.append(np.max(relNutVel))
    ax1=plt.subplot(111)
    pT, = plt.plot(np.arange(len(maxRelNut)), maxRelNut, 'o', color='black')

    plt.ylabel(r'$\nu_t / \nu$')
    plt.xlabel('Sample index')
    figureName = './propagationAnalysis/relativeNut.pdf'
    plt.savefig(figureName)
      
def plotSample_U(scaleV, sampleIdx):

    
    ax1=plt.subplot(111)
    for xPos in xPosList:
        sampleFile = 'pehill_base-tmp_' + str(sampleIdx) + '.0'
        filename = sampleFile + '/postProcessing/sets/' + timeDir + '.000000' + \
                    '/line_x' + str(xPos) + '_UNorm.xy'

        M = np.loadtxt(filename, comments = '%')

        y = M[:, 0]
        x = M[:, 1] * scaleV + float(xPos)
        p1, = plt.plot(x, y, color='#1b9e76',alpha=1.0, lw = 2.0)

    filename = sampleFile + '/' + timeDir + '.000000/nut'
    nutVel = foam.readScalarFromFile(filename)
    relNutVel = nutVel / nu
    plt.title(r'$\nu_t / \nu = $'+str(np.max(relNutVel)))

    p3 = plotBaseline(scaleV, 1, 0, 'U')
    p4 = plotDNS(scaleV, 1, 0, 'U')
    plotDomain()
    plt.axis([-0.5, 11, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+'$ U_x/U_b+x/H$')
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    lg = plt.legend([p1,p3,p4],["samples_"+str(sampleIdx),"baseline","DNS"], bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)
    lg.draw_frame(False)
    matplotlib.rcParams.update({'font.size':12})
    figureName = './propagationAnalysis/Usample_' + str(sampleIdx) + '.pdf'
    plt.savefig(figureName)

def plotSample_Tau(scaleV, sampleIdx):

    plt.clf()
    ax1=plt.subplot(111)
    for xPos in xPosList:
        sampleFile = 'pehill_base-tmp_' + str(sampleIdx) + '.0'
        filename = sampleFile + '/postProcessing/sets/' + timeDir + '.000000' + \
                    '/line_x' + str(xPos) + '_TauNorm.xy'

        M = np.loadtxt(filename, comments = '%')

        y = M[:, 0]
        x = M[:, 1] * scaleV + float(xPos)
        p1, = plt.plot(x, y, color='#1b9e76',alpha=1.0, lw = 2.0)

    filename = sampleFile + '/' + timeDir + '.000000/nut'
    nutVel = foam.readScalarFromFile(filename)
    relNutVel = nutVel / nu
    plt.title(r'$\nu_t / \nu = $'+str(np.max(relNutVel)))
    
    p3 = plotBaseline(scaleV, 1, 0, 'Tau')
    p4 = plotDNS(scaleV, 1, 0, 'Tau')
    plotDomain()
    plt.axis([-0.5, 9.5, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+'$ U_x/U_b+x/H$')
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    lg = plt.legend([p1,p3,p4],["samples_"+str(sampleIdx),"baseline","DNS"], bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)
    lg.draw_frame(False)
    matplotlib.rcParams.update({'font.size':12})
    figureName = './propagationAnalysis/Tauxx_sample_' + str(sampleIdx) + '.pdf'
    plt.savefig(figureName) 


    plt.clf()
    ax1=plt.subplot(111)
    for xPos in xPosList:
        sampleFile = 'pehill_base-tmp_' + str(sampleIdx) + '.0'
        filename = sampleFile + '/postProcessing/sets/' + timeDir + '.000000' + \
                    '/line_x' + str(xPos) + '_TauNorm.xy'

        M = np.loadtxt(filename, comments = '%')

        y = M[:, 0]
        x = M[:, 2] * scaleV + float(xPos)
        p1, = plt.plot(x, y, color='#1b9e76',alpha=1.0, lw = 2.0)

    filename = sampleFile + '/' + timeDir + '.000000/nut'
    nutVel = foam.readScalarFromFile(filename)
    relNutVel = nutVel / nu
    plt.title(r'$\nu_t / \nu = $'+str(np.max(relNutVel)))

    p3 = plotBaseline(scaleV, 2, 0, 'Tau')
    p4 = plotDNS(scaleV, 2, 0, 'Tau')
    plotDomain()
    plt.axis([-0.5, 9.5, 0, 3.05])
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad$'+r'$'+str(scaleV)+'$'+'$ U_x/U_b+x/H$')
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    lg = plt.legend([p1,p3,p4],["samples_"+str(sampleIdx),"baseline","DNS"], bbox_to_anchor=(0.05, 1.01), loc=3, ncol=3)
    lg.draw_frame(False)
    matplotlib.rcParams.update({'font.size':12})
    figureName = './propagationAnalysis/Tauxy_sample_' + str(sampleIdx) + '.pdf'
    plt.savefig(figureName)



def plotBaseline(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = 'pehill_base-Run/postProcessing/sets/'+timeDir +'.000000'+ \
                  '/line_x' + str(xPos) + '_' + para + 'Norm.xy'
        p3 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'darkred', (12,3), para)
    return p3

def plotDNS(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        if para == 'U':
            filename = 'DNS/X/x' + str(xPos) + '.0_U2.xy'
        else:
            filename = 'DNS/X/dnsTau' + str(xPos) 
        p4 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'black', 'None', para)
    return p4

def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
    M = np.loadtxt(filename, comments = '%')
    y = M[:,yCol]
    x = M[:,xCol] * scaleV + float(xPos)
    if dash == 'None':
        p, = plt.plot(x, y, ls, color = cl, lw=2)
    else:
        p, = plt.plot(x, y, ls, color = cl, lw=2, dashes=dash)
    return p

def plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    plt.plot(y,h,'g-')

if __name__ == "__main__":
    
    paramDict = readInputData('plotInfo.in')
    timeDir = paramDict['timeDir']
    caseName = paramDict['caseName']
    #obtainMaxNut()
    xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
    
    
    for sampleIdx in range(100):
        print "sample = ", sampleIdx+1
        scaleV_U = 2.0    
        plotSample_U(scaleV_U, sampleIdx+1)
        scaleV_Tau = 20.0
        plotSample_Tau(scaleV_Tau, sampleIdx+1)
