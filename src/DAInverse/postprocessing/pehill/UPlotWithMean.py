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

if not os.path.exists('figures/'):
    os.mkdir('figures/')

xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
samplePlotList = range(1,61)
#samplePlotList = [2]
#samplePlotList = range(28,30)

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
        xSum = 0
        xMean = 0
        Nsample = 0
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                        Nsample = Nsample + 1
                if plotFlag ==True:
                    filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + para + 'Norm.xy'
                    M = np.loadtxt(filename, comments = '%')
                    y = M[:,yCol]
                    x = M[:,xCol] * scaleV + float(xPos)
                    xSum = xSum + x
                    p1, = plt.plot(x, y, 'grey', alpha = 0.5, lw = 1, mfc = 'none')
        xMean = xSum / Nsample
        #p2, = plt.plot(xMean, y, 'b-',lw=2,dashes=(7,2,2,2,2,2))
        p2, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))
    return p1, p2

def calculateTKE():
    for xPos in xPosList:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + 'TauNorm.xy'
                    M = np.loadtxt(filename, comments = '%')
                    tke = 0.5*(M[:,1]+ M[:,4] + M[:,6])
                    Mfull = np.vstack((M[:,0:7].transpose(),tke))
                    np.savetxt(filename, Mfull.transpose())
        filename = 'sets/X/dnsTau' + str(xPos)       
        M = np.loadtxt(filename, comments = '%')
        tke = 0.5*(M[:,1]+ M[:,4] + M[:,6])
        Mfull = np.vstack((M[:,0:7].transpose(),tke))
        np.savetxt(filename, Mfull.transpose())
        filename = 'pehill_base-observation/postProcessing/sets/0' + \
                  '/line_x' + str(xPos) + '_' + 'TauNorm.xy'       
        M = np.loadtxt(filename, comments = '%')
        tke = 0.5*(M[:,1]+ M[:,4] + M[:,6])
        Mfull = np.vstack((M[:,0:7].transpose(),tke))
        np.savetxt(filename, Mfull.transpose())


def calculateUV():
    for xPos in xPosList:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + 'TauNorm.xy'
                    M = np.loadtxt(filename, comments = '%')
                    UV = (M[:,1]- M[:,4])
                    Mfull = np.vstack((M[:,0:7].transpose(),UV))
                    np.savetxt(filename, Mfull.transpose())
        filename = 'sets/X/dnsTau' + str(xPos)       
        M = np.loadtxt(filename, comments = '%')
        UV = (M[:,1]- M[:,4])
        Mfull = np.vstack((M[:,0:7].transpose(),UV))
        np.savetxt(filename, Mfull.transpose())
        filename = 'pehill_base-observation/postProcessing/sets/0' + \
                  '/line_x' + str(xPos) + '_' + 'TauNorm.xy'       
        M = np.loadtxt(filename, comments = '%')
        UV = (M[:,1]- M[:,4])
        Mfull = np.vstack((M[:,0:7].transpose(),UV))
        np.savetxt(filename, Mfull.transpose())

def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
    M = np.loadtxt(filename, comments = '%')
    y = M[:,yCol]
    x = M[:,xCol] * scaleV + float(xPos)
    if dash == 'None':
        p, = plt.plot(x, y, ls, color = cl, lw=2)
    else:
        p, = plt.plot(x, y, ls, color = cl, lw=2, dashes=dash)
    return p

def plotBaseline(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = 'pehill_base-observation/postProcessing/sets/0' + \
                  '/line_x' + str(xPos) + '_' + para + 'Norm.xy'
        p3 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'darkred', (12,3), para)
    return p3

def plotDNS(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        if para == 'U':
            filename = 'sets/X/x' + str(xPos) + '.0_U2.xy'
        else:
            filename = 'sets/X/dnsTau' + str(xPos)
        #p4 = plotProfile(filename, xPos, scaleV, xCol, yCol, 'k-', (9,2,2,2), para)
        p4 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'black', 'None', para)
    return p4

def plotObs(scaleV):
    filename = 'obsVPlot'
    M = np.loadtxt(filename)
    y = M[:,1]
    x = M[:,3]
    xPos = M[:,0]
    x = x*scaleV/0.028 + xPos
    p5, = plt.plot(x,y,'kx',markersize=8,markeredgewidth=2.5)
    return p5

def plotProfiles(scaleV, xCol, yCol, para = 'U'):
    p1, p2 = plotSamples(scaleV, xCol, yCol, para)
    p3 = plotBaseline(scaleV, xCol, yCol, para)
    p4 = plotDNS(scaleV, xCol, yCol, para)
    return p1, p2, p3, p4

def main(iShow=False):

    #plot Ux
    if plotUx == 'True':
        scaleV = 2
        plt.figure();
        ax1=plt.subplot(111)
        plotDomain()
        p1, p2, p3, p4 = plotProfiles(scaleV, 1, 0)
        p5 = plotObs(scaleV)
        plt.axis([-0.5, 11, 0, 3.05])
        plt.ylabel("$y/H$")
        #plt.xlabel('normalized velocity'+r'$\, (U_x)$')
        plt.xlabel(r'$x/H;\quad  2U_x/U_b+x/H$')
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        #lg = plt.legend([p1,p2,p3,p4,p5],["samples","sample mean","baseline","DNS","observations"],loc = 7)
        #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 7)
        #lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':12})
        #plt.savefig("Ux"+caseName+".pdf")
        plt.savefig("./figures/Ux-pehill-"+caseName+'_'+str(timeDir)+".pdf")

    #plot Tau 
    scaleVList = [15,30,15,15,15,15,10]
    calculateUV()
    if plotTau == 'True':
        TauLabels = [r'$\tau_{xx}$',r'$\tau_{xy}$', r'$\tau_{xz}$', r'$\tau_{yy}$', r'$\tau_{yz}$', r'$\tau_{zz}$', r'$uu - vv$']
        for iter in range(1,8):
            scaleV = scaleVList[iter-1] 
            plt.figure();
            ax1=plt.subplot(111)
            plotDomain()
            p1, p2, p3, p4 = plotProfiles(scaleV, iter, 0, 'Tau')
            plt.axis([-0.5, 9.5, 0, 3.05])
            plt.ylabel("$y/H$")
            plt.xlabel(r'$normalized\, $' + TauLabels[iter-1])
            #plt.xlabel(r'$x/H;\quad  10k/U_b^2+x/H $')
            ax1.set_aspect(1.3)
            fig = plt.gcf()
            fig.set_size_inches(8, 3.7)
            #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 7)
            #lg.draw_frame(False)
            matplotlib.rcParams.update({'font.size':12})
            plt.savefig("./figures/tau" + str(iter) + "_" + caseName+'_'+str(timeDir)+ ".pdf")
            #plt.savefig("tke-" + caseName+ ".pdf")

    if iShow:
        plt.show()

if __name__ == "__main__":
   #calculateTKE()
   main(False)
