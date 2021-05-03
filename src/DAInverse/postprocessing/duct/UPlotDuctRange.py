#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import pdb
from pylab import *
from scipy import interpolate
from utilities import readInputData
from scipy.stats import gaussian_kde

paramDict = readInputData('plotInfo.in')
baseLine = paramDict['baseLine']
perturbLine = paramDict['perturbLine']
benchLine = paramDict['benchLine']
#timeDir = paramDict['timeDir']
outputFlag = paramDict['outputFlag']
outputDir = paramDict['outputDir']
plotUy = paramDict['plotUy']
plotUz = paramDict['plotUz']
plotTau = paramDict['plotTau']
plotDeltaTau = paramDict['plotDeltaTau']
caseName = paramDict['caseName']

if not os.path.exists('figures/'):
    os.mkdir('figures/')

xPosList = [ 1, 2, 3, 4]
#xPosList = [ 2,4,6,8]
allLines = [perturbLine,baseLine,benchLine]
allWidth = [2,2,2]
allDashes = [(12,0),(12,0),(12,0)]
allLabels = ['samples', 'baseline case', 'benchmark data']
xLabel = '$y$'
yLabels = ['$U_x$', '$U_y$']

samplePlotList = range(1,101)
#samplePlotList = [3,4,5,6,9,10]+range(12,52)+range(53,55)+range(60,101) 
#samplePlotList = [9] 

def getBound(x,CDF,confidence):
    for iter in range(0,len(CDF)-1):
        if (CDF[iter] - 1 + confidence) * (CDF[iter+1] -1 + confidence) <= 0:
            xLow = x[iter+1]
        elif (CDF[iter] - confidence) * (CDF[iter+1] - confidence) <= 0:
            xUp = x[iter]
    return [xLow,xUp]

def confRange(x, confidence=0.975):
    xmin = min(x)
    xmax = max(x)
    density = gaussian_kde(x)
    dx = (xmax - xmin)/100
    xr = np.arange(xmin,xmax,dx)
    density = density(xr)
    density /= (dx*density).sum()
    CDF = np.cumsum(density*dx)
    #plt.plot(xr,CDF)
    #plt.show()
    xBound = getBound(xr,CDF,confidence)
    return xBound

def coverRatio(xCover, xTruth, cv):
    cover = 0;
    uncover = 0;
    for iter in range(0,xCover.shape[0]):
        if xTruth[iter] >= xCover[iter,0] and xTruth[iter] <= xCover[iter,1]:
            cover = cover + cv[iter]
        else:
            uncover = uncover + cv[iter]
    return cover/cv.sum()

def plotSamples(timeDir,cl,priorFactor,ap,scaleV, xCol, yCol, para = 'U', tauN = 0):
    covRsum = 0
    for xPos in xPosList:
        xSum = 0
        xMean = 0
        Nsample = 0
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'duct-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                        Nsample = Nsample + 1
                if plotFlag ==True:
                    filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + para + '.xy'
                    M = np.loadtxt(filename, comments = '%')
                    filename_cv = file + '/postProcessing/sets/0/line_x' + str(xPos) + '_cv.xy'
                    M_cv = np.loadtxt(filename_cv, comments = '%')
                    cv = M_cv[:,1]
                    y = M[:,yCol]
                    x = M[:,xCol] * scaleV * priorFactor + float(xPos) * 0.25
                    y = y*2
                    xSum = xSum + x
                    if Nsample == 1:
                        xArray = x
                    else:
                        xArray = np.vstack([xArray,x])

        for iter in range(0,xArray.shape[1]):
            xStates = xArray[:,iter]
            xRange = confRange(xStates)
            if iter ==0:
                xList = xRange
            else:
                xList = np.vstack([xList,xRange])
        filename_dns = 'DNS/Uorg/U' + str(xPos)
        M_dns = np.loadtxt(filename_dns, comments = '%')
        xDNS = M_dns[:,4] * scaleV * priorFactor + float(xPos) * 0.25
        yDNS = M_dns[:,yCol]
        f = interpolate.interp1d(yDNS,xDNS)
        xTruth = f(y/2)
        #xTruth = scipy.interp(y,yDNS,xDNS)
        covR = coverRatio(xList, xTruth, cv)
        covRsum = covRsum + covR
        print 'The cover ratio at x/H=', xPos, ' is ', covR ,'\n'
        plt.plot(xList[:,0],y,color=cl,lw=1)
        p1, = plt.plot(xList[:,1],y,color=cl,lw=1)
        plt.fill_betweenx(y,xList[:,0],xList[:,1],linewidth=0,facecolor = cl, alpha = ap)
    print 'The average cover ratio is ', covRsum/len(xPosList) ,'\n'    
    return p1

def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
    M = np.loadtxt(filename, comments = '%')
    y = M[:,yCol]
    x = M[:,xCol] * scaleV + float(xPos) * 0.25
    y = y*2
    if dash == 'None':
        p, = plt.plot(x, y, ls,lw=3, color=cl, markevery = 3,  mfc = 'none')
    else:
        p, = plt.plot(x, y, ls, color=cl, dashes=dash,lw=2, markevery = 3,  mfc = 'none')
    return p

def plotBaseline(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = 'duct-observation/postProcessing/sets/0/line_x' + str(xPos) + '_' + para + '.xy'
        p3 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'darkred', (12,3),para)
    return p3

def plotDNS(scaleV, xCol, yCol, para = 'U'):
    for xPos in xPosList:
        filename = 'DNS/' + para + str(xPos)
        p4 = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'black', 'None', para)
    return p4

def plotObs(scaleV, xCol, yCol):
    filename = 'obsVPlot'
    M = np.loadtxt(filename)
    y = M[:,yCol]
    x = M[:,xCol]
    x = x*scaleV + 2 * 0.25
    y = y*2
    p5, = plt.plot(x,y,'kx',markersize=8,markeredgewidth=1.5)
    return p5

def plotProfiles(scaleV, xCol, yCol, para = 'U'):
    p1 = plotSamples('10000','darkred',0.3,0.2,scaleV, xCol, yCol, para)
    p2 = plotSamples('200000','darkblue',1,0.5,scaleV, xCol, yCol, para)
    #p3 = plotBaseline(scaleV, xCol, yCol, para)
    #p3 = 0
    return p1, p2

def main(iShow=False):

    #plot Uy
    if plotUy == 'True':
        scaleV = 0.5 
        plt.figure();
        ax1=plt.subplot(111)
        p1, p2 = plotProfiles(scaleV, 2, 0)
        p4 = plotDNS(scaleV, 4, 0)
        #p5 = plotObs(scaleV,4,2)
        plt.axis([0, 1.2, 0, 1])
        plt.ylabel("$z/h$", fontsize=22)
        plt.xlabel(r'$y/h;\quad  y/h+0.5U_y$', fontsize=22)
        plt.xticks([0,0.25,0.5,0.75,1],['0','0.25','0.5','0.75','1'])
        plt.minorticks_on()
        plt.tight_layout()
        #plt.title(r'$U_y$', fontsize=16)
        #ax1.set_aspect(1)
        #rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        #fig.set_size_inches(8, 3.7)
        #lg = plt.legend([p1,p2,p3,p4,p5],["samples","sample mean","baseline","DNS","observation"],loc = 0)
        #lg = plt.legend([p2,p3,p4],["posterior","baseline","DNS"],loc = 0)
        #lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':20})
        #plt.savefig("Uy"+timeDir+caseName+".pdf")
        plt.savefig("Uy-duct-confidence.pdf")

    
    #plot Uz
    if plotUz == 'True':
        scaleV = 0.5 
        plt.figure();
        ax1=plt.subplot(111)
        p1, p2 = plotProfiles(scaleV, 3, 0)
        p4 = plotDNS(scaleV, 3, 0)
        #p5 = plotObs(scaleV,5,2)
        plt.axis([0, 1.2, 0, 1])
        plt.ylabel("$z/h$", fontsize=22)
        plt.xlabel(r'$y/h;\quad  y/h+0.5U_z$', fontsize=22)
        plt.xticks([0,0.25,0.5,0.75,1],['0','0.25','0.5','0.75','1'])
        plt.minorticks_on()
        plt.tight_layout()
        #plt.title(r'$U_z$', fontsize=16)
        #ax1.set_aspect(1)
        #rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        #fig.set_size_inches(8, 3.7)
        #lg = plt.legend([p1,p2,p3,p4,p5],["samples","sample mean","baseline","DNS","observation"],loc = 0)
        #lg = plt.legend([p2,p3,p4],["posterior","baseline","DNS"],loc = 0)
        #lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':20})
        #plt.savefig("Uz"+timeDir+caseName+".pdf")
        plt.savefig("Uz-duct-confidence.pdf")

    #plot Tau 
    scaleVList = [0.01,0.1,0.1,0.05,0.1,0.05,]
    if plotTau == 'True':
        TauLabels = [r'$\tau_{xx}$',r'$\tau_{xy}$', r'$\tau_{xz}$', r'$\tau_{yy}$', r'$\tau_{yz}$', r'$\tau_{zz}$']
        for iter in range(1,7):
            scaleV = scaleVList[iter-1] 
            plt.figure();
            ax1=plt.subplot(111)
            #plotDomain()
            p1, p2, p3 = plotProfiles(scaleV, iter, 0, 'Tau')
            p4 = plotDNS(scaleV, iter+1, 1, 'Tau')
            plt.axis([0, 1.2, 0, 1])
            plt.xlabel("$y$", fontsize=18)
            plt.ylabel("$z$", fontsize=18)
            plt.title(str(scaleV) + TauLabels[iter-1], fontsize=18)
            #ax1.set_aspect(1.3)
            fig = plt.gcf()
            #fig.set_size_inches(8, 3.7)
            lg = plt.legend([p2,p3,p4],["posterior","baseline","DNS"],loc = 0)
            lg.draw_frame(False)
            matplotlib.rcParams.update({'font.size':20})
            plt.savefig("tau" + str(iter) + "_" + timeDir +caseName+ ".pdf")

    if outputFlag == 'True':
        os.system('cp pehill*.pdf ' + outputDir)
        os.system('cp tau*.pdf ' + outputDir)
        os.system('cp deltaTau*.pdf ' + outputDir)
    if iShow:
        plt.show()

if __name__ == "__main__":
   main(True)
