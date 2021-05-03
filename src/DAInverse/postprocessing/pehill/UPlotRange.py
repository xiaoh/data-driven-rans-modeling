#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import hillShape
import pdb
import scipy
from pylab import *
from scipy import interpolate
from utilities import readInputData
from scipy.stats import gaussian_kde

paramDict = readInputData('plotInfo.in')
baseLine = paramDict['baseLine']
perturbLine = paramDict['perturbLine']
benchLine = paramDict['benchLine']
#timeDir = paramDict['timeDir']
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
    dx = (xmax - xmin)/500
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


def plotSamples(timeDir,cl,ap,scaleV, xCol, yCol, para = 'U'):
    covRsum = 0
    for xPos in xPosList:
        xSum = 0
        xMean = 0
        Nsample = 0
        xArray = []
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                        Nsample = Nsample + 1
                if plotFlag ==True:
                    filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + para + 'Norm.xy'
                    filename_cv = file + '/postProcessing/sets/0' + \
                                '/line_x' + str(xPos) + '_cv.xy'
                    M = np.loadtxt(filename, comments = '%')
                    M_cv = np.loadtxt(filename_cv, comments = '%')
                    y = M[:,yCol]
                    x = M[:,xCol] * scaleV + float(xPos)
                    cv =M_cv[:,1]
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
        filename_dns = 'sets/X/x' + str(xPos) + '.0_U2.xy'
        M_dns = np.loadtxt(filename_dns, comments = '%')
        xDNS = M_dns[:,xCol] * scaleV + float(xPos)
        yDNS = M_dns[:,yCol]
        yDNS[-1] = y[0]*(1-1e-6)
        f = interpolate.interp1d(yDNS[::-1],xDNS[::-1])
        xTruth = f(y)
        #xTruth = scipy.interp(y,yDNS,xDNS)
        covR = coverRatio(xList, xTruth, cv)
        covRsum = covRsum + covR
        print 'The cover ratio at x/H=', xPos, ' is ', covR ,'\n'
        plt.plot(xList[:,0],y,color=cl,lw=1)
        p1, = plt.plot(xList[:,1],y,color=cl,lw=1)
        plt.fill_betweenx(y,xList[:,0],xList[:,1],linewidth=0,facecolor = cl, alpha = ap)
        xMean = xSum / Nsample
        #p2, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))
    print 'The average cover ratio is ', covRsum/len(xPosList) ,'\n'
    return p1

def calculateTKE():
    for xPos in xPosList:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill-tmp*'):
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
        filename = 'pehill-observation/postProcessing/sets/0' + \
                  '/line_x' + str(xPos) + '_' + 'TauNorm.xy'       
        M = np.loadtxt(filename, comments = '%')
        tke = 0.5*(M[:,1]+ M[:,4] + M[:,6])
        Mfull = np.vstack((M[:,0:7].transpose(),tke))
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
        filename = 'pehill-observation/postProcessing/sets/0' + \
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
    print 'Prior samples:\n'
    p1 = plotSamples('3000','darkred',0.2,scaleV, xCol, yCol, para)
    print 'Posterior samples:\n'
    p2 = plotSamples('300000','darkblue',0.5,scaleV, xCol, yCol, para)
    #p3 = plotBaseline(scaleV, xCol, yCol, para)
    p4 = plotDNS(scaleV, xCol, yCol, para)
    return p1, p2, p4

def main(iShow=False):

    #plot Ux
    if plotUx == 'True':
        scaleV = 2
        plt.figure();
        ax1=plt.subplot(111)
        plotDomain()
        p1, p2, p4 = plotProfiles(scaleV, 1, 0)
        p1 = Rectangle((0,0),1,1,fc='darkred',alpha=0.2,lw=1)
        p2 = Rectangle((0,0),1,1,fc='darkblue',alpha=0.5,lw=1)
        #p5 = plotObs(scaleV)
        plt.axis([-0.5, 16, 0, 3.05])
        plt.ylabel("$y/H$")
        #plt.xlabel('normalized velocity'+r'$\, (U_x)$')
        plt.xlabel(r'$x/H;\quad  2U_x/U_b+x/H$')
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        lg = plt.legend([p1,p2,p4],["Prior","Posterior","DNS"],loc = 7)
        #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 7)
        lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':12})
        #plt.savefig("Ux"+caseName+".pdf")
        plt.savefig("Ux-pehill-confidence.pdf")

    #plot Tau 
    scaleVList = [15,30,15,15,15,15,10]
    if plotTau == 'True':
        TauLabels = [r'$\tau_{xx}$',r'$\tau_{xy}$', r'$\tau_{xz}$', r'$\tau_{yy}$', r'$\tau_{yz}$', r'$\tau_{zz}$', r'$k$']
        for iter in range(7,8):
            scaleV = scaleVList[iter-1] 
            plt.figure();
            ax1=plt.subplot(111)
            plotDomain()
            p1, p2, p3, p4 = plotProfiles(scaleV, iter, 0, 'Tau')
            plt.axis([-0.5, 9.5, 0, 3.05])
            plt.ylabel("$y/H$")
            #plt.xlabel(r'$normalized\, $' + TauLabels[iter-1])
            plt.xlabel(r'$x/H;\quad  10k/U_b^2+x/H $')
            ax1.set_aspect(1.3)
            fig = plt.gcf()
            fig.set_size_inches(8, 3.7)
            #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 7)
            #lg.draw_frame(False)
            matplotlib.rcParams.update({'font.size':12})
            #plt.savefig("tau" + str(iter) + "_" + caseName+ ".pdf")
            plt.savefig("tke-" + caseName+ ".pdf")

    if iShow:
        plt.show()

if __name__ == "__main__":
   #calculateTKE()
   main(True)
