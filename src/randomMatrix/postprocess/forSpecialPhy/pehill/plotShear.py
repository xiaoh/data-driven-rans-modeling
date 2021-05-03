#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import hillShape
import pdb
import pandas as pd
from scipy.stats import gaussian_kde
from utilities import readInputData
from scipy.interpolate import UnivariateSpline

paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir']
plotUx = paramDict['plotUx']
plotUy = paramDict['plotUy']
plotTau = paramDict['plotTau']
plotDeltaTau = paramDict['plotDeltaTau']
caseName = paramDict['caseName']
Ub = 0.028

if not os.path.exists('figures/'):
    os.mkdir('figures/')

xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
allWidth = [2,2,2]
allDashes = [(12,0),(12,0),(12,0)]
allLabels = ['samples', 'baseline case', 'benchmark data']
xLabel = '$y$'
yLabels = ['$U_x$', '$U_y$']
transparency = 0.3
samplePlotList = range(1,61)
#samplePlotList = [29] 

def detectBubble(xList, tauwList):
    crossPos = [0,0]
    signs = np.sign(tauwList)
    crossList = np.where(np.diff(np.sign(signs)))[0]
    if crossList.shape[0] >= 1:
        if tauwList[0] > 0:
            signs = np.sign(tauwList)
            crossList = np.where(np.diff(np.sign(signs)))[0]
            crossPos[0] = crossList[0]
            crossPos[1] = crossList[1]
        else:
            signs = np.sign(tauwList)
            crossList = np.where(np.diff(np.sign(signs)))[0]
            crossPos[1] = crossList[0]
        signs = np.sign(tauwList)
        crossList = np.where(np.diff(np.sign(signs)))[0]
        iter = -1
        while xList[crossList[iter]] >=7:
            iter = iter-1
        crossPos[1] = crossList[iter]
        #iter = 0
        #while xList[crossList[iter]] < 3 and iter < crossList.shape[0]-2:
        #    iter = iter + 1
        #crossPos[1] = crossList[iter]
    return crossPos[0], crossPos[1]

def getPos(xList, yList, index):
    pos = (xList[index+1]-xList[index])/(yList[index+1]-yList[index])*(0-yList[index])+xList[index]
    if pos < 0:
        pos = 0
    return pos

def ansysBubble(xList, tauwList):
    index1, index2 = detectBubble(xList, tauwList)
    detachPos = getPos(xList, tauwList, index1)
    attachPos = getPos(xList, tauwList, index2)
    bubbleL = attachPos - detachPos
    return bubbleL, attachPos

def plotSamples(scaleV, xCol, yCol):
    ySum = 0
    xMean = 0
    sampleNum = 0
    scatSample = [0,0]

    for file in os.listdir('.'):
          
        if fnmatch.fnmatch(file, 'pehill_base-tmp*'):

            plotFlag = False
            for sample in samplePlotList:
                if fnmatch.fnmatch(file, '*_'+str(sample) + '.0*'):
                    plotFlag = True
                    sampleNum = sampleNum + 1               
            if plotFlag == True:
                filenames = [file + '/postProcessing/surfaces/'+timeDir+'.000000/wallShearStress_bottomWall.raw']
                # Plot
                #pdb.set_trace()
                for i in range(len(filenames)):
                    filename = filenames[i]
                    M = np.loadtxt(filename, comments = '#')
                    M = M[1:M.shape[0]/2,:]
                    y = M[:,yCol]*(-1)*2/Ub/Ub
                    ySum = ySum + y
                    x = M[:,xCol] 
                    xMean = x
                    scat = ansysBubble(x,y)
                    if scat > 0:
                        scatSample = np.vstack((scatSample,scat))
                    p1, = plt.plot(x, y,'-', color='#1b9e76',alpha=transparency, lw = 1, markevery=20,mfc='none')

    scatSample = scatSample[1:,:]
    yMean = ySum/sampleNum
    scatMean = ansysBubble(x,yMean)
    p2, = plt.plot(x, yMean, '--', dashes=(9,2,2,2), color='blue', lw = 2, markevery=20,mfc='none')
    Mmean = np.transpose(np.vstack((x, yMean)))
    np.savetxt('meanShearStress', Mmean)
    return p1, p2, scatSample, scatMean

def plotBaseline(scaleV, xCol, yCol):
    filename = 'pehill_base-Run/postProcessing/surfaces/'+timeDir + \
                '.000000/wallShearStress_bottomWall.raw'
    M = np.loadtxt(filename, comments = '#')
    M = M[1:M.shape[0]/2,:]
    y = M[:,yCol]*(-1)*2/Ub/Ub
    x = M[:,xCol]
    scatBase = ansysBubble(x,y)
    p3, = plt.plot(x, y, '--', dashes=(12,3), color='darkred', lw = 2, markevery=20,mfc='none')
    return p3, scatBase

def plotDNS(scaleV, xCol, yCol):
    filename = 'DNS/X/Cf.csv'
    M = np.loadtxt(filename, comments = '#')
    M = M[1:M.shape[0]/2,:]
    y = M[:,yCol]
    x = M[:,xCol]
    scatDNS = ansysBubble(x[1:],y[1:])
    p4, = plt.plot(x, y, '-', color='black', lw = 2,mfc='none')
    #p4.set_dashes([15,5,2,5])
    return p4, scatDNS

def plotPDF(samples):
    p, x = np.histogram(samples, bins=10)
    x = x[:-1]+(x[1]-x[0])/2
    f = UnivariateSpline(x, p, s=10)
    plt.plot(x,f(x))

def main(iShow=False):
    # plot wall shear stress
    scaleV = 1
    fig = plt.figure(figsize=(8,7));
    gs = matplotlib.gridspec.GridSpec(2,1,height_ratios=[4,1])
    #plt.subplot(2,1,1)
    plt.subplot(gs[0])
    #matplotlib.rc('font',family='Times New Roman')
    p1, p2, scatSample, scatMean = plotSamples(scaleV,0,3);
    p3, scatBase = plotBaseline(scaleV,0,3);
    p4, scatDNS = plotDNS(scaleV,0,3);
    plt.plot([0,9],[0,0],'k--');
    plt.axis([0, 9, -0.02, 0.06])
    plt.locator_params(axis = 'y', nbins = 4)
    plt.ylabel("Wall Shear Stress"+r"$ (\tau_w)$",fontsize=18)
    plt.xlabel(r'$x$',fontsize=18)
    lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS (Breuer et al. 2009)"],loc = 0, prop={'size':12})
    #lg.draw_frame(False)
    matplotlib.rcParams.update({'font.size':18})

    # scatter plot
    #plt.subplot(2,1,2)
    plt.subplot(gs[1])
    #p1, = plt.plot(scatSample[:,1],np.zeros_like(scatSample[:,1]),'o',markeredgecolor='grey', \
    #        markerfacecolor='None',markersize=12, markeredgewidth=2, alpha=0.01,lw=0.8)
    p1, = plt.plot(scatSample[:,1],np.zeros_like(scatSample[:,1]), 'o',markersize=8,\
             markeredgewidth=0, color='#1b9e76',alpha=transparency)
    p3, = plt.plot(scatBase[1],0,'r^',markeredgecolor='darkred', \
            markerfacecolor='None',markersize=8, markeredgewidth=2, mfc='None')
    p2, = plt.plot(scatMean[1],0,'bx',markeredgecolor='b', \
            markerfacecolor='None',markersize=8, markeredgewidth=2, mfc='None')
    p4, = plt.plot(scatDNS[1],0,'ks',markeredgecolor='k', markerfacecolor='None',\
            markersize=8, markeredgewidth=2, mfc='None')
    plt.xlabel(r"$x_{attach}/H$",fontsize=18)
    #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 0)
    #lg.draw_frame(False)
    plt.xlim([0,9])
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    #matplotlib.rcParams.update({'font.size':18})
    fig.tight_layout()
    #fig.savefig("wallShear_"+timeDir+".pdf")
    fig.savefig("./figures/bubble-"+caseName+".pdf")

    np.savetxt('reattch_'+timeDir, scatSample[:,1])

    # pdf plot
    #fig, ax1 = plt.subplots()
    #ax2 = ax1.twinx()
    ##fig = plt.figure()
    #df = pd.DataFrame(scatSample[:,1],columns=['x'])
    #df['x'].hist(ax=ax1, bins=10, color = 'grey')
    #df['x'].plot(ax=ax2, kind='kde', linewidth=3, color='#d95f02')
    #plt.xlabel('reattachment point')
    ##plt.xlim([0,6])
    
    #fig = plt.figure()
    #density = gaussian_kde(scatSample[:,1])
    #x = np.linspace(0,8,200)
    #plt.plot(x,density(x))

    #fig = plt.figure()
    #plotPDF(scatSample[:,1])
    #fig.savefig("reattach_pdf.pdf")

    if iShow:
        plt.show()

if __name__ == "__main__":
   main(False)
