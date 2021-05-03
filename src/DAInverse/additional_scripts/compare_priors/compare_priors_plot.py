
import matplotlib as mpl
import matplotlib.pyplot as plt
import fnmatch
import numpy as np
import os

import hillShape


def plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    p, = plt.plot(y,h,'g-')
    return p



def plotSamples(casedir, casename, timeDir, para, xPosList, samplePlotList, norm=1.0, scaleV=1.0, xCol=0, yCol=1):
    for xPos in xPosList:
        xSum = 0
        xMean = 0
        for case in os.listdir(casedir):
            if fnmatch.fnmatch(case, casename+'-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(case, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = casedir + os.sep + case + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + para + '.xy'
                    M = np.loadtxt(filename, comments = '%')
                    y = M[:,yCol]
                    x = M[:,xCol] * scaleV / norm + float(xPos)
                    xSum = xSum + x
                    p1, = plt.plot(x, y, 'grey', alpha = 0.5, lw = 1, mfc = 'none')
        xMean = xSum / len(samplePlotList)
        p2, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))
    return p1, p2


#casedirs = ['./sens_only/SE_5', './sens_only/GC_5']
casedirs = ['../SE_5', '../GC_5']
casenames2 = ['SE5', 'GC5']
casename = 'pehill_base'
timeDir = '3000'
filename = 'cv_Xi_Eta_k'
xPosList = [1, 4, 8]
Ns = 60; samplePlotList = range(1,Ns+1)
yCol = 0
xCol = [2,3,4]
names = ['Xi', 'Eta', 'k']
scaleV = [1.0, 1.0, 1.0]
norm = [1.0, 1.0, 1.0]

for casedir,casename2 in zip(casedirs, casenames2):
    for i in range(len(xCol)):
        plt.figure();
        ax1=plt.subplot(111)
        p0     = plotDomain()
        p1, p2 = plotSamples(casedir=casedir, casename=casename, timeDir=timeDir, para=filename, xPosList=xPosList, samplePlotList=samplePlotList, \
            norm=norm[i], scaleV=scaleV[i], xCol=xCol[i], yCol=yCol)
        plt.axis([-0.5, 11, 0, 3.05])
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad$'+str(scaleV[i])+'$ F/F_n+x/H$')
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        mpl.rcParams.update({'font.size':12})
        plt.savefig("figures/" + names[i] + "-" + casename + '_' + casename2 + '_t_' + timeDir + ".pdf")

