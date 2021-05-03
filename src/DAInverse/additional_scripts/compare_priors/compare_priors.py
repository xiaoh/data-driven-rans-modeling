
import matplotlib as mpl
import matplotlib.pyplot as plt
import fnmatch
import numpy as np
import os


def plotSamples(casedir, casename, timeDir, para, xPos, samplePlotList, norm=1.0, scaleV=1.0, xCol=0, yCol=1):
    xSum = 0
    xMean = 0
    first = True
    for case in os.listdir(casedir):
        if fnmatch.fnmatch(case, casename+'-tmp*'):
            plotFlag = False
            for sample in samplePlotList:
                if fnmatch.fnmatch(case, '*_' + str(sample) + '.0*'):
                    plotFlag = True
            if plotFlag ==True:
                #filename = casedir + os.sep + case + '/postProcessing/sets/' + timeDir + '.000000' + \
                filename = casedir + os.sep + case + '/postProcessing/sets/' + timeDir + \
                            '/line_x' + str(xPos) + '_' + para + '.xy'
                M = np.loadtxt(filename, comments = '%')
                y = M[:,yCol]
                x = M[:,xCol]
                if first:
                    X = np.expand_dims(np.copy(x),axis=0)
                    first = False
                else:
                    X = np.concatenate((X,np.expand_dims(x,axis=0)),axis=0)
    return y,X


casedirs = ['./SE_5', './GC_5']
casenames2 = ['SE5', 'GC5']
casename = 'pehill_base'
timeDir = '0'
#filename = 'cv_Xi_Eta_k'
filename = 'deltaEta_deltak_deltaXi'
xPosList = [1, 4, 8]
Ns = 60; samplePlotList = range(1,Ns+1)
yCol = 0
#xCol = [2,3,4]
xCol = [1,2,3]
names = ['delXi', 'delEta', 'delk']

for casedir,casename2 in zip(casedirs, casenames2):
    for i in range(len(xCol)):
        for j in range(len(xPosList)):
            y,X = plotSamples(casedir=casedir, casename=casename, timeDir=timeDir, para=filename, xPos=xPosList[j], samplePlotList=samplePlotList, xCol=xCol[i], yCol=yCol)
            mean = np.mean(X,axis=0)
            stddev = np.std(X,axis=0)
            plt.figure()
            for m in range(X.shape[0]):
                plt.plot(X[m,:], y, 'k-', alpha=0.20)
            plt.plot(mean, y, 'b-')
            plt.plot(mean+stddev, y, 'b--')
            plt.plot(mean-stddev, y, 'b--')
            plt.xlabel(names[i])
            plt.ylabel('$y/H$')
            plt.savefig('./figures/' + names[i] + '_' + casename2 + '_x' + str(xPosList[j]) +'.pdf')

