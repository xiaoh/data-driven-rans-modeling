#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
from utilities import readInputData
paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir']

scaleV = 1
#samplePlotList = range(1,101)
samplePlotList = [9] 

def plotSamples(scaleV,iter):
    xSum = 0
    xMean = 0
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, 'duct-tmp*'):
            plotFlag = False
            for sample in samplePlotList:
                if fnmatch.fnmatch(file, '*_' + str(sample) + '.0*'):
                    plotFlag = True
            if plotFlag ==True:
                filename = file + '/postProcessing/sets/' + timeDir + '.000000' + \
                            '/line_x' + str(iter) + '_' + 'Tau' + '.xy'
                M = np.loadtxt(filename, comments = '%')
                y = M[:,0] * 2
                x = (M[:,4]-M[:,6]) * scaleV
                xSum = xSum + x
                p1, = plt.plot(x, y, 'grey', alpha = 0.5, lw = 1, mfc = 'none')
    xMean = xSum / len(samplePlotList)
    p2, = plt.plot(xMean, y, 'b-', lw=2,dashes=(9,2,2,2))
    return p1, p2

def plotBaseline(scaleV, iter):
    filename = 'duct-observation/postProcessing/sets/0/line_x' + str(iter) + '_' + 'Tau' + '.xy'
    M = np.loadtxt(filename, comments = '%')
    y = M[:,0] * 2
    x = (M[:,4]-M[:,6]) * scaleV
    p3, = plt.plot(x, y, '-', color='darkred',dashes=(12,3),lw=2, markevery = 3,  mfc = 'none')
    return p3

def plotDNS(filename,scaleV,iter):
    Mdns = np.loadtxt(filename)
    #x = (Mdns[:,5] - Mdns[:,7])/(Mdns[:,5]+Mdns[:,7])*scaleV + iter*0.125
    x = (Mdns[:,5] - Mdns[:,7])*scaleV
    y = Mdns[:,1] * 2
    p4, = plt.plot(x, y, 'k-',lw=2, mfc = 'none')
    return p4

titleList = ['0.25', '0.5', '0.6', '0.75', '1.0']

xlimLists = [[-0.1,0.1],[-0.05,0.05],[-0.05,0.05],[-0.05,0.05],[-0.05,0.05]]

fig, ax = plt.subplots(nrows=1, ncols=5, sharex=True, sharey=True)
for iter in range(1,6):
    ax[iter-1]=plt.subplot(1,5,iter)
    p1, p2 = plotSamples(scaleV, iter)
    p3 = plotBaseline(scaleV, iter)
    p4 = plotDNS('DNS/Tau'+str(iter), 1, iter)
    plt.ylim([0,1])
    #plt.xlim([-0.1,0.1])
    plt.xlim(xlimLists[iter-1])
    plt.title('$y/h='+titleList[iter-1]+'$',fontsize=12)
    plt.locator_params(axis='x',nbins=3)

for iter in range(1,5):
    plt.setp(ax[iter].get_yticklabels(), visible=False)
plt.tight_layout()
#fig.text(0.5, 0.04, r'$\overline{v^{,2}}- \overline{w^{,2}}$', ha='center', fontsize=16)
fig.text(0.5, 0.04, r'$\tau_{yy}- \tau_{zz}$', ha='center', fontsize=16)
fig.text(0.04, 0.5, "$z/h$", va='center', rotation='vertical', fontsize=18)
#matplotlib.rcParams.update({'font.size':16})
fig.subplots_adjust(bottom=0.12,left=0.1)
#plt.tight_layout()
plt.savefig("vw.pdf")
plt.show()

