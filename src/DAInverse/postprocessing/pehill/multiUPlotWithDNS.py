#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from utilities import readInputData

paramDict = readInputData('plotInfo.in')
baseLine = paramDict['baseLine']
perturbLine = paramDict['perturbLine']
benchLine = paramDict['benchLine']
timeDir = paramDict['timeDir']
outputFlag = paramDict['outputFlag']
outputDir = paramDict['outputDir']
plotUx = paramDict['plotUx']
plotUy = paramDict['plotUy']
plotTau = paramDict['plotTau']
plotDeltaTau = paramDict['plotDeltaTau']

if not os.path.exists('figures/'):
    os.mkdir('figures/')
if not os.path.exists(outputDir):
    os.mkdir(outputDir)

#xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
xPosList = [ 2,4,6,8]
allLines = [perturbLine,baseLine,benchLine]
allWidth = [2,2,2]
allDashes = [(12,0),(12,0),(12,0)]
allLabels = ['samples', 'baseline case', 'benchmark data']
xLabel = '$y$'
yLabels = ['$U_x$', '$U_y$']

samplePlotList = range(1,101)

def plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])

    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)

    plt.plot(y,h,'b-')
    #plt.axis([-1, 10, -1, 3.5])

def plotProfiles(scaleV, xCol, yCol, Para = 'U'):
    legendShown = True
    for xPos in xPosList:
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, 'pehill-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(file, '*_'+str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag == True:
                    if Para == 'U':
                        filenames = [file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_' + Para +'norm.xy',\
                                    file + '/postProcessing/sets/0/line_x' + str(xPos) + '_' + Para + 'norm.xy', \
                                    'sets/X/x' + str(xPos) +'.0_U2.xy']
                    else:
                        filenames = [file + '/postProcessing/sets/'+timeDir+'/line_x' + str(xPos) + '_' + Para +'Norm.xy',\
                                    'pehill-observation/postProcessing/sets/0/line_x' + str(xPos) + '_' + Para + 'Norm.xy', \
                                    'sets/X/dnsTau' + str(xPos)]
                    # Plot
                    for i in range(len(filenames)):
                        lwidth = 1
                        op = 1
                        if i == 0:
                            lwidth = 0.5
                            op = 0.5
                        filename = filenames[i]
                        M = np.loadtxt(filename, comments = '%')
                        y = M[:,yCol]
                        x = M[:,xCol] * scaleV + float(xPos)
                        if legendShown == True:
                            plt.plot(x, y, allLines[i], alpha=op, lw = lwidth, label = allLabels[i], markersize=3, markevery=20,mfc='none')
                        else:
                            plt.plot(x, y, allLines[i], alpha=op, lw = lwidth, markersize=3, markevery=20,mfc='none')
                    legendShown = False

def plotProfile(scaleV, xCol, yCol):
    for xPos in xPosList:
       filenames = ['deltaTau' + str(xPos)]
       # Plot
       for i in range(len(filenames)):
           filename = filenames[i]
           M = np.loadtxt(filename, comments = '%')
           y = M[:,yCol]
           x = M[:,xCol] * scaleV + float(xPos)
           plt.plot(x, y, allLines[i], lw = 1, label = allLabels[i], markevery=1,mfc='none')

def main(iShow=False):

    #plot Ux
    if plotUx == 'True':
        scaleV = 1
        plt.figure();
        ax1=plt.subplot(111)
        plotDomain()
        plotProfiles(scaleV, 1, 0)
        plt.axis([-0.5, 16, 0, 3.05])
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H$;   $U/U_b + x/H$')
        ax1.set_aspect(1.5)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        lg = plt.legend(loc = 7)
        lg.draw_frame(False)
        plt.savefig("Ux"+timeDir+".pdf")

    
    #plot Uy
    if plotUy == 'True':
        scaleV = 2
        plt.figure();
        ax1=plt.subplot(111)
        plotDomain()
        plotProfiles(scaleV, 2, 0)
        plt.axis([-0.5, 16, 0, 3.05])
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H$;   $2 V/U_b + x/H$')
        ax1.set_aspect(1.3)
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        lg = plt.legend(loc = 7)
        lg.draw_frame(False)
        plt.savefig("Uy"+timeDir+".pdf")

    #plot Tau 
    if plotTau == 'True':
        TauLabels = ['$uu$','$uv$', '$uw$', '$vv$', '$vw$', '$ww$']
        for iter in range(1,7):
            scaleV = 5
            plt.figure();
            ax1=plt.subplot(111)
            plotDomain()
            plotProfiles(scaleV, iter, 0, 'Tau')
            plt.axis([-0.5, 15.5, 0, 3.05])
            plt.ylabel("$y/H$")
            plt.xlabel(r'$x/H$;   $5$' + TauLabels[iter-1] + r'$/U_b^2 + x/H$' )
            ax1.set_aspect(1.3)
            fig = plt.gcf()
            fig.set_size_inches(8, 3.7)
            lg = plt.legend(loc = 7)
            lg.draw_frame(False)
            plt.savefig("tau" + str(iter) + "_" + timeDir + ".pdf")

    #plot deltaTau 
    if plotDeltaTau == 'True':
        deltaTauLabels = ['$uu$','$uv$', '$uw$', '$vv$', '$vw$', '$ww$']
        for iter in range(1,7):
            scaleV = 20
            plt.figure();
            ax1=plt.subplot(111)
            plotDomain()
            plotProfile(scaleV, iter, 0)
            plt.axis([-0.5, 10, 0, 3.05])
            plt.ylabel("$y/H$")
            plt.xlabel(r'$x/H$;   $20$' + deltaTauLabels[iter-1] + r'$/U_b^2 + x/H$' )
            ax1.set_aspect(1.3)
            fig = plt.gcf()
            fig.set_size_inches(8, 3.7)
            plt.savefig("deltaTau" + str(iter) + ".pdf")

    if outputFlag == 'True':
        os.system('cp pehill*.pdf ' + outputDir)
        os.system('cp tau*.pdf ' + outputDir)
        os.system('cp deltaTau*.pdf ' + outputDir)
    if iShow:
        plt.show()

if __name__ == "__main__":
   main(True)
