#! /usr/bin/env python
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
from scipy.interpolate import interp1d
import scipy.stats as stats

if not os.path.exists('figures/'):
    os.mkdir('figures/')


# def plotDomain():
#     # Plot the simulation domain
#     y=np.arange(0, 9, 0.01)
#     yext = np.array([9, 9, 0, 0])
#     h=hillShape.profile(y)
#     hext = np.array([1, 3.036, 3.036, 1])
#     y = np.append(y, yext)
#     h = np.append(h, hext)
#     plt.plot(y,h,'g-')

def plotSamples(scaleV, xCol, yCol, para = 'U'):
    dns_all = np.array([])
    bl_all = np.array([])
    po_all = np.array([])
    for xPos in xPosList:
        xSum = 0
        xMean = 0
        for case in os.listdir('.'):
            if fnmatch.fnmatch(case, 'pehill_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(case, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = case + '/postProcessing/sets/' + timeDir + '.000000' + \
                                '/line_x' + str(xPos) + '_' + para + '.xy'
                    M = np.loadtxt(filename, comments = '%')
                    M[:, xCol] = M[:, xCol]/Ub
                    y = M[:,yCol]
                    x = M[:,xCol] * scaleV + float(xPos)
                    xSum = xSum + x
                    #p1, = plt.plot(x, y, '#1b9e76', alpha = 0.5, lw = 1, mfc = 'none')
                    # p1, = plt.plot(x, y, 'grey', alpha = 0.5, lw = 1, mfc = 'none')
        xMean = xSum / len(samplePlotList)
        # p2, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))
#    return xMean

#def plotBaseline(scaleV, xCol, yCol, para = 'U'):
#    for xPos in xPosList:
        filename_bl = 'pehill_base-observation/postProcessing/sets/3000.000000' + \
                  '/line_x' + str(xPos) + '_' + para + '.xy'
        # x = plotProfile(filename_bl, xPos, scaleV, xCol, yCol, '-', 'darkred', (12,3), para)
        bl = np.loadtxt(filename_bl, comments = '%')
    # if para == 'U':
        bl[:, xCol] = bl[:, xCol]/Ub
    # elif para == 'Tau':
        # M[:, xCol] = M[:, xCol]/Ub/Ub
        bl_y = bl[:,yCol]
        bl_x = bl[:,xCol] * scaleV + float(xPos)  
    # return x

# def plotTruth(scaleV, xCol, yCol, para = 'U'):
    # for xPos in xPosList:
        # if para == 'U':
        filename_dns = 'sets/X/x' + str(xPos) + '.0_U2.xy'
        # else:
            # filename = 'sets/X/dnsTau' + str(xPos)
        # x = plotProfile_dns(filename_dns, xPos, scaleV, xCol, yCol, '-', 'black', 'None', para)
        DNS = np.loadtxt(filename_dns, comments = '%')
        DNS[:, xCol] = DNS[:, xCol]
        y_dns = DNS[:,yCol]
        x_dns = DNS[:,xCol] * scaleV + float(xPos)
        f = interp1d(y, x)
        dns_new = f(bl_y)
        dns_all = np.concatenate((dns_all, dns_new))
        bl_all = np.concatenate((bl_all, bl_x))
        po_all = np.concatenate((po_all, xMean))
    # import pdb; pdb.set_trace()
    L2error_posterior = np.linalg.norm(po_all - dns_all)
    L2error_prior = np.linalg.norm(bl_all - dns_all)
    L1error_posterior = np.linalg.norm(po_all - dns_all, ord=1)
    L1error_prior = np.linalg.norm(bl_all - dns_all, ord=1)
    return L2error_posterior, L2error_prior, L1error_posterior, L1error_prior

# def plotObs(scaleV, xCol, para = 'U'):
#     if para == 'U':
#         Uobs = np.loadtxt('pehill_base/observationData/obsVPlot')
#         Uobs[:, xCol+2:] = Uobs[:, xCol+2:]/Ub 
#         xPos = Uobs[:, 0]
#         x = Uobs[:, xCol+2]*scaleV + xPos
#         y = Uobs[:, 1]
#         p5, = plt.plot(x, y, 'kx', markersize=8, markeredgewidth=2.5)
#     elif para == 'Tau':
#         Tauobs = np.loadtxt('pehill_base/observationData/obsTauPlot')
#         Tauobs[:, xCol+2:] = Tauobs[:, xCol+2:]/Ub/Ub 
#         xPos = Tauobs[:, 0]
#         x = Tauobs[:, xCol+2]*scaleV + xPos
#         y = Tauobs[:, 1]
#         p5, = plt.plot(x, y, 'kx', markersize=8, markeredgewidth=2.5)    
#     return p5

# def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
#     M = np.loadtxt(filename, comments = '%')
#     if para == 'U':
#         M[:, xCol] = M[:, xCol]/Ub
#     elif para == 'Tau':
#         M[:, xCol] = M[:, xCol]/Ub/Ub
#     y = M[:,yCol]
#     x = M[:,xCol] * scaleV + float(xPos)    
#     # if dash == 'None':
#     #     p, = plt.plot(x, y, ls, color = cl, lw=2)
#     # else:
#     #     p, = plt.plot(x, y, ls, color = cl, lw=2, dashes=dash)
#     return x

# def plotProfile_dns(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
#     M = np.loadtxt(filename, comments = '%')
#     if para == 'U':
#         M[:, xCol] = M[:, xCol]
#     elif para == 'Tau':
#         M[:, xCol] = M[:, xCol]
#     y = M[:,yCol]
#     x = M[:,xCol] * scaleV + float(xPos)    
#     # if dash == 'None':
#     #     p, = plt.plot(x, y, ls, color = cl, lw=2)
#     # else:
#     #     p, = plt.plot(x, y, ls, color = cl, lw=2, dashes=dash)
#     return x

def plotProfiles(scaleV, xCol, yCol, para = 'U'):
    # p5 = plotObs(scaleV, xCol, para)
    # return p1, p2, p3, p4, p5
    return p1, p2, p4

def main(iShow=False):

    #plot Ux
    if plotUx == 'True':
        scaleV = 2
        # plt.figure();
        # ax1=plt.subplot(111)
        # plotDomain()
        L2error_posterior, L2error_prior, L1error_posterior, L1error_prior = plotSamples(scaleV, 1, 0, 'U')
        # print "Plotting Ux ........."
        error = np.array([L2error_posterior, L2error_prior, L1error_posterior, L1error_prior])
        np.savetxt('data_misfit.dat',error)
    #     plt.axis([-0.5, 11, 0, 3.05])
    #     plt.ylabel("$y/H$")
    #     plt.xlabel(r'$x/H;\quad$'+str(scaleV)+'$ U_x/U_b+x/H$')
    #     ax1.set_aspect(1.3)
    #     rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    #     fig = plt.gcf()
    #     fig.set_size_inches(8, 3.7)
    #     lg = plt.legend([p1,p2,p3,p4,p5],["samples","sample mean","baseline","DNS","observations"],loc = 7,bbox_to_anchor=(1.0, 1.2),ncol=3)
    #     #lg.draw_frame(False)
    #     matplotlib.rcParams.update({'font.size':12})
    #     plt.savefig("figures/Ux-pehill-"+caseName + '_t_' + timeDir + ".pdf")

    # #plot Tau 
    # scaleVList = [20,20,20,20,20,20,20]
    # if plotTau == 'True':
    #     TauLabels = [r'$\tau_{xx}$',r'$\tau_{xy}$', r'$\tau_{xz}$', r'$\tau_{yy}$', r'$\tau_{yz}$', r'$\tau_{zz}$', r'$k$']
    #     for iter in range(1,8):
    #         print "Plotting Tau_" + str(iter) + "........."
    #         scaleV = scaleVList[iter-1] 
    #         plt.figure();
    #         ax1=plt.subplot(111)
    #         plotDomain()
    #         p1, p2, p3, p4, p5 = plotProfiles(scaleV, iter-1, 0, 'Tau')
    #         plt.axis([-0.5, 9.5, 0, 3.05])
    #         plt.ylabel("$y/H$")
    #         plt.xlabel(r'$x/H;\quad  $'+str(scaleVList[iter-1])+TauLabels[iter-1]+r'$/U_b^2+x/H $')
    #         ax1.set_aspect(1.3)
    #         fig = plt.gcf()
    #         fig.set_size_inches(8, 3.7)
    #         #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 7)
    #         #lg.draw_frame(False)
    #         matplotlib.rcParams.update({'font.size':12})
    #         plt.savefig("figures/tau" + str(iter) + "_" + caseName + '_t_' + timeDir + ".pdf")

    # if iShow:
    #     pass

if __name__ == "__main__":
    Ub = 0.028
    paramDict = readInputData('plotInfo.in')
    timeDir = paramDict['timeDir']
    timeDir0 = paramDict['timeDir0']
    plotUx = paramDict['plotUx']
    plotUy = paramDict['plotUy']
    plotTau = paramDict['plotTau']
    plotDeltaTau = paramDict['plotDeltaTau']
    caseName = paramDict['caseName']   
    paramDict_main = readInputData('MainInput.in')
    Ns = int(paramDict_main['Ns'])
    xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
    samplePlotList = range(1,Ns+1)
    main()
    # timeDir = timeDir0
    # main()
