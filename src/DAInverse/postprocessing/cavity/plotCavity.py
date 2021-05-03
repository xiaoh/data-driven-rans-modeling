# Plot cavity field
###############################################################################

#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb

def plotTopU(topUDir, truthDir, figName='figName'):
    """ 

    Arg: 
    topUDir: path of file sampled topline U
    truthDir: path of file sampled true topline U
    
    Regurn: 
    None    
    """
    filenames = [topUDir, truthDir]
    allLines = ['b-o','r-^']
    allWidth = [2,2]
    allDashes = [(12,0),(12,0)]
    allLabels = ['DA case', 'truth']
    xLabel = '$y$'
    yLabels = ['$U_x$', '$U_y$']

    # Plot Ux
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        x = M[:,0]
        Ux = M[:,3]
        plt.plot(x, Ux, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)

    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[0],fontsize = 16);

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)

    lg = plt.legend(loc = 2)
    lg.draw_frame(False)

    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.15)
    plt.savefig(figName+'-topUx.eps')
    # Plot Uy
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        x = M[:,0]
        Uy = M[:,4]
        plt.plot(x, Uy, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)

    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[1],fontsize = 16);

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)

    lg = plt.legend(loc = 1)
    lg.draw_frame(False)

    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.15)
    plt.savefig(figName+'-topUy.eps')

def plotCenterU(centerUDir, truthDir, figName='figName'):
    """ 

    Arg: 
    centerUDir: path of file sampled centerline U
    truthDir: path of file sampled true centerline U
    
    Regurn: 
    None    
    """    
    filenames = [centerUDir, truthDir]
    allLines = ['b-o','r-^']
    allWidth = [2,2]
    allDashes = [(12,0),(12,0)]
    allLabels = ['DA case', 'truth']
    xLabel = '$y$'
    yLabels = ['$U_x$', '$U_y$']

    # Plot Ux
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        y = M[:,1]
        Ux = M[:,3]
        plt.plot(y, Ux, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)

    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[0],fontsize = 16);

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)

    lg = plt.legend(loc = 2)
    lg.draw_frame(False)

    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.15)
    plt.savefig(figName+'-centerUx.eps')
    # Plot Uy
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        y = M[:,1]
        Uy = M[:,4]
        plt.plot(y, Uy, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)

    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[1],fontsize = 16);

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)

    lg = plt.legend(loc = 2)
    lg.draw_frame(False)

    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.15)
    plt.savefig(figName+'-ceterUy.eps') 
