#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb

if not os.path.exists('figures/'):
    os.mkdir('figures/')

xPosList = ['0', '0p5', '1', '2', '3', '4', '5', '6', '7', '8']

for xPos in xPosList:
    filenames = ['./postProcessing/sets/0/line_x' + xPos + '_U.xy', \
                 './postProcessing/sets/3000.000000/line_x' + xPos + '_U.xy', \
                 './postProcessing/sets/6000.000000/line_x' + xPos + '_U.xy', \
                 './postProcessing/sets/9000.000000/line_x' + xPos + '_U.xy' ]
    allLines = ['b-o','r-^','k-*','g-v']
    allWidth = [2,2,2,2]
    allDashes = [(12,0),(12,0),(12,0),(12,0)]
    allLabels = ['t = 0s', 't = 1000s', 't = 4000s', 't = 5000s']
    xLabel = '$y$'
    yLabels = ['$U_x$', '$U_y$']
    
    # Plot Ux
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        y = M[:,0]
        Ux = M[:,1]
        plt.plot(y, Ux, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)
    
    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[0],fontsize = 16);
    
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    
    lg = plt.legend(loc = 0)
    lg.draw_frame(False)
    
    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.15)
    
    plt.savefig("figures/line" + xPos + "_Ux.pdf")
    
    # Plot Uy
    fig = plt.figure()
    plt.hold(True)
    for i in range(len(filenames)):
        filename = filenames[i]
        M = np.loadtxt(filename)
        y = M[:,0]
        Uy = M[:,2]
        plt.plot(y, Uy, allLines[i], lw = 2, label = allLabels[i], markevery=1,mfc='none')
    plt.hold(False)
    
    plt.xlabel(xLabel,fontsize = 16);
    plt.ylabel(yLabels[1],fontsize = 16);
    
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    
    lg = plt.legend(loc = 0)
    lg.draw_frame(False)
    
    plt.gcf().subplots_adjust(bottom = 0.15)
    plt.gcf().subplots_adjust(left = 0.2)
    
    plt.savefig("figures/line" + xPos + "_Uy.pdf")
