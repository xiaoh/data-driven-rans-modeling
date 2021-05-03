#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb

def plotXC(debugDir, TotalDAStep, Npara, Ns, figName='figName'):
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    XCM = np.zeros([TotalDAStep, Npara, Ns])
    for i in np.arange(TotalDAStep):
        XCDir = debugDir + 'XC_' + str(i+1)
        XCEnsemble = np.loadtxt(XCDir)
        XCM[i, :, :] = XCEnsemble
             
    xLabel = '$y$'
    yLabels = 'State'
    # Plot Ux
    fig = plt.figure()
    XCmean = np.zeros(TotalDAStep)
    for j in range(Npara):
        plt.hold(True)
        flag = 1
        XCmean = np.mean(XCM[:, j, :], axis=1) 
        for i in range(Ns):
            if (flag == 1):      
                plt.plot(XCM[:, j, i], 'g--', lw = 0.5, label = 'parameter ensemble', markevery=1, mfc='none')        
            else:
                plt.plot(XCM[:, j, i], 'g--', lw = 0.5, markevery=1, mfc='none')
            flag = 0    
        plt.plot(XCmean, 'ro', lw = 2, label = 'mean of parameter', markevery = 1, mfc = 'none')        
        plt.hold(False)
        plt.xlabel(xLabel,fontsize = 16)
        plt.ylabel('para'+str(j+1),fontsize = 16)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        lg = plt.legend(loc = 1)       
        lg.draw_frame(False)
        #plt.show()        
        plt.savefig(figName+'para-'+str(j+1)+'.eps')
        plt.clf()
        


if __name__ == '__main__':
    debugDir = '../debugData/'
    TotalDAStep = 30.0
    Npara = 6
    Ns = 60
    figName = '../debugData/'
    plotXC(debugDir, TotalDAStep, Npara, Ns, figName)
