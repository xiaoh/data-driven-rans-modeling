# Plot State X with Observation Obs
###############################################################################

#!/usr/bin/python
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb

def plotHXObs(debugDir, DAStep, figName='figName'):
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    HXDir = debugDir + 'HX_' + str(DAStep) + '.txt'
    ObsDir = debugDir + 'Obs_' + str(DAStep) 
    XupdateDir = debugDir + 'updateX_' + str(DAStep)    
    HXEnsemble = np.loadtxt(HXDir)
    ObsEnsemble = np.loadtxt(ObsDir)
    XupdateEnsemble = np.loadtxt(XupdateDir)
    
    
    obsU = ObsEnsemble[0::3, 0]
    obsV = ObsEnsemble[1::3, 0]
    obsW = ObsEnsemble[2::3, 0]
    
    Nstate, Ns = np.shape(HXEnsemble)
    XupdateEnsemble = XupdateEnsemble[0:Nstate, :]
     
    xLabel = '$y$'

    flag = 1
    # Plot Ux
    fig = plt.figure()
    plt.hold(True)
    
    for i in range(Ns):
        HXU = HXEnsemble[0::3, i]
        if(flag == 1):       
            plt.plot(HXU, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else:
            plt.plot(HXU, 'g--', lw = 0.5,  markevery=1, mfc='none')
        flag = 0
    plt.plot(obsU, 'ro', lw = 2, label ='observation', markevery=1, mfc='none')                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('Ux',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-U.eps')

    flag = 1
    fig = plt.figure()
    plt.hold(True)
    for i in range(Ns):
        HXV = HXEnsemble[1::3, i]       
        if(flag == 1):
            plt.plot(HXV, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else:
            plt.plot(HXV, 'g--', lw = 0.5, markevery=1, mfc='none')
        flag = 0
    plt.plot(obsV, 'ro', lw = 2, label ='observation', markevery=1, mfc='none')             
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('Uy',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-V.eps')

    flag = 1
    fig = plt.figure()
    plt.hold(True)
    for i in range(Ns):
        difU = obsU - HXEnsemble[0::3, i]
        if(flag == 1):       
            plt.plot(difU, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else: 
            plt.plot(difU, 'g--', lw = 0.5,  markevery=1, mfc='none')
        flag = 0     
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('diff Ux',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-difU.eps')
    
    flag = 1
    fig = plt.figure()
    plt.hold(True)
    for i in range(Ns):
        difV = obsV - HXEnsemble[1::3, i]
        if (flag == 1):       
            plt.plot(difV, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else:
            plt.plot(difV, 'g--', lw = 0.5, markevery=1, mfc='none')
        flag = 0     
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('diff Uy',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-difV.eps')    
    
    flag = 1
    fig = plt.figure()
    plt.hold(True)
    for i in range(Ns):
        difU = obsU - XupdateEnsemble[0::3, i]       
        if (flag == 1):
            plt.plot(difU, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else: 
            plt.plot(difU, 'g--', lw = 0.5, markevery=1, mfc='none')
        flag = 0    
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('diff Ux',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-difUupdate.eps')

    flag = 1
    fig = plt.figure()
    plt.hold(True)
    for i in range(Ns):
        difV = obsV - XupdateEnsemble[1::3, i]
        if (flag == 1):       
            plt.plot(difV, 'g--', lw = 0.5, label = 'state ensemble', markevery=1, mfc='none')
        else:
            plt.plot(difV, 'g--', lw = 0.5, markevery=1, mfc='none') 
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('diff Uy',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig(figName+'-difVupdate.eps')
    
     
    
if __name__ == '__main__':
    debugDir = '../debugData/'
    finalDAStep = 5.0
    for i in np.arange(finalDAStep):
        figName = '../debugData/DA-'+str(i+1)+'/HXObs-DA_'+str(i+1)
        plotHXObs(debugDir, i+1, figName)
        
