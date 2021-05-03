#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb
from numpy import linalg as LA

def plotIter():
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    
    ''
    xLabel = 'Iterations'
    misfit_L1 = np.loadtxt('misfit_L1.txt')
    misfit_L2 = np.loadtxt('misfit_L2.txt')
    misfit_inf = np.loadtxt('misfit_inf.txt')
    sigmaHX = np.loadtxt('sigmaHX.txt')
    obsSigma = np.loadtxt('obsSigma.txt')
    obsSigmaInf = np.loadtxt('obsSigmaInf.txt')
    
    
    fig = plt.figure()
    plt.hold(True)
    plt.plot(sigmaHX, 'yo-', lw = 2, label ='std of ensemble', markevery=1)
    plt.plot(obsSigma, 'b*-', lw = 2, label ='std of Obs error', markevery=1)                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('standard deviation',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergenceSTD.eps')
        
    
    fig = plt.figure()
    plt.hold(True)
    plt.plot(misfit_L1, 'yo-', lw = 2, label ='misfit L1', markevery=1)
    plt.plot(misfit_L2, 'ro-', lw = 2, label ='misfit L2', markevery=1)
    plt.plot(misfit_inf, 'kx-', lw = 2, label ='misfit inf', markevery=1)
    plt.plot(obsSigma, 'b*-', lw = 2, label ='std of Obs error', markevery=1)                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('misfit',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergence.eps')
    
    fig = plt.figure()
    plt.hold(True)
    plt.plot(misfit_L1, 'yo-', lw = 2, label ='misfit L1', markevery=1)                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('misfit',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergenceL1.eps') 

    fig = plt.figure()
    plt.hold(True)
    plt.plot(misfit_L2, 'ro-', lw = 2, label ='misfit L2', markevery=1)              
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('misfit',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergenceL2.eps')   
    
    fig = plt.figure()
    plt.hold(True)
    plt.plot(misfit_inf, 'kx-', lw = 2, label ='misfit inf', markevery=1)
    plt.plot(obsSigmaInf, 'b*-', lw = 2, label ='infinite norm of Obs error', markevery=1)                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('misfit',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergenceInf.eps')  
    
if __name__ == '__main__':

    plotIter()       
