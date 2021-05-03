#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb
from numpy import linalg as LA
from utilities import readInputData, extractListFromDict


def plotIterPost():
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    
    ''
    paramDict = readInputData('MainInput.in')
    Tend = float(paramDict['Tend'])
    DAInterval = float(paramDict['DAInterval'])
    
    NIter = Tend / DAInterval
    dataDir = './debugData/'
    misfitX_L1 = []; 
    misfitX_L2 = []; 
    misfitX_Linf = [];
    sigmaRX_L2 = [];
    sigmaRX_inf = [];    
    sigmaHX_L2 = [];
     
    paramDict2 = readInputData('forwardModelInput.in')
    caseName = paramDict2['caseName']
    obsSigmaFixedVec = extractListFromDict(paramDict2, 'ObsSigmaFixedVec')
    obsSigmaFixed = np.array([float(pn) for pn in obsSigmaFixedVec])
    
    obsRelCoeffVec = extractListFromDict(paramDict2, 'ObsRelCoeffVec')
    obsRelCoeff = np.array([float(pn) for pn in obsRelCoeffVec])
    
    obsRmaCoeff = float(paramDict2['ObsRmaCoeff'])
    ObsDirReal = caseName + '/observationData/obsVelocity'
    for iterStep in np.arange(NIter):
        
        print iterStep
        HXDir = dataDir + 'HX_' + str(iterStep+1) + '.txt'
        ObsDir = dataDir + 'Obs_' + str(iterStep+1)
        
        
        HXEnsemble = np.loadtxt(HXDir)
        ObsEnsembleV = np.loadtxt(ObsDirReal)
        
        Nobs, Ns = HXEnsemble.shape
        Nnorm = Nobs * Ns 
        ObsEnsemble = ObsEnsembleV
        for i in np.arange(Ns-1):
            ObsEnsemble = np.vstack((ObsEnsemble, ObsEnsembleV))
        ObsEnsemble = ObsEnsemble.T
        misfitTemp = abs(HXEnsemble - ObsEnsemble)
        
        
        misfitNormL1 = np.sum(misfitTemp) / Nnorm
        misfitNormL2 = np.sqrt(np.sum(misfitTemp**2)) / Nnorm
        misfitNormLinf = LA.norm(misfitTemp, np.inf)
        
        misfitX_L1.append(misfitNormL1)
        misfitX_L2.append(misfitNormL2)
        misfitX_Linf.append(misfitNormLinf)
        
        #TODO temp here
        #misfitTemp = np.mean(HXEnsemble, axis=1) - ObsEnsembleV
        #misfitNormL1 = LA.norm(misfitTemp, 1) / Nobs
        #misfitNormL2 = LA.norm(misfitTemp) / Nobs
        #misfitX_L1.append(misfitNormL1)
        #misfitX_L2.append(misfitNormL2)        
        #pdb.set_trace()
        # Obs
        obsU = ObsEnsemble[0::3, 0]
        obsV = ObsEnsemble[1::3, 0]
        obsW = ObsEnsemble[2::3, 0]
        # Observation std
        sigmaObsUVec = obsSigmaFixed[0] * np.ones(Nobs/3) + obsRelCoeff[0] * abs(obsU)
        sigmaObsVVec = obsSigmaFixed[1] * np.ones(Nobs/3) + obsRelCoeff[1] * abs(obsV)
        sigmaObsWVec = obsSigmaFixed[2] * np.ones(Nobs/3) + obsRelCoeff[2] * abs(obsW)
        
        
        # Observation std in R matrix (observation std * coefficient)
        sigmaRUVec = np.sqrt(obsRmaCoeff) * sigmaObsUVec
        sigmaRVVec = np.sqrt(obsRmaCoeff) * sigmaObsVVec
        sigmaRWVec = np.sqrt(obsRmaCoeff) * sigmaObsWVec

        sigmaObs = np.hstack((sigmaRUVec, sigmaRVVec, sigmaRWVec))
        sigmaRNormL2 = LA.norm(sigmaObs) / Nobs
        sigmaRX_L2.append(sigmaRNormL2)
        sigmaRInf = max(sigmaObs)
        sigmaRX_inf.append(sigmaRInf)

        # ensemble std   
        sigmaU = np.std(HXEnsemble[0::3, :], axis=1)
        sigmaV = np.std(HXEnsemble[1::3, :], axis=1)
        sigmaW = np.std(HXEnsemble[2::3, :], axis=1)
        sigmaHX = np.hstack((sigmaU, sigmaV, sigmaW))
        sigmaHXNormL2 = LA.norm(sigmaHX) / Nobs
        sigmaHX_L2.append(sigmaHXNormL2)    
        #pdb.set_trace()
    #pdb.set_trace()
    xLabel = 'Iterations'
    
    
    
    fig = plt.figure()
    plt.hold(True)
    plt.plot(sigmaHX_L2, 'yo-', lw = 2, label ='std of ensemble', markevery=1)
    plt.plot(sigmaRX_L2, 'b*-', lw = 2, label ='std of Obs error', markevery=1)                
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
    plt.plot(misfitX_L1, 'yo-', lw = 2, label ='misfit L1', markevery=1)
    plt.plot(misfitX_L2, 'ro-', lw = 2, label ='misfit L2', markevery=1)
    plt.plot(misfitX_Linf, 'kx-', lw = 2, label ='misfit inf', markevery=1)
    plt.plot(sigmaRX_L2, 'b*-', lw = 2, label ='std of Obs error', markevery=1)                
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
    plt.plot(misfitX_L1, 'yo-', lw = 2, label ='misfit L1', markevery=1)               
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
    plt.plot(misfitX_L2, 'ro-', lw = 2, label ='misfit L2', markevery=1)          
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
    plt.plot(misfitX_Linf, 'kx-', lw = 2, label ='misfit inf', markevery=1)
    plt.plot(sigmaRX_inf, 'b*-', lw = 2, label ='inf of Obs error', markevery=1)                
    plt.hold(False)
    plt.xlabel(xLabel,fontsize = 16)
    plt.ylabel('misfit',fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    lg = plt.legend(loc = 2)
    lg.draw_frame(False)
    plt.savefig('./iterConvergenceInf.eps')  
    
if __name__ == '__main__':

    plotIterPost() 
