#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import pdb
from utilities import readInputData, extractListFromDict

paramDict = readInputData('diagnosisPlot.in')
debugDir = paramDict['debugDir']
finalDAStep = float(paramDict['finalDAStep'])
N_modesXi = int(paramDict['N_modesXi'])
N_modesEta = int(paramDict['N_modesEta'])
N_modesk = int(paramDict['N_modesk'])

# n_omega is number of inferred parameters ( = len(XCVec) )
# n_uloc number of U locations that wanted to be calculated the correlation.

list_ulocStr = extractListFromDict(paramDict, 'list_uloc')
list_uloc = np.array([int(pn) for pn in list_ulocStr])
n_uloc = len(list_ulocStr)
       
n_omegaXi = int(paramDict['n_omegaXi'])
n_omegaEta = int(paramDict['n_omegaEta'])
n_omegak = int(paramDict['n_omegak'])


def plotCorr(debugDir, DAStep, n_uloc, n_omegaXi, N_modesXi, N_modesEta, 
             N_modesk, figName='figName'):
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    HXDir = debugDir + 'HX_' + str(DAStep) + '.txt'
    XCDir = debugDir + 'XC_' + str(DAStep)

    HXEnsemble = np.loadtxt(HXDir)
    XCEnsemble = np.loadtxt(XCDir)
    XCXi = XCEnsemble[0:N_modesXi]
    XCEta = XCEnsemble[N_modesXi:N_modesXi+N_modesEta]
    XCk = XCEnsemble[N_modesXi+N_modesEta:]
    
    ObsQrange = range(n_uloc*3)
    Urange = HXEnsemble[0::3, :]
    Vrange = HXEnsemble[1::3, :]
    #pdb.set_trace()
    # Xi
    f, axarr = plt.subplots(n_uloc, n_omegaXi)
    for i in range(n_uloc):
        for j in range(n_omegaXi):    
            axarr[i, j].scatter(XCXi[j, :], \
                        Urange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$\xi$ coeff $\omega $" + str(j+1))
            axarr[i, j].set_ylabel("Observed U at Loc " + str(list_uloc[i]))
            if j != 0:
                plt.setp([a.get_yticklabels() for a in axarr[:, j]], visible=False)   
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
   
    fig = plt.gcf()
    fig.set_size_inches(18, 14)
    matplotlib.rcParams.update({'font.size':20})
    plt.savefig("./corrParaWithU_Xi_DA"+str(DAStep)+".pdf")
    
    f, axarr = plt.subplots(n_uloc, n_omegaXi)
    for i in range(n_uloc):
        for j in range(n_omegaXi):    
            axarr[i, j].scatter(XCXi[j, :], \
                        Vrange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$\xi$ coeff $\omega $" + str(j+1))
            axarr[i, j].set_ylabel("Observed V at Loc " + str(list_uloc[i]))                            
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
    fig = plt.gcf()
    fig.set_size_inches(18, 14)                
    #plt.show()
    plt.savefig("./corrParaWithV_Xi_DA"+str(DAStep)+".pdf")

    # Eta
    f, axarr = plt.subplots(n_uloc, n_omegaEta)
    for i in range(n_uloc):
        for j in range(n_omegaEta):    
            axarr[i, j].scatter(XCEta[j, :], \
                        Urange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$\eta$ coeff $\omega $" + str(j+1))
            axarr[i, j].set_ylabel("Observed U at Loc " + str(list_uloc[i]))
            if j != 0:
                plt.setp([a.get_yticklabels() for a in axarr[:, j]], visible=False)   
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
   
    fig = plt.gcf()
    fig.set_size_inches(18, 14)
    matplotlib.rcParams.update({'font.size':25})                
    plt.savefig("./corrParaWithU_Eta_DA"+str(DAStep)+".pdf")
    
    f, axarr = plt.subplots(n_uloc, n_omegaEta)
    for i in range(n_uloc):
        for j in range(n_omegaEta):    
            axarr[i, j].scatter(XCEta[j, :], \
                        Vrange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$\eta$ coeff $\omega $" + str(j+1))
            axarr[i, j].set_ylabel("Observed V at Loc " + str(list_uloc[i]))                            
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
    fig = plt.gcf()
    fig.set_size_inches(18, 14)
    matplotlib.rcParams.update({'font.size':30})                
    #plt.show()
    plt.savefig("./corrParaWithV_Eta_DA"+str(DAStep)+".pdf")

    # k
    f, axarr = plt.subplots(n_uloc, n_omegak)
    for i in range(n_uloc):
        for j in range(n_omegak):    
            axarr[i, j].scatter(XCk[j, :], \
                        Urange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$k$ coeff $\omega $" + str(j+1))
            axarr[i, j].set_ylabel("Observed U at Loc " + str(list_uloc[i]))
            if j != 0:
                plt.setp([a.get_yticklabels() for a in axarr[:, j]], visible=False)   
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
   
    fig = plt.gcf()
    fig.set_size_inches(18, 14)                
    plt.savefig("./corrParaWithU_k_DA"+str(DAStep)+".pdf")
    
    f, axarr = plt.subplots(n_uloc, n_omegak)
    for i in range(n_uloc):
        for j in range(n_omegak):    
            axarr[i, j].scatter(XCk[j, :], \
                        Vrange[list_uloc[i], :], alpha=0.5)
            axarr[i, j].set_xlabel(r"$k$ coeff $\omega $" + str(j))
            axarr[i, j].set_ylabel("Observed V at Loc " + str(list_uloc[i]))                            
        if i != n_uloc-1:
            plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
    fig = plt.gcf()
    fig.set_size_inches(18, 14)                
    #plt.show()
    plt.savefig("./corrParaWithV_k_DA"+str(DAStep)+".pdf")



if __name__ == '__main__':
    plotCorr(debugDir, finalDAStep, n_uloc, n_omegaXi, N_modesXi, N_modesEta, 
            N_modesk)    
