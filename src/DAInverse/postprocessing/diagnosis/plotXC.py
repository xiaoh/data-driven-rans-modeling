#!/aoe/bin/python27

# description        :Plot Parameter convergence history

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Feb.24, 2016
# revision           :Mar.13, 2016
####################################################################################################
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import ast
import pdb
from utilities import readInputData, extractListFromDict


def plotXC(debugDir, finalDAStep, Npara, Ns, figName='figName'):
    """ 

    Arg: 
    debugDir: path where HX and Obs stored
    
    Regurn: 
    None    
    """
    XCM = np.zeros([finalDAStep, Npara, Ns])
    for i in np.arange(finalDAStep):
        XCDir = debugDir + 'XC_' + str(i+1)
        XCEnsemble = np.loadtxt(XCDir)
        XCM[i, :, :] = XCEnsemble
             
    xLabel = 'Iteration steps'
    yLabels = 'State'
    # Plot Ux
    fig = plt.figure()
    XCmean = np.zeros(finalDAStep)
    XCEtaTruth = [1.5, 1.0]
    truthString = np.arange(finalDAStep)
    
    
    for j in range(Npara):
        flag = 1
        XCmean = np.mean(XCM[:, j, :], axis=1) 
        for i in range(Ns):
            if (flag == 1):      
                plt.plot(XCM[:, j, i], '#1b9e76', alpha = 0.5, lw = 1, 
                label = 'Samples', markevery=1, mfc='none')        
            else:
                plt.plot(XCM[:, j, i], '#1b9e76',alpha = 0.5, lw = 1,
                markevery=1, mfc='none')
            flag = 0                
        #plt.plot(XCmean, 'r', lw = 2, label = 'Sample Mean', markevery = 1)        
        plt.xlabel(xLabel,fontsize = 18)
        if synFlag:
            plt.plot(truthString, truthOmegaVec[0, j]+np.zeros(finalDAStep), 'k--', label = 'Truth', lw = 2)
        
        if j < Nm_xi:
            plt.ylabel(r'$\omega_{\xi}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)           
            plt.savefig('figures/paraXi-'+str(j+1)+'.pdf')            
            plt.clf()
            
        elif j>= Nm_xi and (j < (Nm_xi + Nm_eta)):
            plt.ylabel(r'$\omega_{\eta}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)           
            plt.savefig('figures/paraEta-'+str(j+1-Nm_xi)+'.pdf') 
            plt.clf()
                    
        elif (j>= (Nm_xi + Nm_eta)) and (j < (Nm_xi + Nm_eta + Nm_k)):           
            plt.ylabel(r'$\omega_{k}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)             
            plt.savefig('figures/paraK-'+str(j+1-Nm_xi-Nm_eta)+'.pdf')
            plt.clf()
                    
        elif (j >= (Nm_xi + Nm_eta + Nm_k)) and (j < (Nm_xi + Nm_eta + Nm_k + Nm_VA)):
            plt.ylabel(r'$\omega_{\varphi_1}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)             
            plt.savefig('figures/paraVA-'+str(j+1-Nm_xi-Nm_eta-Nm_k)+'.pdf')
            plt.clf()        
        elif (j >= (Nm_xi + Nm_eta + Nm_k + Nm_VA)) and (j < (Nm_xi + Nm_eta + Nm_k + Nm_VA + Nm_VB)):
            plt.ylabel(r'$\omega_{\varphi_2}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)             
            plt.savefig('figures/paraVB-'+str(j+1-Nm_xi-Nm_eta-Nm_k-Nm_VA)+'.pdf')        
        else:
            plt.ylabel(r'$\omega_{\varphi_3}^{, ('+str(j+1)+')}$', fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.legend(fontsize = 18)
            lg = plt.legend(loc = 4)             
            plt.savefig('figures/paraVC-'+str(j+1-Nm_xi-Nm_eta-Nm_k-Nm_VA-Nm_VB)+'.pdf')
                      

if __name__ == '__main__':
    
    figName = './'    
    debugDir = './debugData/'
    
    paramDict_plot = readInputData('plotInfo.in')
    paramDict_DA = readInputData('MainInput.in')
    paramDict_model = readInputData('forwardModelInput.in')

    synFlag = ast.literal_eval(paramDict_plot['synFlag'])    
    nModes = int(paramDict_model['Nmode'])
    Nmodelist = np.ones(6) * nModes
    Ns = int(paramDict_DA['Ns'])
    Npara = 6 * nModes
    if synFlag:
        paramDict = readInputData('./'+ 'pehill' +'_truth/psuTruthInfo.in')
        # read coefficients to be inferred (omega for deltaXi, deltaEta, deltak, deltaVA, deltaVB, and deltaVC)
        omegaXi = extractListFromDict(paramDict, 'omegaXiVec')
        omegaXiVec = np.array([[float(pn) for pn in omegaXi]])
        print "Reading true omega vector for delta Xi = ", omegaXiVec
        omegaEta = extractListFromDict(paramDict, 'omegaEtaVec')
        omegaEtaVec = np.array([[float(pn) for pn in omegaEta]])
        print "Reading true omega vector for delta Eta = ", omegaEtaVec
        omegak = extractListFromDict(paramDict, 'omegakVec')
        omegakVec = np.array([[float(pn) for pn in omegak]])
        print "Reading true omega vector for delta k = ", omegakVec
        omegaVA = extractListFromDict(paramDict, 'omegaVAVec')
        omegaVAVec = np.array([[float(pn) for pn in omegaVA]])
        print "Reading true omega vector for delta VA = ", omegaVAVec
        omegaVB = extractListFromDict(paramDict, 'omegaVBVec')
        omegaVBVec = np.array([[float(pn) for pn in omegaVB]])
        print "Reading true omega vector for delta VB = ", omegaVBVec
        omegaVC = extractListFromDict(paramDict, 'omegaVCVec')
        omegaVCVec = np.array([[float(pn) for pn in omegaVC]])
        print "Reading true omega vector for delta VC = ", omegaVCVec        
        truthOmegaVec = np.hstack((omegaXiVec, omegaEtaVec, omegakVec, omegaVAVec, omegaVBVec, omegaVCVec))

    Nm_xi = Nmodelist[0]
    Nm_eta = Nmodelist[1]
    Nm_k = Nmodelist[2]
    Nm_VA = Nmodelist[3]
    Nm_VB = Nmodelist[4]
    Nm_VC = Nmodelist[5]
    
    finalDAStep = float(paramDict_plot['finalEnKFStep'])
    
    plotXC(debugDir, finalDAStep, Npara, Ns, figName)

    
