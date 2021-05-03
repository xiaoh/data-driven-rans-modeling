#!/aoe/bin/python27

# description        :Generate synthetic truth with given KL coefficients (omega).

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's G-roup
# date               :Dec.03, 2015
# revision           :Mar.09, 2016

import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb
from utilities import readInputData, extractListFromDict
import foamFileOperation as foam

NVec = 3
NTensor = 6

class plotEnsemble:
    
    def __init__(self):
        pass
    
     
    def _splitUTauObs(self, XObs):
        """

        :param XU:
        :return:
        """
        UObs = XObs[0 : (NVec * NcellObs)]
        TauObs = XObs[(NVec * NcellObs) :]

        return UObs, TauObs

    def plotPriorVar(self, debugDir, DAStep, obsSigmaFixed, obsRelCoeff, obsRmaCoeff, figName='figName'):
        """ 

        Arg: 
        debugDir: path where HX and Obs stored
        
        Regurn: 
        None    
        """
        baseFolder = problem 
        baseUFile = baseFolder + '/0/U'
        base = foam.readVelocityFromFile(baseUFile)        
        U1_b = base[0::3]
        U2_b = base[1::3]
        U3_b = base[2::3]
        baseTauFile = baseFolder + '/0/Tau'
        base = foam.readTurbStressFromFile(baseTauFile)
        R1_b = base[:, 0]
        R2_b = base[:, 1]
        R3_b = base[:, 2]
        R4_b = base[:, 3]
        R5_b = base[:, 4]
        R6_b = base[:, 5]
        
        HXDir = debugDir + 'HX_' + str(DAStep) + '.txt'
        ObsRDir = debugDir + 'Obs_' + str(DAStep) 
        XupdateDir = debugDir + 'updateX_' + str(DAStep)    
              
        HXEnsemble = np.loadtxt(HXDir)    
        HXEnsemble_update = np.loadtxt(XupdateDir)
        if TauOnFlag:
            HXEnsemble_update = HXEnsemble_update[0:(NVec*Ncell+NTensor*Ncell), :]
        else:
            HXEnsemble_update = HXEnsemble_update[0:NVec*Ncell, :]
        ObsEnsembleR = np.loadtxt(ObsRDir)
        ObsR = ObsEnsembleR[:, 0]
        
        UHXEnsemble = np.zeros([Ns, (NVec * NcellObs)])    
        TauHXEnsemble = np.zeros([Ns, (NTensor * NcellObs)])

        UHXEnsemble_update = np.zeros([Ns, (NVec * NcellObs)])    
        TauHXEnsemble_update = np.zeros([Ns, (NTensor * NcellObs)])
        
        if TauOnFlag:
            UObs, TauObs = self._splitUTauObs(ObsR)
            for icase in range(Ns):
                UHXEnsemble[icase, :], TauHXEnsemble[icase, :] = self._splitUTauObs(HXEnsemble[:, icase])
                if NcellObs == Ncell:
                    UHXEnsemble_update[icase, :], TauHXEnsemble_update[icase, :] = self._splitUTauObs(HXEnsemble_update[:, icase])
            
            if NcellObs == Ncell:
                UHX1_u = UHXEnsemble_update[:, 0::3]
                UHX2_u = UHXEnsemble_update[:, 1::3]
                UHX3_u = UHXEnsemble_update[:, 2::3]
                TauHX1_u = TauHXEnsemble_update[:, 0::6]
                TauHX2_u = TauHXEnsemble_update[:, 1::6]
                TauHX3_u = TauHXEnsemble_update[:, 2::6]
                TauHX4_u = TauHXEnsemble_update[:, 3::6]
                TauHX5_u = TauHXEnsemble_update[:, 4::6]
                TauHX6_u = TauHXEnsemble_update[:, 5::6]                                    
            
            UHX1 = UHXEnsemble[:, 0::3]
            UHX2 = UHXEnsemble[:, 1::3]
            UHX3 = UHXEnsemble[:, 2::3]
            TauHX1 = TauHXEnsemble[:, 0::6]
            TauHX2 = TauHXEnsemble[:, 1::6]
            TauHX3 = TauHXEnsemble[:, 2::6]
            TauHX4 = TauHXEnsemble[:, 3::6]
            TauHX5 = TauHXEnsemble[:, 4::6]
            TauHX6 = TauHXEnsemble[:, 5::6]
            
            obsU1 = UObs[0::3]
            obsU2 = UObs[1::3]
            obsU3 = UObs[2::3]            
            obsTau1 = TauObs[0::6]
            obsTau2 = TauObs[1::6]
            obsTau3 = TauObs[2::6]
            obsTau4 = TauObs[3::6]
            obsTau5 = TauObs[4::6]
            obsTau6 = TauObs[5::6]
            
        else:
            UHXEnsemble = HXEnsemble 
            UHX1 = UHXEnsemble[0::3, :].T
            UHX2 = UHXEnsemble[1::3, :].T
            UHX3 = UHXEnsemble[2::3, :].T     
            if NcellObs == Ncell:
                UHX1_u = HXEnsemble_update[0::3, :].T
                UHX2_u = HXEnsemble_update[1::3, :].T
                UHX3_u = HXEnsemble_update[2::3, :].T   
            obsU1 = ObsR[0::3]
            obsU2 = ObsR[1::3]
            obsU3 = ObsR[2::3]
   
        #pdb.set_trace() 
        xLabel = 'cell idx'               
        flag = 1
        
        # Plot U1
        fig = plt.figure()
        plt.clf()
        plt.hold(True)
        #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
        for i in range(Ns):
            if(flag == 1):       
                plt.plot(UHX1[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                flag = 0
            else:
                plt.plot(UHX1[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
            flag = 0                
        flag = 1
        if NcellObs == Ncell:        
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(UHX1_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(UHX1_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
                flag = 0
        plt.plot(obsU1, '-k*', lw = 2, label ='observation', markevery=1)
        plt.plot(U1_b, '-ro', lw = 2, label ='baseline', markevery=1)  
        plt.hold(False)
        plt.xlabel(xLabel,fontsize = 16)
        plt.ylabel('Ux',fontsize = 16)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        lg = plt.legend(loc = 2)
        lg.draw_frame(False)
        plt.savefig(figName+'-U1.eps')

        # Plot U1
        fig = plt.figure()
        plt.clf()        
        plt.hold(True)
        #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)  
        for i in range(Ns):
            if(flag == 1):       
                plt.plot(UHX2[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                flag = 0
            else:
                plt.plot(UHX2[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
            flag = 0                
        flag = 1
        if NcellObs == Ncell:        
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(UHX2_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(UHX2_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
        plt.plot(obsU2, '-k*', lw = 2, label ='observation', markevery=1)  
        plt.plot(U2_b, '-ro', lw = 2, label ='baseline', markevery=1)
        plt.hold(False)
        plt.xlabel(xLabel,fontsize = 16)
        plt.ylabel('Ux',fontsize = 16)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        lg = plt.legend(loc = 2)
        lg.draw_frame(False)
        plt.savefig(figName+'-U2.eps')

        # Plot U3
        fig = plt.figure()
        plt.clf()        
        plt.hold(True)
        #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)   
        for i in range(Ns):
            if(flag == 1):       
                plt.plot(UHX3[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                flag = 0
            else:
                plt.plot(UHX3[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
            flag = 0                
        flag = 1
        if NcellObs == Ncell:        
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(UHX3_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(UHX3_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
        plt.plot(obsU3, '-k*', lw = 2, label ='observation', markevery=1) 
        plt.plot(U3_b, '-ro', lw = 2, label ='baseline', markevery=1)
        plt.hold(False)
        plt.xlabel(xLabel,fontsize = 16)
        plt.ylabel('Ux',fontsize = 16)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        lg = plt.legend(loc = 2)
        lg.draw_frame(False)
        plt.savefig(figName+'-U3.eps')   

        if TauOnFlag:
            # Plot Tau1
            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
               
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX1[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX1[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0
                #pdb.set_trace()                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX1_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX1_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
            plt.plot(obsTau1, '-k*', lw = 2, label ='observation', markevery=1) 
            plt.plot(R1_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R1',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R1.eps')     
              

            # Plot Tau2
            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
             
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX2[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX2[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX2_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX2_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
            plt.plot(obsTau2, '-k*', lw = 2, label ='observation', markevery=1)   
            plt.plot(R2_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R2',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R2.eps')  


            # Plot Tau3

            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
            
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX3[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX3[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX3_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX3_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')           
            plt.plot(obsTau3, '-k*', lw = 2, label ='observation', markevery=1)    
            plt.plot(R3_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R3',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R3.eps')

            # Plot Tau4
            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
               
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX4[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX4[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX4_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX4_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')            
            plt.plot(obsTau4, '-k*', lw = 2, label ='observation', markevery=1) 
            plt.plot(R4_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R4',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R4.eps')

            # Plot Tau5
            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            #plt.plot(obsU, '-ro', lw = 2, label ='truth', markevery=1)
               
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX5[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX5[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX5_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX5_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
            plt.plot(obsTau5, 'g*', lw = 2, label ='observation', markevery=1) 
            plt.plot(R5_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R4',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R5.eps')
         

            # Plot Tau6
            fig = plt.figure()
            plt.clf()
            plt.hold(True)
            for i in range(Ns):
                if(flag == 1):       
                    plt.plot(TauHX6[i, :], '-y.', lw = 0.5, alpha=0.2, label = 'state ensemble', markevery=1, mfc='none')
                    flag = 0
                else:
                    plt.plot(TauHX6[i, :], '-y.', lw = 0.5, alpha=0.2,  markevery=1, mfc='none')
                flag = 0                
            flag = 1
            if NcellObs == Ncell:        
                for i in range(Ns):
                    if(flag == 1):       
                        plt.plot(TauHX6_u[i, :], '-g.', lw = 0.5, alpha=0.1, label = 'state ensemble (update)', markevery=1, mfc='none')
                        flag = 0
                    else:
                        plt.plot(TauHX6_u[i, :], '-g.', lw = 0.5, alpha=0.1,  markevery=1, mfc='none')
            plt.plot(obsTau6, 'g*', lw = 2, label ='observation', markevery=1) 
            plt.plot(R6_b, '-ro', lw = 2, label ='baseline', markevery=1)
            plt.hold(False)
            plt.xlabel(xLabel,fontsize = 16)
            plt.ylabel('R4',fontsize = 16)
            plt.xticks(fontsize = 16)
            plt.yticks(fontsize = 16)
            lg = plt.legend(loc = 2)
            lg.draw_frame(False)
            plt.savefig(figName+'-R6.eps')
    
if __name__ == '__main__':
    debugDir = './debugData/'
    
    paramDict_mode = readInputData('forwardModelInput.in')
    paramDict_plot = readInputData('plotInfo.in')
    paramDict_main = readInputData('MainInput.in')
    
    Ns = int(paramDict_main['Ns'])
        
    EnKFStep = float(paramDict_plot['EnKFStep'])
    EnKFStep = str(EnKFStep)
    
    TauOnFlag = eval(paramDict_mode['TauOnFlag'])
    NcellObs = int(paramDict_mode['NcellObs'])
    Ncell = int(paramDict_mode['Ncell'])
    problem = paramDict_mode['caseName']
  
    # observation er-ror
    obsSigmaFixedVec = extractListFromDict(paramDict_mode, 'ObsSigmaFixedVec')    
    obsSigmaFixed = np.array([float(pn) for pn in obsSigmaFixedVec])    
    obsRelCoeffVec = extractListFromDict(paramDict_mode, 'ObsRelCoeffVec')
    obsRelCoeff = np.array([float(pn) for pn in obsRelCoeffVec])    
    obsRmaCoeff = float(paramDict_mode['ObsRmaCoeff'])    
            
    figName = 'figures/HXObs-DA_'+str(EnKFStep)
    # plot ensemble variance at observed locations
    plotE = plotEnsemble()
    plotE.plotPriorVar(debugDir, EnKFStep, obsSigmaFixed, obsRelCoeff, obsRmaCoeff, figName)   
        
