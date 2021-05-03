#!/aoe/bin/python27

# description        :Contain all functions for postprocessing the result of Physics-based UQ

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Feb.24, 2016
# revision           :Feb.24, 2016
####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.stats as ss
# system, file operation
import pdb
import sys
import ast
import os
import os.path as ospt
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
from scipy import interpolate
# Import local modules
import ReynoldsStressRF as ReRF
from utilities import readInputData, extractListFromDict

def getPhyComponents(EnKFStep):
    """
    Output physical components (c1, c2, c3, xi, eta, k, XC1, XC2, VA, VB, VC)
    """
    c1Fields = np.zeros([Ns, Ncell])
    c2Fields = np.zeros([Ns, Ncell])
    c3Fields = np.zeros([Ns, Ncell])
    
    XC1Fields = np.zeros([Ns, Ncell])
    XC2Fields = np.zeros([Ns, Ncell])

    XiFields = np.zeros([Ns, Ncell])
    EtaFields = np.zeros([Ns, Ncell])
    kFields = np.zeros([Ns, Ncell])
    VAFields = np.zeros([Ns, Ncell])
    VBFields = np.zeros([Ns, Ncell])
    VCFields = np.zeros([Ns, Ncell])     
       

    sampleIndex = np.linspace(0, Ns-1, Ns)
    patch = 10
    DAstep = "%.1f"%EnKFStep
    for idx in sampleIndex:    
        if np.remainder(idx, patch) == 0.0:
            print "processing " + str(idx+1) + '/' + str(Ns) + ' samples'
        idx = int(idx)
        TauNew = np.loadtxt(dataFolder+'/DA-'+DAstep+'/Tau/Tau_sample-'+str(idx))
        k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(TauNew) # collapse time = 1.005s (3000 cells)
        X = mapTau._C2X(C)  # collapse time = 0.02s (3000 cells)
        RS = mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
        VA, VB, VC = mapTau.getThetaVABC(TauNew) # collapse time = 1.005s (3000 cells)

        c1Fields[idx, :] = C[:, 0]
        c2Fields[idx, :] = C[:, 1]
        c3Fields[idx, :] = C[:, 2]
        
        XC1Fields[idx, :] = X[:, 0]
        XC2Fields[idx, :] = X[:, 1]
        
        XiFields[idx, :] = RS[:, 0]
        EtaFields[idx, :] = RS[:, 1]
        kFields[idx, :] = k[:, 0]
        VAFields[idx, :] = VA[:, 0]
        VBFields[idx, :] = VB[:, 0]
        VCFields[idx, :] = VC[:, 0]    
    
    # baseline
    k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(mapTau.tau)
    X = mapTau._C2X(C)
    RS = mapTau._phys2Natural(X)
    VA, VB, VC = mapTau.getThetaVABC(mapTau.tau)
    
    c1Field_base = C[:, 0]
    c2Field_base = C[:, 1]
    c3Field_base = C[:, 2]
    
    XC1Field_base = X[:, 0]
    XC2Field_base = X[:, 1]
    
    XiField_base = RS[:, 0]
    EtaField_base = RS[:, 1]
    kField_base = k[:, 0]
        
    VAField_base = VA[:, 0]
    VBField_base = VB[:, 0]
    VCField_base = VC[:, 0]
    
    # output all components
    np.savetxt(dataFolder+'/DA-'+DAstep+'/XC1_s', XC1Fields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/XC2_s', XC2Fields)        
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c1_s', c1Fields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c2_s', c2Fields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c3_s', c3Fields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/Xi_s', XiFields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/Eta_s', EtaFields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/TKE_s', kFields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VA_s', VAFields)                        
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VB_s', VBFields)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VC_s', VCFields)
    
    
    np.savetxt(dataFolder+'/DA-'+DAstep+'/XC1_base', XC1Field_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/XC2_base', XC2Field_base)        
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c1_base', c1Field_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c2_base', c2Field_base)        
    np.savetxt(dataFolder+'/DA-'+DAstep+'/c3_base', c3Field_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/Xi_base', XiField_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/Eta_base', EtaField_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/TKE_base', kField_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VA_base', VAField_base)                        
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VB_base', VBField_base)
    np.savetxt(dataFolder+'/DA-'+DAstep+'/VC_base', VCField_base)

       
if __name__ == '__main__':        
    dataFolder = './debugData'
    mainInputFile = './MainInput.in'
    modelInputFile = './forwardModelInput.in'
    plotInputFile = './plotInfo.in'
    paramDict_main = readInputData(mainInputFile)
    paramDict_model = readInputData(modelInputFile)
    paramDict_plot = readInputData(plotInputFile)
    
    Ns = int(paramDict_main['Ns'])
    DAInterval = float(paramDict_main['DAInterval'])
    
    problem = paramDict_model['caseName']
    Ncell = int(paramDict_model['Ncell'])
        
    EnKFStep = int(paramDict_plot['EnKFStep'])
    timeStep = DAInterval * EnKFStep
    
    baseCaseDir = ospt.join(os.getcwd(), problem, '0/') 
    
    mapTau = ReRF.ReynoldsStressRF(baseCaseDir, 'Tau', Ncell, 1, True)
    
    getPhyComponents(EnKFStep)
    


