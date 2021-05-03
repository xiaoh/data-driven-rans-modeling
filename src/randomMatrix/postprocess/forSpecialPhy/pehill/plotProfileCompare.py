#!/aoe/bin/python27

# description        :Contain all functions for postprocessing the result of maxEnt

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.03, 2015
# revision           :Dec.05, 2015
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
# plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
# Import local modules
import hillShape
from utilities import readInputData

lineWidth1 = 2.0
lineWidth2 = 3.0
colorRMT = 'blue'
colorPhy = 'magenta'
colorMean = 'black'
colorBase = 'darkred'

labelDict1 = {
            'c1_SField':r'$S(C_1)$', 'c2_SField':r'$S(C_2$)', 'c3_SField':r'$S(C_3)$',
            'Xi_SField':r'$S(\xi)$', 'Eta_SField':r'$S(\eta)$', 'TKE_SField':r'$S(k/U_b^2)$',
            'VA_SField':r'$S(\varphi_1)$', 'VB_SField':r'$S(\varphi_2)$', 'VC_SField':r'$S(\varphi_3)$',
            'deltaXi_SField':r'$S(\Delta\xi)$', 'deltaEta_SField':r'$S(\Delta\eta)$', 'deltaK_SField':r'$S(\Delta$ln$k)$',
            'deltaVA_SField':r'$S(\Delta\varphi_1)$', 'deltaVB_SField':r'$S(\Delta\varphi_2)$', 
            'deltaVC_SField':r'$S(\Delta\varphi_3)$', 
            }

labelDict2 = {
            'c1_SField':r'$S(C_1)^{(RMT)} - S(C_1)^{(Phy)}$', 'c2_SField':r'$S(C_2)^{(RMT)} - S(C_2)^{(Phy)}$', 
            'c3_SField':r'$S(C_3)^{(RMT)} - S(C_3)^{(Phy)}$', 'Xi_SField':r'$S(\xi)^{(RMT)} - S(\xi)^{(Phy)}$',
            'Eta_SField':r'$S(\eta)^{(RMT)} - S(\eta)^{(Phy)}$', 'TKE_SField':r'$S(k/U_b^2)^{(RMT)} - S(k/U_b^2)^{(Phy)}$',
            'VA_SField':r'$S(\varphi_1)^{(RMT)} - S(\varphi_1)^{(Phy)}$', 'VB_SField':r'$S(\varphi_2)^{(RMT)} - S(\varphi_2)^{(Phy)}$', 
            'VC_SField':r'$S(\varphi_3)^{(RMT)} - S(\varphi_2)^{(Phy)}$',
            'deltaXi_SField':r'$S(\Delta\xi)^{(RMT)} - S(\Delta\xi)^{(Phy)}$', 'deltaEta_SField':r'$S(\Delta\eta^{(RMT)} - S(\Delta\eta)^{(Phy)}$', 
            'deltaK_SField':r'$S(\Delta$ln$k)^{(RMT)} - S(\Delta$ln$k)^{(Phy)}$',
            'deltaVA_SField':r'$S(\Delta\varphi_1)^{(RMT)} - S(\Delta\varphi_1)^{(Phy)}$', 
            'deltaVB_SField':r'$S(\Delta\varphi_2)^{(RMT)} - S(\Delta\varphi_2)^{(Phy)}$', 
            'deltaVC_SField':r'$S(\Delta\varphi_3)^{(RMT)} - S(\Delta\varphi_3)^{(Phy)}$', 
            }

labelDict3 = {
            'c1_KLDistField':r'$KL(C_1)$', 'c2_KLDistField':r'$KL(C_2$)', 'c3_KLDistField':r'$KL(C_3)$',
            'Xi_KLDistField':r'$KL(\xi)$', 'Eta_KLDistField':r'$KL(\eta)$', 'TKE_KLDistField':r'$KL(k/U_b^2)$',
            'VA_KLDistField':r'$KL(\varphi_1)$', 'VB_KLDistField':r'$KL(\varphi_2)$', 'VC_KLDistField':r'$KL(\varphi_3)$',
            'deltaXi_KLDistField':r'$KL(\Delta\xi)$', 'deltaEta_KLDistField':r'$KL(\Delta\eta)$', 'deltaK_KLDistField':r'$KL(\Delta$ln$k)$',
            'deltaVA_KLDistField':r'$KL(\Delta\varphi_1)$', 'deltaVB_KLDistField':r'$KL(\Delta\varphi_2)$', 'deltaVC_KLDistField':r'$KL(\Delta\varphi_3)$', 
            }            

def plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV):
    plt.close()
    plt.figure()
    ax1=plt.subplot(111)
    base, p_RMT = _plotProfile_AtMultiline(xPosList, sampleSet_RMT, '0', componentName, scaleV, colorRMT)
    base, p_Phy = _plotProfile_AtMultiline(xPosList, sampleSet_Phy, '0', componentName, scaleV, colorPhy, True)      
    plt.axis([-0.5, 9, 0, 3.05])
    
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':12})
    _plotDomain()
    
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'$'+labelDict1[componentName]+r'$+x/H $')
    
    ax1.legend([base, p_RMT, p_Phy], ["Zero", "RMT", "Phy"], prop={'size':12}, numpoints=1,
               bbox_to_anchor=(0.22, 1.01), loc=3, ncol=3)

    
    fname = figurefolder + 'profile_'+componentName+'_compare_'+ caseName +'.pdf'
    plt.savefig(fname) 
    
def plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV):
    """
    
    """
    
    plt.close()
    plt.figure()
    ax1=plt.subplot(111)
    base, p = _plotProfile_AtMultilineDiff(xPosList, sampleSet_RMT, sampleSet_Phy, '0', componentName, scaleV, colorBase, True)    
    plt.axis([-0.5, 9, 0, 3.05])
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':12})
    ax1.legend([base, p], ["Zero", "Entropy Difference"], prop={'size':12}, numpoints=1,
               bbox_to_anchor=(0.22, 1.01), loc=3, ncol=2)    
    _plotDomain()
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'($'+labelDict2[componentName]+r'$)+x/H $')    
    fname = figurefolder + 'profile_'+componentName+'diff_compare_'+ caseName +'.pdf'
    plt.savefig(fname)
    
def plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV):
    """
    
    """ 
    plt.close()
    plt.figure()
    ax1=plt.subplot(111)
    base, p = _plotProfile_AtMultiline(xPosList, sampleSet_RMT, '0', componentName, scaleV, colorPhy, True)   
    plt.axis([-0.5, 10, 0, 3.05])    
    ax1.set_aspect(1.3)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    matplotlib.rcParams.update({'font.size':12})
    _plotDomain()
    ax1.legend([base, p], ["Zero", "KL Divergence"], prop={'size':12}, numpoints=1,
               bbox_to_anchor=(0.22, 1.01), loc=3, ncol=2) 
    plt.ylabel("$y/H$")
    plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'($'+labelDict3[componentName]+r'$)+x/H $')    
  
    fname = figurefolder + 'profile_'+componentName+'_KLDiverg_'+ caseName +'.pdf'
    plt.savefig(fname)
                 
def _plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    plt.plot(y, h, 'g-', lw=lineWidth2)

def _plotProfile_AtMultilineDiff(xPosList, sampleSet1, sampleSet2, timeStep, componentName, scaleV, colorS, baseFlag=False):
    for xPos in xPosList:
        fileName1 = sampleSet1 + timeStep + '/line_x' + str(xPos) + '_' + componentName + '.xy'
        fileName2 = sampleSet2 + timeStep + '/line_x' + str(xPos) + '_' + componentName + '.xy'
        M1 = np.loadtxt(fileName1, comments = '%')
        M2 = np.loadtxt(fileName2, comments = '%')
        yCoor = M1[:, 0]
        M1mM2 = (M1[:, 1] - M2[:, 1]) * scaleV + float(xPos)
        if baseFlag:
            base, = plt.plot(M1[:, 1]*0.0+float(xPos), yCoor , lw=lineWidth2, color=colorMean, ls='--')
        else:
            base = None
        
        pCom, = plt.plot(M1mM2, yCoor, lw=lineWidth1, color=colorS)
    return base, pCom
            
def _plotProfile_AtMultiline(xPosList, sampleSet, timeStep, componentName, scaleV, colorS, baseFlag=False):
    for xPos in xPosList:
        fileName = sampleSet + timeStep + '/line_x' + str(xPos) + '_' + componentName + '.xy'
        M = np.loadtxt(fileName, comments = '%')
        yCoor = M[:, 0]
        component = M[:, 1] * scaleV + float(xPos)
        if baseFlag:
            base, = plt.plot(M[:, 1]*0.0+float(xPos), yCoor , lw=lineWidth2, color=colorMean, ls='--')
        else:
            base = None
        pCom, = plt.plot(component, yCoor, lw=lineWidth1, color=colorS)
    
    return base, pCom


def obtainRij():
    """
    plot R_ij statistics at one pt
    """
    Ub = 0.028
    R_11_rmt = np.zeros([self.nSample, self.nCell])
    R_12_rmt = np.zeros([self.nSample, self.nCell])
    R_13_rmt = np.zeros([self.nSample, self.nCell])
    R_22_rmt = np.zeros([self.nSample, self.nCell])
    R_23_rmt = np.zeros([self.nSample, self.nCell])
    R_33_rmt = np.zeros([self.nSample, self.nCell])

    R_11_phy = np.zeros([self.nSample, self.nCell])
    R_12_phy = np.zeros([self.nSample, self.nCell])
    R_13_phy = np.zeros([self.nSample, self.nCell])
    R_22_phy = np.zeros([self.nSample, self.nCell])
    R_23_phy = np.zeros([self.nSample, self.nCell])
    R_33_phy = np.zeros([self.nSample, self.nCell])
    
    for isample in range(self.nSample):
        
        R_rmt_i = np.loadtxt(resultDir_RMT+'R_samples/'+'R_s'+str(isample))
        R_phy_i = np.loadtxt(resultDir_Phy+'Tau/Tau_init_sample-'+str(isample))
        
        R_11_rmt[isample, :] = R_rmt_i[:, 0]/Ub/Ub
        R_12_rmt[isample, :] = R_rmt_i[:, 1]/Ub/Ub
        R_13_rmt[isample, :] = R_rmt_i[:, 2]/Ub/Ub
        R_22_rmt[isample, :] = R_rmt_i[:, 3]/Ub/Ub
        R_23_rmt[isample, :] = R_rmt_i[:, 4]/Ub/Ub
        R_33_rmt[isample, :] = R_rmt_i[:, 5]/Ub/Ub
    
        R_11_phy[isample, :] = R_phy_i[:, 0]/Ub/Ub
        R_12_phy[isample, :] = R_phy_i[:, 1]/Ub/Ub
        R_13_phy[isample, :] = R_phy_i[:, 2]/Ub/Ub
        R_22_phy[isample, :] = R_phy_i[:, 3]/Ub/Ub
        R_23_phy[isample, :] = R_phy_i[:, 4]/Ub/Ub
        R_33_phy[isample, :] = R_phy_i[:, 5]/Ub/Ub

    return R_11_rmt, R_12_rmt, R_13_rmt, R_22_rmt, R_23_rmt, R_33_rmt, \
           R_11_phy, R_12_phy, R_13_phy, R_22_phy, R_23_phy, R_33_phy
    
        
def main(iShow = False):
    xPosList = [ 1, 2, 3, 4, 5, 6, 7, 8]
    
    if allEntFlag or EntCFlag:
        namePool = ['c1_SField', 'c2_SField', 'c3_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
        namePool = ['c1_KLDistField', 'c2_KLDistField', 'c3_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)
            
    if allEntFlag or EntXiEtaFlag:
        namePool = ['Xi_SField', 'Eta_SField']
        for componentName in namePool:        
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)    
        namePool = ['Xi_KLDistField', 'Eta_KLDistField']
        for componentName in namePool:        
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)
                        
    if allEntFlag or EntTKEFlag:
        namePool = ['TKE_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
        namePool = ['TKE_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)    
                
    if allEntFlag or EntVFlag:
        namePool = ['VA_SField', 'VB_SField', 'VC_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
        namePool = ['VA_KLDistField', 'VB_KLDistField', 'VC_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)    
            
    if allEntFlag or EntdeltaXiEtaFlag:
        namePool = ['deltaXi_SField', 'deltaEta_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
        namePool = ['deltaXi_KLDistField', 'deltaEta_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)                
        
    if allEntFlag or EntdeltaXiEtaFlag:
        namePool = ['deltaK_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)
        namePool = ['deltaK_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)
                 
    if allEntFlag or EntdeltaVFlag:
        namePool = ['deltaVA_SField', 'deltaVB_SField', 'deltaVC_SField']
        for componentName in namePool:
            print "Plotting Entropy field for " + componentName + " ..."
            scaleV = 1
            plotProfile_EntropyDiff(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)            
            scaleV = 0.1
            plotProfile_Entropy(sampleSet_RMT, sampleSet_Phy, xPosList, componentName, scaleV)    
        namePool = ['deltaVA_KLDistField', 'deltaVB_KLDistField', 'deltaVC_KLDistField']
        for componentName in namePool:
            scaleV = 0.3
            plotProfile_KLDiverg(sampleSet_RMT, xPosList, componentName, scaleV)
    
              
    if iShow:
        plt.show()
            

if __name__ == "__main__":
    mainInputFile = './mainInput.in'
    plotInputFile = './plotInfo.in'
    paramDict = readInputData(plotInputFile)
    paramDict_main = readInputData(mainInputFile)
    figurefolder = './figures/'
    
    timeDir = paramDict['timeDir']
    caseName = paramDict['caseName']
    case = paramDict_main['caseName']
    caseRMT = '.'
    casePhy = paramDict['case_compared']     
    sampleSet_RMT = caseRMT+'/'+case+'_sample/postProcessing/sets/'
    sampleSet_Phy = casePhy+'/'+case+'_sample/postProcessing/sets/'
    
    # Plot control
    plotUx = paramDict['plotUx']
    plotUy = paramDict['plotUy']
    plotTau = paramDict['plotTau']
    plotDeltaTau = paramDict['plotDeltaTau']
    
    allEntFlag = ast.literal_eval(paramDict['scalarEntFlag'])
    EntCFlag = ast.literal_eval(paramDict['EntCFlag'])
    EntXiEtaFlag = ast.literal_eval(paramDict['EntXiEtaFlag'])
    EntTKEFlag = ast.literal_eval(paramDict['EntTKEFlag'])
    EntVFlag = ast.literal_eval(paramDict['EntVFlag'])
    EntdeltaXiEtaFlag = ast.literal_eval(paramDict['EntdeltaXiEtaFlag'])
    EntdeltaKFlag = ast.literal_eval(paramDict['EntdeltaKFlag'])    
    EntdeltaVFlag = ast.literal_eval(paramDict['EntdeltaVFlag'])
        
    if not os.path.exists('figures/'):
        os.mkdir('figures/')
    
    main(False)

