#!/aoe/bin/python27

# description        :comparing profiles btween Phy and RMT

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :May.03, 2016
# revision           :May.05, 2016
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
import fnmatch
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
import postFuns as pF


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
            
labelDict_tau = {
                    'R11':r'$R_{11}$', 'R12':r'$R_{12}$', 'R13':r'$R_{13}$', 'R22':r'$R_{22}$', 'R23':r'$R_{23}$','R33':r'$R_{33}$' 
                }            

xVec = [1., 2., 3., 4., 5., 6., 7., 8.]
y0 = 2
z0 = 0.05
Ub = 0.028
def getProfiles_MeshIdx(caseMesh, xVec, y0, z0):
    """
    To extract profiles on give xVec = [x0, x1, ..., xn]
    
    Arg:
    """
    
    meshIdx_dict = {} # define a dictionary for store profiles at given x loc
    meshYcoord_dict = {}
    for x0 in xVec:   
        hCoord, vCoord = pF.getVerandHorCoord(x0, y0, z0, caseMesh)
        idx_v = pF.coord2cellIdx(vCoord, caseMesh)    
        idx_locs = idx_v[:, 1].tolist()              
        meshIdx_dict[str(x0)] = idx_locs
        meshYcoord_dict[str(x0)] = vCoord
    
    return meshIdx_dict, meshYcoord_dict

def _plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])
    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)
    plt.plot(y, h, 'g-', lw=lineWidth2)

def plotProfile(filename, xPos, scaleV, xCol, yCol, ls, cl, dash, para = 'U'):
    M = np.loadtxt(filename, comments = '%')
    y = M[:,yCol]
    if para == 'Tau':
        #scaleV = -scaleV
        pass
    x = M[:,xCol] * scaleV + float(xPos)
    if dash == 'None':
        p, = plt.plot(x, y, ls, color = cl, lw=2)
    else:
        p, = plt.plot(x, y, ls, color = cl, lw=2, dashes=dash)
    return p

def ansysBubble(xList, tauwList):
    """
    calculate bubble length
    """
    index1, index2 = detectBubble(xList, tauwList)
    detachPos = getPos(xList, tauwList, index1)
    attachPos = getPos(xList, tauwList, index2)
    bubbleL = attachPos - detachPos
    return bubbleL, attachPos

def detectBubble(xList, tauwList):
    crossPos = [0,0]
    signs = np.sign(tauwList)
    crossList = np.where(np.diff(np.sign(signs)))[0]
    if crossList.shape[0] >= 1:
        if tauwList[0] > 0:
            signs = np.sign(tauwList)
            crossList = np.where(np.diff(np.sign(signs)))[0]
            crossPos[0] = crossList[0]
            crossPos[1] = crossList[1]
        else:
            signs = np.sign(tauwList)
            crossList = np.where(np.diff(np.sign(signs)))[0]
            crossPos[1] = crossList[0]
        signs = np.sign(tauwList)
        crossList = np.where(np.diff(np.sign(signs)))[0]
        iter = -1
        while xList[crossList[iter]] >=7:
            iter = iter-1
        crossPos[1] = crossList[iter]
        #iter = 0
        #while xList[crossList[iter]] < 3 and iter < crossList.shape[0]-2:
        #    iter = iter + 1
        #crossPos[1] = crossList[iter]
    return crossPos[0], crossPos[1]

def getPos(xList, yList, index):
    pos = (xList[index+1]-xList[index])/(yList[index+1]-yList[index])*(0-yList[index])+xList[index]
    if pos < 0:
        pos = 0
    return pos

    
class postPlotCompare:

    """

    Arg:


    """

    def __init__(self, mainInput):
        paramDict = readInputData(mainInputFile)
        self.case = paramDict['caseName']
        self.nSample = int(paramDict['nSample'])
        self.nSample_propage = int(paramDict['nSample_propagate'])
        self.nCell = int(paramDict['nCell_cfd'])
        # get cell idx to be analyzed 
        meshcase = self.case+'_mesh'
        
        self.MeshIdx, self.meshYcoord = getProfiles_MeshIdx(meshcase, xVec, y0, z0)
        
        #self.R11_rmt, self.R12_rmt, self.R13_rmt, self.R22_rmt, self.R23_rmt, self.R33_rmt,\
        #self.R11_phy, self.R12_phy, self.R13_phy, self.R22_phy, self.R23_phy, self.R33_phy = self.obtainRij()

    
    def plotTauSamples(self, tau_ij_rmt, tau_ij_phy, scaleV, componentName):
        """
        plot tau component samples comparison
        """
        plt.figure()
        ax1=plt.subplot(111)
        
        for xPos in xVec:
            vCoord = self.meshYcoord[str(xPos)]
            yy = vCoord[:, 1]
            
            tau_ij_xVec_rmt = tau_ij_rmt[:, self.MeshIdx[str(xPos)]]
            tau_ij_xVec_phy = tau_ij_phy[:, self.MeshIdx[str(xPos)]]        
            
            xx_rmt = scaleV * tau_ij_xVec_rmt + xPos
            xx_phy = scaleV * tau_ij_xVec_phy + xPos
            
            for i in range(self.nSample):
                ps_rmt, = plt.plot(xx_rmt[i, :], yy, color='green', alpha=0.2, lw = 2.0, mfc = 'none')            
                ps_phy, = plt.plot(xx_phy[i, :], yy, color='red', alpha=0.2, lw = 2.0, mfc = 'none')
            
            tau_ij_xVecMean_rmt = np.mean(tau_ij_xVec_rmt, axis=0)
            tau_ij_xVecMean_phy = np.mean(tau_ij_xVec_phy, axis=0)
            
            xMean_rmt = scaleV * tau_ij_xVecMean_rmt + xPos
            xMean_phy = scaleV * tau_ij_xVecMean_phy + xPos
            
            pMean_rmt, = plt.plot(xMean_rmt, yy, '-',color ='blue',lw=2, dashes=(9,2,2,2))        
            pMean_phy, = plt.plot(xMean_phy, yy, '-',color ='black',lw=2, dashes=(9,2,2,2))
        
        
            
        plt.axis([-0.5, 9, 0, 3.05])
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        matplotlib.rcParams.update({'font.size':12})
        ax1.legend([ps_rmt, ps_phy, pMean_rmt, pMean_phy], ["RMT", "Phy"], prop={'size':12}, numpoints=1,
                   bbox_to_anchor=(0.22, 1.01), loc=3, ncol=2)    
        _plotDomain()
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'($'+labelDict_tau[componentName]+r'$)+x/H $')    
        fname = figurefolder + 'profile_'+componentName+'_compare_'+ caseName +'.pdf'
        plt.savefig(fname)     


    def plotTauArea(self, tau_ij_rmt, tau_ij_phy, scaleV, componentName):
        """
        plot tau component samples comparison
        """
        idx_tau = {
                    'R11':1, 'R12':2, 'R13':3, 'R22':4, 'R23':5,'R33':6 
                }         
        plt.figure()
        ax1=plt.subplot(111)
        
        for xPos in xVec:
            vCoord = self.meshYcoord[str(xPos)]
            yy = vCoord[:, 1]
            
            tau_ij_xVec_rmt = tau_ij_rmt[:, self.MeshIdx[str(xPos)]]
            tau_ij_xVec_phy = tau_ij_phy[:, self.MeshIdx[str(xPos)]]        
                        
            tau_ij_std_rmt = np.std(tau_ij_xVec_rmt, axis=0)
            tau_ij_std_phy = np.std(tau_ij_xVec_phy, axis=0)
            
            tau_ij_xVecMean_rmt = np.mean(tau_ij_xVec_rmt, axis=0)
            tau_ij_xVecMean_phy = np.mean(tau_ij_xVec_phy, axis=0)
            
            xUp_rmt = scaleV * (tau_ij_xVecMean_rmt + 2 * tau_ij_std_rmt) + xPos
            xUp_phy = scaleV * (tau_ij_xVecMean_phy + 2 * tau_ij_std_phy) + xPos

            xlow_rmt = scaleV * (tau_ij_xVecMean_rmt - 2 * tau_ij_std_rmt) + xPos
            xlow_phy = scaleV * (tau_ij_xVecMean_phy - 2 * tau_ij_std_phy) + xPos            
            
            xMean_rmt = scaleV * tau_ij_xVecMean_rmt + xPos
            xMean_phy = scaleV * tau_ij_xVecMean_phy + xPos
                                  
            pUp_rmt = plt.plot(xUp_rmt, yy, '-',color ='darkblue',lw=0.5)
            plow_rmt = plt.plot(xlow_rmt, yy, '-',color ='darkblue',lw=0.5)
            plt.fill_betweenx(yy, xlow_rmt, xUp_rmt, linewidth=0, facecolor = 'darkblue', alpha=0.3)
            
            pUp_phy = plt.plot(xUp_phy, yy, '-',color ='darkred',lw=0.5)
            plow_phy = plt.plot(xlow_phy, yy, '-',color ='darkred',lw=0.5)
            plt.fill_betweenx(yy, xlow_phy, xUp_phy, linewidth=0, facecolor = 'darkred', alpha=0.3)
            
            pMean_rmt, = plt.plot(xMean_rmt, yy, '-',color ='darkblue',lw=2, dashes=(9,2,2,2))        
            pMean_phy, = plt.plot(xMean_phy, yy, '-',color ='darkred',lw=2, dashes=(2,2,1))
        
        pDNS, scatDNS = self.plotDNS(scaleV, idx_tau[componentName], 0, 'Tau')
        
        pRMT = plt.Rectangle((0, 0), 1, 1, fc='darkblue', alpha=0.3)
        pPhy = plt.Rectangle((0, 0), 1, 1, fc='darkred', alpha=0.3)
            
        plt.axis([-0.5, 9, 0, 3.05])
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        matplotlib.rcParams.update({'font.size':12})
        ax1.legend([pMean_phy, pMean_rmt, pRMT, pPhy], ["RMT Mean", "RMT Phy", "RMT 2-$\sigma$", "Phy 2-$\sigma$"], prop={'size':12}, numpoints=1,
                   bbox_to_anchor=(0.02, 1.01), loc=3, ncol=4)   
        _plotDomain()
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'($'+labelDict_tau[componentName]+r'$)+x/H $')    
        fname = figurefolder + 'profile_'+componentName+'_compareArea_'+ caseName +'.pdf'
        plt.savefig(fname)

    def plotUArea(self, scaleV=1.0):
        """
        plot U samples comparison (RMT and Phy)
        """
        plt.figure()
        ax1=plt.subplot(111)
        _plotDomain()               
        
        for xPos in xVec:
            i = 0
            filename = 'pehill_base-tmp_1.0/postProcessing/sets/3000.000000' + '/line_x' + str(int(xPos)) + '_UNorm.xy'
            M = np.loadtxt(filename, comments = '%')
            nVcell = M.shape[0]
            yy = M[:, 0]
            M_RMT = np.zeros([self.nSample_propage, nVcell])
            for f in os.listdir('.'):
                if fnmatch.fnmatch(f, 'pehill_base-tmp*'):                    
                    filename = f + '/postProcessing/sets/3000.000000' + '/line_x' + str(int(xPos)) + '_UNorm.xy'
                    M = np.loadtxt(filename, comments = '%')
                    M_RMT[i, :] = M[:, 1]                    
                    #ps_rmt, = plt.plot(M[:, 1] + xPos, yy, color='green', alpha=0.2, lw = 2.0, mfc = 'none')
                    i = i + 1

            Mean_rmt = np.mean(M_RMT, axis=0)
            std_rmt = np.std(M_RMT, axis=0)        
            
            xMean_rmt = scaleV * Mean_rmt + xPos
            xUp_rmt = scaleV * (Mean_rmt + 2 * std_rmt) + xPos
            xLow_rmt = scaleV * (Mean_rmt - 2 * std_rmt) + xPos
            
            pMean_rmt, = plt.plot(xMean_rmt, yy, '-',color ='darkblue',lw=2, dashes=(9,2,2,2))   
            pUp_rmt, = plt.plot(xUp_rmt, yy, '-',color ='darkblue',lw=0.5)
            plow_rmt, = plt.plot(xLow_rmt, yy, '-',color ='darkblue',lw=0.5)
            plt.fill_betweenx(yy, xLow_rmt, xUp_rmt, linewidth=0, facecolor = 'darkblue', alpha=0.3)
            
            
            i = 0
            M_Phy = np.zeros([self.nSample_propage, nVcell])
            for f in os.listdir(case_compared):
                if fnmatch.fnmatch(f, 'pehill_base-tmp*'):                   
                    filename = case_compared + '/' +  f + '/postProcessing/sets/3000.000000' + '/line_x' + str(int(xPos)) + '_U.xy'
                    M = np.loadtxt(filename, comments = '%')
                    M_Phy[i, :] = M[:, 1]/Ub                    
                    #ps_rmt, = plt.plot(M[:, 1] + xPos, yy, color='green', alpha=0.2, lw = 2.0, mfc = 'none')
                    i = i + 1

            Mean_phy = np.mean(M_Phy, axis=0)
            std_phy = np.std(M_Phy, axis=0)        
            
            xMean_phy = scaleV * Mean_phy + xPos
            xUp_phy = scaleV * (Mean_phy + 2 * std_phy) + xPos
            xLow_phy = scaleV * (Mean_phy - 2 * std_phy) + xPos
            
            pMean_phy, = plt.plot(xMean_phy, yy, '-',color ='darkred',lw=2, dashes=(2,2,1))   
            pUp_phy, = plt.plot(xUp_phy, yy, '-',color ='darkred',lw=0.5)
            plow_phy, = plt.plot(xLow_phy, yy, '-',color ='darkred',lw=0.5)
            plt.fill_betweenx(yy, xLow_phy, xUp_phy, linewidth=0, facecolor = 'darkred', alpha=0.3)
            
                    
        pDNS, scatDNS = self.plotDNS(scaleV, 1, 0, 'U')                       
        #pRANS, scatBase = self.plotBaseline(scaleV, 1, 0, 'U')     
        
        pRMT = plt.Rectangle((0, 0), 1, 1, fc='darkblue', alpha=0.3)
        pPhy = plt.Rectangle((0, 0), 1, 1, fc='darkred', alpha=0.3)
            
        plt.axis([-0.5, 11, 0, 3.05])
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)
        matplotlib.rcParams.update({'font.size':12})
        ax1.legend([pMean_phy, pMean_rmt, pRMT, pPhy], ["RMT Mean", "RMT Phy", "RMT 2-$\sigma$", "Phy 2-$\sigma$"], prop={'size':12}, numpoints=1,
                   bbox_to_anchor=(0.02, 1.01), loc=3, ncol=4)    

        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad  $'+r'$'+str(scaleV)+'(U_x/U_b)+x/H $')
            
        fname = figurefolder + 'profile_U_compareArea_'+ caseName +'.pdf'
        plt.savefig(fname)        


    def plotWallShearArea(self, scaleV=1.0):
        """
        plot wallShear comparison (RMT and Phy)
        """
        fig = plt.figure(figsize=(8,7));
        gs = matplotlib.gridspec.GridSpec(2,1,height_ratios=[4,1])
        #plt.subplot(2,1,1)
        plt.subplot(gs[0])
        
        
        ySum = 0
        scatSample = [0,0]
        filename = './pehill_base-tmp_1.0/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
        M = np.loadtxt(filename, comments = '#')
        M = M[1:M.shape[0]/2,:]
        nXCell = M.shape[0]
        shearSamples = np.zeros([self.nSample_propage, nXCell])
        i = 0
        for f in os.listdir('.'):              
            if fnmatch.fnmatch(f, 'pehill_base-tmp*'):
                filename = f + '/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
                M = np.loadtxt(filename, comments = '#')
                M = M[1:M.shape[0]/2,:]
                y = M[:,3]*(-1)*2/Ub/Ub
                ySum = ySum + y
                x = M[:,0]                 
                scat = ansysBubble(x,y)
                shearSamples[i, :] = y
                if scat > 0:
                    scatSample = np.vstack((scatSample,scat))
                #p1, = plt.plot(x, y,'-', color='#1b9e76',alpha=0.3, lw = 1, markevery=20,mfc='none')
                i = i + 1
                
        shearSamplesMean = np.mean(shearSamples, axis=0)
        shearSamplesStd = np.std(shearSamples, axis=0)
        shearSamplesUp = shearSamplesMean + 2*shearSamplesStd
        shearSamplesLow = shearSamplesMean - 2*shearSamplesStd
        
        p2, = plt.plot(x, shearSamplesMean, '--', dashes=(9,2,2,2), color='darkblue', lw = 2, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesUp, '-', dashes=(9,2,2,2), color='darkblue', lw = 0.5, alpha=1.0, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesLow, '-', dashes=(9,2,2,2), color='darkblue', lw = 0.5, alpha=1.0, markevery=20,mfc='none')
        plt.fill_between(x, shearSamplesLow, shearSamplesUp, linewidth=0, facecolor = 'darkblue', alpha=0.3)
        scatSample = scatSample[1:,:]
        scatMean = ansysBubble(x,shearSamplesMean)
        
        Mmean = np.transpose(np.vstack((x, shearSamplesMean)))
        
        
        np.savetxt('meanShearStress', Mmean)                               

        ySum = 0
        scatSample_phy = [0,0]
        filename = case_compared +'/pehill_base-tmp_1.0/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
        M = np.loadtxt(filename, comments = '#')
        M = M[1:M.shape[0]/2,:]
        nXCell = M.shape[0]
        shearSamples_phy = np.zeros([self.nSample_propage, nXCell])
        i = 0
        for f in os.listdir(case_compared):              
            if fnmatch.fnmatch(f, 'pehill_base-tmp*'):
                filename = case_compared + '/' +  f  + '/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
                M = np.loadtxt(filename, comments = '#')
                M = M[1:M.shape[0]/2,:]
                y = M[:,3]*(-1)*2/Ub/Ub
                x = M[:,0]                 
                scat_phy = ansysBubble(x,y)
                shearSamples_phy[i, :] = y
                if scat_phy > 0:
                    scatSample_phy = np.vstack((scatSample_phy,scat_phy))
                #p1, = plt.plot(x, y,'-', color='red',alpha=0.3, lw = 1, markevery=20,mfc='none')
                i = i + 1
                
        shearSamplesMean_phy = np.mean(shearSamples_phy, axis=0)
        shearSamplesStd_phy = np.std(shearSamples_phy, axis=0)
        shearSamplesUp_phy = shearSamplesMean_phy + 2*shearSamplesStd_phy
        shearSamplesLow_phy = shearSamplesMean_phy - 2*shearSamplesStd_phy
        
        p2, = plt.plot(x, shearSamplesMean_phy, '--', dashes=(9,2,2,2), color='darkred', lw = 2, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesUp_phy, '-', dashes=(9,2,2,2), color='darkred', lw = 0.5, alpha=1.0, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesLow_phy, '-', dashes=(9,2,2,2), color='darkred', lw = 0.5, alpha=1.0, markevery=20,mfc='none')
        plt.fill_between(x, shearSamplesLow_phy, shearSamplesUp_phy, linewidth=0, facecolor = 'darkred', alpha=0.3)
        scatSample_phy = scatSample_phy[1:,:]
        scatMean_phy = ansysBubble(x,shearSamplesMean_phy)
        
        Mmean_phy = np.transpose(np.vstack((x, shearSamplesMean_phy)))
        
        #plot base
        pBase, scatBase = self.plotBaseline(scaleV,0,3, 'wallShear');
        pDNS, scatDNS = self.plotDNS(scaleV,0,3, 'wallShear');

        plt.plot([0,9],[0,0],'k--');
        plt.axis([0, 9, -0.02, 0.06])
        plt.locator_params(axis = 'y', nbins = 4)
        plt.ylabel("Wall Shear Stress"+r"$ (\tau_w)$",fontsize=18)
        plt.xlabel(r'$x$',fontsize=18)
        pRMT = plt.Rectangle((0, 0), 1, 1, fc='darkblue', alpha=0.3)
        pPhy = plt.Rectangle((0, 0), 1, 1, fc='darkred', alpha=0.3)
        
        lg = plt.legend([pRMT, pPhy, pBase,pDNS],["RMT", "Phy", "baseline","DNS (Breuer et al. 2009)"],loc = 0, prop={'size':12})
        #lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':18})

        plt.subplot(gs[1])
        p1, = plt.plot(scatSample[:,1],np.zeros_like(scatSample[:,1]), 'o',markersize=8,\
             markeredgewidth=0, color='darkblue',alpha=0.05)
        p1, = plt.plot(scatSample_phy[:,1],np.zeros_like(scatSample_phy[:,1]), 'o',markersize=8,\
             markeredgewidth=0, color='darkred',alpha=0.05)
        p3, = plt.plot(scatBase[1],0,'r^',markeredgecolor='darkgreen', \
                markerfacecolor='None',markersize=8, markeredgewidth=2, mfc='None')
        
        p2, = plt.plot(scatMean[1],0,'bx',markeredgecolor='darkblue', \
                markerfacecolor='None',markersize=8, markeredgewidth=2, mfc='None')
        p5, = plt.plot(scatMean_phy[1],0,'ro',markeredgecolor='darkred', \
                markerfacecolor='None',markersize=8, markeredgewidth=2, mfc='None')
                
        p4, = plt.plot(scatDNS[1],0,'ks',markeredgecolor='k', markerfacecolor='None',\
                markersize=8, markeredgewidth=2, mfc='None')


        plt.xlabel(r"$x_{attach}/H$",fontsize=18)
        #lg = plt.legend([p1,p2,p3,p4],["samples","sample mean","baseline","DNS"],loc = 0)
        #lg.draw_frame(False)
        plt.xlim([0,9])
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        #matplotlib.rcParams.update({'font.size':18})
        fig.tight_layout()

        fig.savefig("./figures/bubble-"+caseName+".pdf")

    def plotWallShearArea_clean(self, scaleV=1.0):
        """
        plot wallShear comparison (RMT and Phy)
        """
        fig = plt.figure()
        plt.subplot(111)
               
        ySum = 0
        scatSample = [0,0]
        filename = './pehill_base-tmp_1.0/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
        M = np.loadtxt(filename, comments = '#')
        M = M[1:M.shape[0]/2,:]
        nXCell = M.shape[0]
        shearSamples = np.zeros([self.nSample_propage, nXCell])
        i = 0
        for f in os.listdir('.'):              
            if fnmatch.fnmatch(f, 'pehill_base-tmp*'):
                filename = f + '/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
                M = np.loadtxt(filename, comments = '#')
                M = M[1:M.shape[0]/2,:]
                y = M[:,3]*(-1)*2/Ub/Ub
                ySum = ySum + y
                x = M[:,0]                 
                scat = ansysBubble(x,y)
                shearSamples[i, :] = y
                if scat > 0:
                    scatSample = np.vstack((scatSample,scat))
                #p1, = plt.plot(x, y,'-', color='#1b9e76',alpha=0.3, lw = 1, markevery=20,mfc='none')
                i = i + 1
                
        shearSamplesMean = np.mean(shearSamples, axis=0)
        shearSamplesStd = np.std(shearSamples, axis=0)
        shearSamplesUp = shearSamplesMean + 2*shearSamplesStd
        shearSamplesLow = shearSamplesMean - 2*shearSamplesStd
        
        p1, = plt.plot(x, shearSamplesMean, '--', dashes=(9,2,2,2), color='darkblue', lw = 2, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesUp, '-', dashes=(9,2,2,2), color='darkblue', lw = 0.5, alpha=1.0, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesLow, '-', dashes=(9,2,2,2), color='darkblue', lw = 0.5, alpha=1.0, markevery=20,mfc='none')
        plt.fill_between(x, shearSamplesLow, shearSamplesUp, linewidth=0, facecolor = 'darkblue', alpha=0.3)
        scatSample = scatSample[1:,:]
        scatMean = ansysBubble(x,shearSamplesMean)
        
        Mmean = np.transpose(np.vstack((x, shearSamplesMean)))
        
        
        np.savetxt('meanShearStress', Mmean)                               

        ySum = 0
        scatSample_phy = [0,0]
        filename = case_compared +'/pehill_base-tmp_1.0/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
        M = np.loadtxt(filename, comments = '#')
        M = M[1:M.shape[0]/2,:]
        nXCell = M.shape[0]
        shearSamples_phy = np.zeros([self.nSample_propage, nXCell])
        i = 0
        for f in os.listdir(case_compared):              
            if fnmatch.fnmatch(f, 'pehill_base-tmp*'):
                filename = case_compared + '/' +  f  + '/postProcessing/surfaces/3000.000000/wallShearStress_bottomWall.raw'
                M = np.loadtxt(filename, comments = '#')
                M = M[1:M.shape[0]/2,:]
                y = M[:,3]*(-1)*2/Ub/Ub
                x = M[:,0]                 
                scat_phy = ansysBubble(x,y)
                shearSamples_phy[i, :] = y
                if scat_phy > 0:
                    scatSample_phy = np.vstack((scatSample_phy,scat_phy))
                #p1, = plt.plot(x, y,'-', color='red',alpha=0.3, lw = 1, markevery=20,mfc='none')
                i = i + 1
                
        shearSamplesMean_phy = np.mean(shearSamples_phy, axis=0)
        shearSamplesStd_phy = np.std(shearSamples_phy, axis=0)
        shearSamplesUp_phy = shearSamplesMean_phy + 2*shearSamplesStd_phy
        shearSamplesLow_phy = shearSamplesMean_phy - 2*shearSamplesStd_phy
        
        p2, = plt.plot(x, shearSamplesMean_phy, '-', dashes=(2,2,1), color='darkred', lw = 2, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesUp_phy, '-', dashes=(9,2,2,2), color='darkred', lw = 0.5, alpha=1.0, markevery=20,mfc='none') 
        plt.plot(x, shearSamplesLow_phy, '-', dashes=(9,2,2,2), color='darkred', lw = 0.5, alpha=1.0, markevery=20,mfc='none')
        plt.fill_between(x, shearSamplesLow_phy, shearSamplesUp_phy, linewidth=0, facecolor = 'darkred', alpha=0.3)
        scatSample_phy = scatSample_phy[1:,:]
        scatMean_phy = ansysBubble(x,shearSamplesMean_phy)
        
        Mmean_phy = np.transpose(np.vstack((x, shearSamplesMean_phy)))
        
        #plot base
        #pBase, scatBase = self.plotBaseline(scaleV,0,3, 'wallShear');
        pDNS, scatDNS = self.plotDNS(scaleV,0,3, 'wallShear');

        plt.plot([0,9],[0,0],'k--');
        plt.axis([0, 9, -0.02, 0.06])
        plt.locator_params(axis = 'y', nbins = 4)
        plt.ylabel("Wall Shear Stress"+r"$ (\tau_w)$",fontsize=18)
        plt.xlabel(r'$x/H$',fontsize=18)
        pRMT = plt.Rectangle((0, 0), 1, 1, fc='darkblue', alpha=0.3)
        pPhy = plt.Rectangle((0, 0), 1, 1, fc='darkred', alpha=0.3)
        
        lg = plt.legend([pRMT, pPhy, p1, p2, pDNS],["RMT 2-$\sigma$", "Phy 2-$\sigma$", "RMT Mean", "Phy Mean", "DNS (Breuer et al. 2009)"],loc = 0, prop={'size':12})
        #lg.draw_frame(False)
        matplotlib.rcParams.update({'font.size':18})
        fig.tight_layout()
        fig.savefig("./figures/bubble-"+caseName+".pdf")
        
        
    def plotBaseline(self, scaleV, xCol, yCol, para = 'U'):
        """
        plot baseline results
        """
        
        if para == 'wallShear':
            filename = 'pehill_base-Run/postProcessing/surfaces/'+timeDir + \
                        '.000000/wallShearStress_bottomWall.raw'
            M = np.loadtxt(filename, comments = '#')
            M = M[1:M.shape[0]/2,:]
            y = M[:,yCol]*(-1)*2/Ub/Ub
            x = M[:,xCol]
            scatBase = ansysBubble(x,y)
            pBase, = plt.plot(x, y, '--', dashes=(12,3), color='darkgreen', lw = 2, markevery=20,mfc='none')        
        else:
            scatBase = 0
            for xPos in xVec:
                filename = 'pehill_base-Run/postProcessing/sets/3000.000000'+ \
                          '/line_x' + str(int(xPos)) + '_' + para + 'Norm.xy'
                pBase = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'darkgreen', (12,3), para)
        
        return pBase, scatBase    
    
    def plotDNS(self, scaleV, xCol, yCol, para = 'U'):
        """
        plot DNS results
        """
        if para == 'wallShear':
            filename = 'DNS/X/Cf.csv'
            M = np.loadtxt(filename, comments = '#')
            M = M[1:M.shape[0]/2,:]
            y = M[:,yCol]
            x = M[:,xCol]
            scatDNS = ansysBubble(x[1:],y[1:])
            pDNS, = plt.plot(x, y, '-', color='black', lw = 2,mfc='none')          
        else:
            scatDNS = 0
            for xPos in xVec:
                if para == 'U':
                    filename = 'DNS/X/x' + str(int(xPos)) + '.0_U2.xy'
                                                  
                else:
                    filename = 'DNS/X/dnsTau' + str(int(xPos)) 
                pDNS = plotProfile(filename, xPos, scaleV, xCol, yCol, '-', 'black', 'None', para)
        
        return pDNS, scatDNS            


    def obtainRij(self):
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
    """
    Main plot files
    """
    posPlt = postPlotCompare(mainInputFile)
    if compareTau:
        for R_ij in ['R11', 'R12', 'R22', 'R33']:
            Rcom_rmt = eval('posPlt.' + R_ij + '_rmt')
            Rcom_phy = eval('posPlt.' + R_ij + '_phy')
            posPlt.plotTauArea(Rcom_rmt, Rcom_phy, 20.0, R_ij)

    if compareU:
        posPlt.plotUArea()
    
    if compareWallShear:
        #posPlt.plotWallShearArea()
        posPlt.plotWallShearArea_clean()    




            

if __name__ == "__main__":
    mainInputFile = './mainInput.in'
    plotInputFile = './plotInfo.in'
    paramDict = readInputData(plotInputFile)
    
    case_compared = paramDict['case_compared'] 
    timeDir = paramDict['timeDir']
    caseName = paramDict['caseName']
    
    paramDict_main = readInputData(mainInputFile)    
    case = paramDict_main['caseName']

    resultDir_RMT = 'resultData/'
    resultDir_Phy = case_compared + '/debugData/init/' 

    figurefolder = './figures/'        
    if not os.path.exists('figures/'):
        os.mkdir('figures/')    
    
    compareTau = ast.literal_eval(paramDict['compareTau'])
    compareU = ast.literal_eval(paramDict['compareU'])
    compareWallShear = ast.literal_eval(paramDict['compareWallShear'])
    
    
    main(False)

