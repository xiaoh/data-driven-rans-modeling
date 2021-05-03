#!/aoe/bin/python27

# description        :Plot profiles in Pehill domain

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Feb.24, 2016
# revision           :Mar.13, 2016
####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.stats as ss
# system, file operation
import pdb
import sys
import ast
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
from scipy import interpolate
# Import local modules
import postFuns as pF
import hillShape
from utilities import readInputData, extractListFromDict
from ReynoldsStressRF import ReynoldsStressRF

Ub = 0.028
transparency = 0.1
lineWidth1 = 2.0
lineWidth2 = 3.0
dashType1 = (1, 2, 2)
dashType2 = (9, 2, 2, 2)
colorSample = 'blue'
colorMean = 'black'
colorBase = 'darkred'
colorGaussian = 'darkgreen'
figurefolder = './figures/'

labelDict = {'Xi':r'$\xi$', 'Eta':r'$\eta$', 'TKE':r'$k/{U_b}^2$', 'VA':r'$\phi_1$', 'VB':r'$\phi_2$', 'VC':r'$\phi_3$'}
class postProfile:

    """
    Arg:


    """

    def __init__(self):
        
        DAstep = "%.1f"%EnKFStep
        meshcase = problem +'_mesh'
        self.resultDir = dataFolder+'/DA-'+DAstep + '/'
        self.ncell_vline = 31
        self.nvline = 8
        self.cellidx_v = np.zeros((self.ncell_vline, self.nvline))
        self.x_array = np.array([1., 2., 3., 4., 5., 6., 7., 8.])        
        
        y = 2.0
        z = 0.05
        
        i = 0        
        self.y_sample = np.zeros((self.ncell_vline, self.nvline))
        for x in self.x_array:
            hCoord, vCoord = pF.getVerandHorCoord(x, y, z, meshcase)
            self.y_sample[:, i] = vCoord[:, 1]
            idx_v = pF.coord2cellIdx(vCoord, meshcase)
            self.cellidx_v[:, i] = idx_v[:, 1]
            i = i + 1
        if synFlag:
            self.XC1_dns, self.XC2_dns, self.Xi_dns, self.Eta_dns, self.TKE_dns, \
            self.VA_dns, self.VB_dns, self.VC_dns, self.XC1_rans, self.XC2_rans, \
            self.Xi_rans, self.Eta_rans, self.TKE_rans, self.VA_rans, self.VB_rans, \
            self.VC_rans= self._readDNS_Syn()
        else:
            # Get DNS data
            self.XC1_dns, self.XC2_dns, self.Xi_dns, self.Eta_dns, self.TKE_dns, \
            self.VA_dns, self.VB_dns, self.VC_dns, self.XC1_rans, self.XC2_rans, \
            self.Xi_rans, self.Eta_rans, self.TKE_rans, self.VA_rans, self.VB_rans, \
            self.VC_rans= self._readDNS()

        
    def mainPlot_component(self, QoI, scaleV):
        
        plt.figure();
        ax1=plt.subplot(111)
        self._plotDomain()
        p_sample = self._plotProfileSamples_component(QoI, scaleV)
        p_dns = self._plotDNS(QoI, scaleV)
        p_rans = self._plotRANS(QoI, scaleV) 
        plt.legend([p_rans, p_dns], ["Baseline", "Truth"], prop={'size':12}, numpoints=1,
                    bbox_to_anchor=(0.22, 1.01), loc=3, ncol=2)
        plt.ylabel("$y/H$")
        plt.xlabel(r'$x/H;\quad$ '+ str(scaleV)+labelDict[QoI] + ' $+x/H$') 
        plt.axis([-0.5, 11, 0, 3.05])
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)        
        mp.rcParams.update({'font.size':12})        
        figureName = figurefolder + 'profile_' + QoI +'_DA'+str(EnKFStep)+'_' + caseName + '.pdf'
        plt.savefig(figureName)

    def mainPlot_delta(self, QoI, scaleV):
        
        plt.figure();
        ax1=plt.subplot(111)
        self._plotDomain()
        p_sample = self._plotProfileSamples_delta(QoI, scaleV)
        p_deltaDns = self._plotDeltaRANS(QoI, scaleV)  
        plt.axis([-0.5, 11, 0, 3.05])
        ax1.set_aspect(1.3)
        rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
        fig = plt.gcf()
        fig.set_size_inches(8, 3.7)        
        mp.rcParams.update({'font.size':12})        
        figureName = figurefolder + 'profile_' + QoI +'_DA'+str(EnKFStep)+'_' + caseName + '.pdf'
        plt.savefig(figureName)

                
    def _plotProfileSamples_component(self, QoI, scaleV):
        
        QoIField = np.loadtxt(self.resultDir+QoI+'Field')
        if QoI == 'TKE':
            QoIField = QoIField/Ub/Ub                                
        for isample in np.arange(nSample):
            QoIplot = np.zeros(self.cellidx_v.shape)            
            for i in np.arange(len(self.x_array)):
                idx_locs = self.cellidx_v[:, i]
                idx_locs = idx_locs.tolist()
                QoIplot[:, i] = QoIField[isample, idx_locs]            
            ps = self._plotProfile(self.y_sample, QoIplot, scaleV, '#1b9e76', transparency, lineWidth1)
        
        return ps     

    def _plotProfileSamples_delta(self, QoI, scaleV):
                
        QoIField = np.loadtxt(self.resultDir+QoI+'Field')
        if QoI == 'TKE':
            QoIField = QoIField/Ub/Ub  
        for isample in np.arange(nSample):
                          
            QoIplot = np.zeros(self.cellidx_v.shape)
            for i in np.arange(len(self.x_array)):
                idx_locs = self.cellidx_v[:, i]
                idx_locs = idx_locs.tolist()
                QoIplot[:, i] = QoIField[idx_locs]
            pdb.set_trace()
            ps = self._plotProfile(self.y_sample, QoIplot, scaleV, '#1b9e76', transparency, lineWidth1)
        
        return ps

   
    def _plotRANS(self, QoI, scaleV):
        QoIDict = {'Xi':'Xi', 'Eta':'Eta', 'TKE':'TKE', 'VA':'VA', 'VB':'VB', 'VC':'VC'}
        QoIField = np.loadtxt(dataFolder+'/init/'+QoIDict[QoI]+'_base')
        QoIplot = np.zeros(self.cellidx_v.shape)
        if QoI == 'TKE':
            QoIField = QoIField/Ub/Ub
        for i in np.arange(len(self.x_array)):
            idx_locs = self.cellidx_v[:, i]
            idx_locs = idx_locs.tolist()
            QoIplot[:, i] = QoIField[idx_locs]
        pa = self._plotProfile(self.y_sample, QoIplot, scaleV, 'red', 1, lineWidth1)        
        return pa

    def _plotDNS(self, QoI, scaleV):
        QoIplot = eval('self.'+QoI+'_dns')
        if QoI == 'TKE':
            QoIplot = QoIplot/Ub/Ub        
        pd = self._plotProfile(self.y_sample, QoIplot, scaleV, 'black', 1, 2.0)
        return pd

    def _plotDeltaRANS(self, QoI, scaleV):
        QoIDict = {'deltaXi':'Xi', 'deltaEta':'Eta', 'deltaTKE':'TKE', \
        'deltaVA':'VA', 'deltaVB':'VB', 'deltaVC':'VC'}

        dns = eval('self.'+QoIDict[QoI]+'_dns')
        rans = eval('self.'+QoIDict[QoI]+'_rans')
        if QoI == 'deltaXi' or QoI == 'deltaEta':
            delta_dns2rans = dns - rans
        elif QoI == 'deltaTKE':
            delta_dns2rans = np.log2(dns) - np.log2(rans)
        elif QoI == 'deltaVA' or QoI == 'deltaVB' or QoI == 'deltaVC':
            delta_dns2rans = (dns -rans)*180.0/np.pi
            for i in np.arange(self.ncell_vline):
                for j in np.arange(self.nvline):
                    if (delta_dns2rans[i, j] > 90.0):
                        delta_dns2rans[i, j] = 180.0 - delta_dns2rans[i, j]
                    elif (delta_dns2rans[i, j] < -90.0):
                        delta_dns2rans[i, j] = 180.0 + delta_dns2rans[i, j]  
            delta_dns2rans = delta_dns2rans*np.pi/180.0
        #TODO
        delta_dns2rans = abs(delta_dns2rans)
        np.savetxt('./debugData/dns2rans_'+QoI, delta_dns2rans)
        np.savetxt('./debugData/dns2rans_yline', self.y_sample)
        pd = self._plotProfile(self.y_sample, delta_dns2rans, scaleV, 'black', 1, 2.0)
        return pd
    
    def _plotProfile(self, vCoord, QoI, scaleV, Mcolor, transpa, linew):
        i = 0
        for xPos in self.x_array:
            x = scaleV * QoI[:, i]   + float(xPos)
            y = vCoord[:, i]
            xzero = np.zeros(x.shape) + float(xPos)
            p1, = plt.plot(x, y, color=Mcolor, alpha=transpa, lw = linew, mfc = 'none')
            p2, = plt.plot(xzero, y, '--', color='black', lw = 1.5)
            i = i + 1
        return p1          

    def _plotDomain(self):
        # Plot the simulation domain
        y=np.arange(0, 9, 0.01)
        yext = np.array([9, 9, 0, 0])
        h=hillShape.profile(y)
        hext = np.array([1, 3.036, 3.036, 1])
        y = np.append(y, yext)
        h = np.append(h, hext)
        plt.plot(y,h,'g-')

    def _readDNS_Syn(self):
        """
        get Synthetic DNS components
        """
        synFolder = 'pehill_truth'
        XC1 = np.loadtxt(synFolder + '/XC1_new.dat')
        XC2 = np.loadtxt(synFolder + '/XC2_new.dat')
        Xi = np.loadtxt(synFolder + '/Xi_new.dat')
        Eta = np.loadtxt(synFolder + '/Eta_new.dat')
        k = np.loadtxt(synFolder + '/k_new.dat')
        VA = np.loadtxt(synFolder + '/VA_new.dat')
        VB = np.loadtxt(synFolder + '/VB_new.dat')
        VC = np.loadtxt(synFolder + '/VC_new.dat')

        XC1_base = np.loadtxt(synFolder + '/XC1_base.dat')
        XC2_base = np.loadtxt(synFolder + '/XC2_base.dat')
        Xi_base = np.loadtxt(synFolder + '/Xi_base.dat')
        Eta_base = np.loadtxt(synFolder + '/Eta_base.dat')
        k_base = np.loadtxt(synFolder + '/k_base.dat')
        VA_base = np.loadtxt(synFolder + '/VA_base.dat')
        VB_base = np.loadtxt(synFolder + '/VB_base.dat')
        VC_base = np.loadtxt(synFolder + '/VC_base.dat')
        
        comPool = ['XC1', 'XC2', 'Xi', 'Eta', 'k', 'VA', 'VB', 'VC' ]
        for com in comPool:
            QoIplot = np.zeros(self.cellidx_v.shape)                    
            for i in np.arange(len(self.x_array)):
                idx_locs = self.cellidx_v[:, i]
                idx_locs = idx_locs.tolist()
                exec('QoIplot[:, i] = ' + com + '_base[idx_locs]') 
            exec(com+'_rans = QoIplot')
            
            QoIplot = np.zeros(self.cellidx_v.shape)                    
            for i in np.arange(len(self.x_array)):
                idx_locs = self.cellidx_v[:, i]
                idx_locs = idx_locs.tolist()
                exec('QoIplot[:, i] = ' + com + '[idx_locs]')
            exec(com+'_dns = QoIplot')       

        return XC1_dns, XC2_dns, Xi_dns, Eta_dns, k_dns, VA_dns, VB_dns, VC_dns, \
               XC1_rans, XC2_rans, Xi_rans, Eta_rans, k_rans, VA_rans, VB_rans, VC_rans
        
        
    def _readDNS(self):
        """
        get DNS component (c1, c2, c3, XC1, XC2, Xi, Eta, TKE, VA, VB, VC)
        """
        x_array = np.array([1., 2., 3., 4., 5., 6., 7., 8.])
        XC1_dns = np.zeros((self.ncell_vline, self.nvline))
        XC2_dns = np.zeros((self.ncell_vline, self.nvline))
        Xi_dns =  np.zeros((self.ncell_vline, self.nvline))
        Eta_dns = np.zeros((self.ncell_vline, self.nvline))
        k_dns = np.zeros((self.ncell_vline, self.nvline))
        VA_dns = np.zeros((self.ncell_vline, self.nvline))
        VB_dns = np.zeros((self.ncell_vline, self.nvline))
        VC_dns = np.zeros((self.ncell_vline, self.nvline))
        
        XC1_rans = np.zeros((self.ncell_vline, self.nvline))
        XC2_rans = np.zeros((self.ncell_vline, self.nvline))
        Xi_rans =  np.zeros((self.ncell_vline, self.nvline))
        Eta_rans = np.zeros((self.ncell_vline, self.nvline))
        k_rans = np.zeros((self.ncell_vline, self.nvline))
        VA_rans = np.zeros((self.ncell_vline, self.nvline))
        VB_rans = np.zeros((self.ncell_vline, self.nvline))
        VC_rans = np.zeros((self.ncell_vline, self.nvline))        
        i = 0
        for x in x_array:
            tau_vline = np.loadtxt('sets/X/ransTauFine'+str(int(x)))
            y_vline = tau_vline[:, 0]
            tau_dns_vline = tau_vline[:, 1:7]
           
            tau_rans_vline = np.loadtxt('sets/X/line_x'+str(int(x))+'_Tau.xy')
            y_vline_rans = tau_rans_vline[:, 0]
            tau_rans_vline = tau_rans_vline[:, 1:]
                    
            mapTau = ReynoldsStressRF('None', tau_rans_vline, 400, 1, 'True')
            k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(tau_dns_vline)
            X = mapTau._C2X(C)
            RS = mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
            VA, VB, VC = mapTau.getThetaVABC(tau_dns_vline) # collapse time = 1.005s (3000 cells)        
            
            # interpolate to coarse grid (400 to 31)
            XC1_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, X[:, 0])
            XC2_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, X[:, 1])
            Xi_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, RS[:, 0])
            Eta_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, RS[:, 1])
            k_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, k[:, 0])
            VA_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, VA[:, 0])
            VB_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, VB[:, 0])
            VC_dns[:, i] = np.interp(self.y_sample[:, i], y_vline, VC[:, 0])                         

            k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(tau_rans_vline)
            X = mapTau._C2X(C)
            RS = mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
            VA, VB, VC = mapTau.getThetaVABC(tau_rans_vline) # collapse time = 1.005s (3000 cells)        
            
            # interpolate to coarse grid (400 to 31)
            XC1_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, X[:, 0])
            XC2_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, X[:, 1])
            Xi_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, RS[:, 0])
            Eta_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, RS[:, 1])
            k_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, k[:, 0])
            VA_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, VA[:, 0])
            VB_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, VB[:, 0])
            VC_rans[:, i] = np.interp(self.y_sample[:, i], y_vline, VC[:, 0]) 
          
            i = i + 1
        return XC1_dns, XC2_dns, Xi_dns, Eta_dns, k_dns, VA_dns, VB_dns, VC_dns, \
               XC1_rans, XC2_rans, Xi_rans, Eta_rans, k_rans, VA_rans, VB_rans, VC_rans
    
if __name__ == "__main__":
    dataFolder = './debugData'
    figurefolder = './figures/'
    mainInputFile = './MainInput.in'
    modelInputFile = './forwardModelInput.in'
    plotInputFile = './plotInfo.in'
    paramDict_main = readInputData(mainInputFile)
    paramDict_model = readInputData(modelInputFile)
    paramDict_plot = readInputData(plotInputFile)
    
    nSample = int(paramDict_main['Ns'])    
    problem = paramDict_model['caseName']
    caseName = paramDict_plot['caseName']
    synFlag = ast.literal_eval(paramDict_plot['synFlag'])
    EnKFStep = int(paramDict_plot['EnKFStep'])    
    
    XiEtaFlag = ast.literal_eval(paramDict_plot['XiEtaFlag'])
    KFlag = ast.literal_eval(paramDict_plot['KFlag'])
    VFlag = ast.literal_eval(paramDict_plot['VFlag'])
    deltaXiEtaFlag = ast.literal_eval(paramDict_plot['deltaXiEtaFlag'])
    deltaTKEFlag = ast.literal_eval(paramDict_plot['deltaTKEFlag'])
    deltaVFlag = ast.literal_eval(paramDict_plot['deltaVFlag'])
    ps = postProfile()
    if XiEtaFlag:
        print 'process profile: Xi and Eta'
        ps.mainPlot_component('Xi', 0.8)
        ps.mainPlot_component('Eta', 0.8)     
    if KFlag:
        print 'process profile: Log2K'
        ps.mainPlot_component('TKE', 10)
    if VFlag:
        print 'process profile: VA, VB, VC'
        ps.mainPlot_component('VA', 1)
        ps.mainPlot_component('VB', 1)
        ps.mainPlot_component('VC', 1)
    if deltaXiEtaFlag:
        print 'process profile: deltaXi and deltaEta'
        ps.mainPlot_delta('deltaXi', 0.5)
        ps.mainPlot_delta('deltaEta', 0.5)
    if deltaTKEFlag:
        print 'process profile: deltaTKE'
        ps.mainPlot_delta('deltaTKE', 0.5) 
    if deltaVFlag:
        print 'process profile: deltaVA, deltaVB, deltaVC'
        ps.mainPlot_delta('deltaVA', 0.5)
        ps.mainPlot_delta('deltaVB', 0.5)
        ps.mainPlot_delta('deltaVC', 0.5)        
                        
