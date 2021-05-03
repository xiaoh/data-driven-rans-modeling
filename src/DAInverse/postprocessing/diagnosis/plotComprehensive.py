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

transparency = 0.2
lineWidth1 = 2.0
lineWidth2 = 3.0
dashType1 = (1, 2, 2)
dashType2 = (9, 2, 2, 2)
colorSample = 'blue'
colorMean = 'black'
colorBase = 'darkred'
colorGaussian = 'darkgreen'
figurefolder = './figures/'

class postPlot:

    """
    Arg:


    """

    def __init__(self, onePtCoord):
        
        self.onePtCoord = onePtCoord # one pt coordinate
        DAstep = "%.1f"%EnKFStep
        self.resultDir = dataFolder+'/DA-'+DAstep + '/'
        self.XC1Fields = np.loadtxt(self.resultDir+'XC1_s')
        self.XC2Fields = np.loadtxt(self.resultDir+'XC2_s')
        self.XC1Field_base = np.loadtxt(dataFolder+'/init/'+'XC1_base')
        self.XC2Field_base = np.loadtxt(dataFolder+'/init/'+'XC2_base')
        meshcase = problem +'_mesh'
        self.cellidx = pF.coord2cellIdx(onePtCoord, meshcase)
        hCoord, vCoord = pF.getVerandHorCoord(onePtCoord[0, 0], onePtCoord[0, 1], onePtCoord[0, 2], meshcase)
        self.cellidx_h = pF.coord2cellIdx(hCoord, meshcase)
        self.cellidx_v = pF.coord2cellIdx(vCoord, meshcase)
        if synFlag:
            self.C_dns, self.XC_dns, self.k_dns, self.VA_dns, self.VB_dns, self.VC_dns, \
            self.C_dnsV, self.XC_dnsV, self.k_dnsV, self.VA_dnsV, self.VB_dnsV, self.VC_dnsV, self.XC_dnsH = self._readDNS_Syn()            
        else:
            self.C_dns, self.XC_dns, self.k_dns, self.VA_dns, self.VB_dns, self.VC_dns, \
            self.C_dnsV, self.XC_dnsV, self.k_dnsV, self.VA_dnsV, self.VB_dnsV, self.VC_dnsV= self._readDNS()
                                     
    def plotTau1LocSamplesOnBay(self, caseName):    
        """
        plot ensemble of one location of Tau fields on Baycentric triangle
        :return:
        """
        plt.clf()
        XC1meanField, XC1varField, XC1stdField = pF.statEva_scalarField(self.XC1Fields)
        XC2meanField, XC2varField, XC2stdField = pF.statEva_scalarField(self.XC2Fields)
        idx_loc = int(self.cellidx[1])
        asample, = plt.plot(self.XC1Fields[:, idx_loc], self.XC2Fields[:, idx_loc], 'o', markersize=6, color='#1b9e76',alpha=transparency)
        abase, = plt.plot(self.XC1Field_base[idx_loc], self.XC2Field_base[idx_loc], 'o',markersize=10, c='red')
        amean, = plt.plot(XC1meanField[idx_loc], XC2meanField[idx_loc], 'o',markersize=10, c='yellow', alpha=1)
        
        aDNS, = plt.plot(self.XC_dns[0], self.XC_dns[1], 'D',markersize=10, c='blue') 
        if legendFlag == True:    
            plt.legend([asample, amean, abase, aDNS],["Samples", "Sample Mean", "Baseline", "Benchmark"],prop={'size':16},numpoints=1, bbox_to_anchor=(1.1, 1.1), loc=1) 
                    
        self._plotBaryTriangle()
        if legendFlag == True:
            subfig = plt.axes([0.1, 0.7, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)        
        fname = figurefolder + 'bay_idx_'+str(idx_loc) +'_DA'+str(EnKFStep)+'_' + caseName +'.pdf'
        plt.savefig(fname) 

    def plotHlineonBay(self, nSample, caseName):
        """

        :param nSample:
        :return:
        """
        idx_loc = self.cellidx[1]
        cellidxs = self.cellidx_h
        #pdb.set_trace()
        figureName = 'idx_'+str(idx_loc)+'_hline_onBay' +'_DA'+str(EnKFStep)+'_' + caseName
        self.plotMultiLocsOnBay(cellidxs, nSample, figureName)

    def plotVlineonBay(self, nSample, caseName):
        """

        :param nSample:
        :return:
        """
        idx_loc = self.cellidx[1]
        cellidxs = self.cellidx_v

        figureName = 'idx_'+str(idx_loc)+'_vline_onBay_'+'_DA'+str(EnKFStep)+'_' +caseName
        self.plotMultiLocsOnBay(cellidxs, nSample, figureName, True)

    def plotMultiLocsOnBay(self, cellidxs, nSample, figureName, vline=False):
        """

        :param cellidxs:
        :return:
        """
        plt.clf()
        idx_locs = cellidxs[:, 1].tolist()
        for isample in np.arange(nSample):
            asample, = plt.plot(self.XC1Fields[isample, idx_locs], self.XC2Fields[isample, idx_locs], '-o', markersize=5,alpha=transparency)
            plt.plot(self.XC1Fields[isample, idx_locs[0]], self.XC2Fields[isample, idx_locs[0]], '*', markersize=8, color='k')
            plt.plot(self.XC1Fields[isample, idx_locs[-1]], self.XC2Fields[isample, idx_locs[-1]], '^', markersize=8, color='blue')
        am, = plt.plot(self.XC1Field_base[idx_locs], self.XC2Field_base[idx_locs], '-o',markersize=5, c='red')
        if vline == True:
            aDns, = plt.plot(self.XC_dnsV[:,0], self.XC_dnsV[:,1], '-s',markersize=5, c='blue')
            plt.plot(self.XC_dnsV[0,0], self.XC_dnsV[0,1], '*', markersize=8, color='red')
            plt.plot(self.XC_dnsV[-1,0], self.XC_dnsV[-1,1], '^', markersize=8, color='red')
        elif synFlag:
            aDns, = plt.plot(self.XC_dnsH[:,0], self.XC_dnsH[:,1], '-s',markersize=5, c='blue')
            plt.plot(self.XC_dnsH[0,0], self.XC_dnsH[0,1], '*', markersize=8, color='red')
            plt.plot(self.XC_dnsH[-1,0], self.XC_dnsH[-1,1], '^', markersize=8, color='red')
        
        #plt.legend([am, asample],["Baseline", "Samples"],prop={'size':10},numpoints=1)
        self._plotBaryTriangle()
        fname = figurefolder + figureName +'.pdf'
        if legendFlag == True:
            subfig = plt.axes([0.1, 0.7, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)   
        plt.savefig(fname)
        
    def _plotDomain(self, figHand):
        """
        plot domain shape of periodic hills
        """
        y=np.arange(0, 9, 0.01)
        yext = np.array([9, 9, 0, 0])
        h = hillShape.profile(y)
        hext = np.array([1, 3.036, 3.036, 1])
        y = np.append(y, yext)
        h = np.append(h, hext)
        ax, = plt.plot(y, h, '-', color='darkgreen', lw=2.0)
        figHand.plot(self.onePtCoord[0,0], self.onePtCoord[0, 1], 'o', markersize=8, color='k')
        #plt.text(self.onePtCoord[0,0]+2.0, self.onePtCoord[0,1]+1.0, r'(x='+str(self.onePtCoord[0,0])+', y='+str(self.onePtCoord[0,1])+')',
        #        horizontalalignment='center', verticalalignment='center', fontsize=16)
        gap = 0.01
        pylab.xticks([])
        pylab.yticks([])
        plt.axis('off')
        pylab.xlim([0-gap, 9.0+gap])
        pylab.ylim([0-0.5, 3.0+0.5])
        figHand.set_aspect('equal')
        
    def _plotBaryTriangle(self):
        """
        plot Barycentric triangle with annotated texts
        """

        plt.plot([0,1,0.5,0.5,0], [0,0,3**0.5/2.0,3**0.5/2.0,0], 'k-', lw=2.0)
        plt.text(-0.0, -0.05, r'2-Comp', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(1.0, -0.05, r'1-Comp', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(0.5, 0.5*np.sqrt(3)+0.05, r'3-Comp', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(0.2, 0.25*np.sqrt(3), r'$C_1$', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(0.8, 0.25*np.sqrt(3), r'$C_2$', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(0.5, -0.05, r'$C_3$', horizontalalignment='center', verticalalignment='center', fontsize=20)
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        frame.axes.get_yaxis().set_ticks([])
        plt.axis('off')
        gap = 0.01
        frame.set_xlim([0-gap, 1.0+gap])
        frame.set_ylim([0-gap, 0.5*np.sqrt(3)+gap])


    def _readDNS(self):
        """
        get DNS component (c1, c2, c3, XC1, XC2, Xi, Eta, k, VA, VB, VC)
        """
        x = self.onePtCoord[0, 0]; y = self.onePtCoord[0, 1];
        tau_vline = np.loadtxt('sets/X/ransTauFine'+str(int(x)))
        y_vline = tau_vline[:, 0]
        tau_dns_vline = tau_vline[:, 1:7]
        
        idxs = np.where(abs(y_vline-y)<5e-3)
        tau_dns = tau_dns_vline[idxs[0][0], :]
        
        tau_rans_vline = np.loadtxt('sets/X/line_x'+str(int(x))+'_Tau.xy')
        y_vline_rans = tau_rans_vline[:, 0]
        tau_rans_vline = tau_rans_vline[:, 1:]
        tau_rans = tau_rans_vline[idxs[0][0], :]
        
        mapTau = ReynoldsStressRF('None', tau_rans_vline, 31, 1, 'True')
        k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(tau_dns_vline)
        X = mapTau._C2X(C)
        RS = mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
        VA, VB, VC = mapTau.getThetaVABC(tau_dns_vline) # collapse time = 1.005s (3000 cells)
        
        C_dns = C[idxs[0][0], :]
        XC_dns = X[idxs[0][0], :]
        k_dns = k[idxs[0][0], 0]        
        VA_dns = VA[idxs[0][0], 0]
        VB_dns = VB[idxs[0][0], 0]
        VC_dns = VC[idxs[0][0], 0]

        return C_dns, XC_dns, k_dns, VA_dns, VB_dns, VC_dns, C, X, k, VA, VB, VC

    def _readDNS_Syn(self):
        """
        get Synthetic DNS components
        """
        idx_loc = self.cellidx[1]
        cellidxs = self.cellidx_v
        idx_locs = cellidxs[:, 1].tolist()
        cellidxs_H = self.cellidx_h
        idx_locs_H = cellidxs_H[:, 1].tolist()
        
        synFolder = 'pehill_truth'
        XC1 = np.loadtxt(synFolder + '/XC1_new.dat')
        XC2 = np.loadtxt(synFolder + '/XC2_new.dat')
        Xi = np.loadtxt(synFolder + '/Xi_new.dat')
        Eta = np.loadtxt(synFolder + '/Eta_new.dat')
        k = np.loadtxt(synFolder + '/k_new.dat')
        VA = np.loadtxt(synFolder + '/VA_new.dat')
        VB = np.loadtxt(synFolder + '/VB_new.dat')
        VC = np.loadtxt(synFolder + '/VC_new.dat')

        C_dns = []
        C = []
        
        XC_dns = [XC1[idx_loc], XC2[idx_loc]]
        X = np.vstack((XC1[idx_locs], XC2[idx_locs])).T
        X_H = np.vstack((XC1[idx_locs_H], XC2[idx_locs_H])).T        
        k_dns = k[idx_loc]
        k = k[idx_locs]
        
        VA_dns = VA[idx_loc]
        VA_ = VA[idx_locs]

        VB_dns = VB[idx_loc]
        VB_ = VB[idx_locs]

        VC_dns = VC[idx_loc]
        VC_ = VC[idx_locs]                

        return C_dns, XC_dns, k_dns, VA_dns, VB_dns, VC_dns, C, X, k, VA, VB, VC, X_H

                        
if __name__ == '__main__':        
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
    plotAll = ast.literal_eval(paramDict_plot['plotAll'])
    plotPtBary = ast.literal_eval(paramDict_plot['plotPtBary'])
    plotLineBary = ast.literal_eval(paramDict_plot['plotLineBary'])
    legendFlag = ast.literal_eval(paramDict_plot['legendFlag'])
    onePtCoord = extractListFromDict(paramDict_plot, 'onePtCoord')
    onePtCoord = np.array([[float(pn) for pn in onePtCoord]])        
    
    post = postPlot(onePtCoord)
    
    if plotAll or plotPtBary:
        post.plotTau1LocSamplesOnBay(caseName)
    
    if plotAll or plotLineBary:
        post.plotHlineonBay(10, caseName)
        post.plotVlineonBay(10, caseName)
    
    
    
