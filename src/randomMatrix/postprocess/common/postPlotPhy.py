#!/aoe/bin/python27

# description        :Contain all functions for postprocessing the result of Physics-based UQ

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
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
from scipy import interpolate
# Import local modules
import postFuns as pF
from utilities import readInputData, extractListFromDict, replace
import hillShape
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

class postPlotEnKF:

    """

    Arg:


    """

    def __init__(self, mainInput, onePtCoord):
        paramDict = readInputData(mainInputFile)
        self.nSample = int(paramDict['Ns'])
        self.resultDir = 'debugData/init/'
        self.onePtCoord = onePtCoord
        self.caseName = 'pehill_base'
        self.C1Fields = np.loadtxt(self.resultDir+'XC1_s')
        self.C2Fields = np.loadtxt(self.resultDir+'XC2_s')
        self.C1Field_base = np.loadtxt(self.resultDir+'XC1_base')
        self.C2Field_base = np.loadtxt(self.resultDir+'XC2_base')
        meshcase = self.caseName+'_mesh'
        self.cellidx = pF.coord2cellIdx(onePtCoord, meshcase)
        hCoord, vCoord = pF.getVerandHorCoord(onePtCoord[0, 0], onePtCoord[0, 1], onePtCoord[0, 2], meshcase)
        self.cellidx_h = pF.coord2cellIdx(hCoord, meshcase)
        self.cellidx_v = pF.coord2cellIdx(vCoord, meshcase)
        # get DNS component
        #self.C_dns, self.XC_dns, self.k_dns, self.VA_dns, self.VB_dns, self.VC_dns = self._readDNS()
    # Plot one location statistics
    def plotTau1LocSamplesOnBay(self, caseName):
        """
        plot ensemble of one location of Tau fields on Baycentric triangle
        :return:
        """
        plt.clf()
        C1meanField, C1varField, C1stdField = pF.statEva_scalarField(self.C1Fields)
        C2meanField, C2varField, C2stdField = pF.statEva_scalarField(self.C2Fields)
        idx_loc = int(self.cellidx[1])
        asample, = plt.plot(self.C1Fields[:, idx_loc], self.C2Fields[:, idx_loc], 'o', markersize=6, color='#1b9e76',alpha=transparency)
        abase, = plt.plot(self.C1Field_base[idx_loc], self.C2Field_base[idx_loc], 'o',markersize=10, c='red')
        amean, = plt.plot(C1meanField[idx_loc], C2meanField[idx_loc], 'o',markersize=10, c='yellow', alpha=1)
        
        xDns = 0.5302275; yDns = 0.54465
        #xDns = 0.16; yDns = 0.00001
        aDNS, = plt.plot(xDns, yDns, 'D',markersize=10, c='blue')        
        if legendFlag == 'True':
            plt.annotate(r'Benchmark', xy=(xDns, yDns), xycoords='data',
                     xytext=(150, -100), textcoords='offset points', size=20,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.1'))   
        
        if compareDNSFlag == 'True':
            aDNS, = plt.plot(self.XC_dns[0], self.XC_dns[1], 'D',markersize=10, c='blue')
            if legendFlag == 'True':
                plt.annotate(r'Benchmark', xy=(self.XC_dns[0], self.XC_dns[1]), xycoords='data',
                         xytext=(150, 100), textcoords='offset points', size=20,
                         arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.1'))            
        if legendFlag == 'True':    
            plt.legend([asample, amean],["Samples", "Sample Mean"],prop={'size':18},numpoints=1, bbox_to_anchor=(1.1, 1.1), loc=1)
            plt.annotate(r'Baseline', xy=(self.C1Field_base[idx_loc], self.C2Field_base[idx_loc]), xycoords='data',
                     xytext=(150, -100), textcoords='offset points', size=20,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.1'))
                     
        self._plotBaryTriangle()
        if legendFlag == 'True':
            subfig = plt.axes([0.1, 0.7, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)

        fname = figurefolder + 'bay_idx_'+str(idx_loc)+'_Phy_'+ caseName +'.pdf'
        plt.savefig(fname)        


    def plotTau1LocContourOnBay(self, caseName):
        """

        :return:
        """
        plt.clf()
        fig = plt.figure()
        nbins = 100
        line_colors = ('BlueViolet', 'Crimson', 'ForestGreen',
        'Indigo', 'Tomato', 'Maroon')
        idx_loc = int(self.cellidx[1])
        scalarField1 = self.C1Fields[:, idx_loc]
        scalarField2 = self.C2Fields[:, idx_loc]
        density2D = pF.statEva_2scalarFields(scalarField1, scalarField2)
        xPlot = np.linspace(np.min(scalarField1), np.max(scalarField1), nbins)
        yPlot = np.linspace(np.min(scalarField2), np.max(scalarField2), nbins)
        XPlot, YPlot = np.meshgrid(xPlot, yPlot)
        position = np.vstack([XPlot.ravel(), YPlot.ravel()])
        density2DMap = np.reshape(density2D(position).T, XPlot.shape)
        cs = plt.contour(XPlot, YPlot, density2DMap, 10, lw=3.0)
        #plt.clabel(cs,inline=1,fontsize=10)
        self._plotBaryTriangle()
        # subfig = plt.axes([0.15, 0.65, 0.3, 0.25]) # first 2 are location, second two are shape
        # self._plotDomain(subfig)
        fname = figurefolder + 'bayContour_idx_'+str(idx_loc)+'_Phy_'+ caseName +'.pdf'
        plt.savefig(fname)

    def plotScalarStat_1Loc(self, scalarFields, scalarName, caseName, CDFFlag='False',
                            xRange='None', compareGaussian='True', shadedFlag='False',
                            kde_bw=0.3):
        """
        plot PDF and CDF of scalar at one location
        
        :Arg:
            scalarFields:   scalar field to be analyzed
            scalarName:     name of component (string)
            caseName:       name of case
         
        :return:
        """
        plt.clf()
        idx_loc = int(self.cellidx[1])
        meanField, varField, stdField = pF.statEva_scalarField(scalarFields)
        mp.rc('xtick', labelsize=15)
        mp.rc('ytick', labelsize=15)
        sample_1Loc = scalarFields[:, idx_loc]
        # do statistics
        # TODO: bw_method for angle (degree), the bandwidth is too small, to be checked
        #pdf = ss.kde.gaussian_kde(sample_1Loc, bw_method=0.2/sample_1Loc.std(ddof=1))
        pdf = ss.kde.gaussian_kde(sample_1Loc, bw_method=kde_bw)
        scalarMax = np.max(sample_1Loc)
        scalarMin = np.min(sample_1Loc)
        if xRange == 'None':
            scalarMax = np.max(sample_1Loc)+0.1*abs(np.max(sample_1Loc))
            scalarMin = np.min(sample_1Loc)-0.1*abs(np.max(sample_1Loc))
        n_bins = 1000
        span = (scalarMax - scalarMin) / float(n_bins)
        bins = np.linspace(scalarMin, scalarMax, n_bins)
        cdf = np.cumsum(pdf(bins)*span)
        # analytical Gaussian PDF and CDF
        pdf_analytical = ss.norm.pdf(bins, loc=meanField[idx_loc], scale=stdField[idx_loc])
        cdf_analytical = ss.norm.cdf(bins, loc=meanField[idx_loc], scale=stdField[idx_loc])
        # plot figure
        plt.figure(1)
        plt.clf()
        sp, = plt.plot(bins, pdf(bins), lw=lineWidth2, color=colorSample)
        at = plt.axvline(x=meanField[idx_loc], lw=lineWidth2, color=colorMean, ls='--')
        if shadedFlag == 'True':
            xbegin = meanField[idx_loc]-2*stdField[idx_loc]
            xend = meanField[idx_loc]+2*stdField[idx_loc]
            xx = np.linspace(xbegin, xend, 500)
            density_be = pdf(xx)
            density_an_be = ss.norm.pdf(xx, loc=meanField[idx_loc], scale=stdField[idx_loc])        
            self._fillPDFInside(xx, density_be, 'darkred')
        
        # determine if compared to Gaussian
        if compareGaussian == 'True':
            spa, = plt.plot(bins, pdf_analytical, lw=lineWidth2, color=colorGaussian, ls='--', dashes=dashType1)
            if shadedFlag == 'True':
                self._fillPDFInside(xx, density_an_be, 'darkblue', 'black', 0.3, '\\\\')        
        
        if scalarName in ['c1', 'c2', 'c3', 'XC1', 'XC2', 'Xi', 'Eta', 'VA', 'VB', 'VC']:
            scalar_base = np.loadtxt(self.resultDir+scalarName+'_base')
            aBase = plt.axvline(x=scalar_base[idx_loc], lw=lineWidth2, color=colorBase, ls='--', dashes=dashType2)
        
        frame = plt.gca()
        frame.set_ylim([0.0, plt.ylim()[1]])        
        if xRange != 'None':
            frame.set_xlim([xRange[0], xRange[1]])
            frame.xaxis.set_ticks(np.arange(xRange[0], xRange[1], xRange[2]))
            frame.set_ylim([0.0, plt.ylim()[1]])
            
        if legendFlag == 'True':
            if scalarName in ['c1', 'c2', 'c3', 'XC1', 'XC2', 'Xi', 'Eta', 'VA', 'VB', 'VC']:    
                plt.legend([sp, at, aBase],[r'Sample', r'Sample Mean', r'Baseline'],prop={'size':15},numpoints=1)
            else:
                plt.legend([sp, at],[r'Sample', r'Sample Mean'],prop={'size':15},numpoints=1) 
            if compareGaussian == 'True':                
                if scalarName in ['c1', 'c2', 'c3', 'XC1', 'XC2', 'Xi', 'Eta', 'VA', 'VB', 'VC']:    
                    plt.legend([sp, spa, at, aBase],[r'Sample', r'Gaussian', r'Sample Mean', r'Baseline'],prop={'size':15},numpoints=1)
                else:
                    plt.legend([sp, spa, at],[r'Sample', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1) 
        # x and y label

        labelDict = {
                    'c1':r'$C_1$', 'c2':r'$C_2$', 'c3':r'$C_3$',
                    'Xi':r'$\xi$', 'Eta':r'$\eta$', 'TKE':r'$k/U_b^2$',
                    'VA':r'$\varphi_1$ (degree)', 'VB':r'$\varphi_2$ (degree)', 'VC':r'$\varphi_3$ (degree)',
                    'deltaXi':r'$\Delta\xi$', 'deltaEta':r'$\Delta\eta$', 'deltaK':r'$\Delta$ln$k$',
                    'deltaVA':r'$\Delta\varphi_1$ (degree)', 'deltaVB':r'$\Delta\varphi_2$ (degree)', 'deltaVC':r'$\Delta\varphi_3$ (degree)', 
                    }
        plt.xlabel(labelDict[scalarName], fontsize=20)
        plt.ylabel(r'PDF', fontsize=20)        
        
        fname = figurefolder + scalarName+'_pdf_'+str(idx_loc)+'_Phy_'+ caseName +'.pdf'
        plt.savefig(fname)
    
    def log2kTologek(self, deltaLog2k):
        """
        convert log2k to logek
        """
        kRatiokbase = 2**deltaLog2k
        deltaLogek = np.log(kRatiokbase)
        return deltaLogek

########################################## Private functions ##########################################
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
        pass
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
        pdb.set_trace()
        return C_dns, XC_dns, k_dns, VA_dns, VB_dns, VC_dns
        
        
                
if __name__ == '__main__':
    mainInputFile = './MainInput.in'
    plotInputFile = './plotInfo.in'
    resultDir = 'debugData/init/'
    paramDict = readInputData(plotInputFile)
    caseName = paramDict['caseName']
    # parse plot control
    plotAllTauComponents = ast.literal_eval(paramDict['plotAllTauComponents'])
    plotPtBary = ast.literal_eval(paramDict['plotPtBary'])
    plotVlineBary = ast.literal_eval(paramDict['plotVlineBary'])
    plotHlineBary = ast.literal_eval(paramDict['plotHlineBary'])
    plotAllDeltaComponents = ast.literal_eval(paramDict['plotAllDeltaComponents'])
    plotAllNonDeltaComponents = ast.literal_eval(paramDict['plotAllNonDeltaComponents'])
    plotC = ast.literal_eval(paramDict['plotC'])
    plotXiEta = ast.literal_eval(paramDict['plotXiEta'])
    plotTKE = ast.literal_eval(paramDict['plotTKE'])
    plotTheta = ast.literal_eval(paramDict['plotTheta'])
    plotDeltaXiEta = ast.literal_eval(paramDict['plotDeltaXiEta'])
    plotDeltaLnK = ast.literal_eval(paramDict['plotDeltaLnK'])
    plotDeltaTheta = ast.literal_eval(paramDict['plotDeltaTheta'])
    plotDeltaC = ast.literal_eval(paramDict['plotDeltaC'])
    # parse plot coefficient
    ns_VH = int(paramDict['ns_VH'])
    onePtCoord = extractListFromDict(paramDict, 'onePtCoord')
    onePtCoord = np.array([[float(pn) for pn in onePtCoord]])
    pltPost = postPlotEnKF(mainInputFile, onePtCoord)
    compareGaussianFlag = paramDict['compareGaussian']
    compareDNSFlag = paramDict['compareDNSFlag']
    legendFlag = paramDict['legendFlag']
    shadedFlag = paramDict['shadedFlag']
    CDFFlag = paramDict['CDFFlag']
    
    if plotAllTauComponents or plotPtBary:
        # Scattering plot of Tau at one location on Barycentric triangle        
        pltPost.plotTau1LocSamplesOnBay(caseName)
        # Contour plot of Tau at one location on Barycentric triangle
        pltPost.plotTau1LocContourOnBay(caseName)
    
    if plotAllTauComponents or plotAllNonDeltaComponents or plotC:
        xRange = [0, 1.1, 0.2] # x- range for c1, c2, and c3 [min, max, interval] 
        c1Fields = np.loadtxt(resultDir+'c1_s')
        print 'processing the component: c1'
        pltPost.plotScalarStat_1Loc(c1Fields, 'c1', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
        c2Fields = np.loadtxt(resultDir+'c2_s')
        print 'processing the component: c2'
        pltPost.plotScalarStat_1Loc(c2Fields, 'c2', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
        c3Fields = np.loadtxt(resultDir+'c3_s')
        print 'processing the component: c3'
        pltPost.plotScalarStat_1Loc(c3Fields, 'c3', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)               
    
    if plotAllTauComponents or plotAllNonDeltaComponents or plotXiEta:
        xRange = 'None' # x- range for Xi and Eta [min, max, interval]
        XiFields = np.loadtxt(resultDir+'Xi_s')
        print 'processing the component: Xi'
        pltPost.plotScalarStat_1Loc(XiFields, 'Xi', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
        EtaFields = np.loadtxt(resultDir+'Eta_s')
        print 'processing the component: Eta'
        pltPost.plotScalarStat_1Loc(EtaFields, 'Eta', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
    
    if plotAllTauComponents or plotAllNonDeltaComponents or plotTKE:
        xRange = 'None' # x- range for c1, c2, and c3 [min, max, interval]
        TKEFields = np.loadtxt(resultDir+'TKE_s')
        # normalization
        Ub = 0.028
        TKEFields = TKEFields / Ub / Ub
        print 'processing the component: TKE'
        pltPost.plotScalarStat_1Loc(TKEFields, 'TKE', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
    
    if plotAllTauComponents or plotAllNonDeltaComponents or plotTheta:
        xRange = [0, 180, 30] # x- range for c1, c2, and c3 [min, max, interval] 
        VAFields = np.loadtxt(resultDir+'VA_s'); VAFields = 180.0*VAFields/np.pi
        VBFields = np.loadtxt(resultDir+'VB_s'); VBFields = 180.0*VBFields/np.pi
        VCFields = np.loadtxt(resultDir+'VC_s'); VCFields = 180.0*VCFields/np.pi
        print 'processing the component: VA'
        pltPost.plotScalarStat_1Loc(VAFields, 'VA', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
        print 'processing the component: VB'
        pltPost.plotScalarStat_1Loc(VBFields, 'VB', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag) 
        print 'processing the component: VC'
        pltPost.plotScalarStat_1Loc(VCFields, 'VC', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
    
    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaXiEta:
        xRange = 'None' # x- range for deltaXi and deltaEta [min, max, interval]
        deltaXiFields = np.loadtxt(resultDir+'deltaXi_s')
        print 'processing the component: deltaXi'
        pltPost.plotScalarStat_1Loc(deltaXiFields, 'deltaXi', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
        deltaEtaFields = np.loadtxt(resultDir+'deltaEta_s')
        print 'processing the component: deltaEta'
        pltPost.plotScalarStat_1Loc(deltaEtaFields, 'deltaEta', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)
    
    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaLnK:
        xRange = 'None' # x- range for deltaLnK [min, max, interval]
        kde_bw = 0.1
        deltaKFields = np.loadtxt(resultDir+'deltaK_s')
        #deltaKFields = pltPost.log2kTologek(deltaKFields)
        print 'processing the component: deltaLnK'
        pltPost.plotScalarStat_1Loc(deltaKFields, 'deltaK', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag, kde_bw)
    
    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaTheta:
        xRange = [-90, 91, 30] # x- range for deltaLnK [min, max, interval]
        deltaVAFields = np.loadtxt(resultDir+'deltaVA_s'); deltaVAFields = 180.0*deltaVAFields/np.pi
        deltaVBFields = np.loadtxt(resultDir+'deltaVB_s'); deltaVBFields = 180.0*deltaVBFields/np.pi
        deltaVCFields = np.loadtxt(resultDir+'deltaVC_s'); deltaVCFields = 180.0*deltaVCFields/np.pi
        print 'processing the component: deltaVA'
        pltPost.plotScalarStat_1Loc(deltaVAFields, 'deltaVA', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)        
        print 'processing the component: deltaVB'
        pltPost.plotScalarStat_1Loc(deltaVBFields, 'deltaVB', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag) 
        print 'processing the component: deltaVC'
        pltPost.plotScalarStat_1Loc(deltaVCFields, 'deltaVC', caseName, CDFFlag, xRange, compareGaussianFlag, shadedFlag)                                       
