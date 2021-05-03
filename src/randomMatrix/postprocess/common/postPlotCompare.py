#!/aoe/bin/python27

# description        :comparing physics-based and RMT-based results

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Jan.12, 2016
# revision           :Jan.12, 2016
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
transparency1 = 0.25
lineWidth1 = 2.0
lineWidth2 = 3.0
dashType1 = (1, 2, 2)
dashType2 = (9, 2, 2, 2)
colorSample_RMT = 'blue'
colorSample_Phy = 'magenta'
colorMean = 'black'
colorBase = 'darkred'
colorGaussian = 'darkgreen'
figurefolder = './figures/'

class postPlotCompare:

    """

    Arg:


    """

    def __init__(self, mainInput, onePtCoord):
        paramDict = readInputData(mainInputFile)
        self.case = paramDict['caseName']
        self.nSample = int(paramDict['nSample'])
        self.onePtCoord = onePtCoord
        self.nCell = int(paramDict['nCell_cfd'])
        # get cell idx to be analyzed 
        meshcase = self.case+'_mesh'
        self.cellidx = pF.coord2cellIdx(onePtCoord, meshcase)
        hCoord, vCoord = pF.getVerandHorCoord(onePtCoord[0, 0], onePtCoord[0, 1], onePtCoord[0, 2], meshcase)
        self.cellidx_h = pF.coord2cellIdx(hCoord, meshcase)
        self.cellidx_v = pF.coord2cellIdx(vCoord, meshcase)

        # load XC1 XC2 for RMT
        self.XC1Fields_RMT = np.loadtxt(resultDir_RMT + RComponentDir + 'XC1_s')
        self.XC2Fields_RMT = np.loadtxt(resultDir_RMT + RComponentDir + 'XC2_s')        
        # load XC1 XC2 for Phy
        self.XC1Fields_Phy = np.loadtxt(resultDir_Phy + 'XC1_s')
        self.XC2Fields_Phy = np.loadtxt(resultDir_Phy + 'XC2_s')        
        # load XC1 XC2 for baseline
        self.XC1Field_base = np.loadtxt(resultDir_RMT + RComponentDir + 'XC1_base')
        self.XC2Field_base = np.loadtxt(resultDir_RMT + RComponentDir + 'XC2_base')
        # load natural coordinates
        self.XiFields_RMT = np.loadtxt(resultDir_RMT + RComponentDir + 'Xi_s')
        self.EtaFields_RMT = np.loadtxt(resultDir_RMT + RComponentDir + 'Eta_s') 
        self.XiFields_Phy = np.loadtxt(resultDir_Phy + 'Xi_s')
        self.EtaFields_Phy = np.loadtxt(resultDir_Phy + 'Eta_s')
           
    # Plot one location statistics
    def plotTau1LocSamplesOnBay(self, caseName):
        """
        scattering plot Tau at one location on Baycentric triangle
        
        """
        plt.clf()
        # statistical calculation for XC1 and XC2 for RMT
        XC1meanField_RMT, XC1varField_RMT, XC1stdField_RMT = pF.statEva_scalarField(self.XC1Fields_RMT)
        XC2meanField_RMT, XC2varField_RMT, XC2stdField_RMT = pF.statEva_scalarField(self.XC2Fields_RMT)
        # statistical calculation for XC1 and XC2 for Phy
        XC1meanField_Phy, XC1varField_Phy, XC1stdField_Phy = pF.statEva_scalarField(self.XC1Fields_Phy)
        XC2meanField_Phy, XC2varField_Phy, XC2stdField_Phy = pF.statEva_scalarField(self.XC2Fields_Phy)
        # get location cell index        
        idx_loc = int(self.cellidx[1])
        # plot samples for RMT and Phy
        asample_RMT, = plt.plot(self.XC1Fields_RMT[:, idx_loc], self.XC2Fields_RMT[:, idx_loc], 
                            'o', markersize=6, color='#1b9e76',alpha=transparency1)
        asample_Phy, = plt.plot(self.XC1Fields_Phy[:, idx_loc], self.XC2Fields_Phy[:, idx_loc], 
                            'D', markersize=6, color='red',alpha=transparency)
        # plot baseline                     
        abase, = plt.plot(self.XC1Field_base[idx_loc], self.XC2Field_base[idx_loc], 'o',markersize=10, c='red') 
        # plot sample mean for RMT and Phy        
        amean_RMT, = plt.plot(XC1meanField_RMT[idx_loc], XC2meanField_RMT[idx_loc], 'o',markersize=10, c='yellow', alpha=1)
        amean_Phy, = plt.plot(XC1meanField_Phy[idx_loc], XC2meanField_Phy[idx_loc], 'D',markersize=10, c='blue', alpha=1)               
        
        if legendFlag == 'True':    
            plt.legend([asample_RMT, asample_Phy, amean_RMT, amean_Phy], ["Samples (RMT)", "Samples (Phy)", 
                        "Sample Mean (RMT)", "Sample Mean (Phy)"], prop={'size':16}, numpoints=1, 
                        bbox_to_anchor=(1.12, 1.12), loc=1)
            plt.annotate(r'Baseline', xy=(self.XC1Field_base[idx_loc], self.XC2Field_base[idx_loc]), xycoords='data',
                     xytext=(150, 40), textcoords='offset points', size=20,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.2'))
        # plot Barycentric triangle
        self._plotBaryTriangle()
        if legendFlag == 'True':
            subfig = plt.axes([0.1, 0.7, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)
        fname = figurefolder + 'bay_idx_'+str(idx_loc)+'_compare_'+ caseName +'.pdf'
        plt.savefig(fname)                                 

    def plotTau1LocContourOnBay(self, caseName):
        """
        contour plot Tau at one location on Baycentric triangle

        """
        plt.close("all")
        fig = plt.figure()
        nbins = 100
        idx_loc = int(self.cellidx[1])
        # Plot RMT        
        scalarField1_RMT = self.XC1Fields_RMT[:, idx_loc]
        scalarField2_RMT = self.XC2Fields_RMT[:, idx_loc]
        
        density2D_RMT = pF.statEva_2scalarFields(scalarField1_RMT, scalarField2_RMT)
        xPlot_RMT = np.linspace(np.min(scalarField1_RMT), np.max(scalarField1_RMT), nbins)
        yPlot_RMT = np.linspace(np.min(scalarField2_RMT), np.max(scalarField2_RMT), nbins)
        XPlot_RMT, YPlot_RMT = np.meshgrid(xPlot_RMT, yPlot_RMT)
        position_RMT = np.vstack([XPlot_RMT.ravel(), YPlot_RMT.ravel()])
        density2DMap_RMT = np.reshape(density2D_RMT(position_RMT).T, XPlot_RMT.shape)
        cs_RMT = plt.contour(XPlot_RMT, YPlot_RMT, density2DMap_RMT, 10, lw=3.0, cmap='Greys')
        #plt.clabel(cs,inline=1,fontsize=10)
        # Plot Phy
        scalarField1_Phy = self.XC1Fields_Phy[:, idx_loc]
        scalarField2_Phy = self.XC2Fields_Phy[:, idx_loc]
        
        density2D_Phy = pF.statEva_2scalarFields(scalarField1_Phy, scalarField2_Phy)
        xPlot_Phy = np.linspace(np.min(scalarField1_Phy), np.max(scalarField1_Phy), nbins)
        yPlot_Phy = np.linspace(np.min(scalarField2_Phy), np.max(scalarField2_Phy), nbins)
        XPlot_Phy, YPlot_Phy = np.meshgrid(xPlot_Phy, yPlot_Phy)
        position_Phy = np.vstack([XPlot_Phy.ravel(), YPlot_Phy.ravel()])
        density2DMap_Phy = np.reshape(density2D_Phy(position_Phy).T, XPlot_Phy.shape)
        cs_Phy = plt.contour(XPlot_Phy, YPlot_Phy, density2DMap_Phy, 10, lw=3.0, cmap='GnBu')        
        
        self._plotBaryTriangle()
        
        fname = figurefolder + 'bayContour_idx_'+str(idx_loc)+'_compare_'+ caseName +'.pdf'
        plt.savefig(fname)

    def plot1LocSampleOnNatural(self, caseName):
        """
        scattering plot Tau at one location on Natural Coordinate
        
        """ 
        pass
    
    def plotTau1LocContourOnNatural(self, caseName):
        """
        scattering plot Tau at one location on Natural Coordinate
        
        """    
        plt.close("all")
        fig = plt.figure()
        nbins = 100
        idx_loc = int(self.cellidx[1])
        bw_method = 0.4
        # Plot RMT        
        scalarField1_RMT = self.XiFields_RMT[:, idx_loc]
        scalarField2_RMT = self.EtaFields_RMT[:, idx_loc]
        density2D_RMT = pF.statEva_2scalarFields(scalarField1_RMT, scalarField2_RMT, bw_method)
        xPlot_RMT = np.linspace(np.min(scalarField1_RMT), np.max(scalarField1_RMT), nbins)
        yPlot_RMT = np.linspace(np.min(scalarField2_RMT), np.max(scalarField2_RMT), nbins)
        XPlot_RMT, YPlot_RMT = np.meshgrid(xPlot_RMT, yPlot_RMT)
        position_RMT = np.vstack([XPlot_RMT.ravel(), YPlot_RMT.ravel()])
        density2DMap_RMT = np.reshape(density2D_RMT(position_RMT).T, XPlot_RMT.shape)
        cs_RMT = plt.contour(XPlot_RMT, YPlot_RMT, density2DMap_RMT, 10, lw=4.0, colors='Blue')
        #cs_RMT = plt.contour(XPlot_RMT, YPlot_RMT, density2DMap_RMT, 15, lw=3.0, cmap='Blues')
        #plt.clabel(cs,inline=1,fontsize=10)
        # Plot Phy
        scalarField1_Phy = self.XiFields_Phy[:, idx_loc]
        scalarField2_Phy = self.EtaFields_Phy[:, idx_loc]
        
        density2D_Phy = pF.statEva_2scalarFields(scalarField1_Phy, scalarField2_Phy, bw_method)
        xPlot_Phy = np.linspace(np.min(scalarField1_Phy), np.max(scalarField1_Phy), nbins)
        yPlot_Phy = np.linspace(np.min(scalarField2_Phy), np.max(scalarField2_Phy), nbins)
        XPlot_Phy, YPlot_Phy = np.meshgrid(xPlot_Phy, yPlot_Phy)
        position_Phy = np.vstack([XPlot_Phy.ravel(), YPlot_Phy.ravel()])
        density2DMap_Phy = np.reshape(density2D_Phy(position_Phy).T, XPlot_Phy.shape)
        cs_Phy = plt.contour(XPlot_Phy, YPlot_Phy, density2DMap_Phy, 10, lw=4.0, colors='red', linestyles='dashed')
        #cs_Phy = plt.contour(XPlot_Phy, YPlot_Phy, density2DMap_Phy, 15, lw=3.0, cmap='Reds')
        # plot natural coordinate
        self._plotNatural()        
        lines = [cs_RMT.collections[0], cs_Phy.collections[0]]
        labels = [r'Random Matrix', r'Phyiscs']
        if legendFlag == 'True':
            plt.legend(lines, labels, prop={'size':15},numpoints=1, loc=3)        
        fname = figurefolder + 'naturalContour_idx_'+str(idx_loc)+'_compare_'+ caseName +'.pdf'
        plt.savefig(fname)
                           
    def plotScalarStat_1Loc(self, scalarFields_RMT, scalarFields_Phy, scalarName, caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='True', 
                            shadedFlag='False', kde_bw=0.3):
        """
        plot PDF and CDF of scalar at one location
        
        :Arg:
            scalarFields_RMT:   scalar field from RMT to be analyzed
            scalarFields_Phy:   scalar field from Phy to be analyzed
            scalarName:         name of component (string)
            caseName:           name of case
         
        :return:
        """
        plt.clf()
        idx_loc = int(self.cellidx[1])
        meanField_RMT, varField_RMT, stdField_RMT = pF.statEva_scalarField(scalarFields_RMT)
        meanField_Phy, varField_Phy, stdField_Phy = pF.statEva_scalarField(scalarFields_Phy)
        mp.rc('xtick', labelsize=15)
        mp.rc('ytick', labelsize=15)
        
        # do statistics for RMT
        sample_1Loc_RMT = scalarFields_RMT[:, idx_loc]
        pdf_RMT = ss.kde.gaussian_kde(sample_1Loc_RMT, bw_method=kde_bw)
        scalarMax_RMT = np.max(sample_1Loc_RMT)
        scalarMin_RMT = np.min(sample_1Loc_RMT)
        if xRange == 'None':
            scalarMax_RMT = np.max(sample_1Loc_RMT)+0.1*abs(np.max(sample_1Loc_RMT))
            scalarMin_RMT = np.min(sample_1Loc_RMT)-0.1*abs(np.max(sample_1Loc_RMT))        
        n_bins = 1000
        span_RMT = (scalarMax_RMT - scalarMin_RMT) / float(n_bins)
        bins_RMT = np.linspace(scalarMin_RMT, scalarMax_RMT, n_bins)
        cdf_RMT = np.cumsum(pdf_RMT(bins_RMT)*span_RMT)
        pdf_analytical_RMT = ss.norm.pdf(bins_RMT, loc=meanField_RMT[idx_loc], scale=stdField_RMT[idx_loc])
        cdf_analytical_RMT = ss.norm.cdf(bins_RMT, loc=meanField_RMT[idx_loc], scale=stdField_RMT[idx_loc])        
        # do statistics for RMT
        sample_1Loc_Phy = scalarFields_Phy[:, idx_loc]
        pdf_Phy = ss.kde.gaussian_kde(sample_1Loc_Phy, bw_method=kde_bw)
        scalarMax_Phy = np.max(sample_1Loc_Phy)
        scalarMin_Phy = np.min(sample_1Loc_Phy)
        if xRange == 'None':
            scalarMax_Phy = np.max(sample_1Loc_Phy)+0.1*abs(np.max(sample_1Loc_Phy))
            scalarMin_Phy = np.min(sample_1Loc_Phy)-0.1*abs(np.max(sample_1Loc_Phy))        
        n_bins = 1000
        span_Phy = (scalarMax_Phy - scalarMin_Phy) / float(n_bins)
        bins_Phy = np.linspace(scalarMin_Phy, scalarMax_Phy, n_bins)
        cdf_Phy = np.cumsum(pdf_Phy(bins_Phy)*span_Phy)
        pdf_analytical_Phy = ss.norm.pdf(bins_Phy, loc=meanField_Phy[idx_loc], scale=stdField_Phy[idx_loc])
        cdf_analytical_Phy = ss.norm.cdf(bins_Phy, loc=meanField_Phy[idx_loc], scale=stdField_Phy[idx_loc])
        
        # K-L divergence
        S_RMT = ss.entropy(pdf_RMT(bins_RMT))
        S_Phy = ss.entropy(pdf_Phy(bins_Phy))
        
        plt.figure(1)
        plt.clf()
        sp_RMT, = plt.plot(bins_RMT, pdf_RMT(bins_RMT), lw=lineWidth2, color=colorSample_RMT)
        at_RMT = plt.axvline(x=meanField_RMT[idx_loc], lw=lineWidth2, color=colorSample_RMT, ls='--')

        sp_Phy, = plt.plot(bins_Phy, pdf_Phy(bins_Phy), lw=lineWidth2, color=colorSample_Phy)
        at_Phy = plt.axvline(x=meanField_Phy[idx_loc], lw=lineWidth2, color=colorSample_Phy, ls='--')                

        if scalarName in ['c1', 'c2', 'c3', 'XC1', 'XC2', 'Xi', 'Eta', 'TKE', 'VA', 'VB', 'VC']:
            scalar_base = np.loadtxt(resultDir_RMT + RComponentDir + scalarName+'_base')
            if scalarName == 'TKE':
                scalar_base = scalar_base / Ub / Ub
            aBase = plt.axvline(x=scalar_base[idx_loc], lw=lineWidth2, color=colorBase, ls='--', dashes=dashType2)

        frame = plt.gca()
        frame.set_ylim([0.0, plt.ylim()[1]])        
        if xRange != 'None':
            frame.set_xlim([xRange[0], xRange[1]])
            frame.xaxis.set_ticks(np.arange(xRange[0], xRange[1], xRange[2]))
            frame.set_ylim([0.0, plt.ylim()[1]])

        if legendFlag == 'True':

            if scalarName in ['c1', 'c2', 'c3', 'XC1', 'XC2', 'Xi', 'Eta', 'VA', 'VB', 'VC']:    
                plt.legend([sp_RMT, sp_Phy, at_RMT, at_Phy, aBase],
                            [r'Sample RM', r'Sample Phy', r'Mean RM', r'Mean Phy', r'Baseline'],
                            prop={'size':16},numpoints=1)
            else:
                plt.legend([sp_RMT, sp_Phy, at_RMT, at_Phy],
                            [r'Sample RM', r'Sample Phy', r'Mean RM', r'Mean Phy'],
                            prop={'size':16},numpoints=2) 

            subfig = plt.axes([0.15, 0.6, 0.3, 0.25])
            #subfig = plt.axes([0.58, 0.15, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)
                    
        # x and y label
        labelDict = {
                    'c1':r'$C_1$', 'c2':r'$C_2$', 'c3':r'$C_3$',
                    'Xi':r'$\xi$', 'Eta':r'$\eta$', 'TKE':r'$k/U_b^2$',
                    'VA':r'$\varphi_1$ (degree)', 'VB':r'$\varphi_2$ (degree)', 'VC':r'$\varphi_3$ (degree)',
                    'deltaXi':r'$\Delta\xi$', 'deltaEta':r'$\Delta\eta$', 'deltaK':r'$\Delta$ln$k$',
                    'deltaVA':r'$\Delta\varphi_1$ (degree)', 'deltaVB':r'$\Delta\varphi_2$ (degree)', 'deltaVC':r'$\Delta\varphi_3$ (degree)', 
                    'R11':r'$\tau_{xx}$', 'R12':r'$\tau_{xy}$', 'R13':r'$\tau_{xz}$', 'R22':r'$\tau_{yy}$', 'R23':r'$\tau_{yz}$','R33':r'$\tau_{zz}$' 
                    }
        plt.xlabel(labelDict[scalarName], fontsize=20)
        plt.ylabel(r'PDF', fontsize=20)
        fname = figurefolder + scalarName+'_pdf_'+str(idx_loc)+'_compare_'+ caseName +'.pdf'
        plt.savefig(fname)

    def log2kTologek(self, deltaLog2k):
        """
        convert log2k to logek
        """
        kRatiokbase = 2**deltaLog2k
        deltaLogek = np.log(kRatiokbase)
        return deltaLogek                                                              


    def plotRij_atOnePt(self, caseName):
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
        
                   
        self.plotScalarStat_1Loc(R_11_rmt, R_11_phy, 'R11', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)
           
        self.plotScalarStat_1Loc(R_12_rmt, R_12_phy, 'R12', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)
        self.plotScalarStat_1Loc(R_13_rmt, R_13_phy, 'R13', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)
        self.plotScalarStat_1Loc(R_22_rmt, R_22_phy, 'R22', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)
        self.plotScalarStat_1Loc(R_23_rmt, R_23_phy, 'R23', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)                            
        self.plotScalarStat_1Loc(R_33_rmt, R_33_phy, 'R33', caseName, 
                            CDFFlag='False', xRange='None', compareGaussian='False', 
                            shadedFlag='True', kde_bw=0.3)                                                        

        


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

    def _plotNatural(self):
        """
        plot Barycentric triangle with annotated texts
        """

        plt.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], 'k-', lw=2.0)
        plt.text(-1, -1.1, r'2', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(1, -1.1, r'1', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(1, 1.1, r'4', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.text(-1, 1.1, r'3', horizontalalignment='center', verticalalignment='center', fontsize=20)
        plt.xlabel(r'$\xi$', fontsize=20)
        plt.ylabel(r'$\eta$', fontsize=20)        
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        frame.axes.get_yaxis().set_ticks([])
        #plt.axis('off')
        gap = 0.01
        frame.set_xlim([-1.0-gap, 1.0+gap])
        frame.set_ylim([-1.0-gap, 1.0+gap])        
                             
if __name__ == '__main__':
    mainInputFile = './mainInput.in'
    plotInputFile = './plotInfo.in'
    paramDict = readInputData(plotInputFile)
    case_compared = paramDict['case_compared'] 
    caseName = paramDict['caseName']
    resultDir_RMT = 'resultData/'
    resultDir_Phy = case_compared + '/debugData/init/'
    RComponentDir = 'RComponent_samples/'
    deltaRComDir = 'deltaRComponent_samples/'
    
    # parse plot control
    plotAllTauComponents = ast.literal_eval(paramDict['plotAllTauComponents'])
    plotPtBary = ast.literal_eval(paramDict['plotPtBary'])
    plotPtNatural = ast.literal_eval(paramDict['plotPtNatural'])
    plotVlineBary = ast.literal_eval(paramDict['plotVlineBary'])
    plotHlineBary = ast.literal_eval(paramDict['plotHlineBary'])
    plotAllDeltaComponents = ast.literal_eval(paramDict['plotAllDeltaComponents'])
    plotAllNonDeltaComponents = ast.literal_eval(paramDict['plotAllNonDeltaComponents'])
    plotTau_ij = ast.literal_eval(paramDict['plotTau_ij'])
    plotC = ast.literal_eval(paramDict['plotC'])
    plotXiEta = ast.literal_eval(paramDict['plotXiEta'])
    plotTKE = ast.literal_eval(paramDict['plotTKE'])
    plotTheta = ast.literal_eval(paramDict['plotTheta'])
    plotDeltaXiEta = ast.literal_eval(paramDict['plotDeltaXiEta'])
    plotDeltaLnK = ast.literal_eval(paramDict['plotDeltaLnK'])
    plotDeltaTheta = ast.literal_eval(paramDict['plotDeltaTheta'])
    plotDeltaC = ast.literal_eval(paramDict['plotDeltaC'])
    scalarEntFlag = ast.literal_eval(paramDict['scalarEntFlag'])
    # parse plot coefficient
    ns_VH = int(paramDict['ns_VH'])
    onePtCoord = extractListFromDict(paramDict, 'onePtCoord')
    onePtCoord = np.array([[float(pn) for pn in onePtCoord]])
    compareGaussianFlag = paramDict['compareGaussian']
    legendFlag = paramDict['legendFlag']
    shadedFlag = paramDict['shadedFlag']
    CDFFlag = paramDict['CDFFlag']
    
    print "Initialize plotting class, readinng related file ..."
    pltPost = postPlotCompare(mainInputFile, onePtCoord)

    if plotAllTauComponents or plotPtBary:
        # Scattering plot of Tau at one location on Barycentric triangle 
        print "plotting Barycentric scattering ..."       
        pltPost.plotTau1LocSamplesOnBay(caseName)
        # Contour plot of Tau at one location on Barycentric triangle 
        print "plotting Barycentric contour ..."
        pltPost.plotTau1LocContourOnBay(caseName)
    
        
    if plotAllTauComponents or plotAllNonDeltaComponents or plotC:
        xRange = [0, 1.1, 0.2] # x- range for c1, c2, and c3 [min, max, interval] 
        print 'processing the component: c1'        
        c1Fields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'c1_s')
        c1Fields_Phy = np.loadtxt(resultDir_Phy+'c1_s')

        pltPost.plotScalarStat_1Loc(c1Fields_RMT, c1Fields_Phy, 'c1', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
                                    
        print 'processing the component: c2'
        c2Fields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'c2_s')
        c2Fields_Phy = np.loadtxt(resultDir_Phy+'c2_s')

        pltPost.plotScalarStat_1Loc(c2Fields_RMT, c2Fields_Phy, 'c2', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
        
        print 'processing the component: c3'                                    
        c3Fields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'c3_s')
        c3Fields_Phy = np.loadtxt(resultDir_Phy+'c3_s')
        pltPost.plotScalarStat_1Loc(c3Fields_RMT, c3Fields_Phy, 'c3', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)                                    
    
    if plotAllTauComponents or plotTau_ij:
        print "processing the components: Tau_ij"
        pltPost.plotRij_atOnePt(caseName)
                        
    if plotAllTauComponents or plotAllNonDeltaComponents or plotXiEta:
        xRange = [-1.0, 1.05, 0.2] # x- range for Xi and Eta [min, max, interval]
        print 'processing the component: Xi'
        XiFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'Xi_s')
        XiFields_Phy = np.loadtxt(resultDir_Phy+'Xi_s')
        pltPost.plotScalarStat_1Loc(XiFields_RMT, XiFields_Phy, 'Xi', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
        
        print 'processing the component: Eta'        
        EtaFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'Eta_s')
        EtaFields_Phy = np.loadtxt(resultDir_Phy+'Eta_s')
        pltPost.plotScalarStat_1Loc(EtaFields_RMT, EtaFields_Phy, 'Eta', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
    if plotAllTauComponents or plotAllNonDeltaComponents or plotTKE:
        xRange = 'None' # x- range for TKE [min, max, interval]
        print 'processing the component: TKE'
        TKEFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'TKE_s')
        TKEFields_Phy = np.loadtxt(resultDir_Phy+'TKE_s')
        # normalization
        Ub = 0.028
        TKEFields_RMT = TKEFields_RMT / Ub / Ub
        TKEFields_Phy = TKEFields_Phy / Ub / Ub
        pltPost.plotScalarStat_1Loc(TKEFields_RMT, TKEFields_Phy, 'TKE', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)

    if plotAllTauComponents or plotAllNonDeltaComponents or plotTheta:
        xRange = [0, 180, 30] # x- range for c1, c2, and c3 [min, max, interval]         
        print 'processing the component: VA'
        VAFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'VA_s'); VAFields_RMT = 180.0*VAFields_RMT/np.pi
        VAFields_Phy = np.loadtxt(resultDir_Phy+'VA_s'); VAFields_Phy = 180.0*VAFields_Phy/np.pi                
        pltPost.plotScalarStat_1Loc(VAFields_RMT, VAFields_Phy, 'VA', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
        print 'processing the component: VB'
        VBFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'VB_s'); VBFields_RMT = 180.0*VBFields_RMT/np.pi        
        VBFields_Phy = np.loadtxt(resultDir_Phy+'VB_s'); VBFields_Phy = 180.0*VBFields_Phy/np.pi        
        pltPost.plotScalarStat_1Loc(VBFields_RMT, VBFields_Phy, 'VB', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag) 
        print 'processing the component: VC'
        VCFields_RMT = np.loadtxt(resultDir_RMT+'RComponent_samples/'+'VC_s'); VCFields_RMT = 180.0*VCFields_RMT/np.pi
        VCFields_Phy = np.loadtxt(resultDir_Phy+'VC_s'); VCFields_Phy = 180.0*VCFields_Phy/np.pi        
        pltPost.plotScalarStat_1Loc(VCFields_RMT, VCFields_Phy, 'VC', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)                                    

    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaXiEta:
        xRange = 'None' # x- range for deltaXi and deltaEta [min, max, interval]
        print 'processing the component: deltaXi'
        deltaXiFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaXi_s')
        deltaXiFields_Phy = np.loadtxt(resultDir_Phy+'deltaXi_s')        
        pltPost.plotScalarStat_1Loc(deltaXiFields_RMT, deltaXiFields_Phy, 'deltaXi', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
                                
        print 'processing the component: deltaEta'
        deltaEtaFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaEta_s')
        deltaEtaFields_Phy = np.loadtxt(resultDir_Phy+'deltaEta_s')        
        pltPost.plotScalarStat_1Loc(deltaEtaFields_RMT, deltaEtaFields_Phy, 'deltaEta', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)

    if plotAllTauComponents or plotPtNatural:
        # Contour plot of Tau at one location on Natural coordinate
        print "plotting natural coord contour ..."
        pltPost.plotTau1LocContourOnNatural(caseName)

    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaLnK:
        xRange = 'None' # x- range for deltaLnK [min, max, interval]
        kde_bw = 0.3
        print 'processing the component: deltaLnK'
        deltaKFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaK_s')
        deltaKFields_RMT = pltPost.log2kTologek(deltaKFields_RMT)
        deltaKFields_Phy = np.loadtxt(resultDir_Phy+'deltaK_s')
        deltaKFields_Phy = pltPost.log2kTologek(deltaKFields_Phy)
        pltPost.plotScalarStat_1Loc(deltaKFields_RMT, deltaKFields_Phy, 'deltaK', caseName, 
                                    CDFFlag, xRange, compareGaussianFlag, shadedFlag, kde_bw)    
                                    
    if plotAllTauComponents or plotAllDeltaComponents or plotDeltaTheta:
        xRange = [-90, 91, 30] # x- range for c1, c2, and c3 [min, max, interval]         
        print 'processing the component: deltaVA'
        deltaVAFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaVA_s'); deltaVAFields_RMT = 180.0*deltaVAFields_RMT/np.pi
        deltaVAFields_Phy = np.loadtxt(resultDir_Phy+'deltaVA_s'); deltaVAFields_Phy = 180.0*deltaVAFields_Phy/np.pi                
        pltPost.plotScalarStat_1Loc(deltaVAFields_RMT, deltaVAFields_Phy, 'deltaVA', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)
        print 'processing the component: deltaVB'
        deltaVBFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaVB_s'); deltaVBFields_RMT = 180.0*deltaVBFields_RMT/np.pi
        deltaVBFields_Phy = np.loadtxt(resultDir_Phy+'deltaVB_s'); deltaVBFields_Phy = 180.0*deltaVBFields_Phy/np.pi                
        pltPost.plotScalarStat_1Loc(deltaVBFields_RMT, deltaVBFields_Phy, 'deltaVB', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag) 
        print 'processing the component: deltaVC'
        deltaVCFields_RMT = np.loadtxt(resultDir_RMT+'deltaRComponent_samples/'+'deltaVC_s'); deltaVCFields_RMT = 180.0*deltaVCFields_RMT/np.pi
        deltaVCFields_Phy = np.loadtxt(resultDir_Phy+'deltaVC_s'); deltaVCFields_Phy = 180.0*deltaVCFields_Phy/np.pi                
        pltPost.plotScalarStat_1Loc(deltaVCFields_RMT, deltaVCFields_Phy, 'deltaVC', caseName, CDFFlag, 
                                    xRange, compareGaussianFlag, shadedFlag)                                                                                                 
