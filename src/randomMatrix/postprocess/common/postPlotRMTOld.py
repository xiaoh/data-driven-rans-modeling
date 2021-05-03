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
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
# Import local modules
import postFuns as pF
from utilities import readInputData, extractListFromDict, replace
import hillShape

transparency = 0.2
figurefolder = './figures/'

class postPlot:
    """

    Arg:


    """

    def __init__(self, mainInput, onePtCoord):
        paramDict = readInputData(mainInputFile)
        self.caseName = paramDict['caseName']
        self.nSample = int(paramDict['nSample'])
        resultDir = 'resultData/'
        RComponentDir = 'RComponent_samples/'
        self.onePtCoord = onePtCoord
        ## load C1 C2 related data
        #pdb.set_trace()
        self.C1Fields = np.loadtxt(resultDir+RComponentDir+'C1_s')
        self.C2Fields = np.loadtxt(resultDir+RComponentDir+'C2_s')
        self.C1Field_base = np.loadtxt(resultDir+RComponentDir+'C1_base')
        self.C2Field_base = np.loadtxt(resultDir+RComponentDir+'C2_base')
        meshcase = self.caseName+'_mesh'
        self.cellidx = pF.coord2cellIdx(onePtCoord, meshcase)
        hCoord, vCoord = pF.getVerandHorCoord(onePtCoord[0, 0], onePtCoord[0, 1], onePtCoord[0, 2], meshcase)
        self.cellidx_h = pF.coord2cellIdx(hCoord, meshcase)
        self.cellidx_v = pF.coord2cellIdx(vCoord, meshcase)

    # Plot one location statistics
    def plotTau1LocSamplesOnBay(self, caseName, legendFlag='False'):
        """
        scattering plot Tau at one location on Baycentric triangle
        
        """
        plt.clf()
        C1meanField, C1varField, C1stdField = pF.statEva_scalarField(self.C1Fields)
        C2meanField, C2varField, C2stdField = pF.statEva_scalarField(self.C2Fields)
        idx_loc = int(self.cellidx[1])
        asample, = plt.plot(self.C1Fields[:, idx_loc], self.C2Fields[:, idx_loc], 'o', markersize=6, color='#1b9e76',alpha=transparency)
        abase, = plt.plot(self.C1Field_base[idx_loc], self.C2Field_base[idx_loc], 'o',markersize=10, c='red')
        amean, = plt.plot(C1meanField[idx_loc], C2meanField[idx_loc], 'o',markersize=10, c='yellow', alpha=1)
        if legendFlag == 'True':    
            plt.legend([asample, amean],["Samples", "Sample Mean"],prop={'size':18},numpoints=1, bbox_to_anchor=(1.1, 1.1), loc=1)
            plt.annotate(r'Baseline', xy=(self.C1Field_base[idx_loc], self.C2Field_base[idx_loc]), xycoords='data',
                     xytext=(150, 40), textcoords='offset points', size=20,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.2'))
        
        # TODO:To be delete
        #xDns = 0.16; yDns = 0.00001
        xDns = 0.5302275; yDns = 0.54465
        plt.plot(xDns, yDns, 's', markersize=10, color='blue',alpha=1)
        if legendFlag == 'True':    
            plt.annotate(r'Benchmark', xy=(xDns, yDns), xycoords='data',
                     xytext=(150, 40), textcoords='offset points', size=20,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=5, lw=1.5, connectionstyle='arc3,rad=0.2'))    
        
        
        self._plotBaryTriangle()
        if legendFlag == 'True':
            subfig = plt.axes([0.1, 0.7, 0.3, 0.25]) # first 2 are location, second two are shape
            self._plotDomain(subfig)
        
    
        fname = figurefolder + 'bay_idx_'+str(idx_loc)+'_'+ caseName +'.pdf'
        plt.savefig(fname)

    def plotTau1LocContourOnBay(self, caseName):
        """
        contour plot Tau at one location on Baycentric triangle

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
        fname = figurefolder + 'bayContour_idx_'+str(idx_loc)+'_'+ caseName +'.pdf'
        plt.savefig(fname)

    def plotHlineonBay(self, nSample, caseName):
        """

        :param nSample:
        :return:
        """
        idx_loc = self.cellidx[1]
        cellidxs = self.cellidx_h
        #pdb.set_trace()
        figureName = 'idx_'+str(idx_loc)+'_hline_onBay' + caseName
        self.plotMultiLocsOnBay(cellidxs, nSample, figureName)

    def plotVlineonBay(self, nSample, caseName):
        """

        :param nSample:
        :return:
        """
        idx_loc = self.cellidx[1]
        cellidxs = self.cellidx_v

        figureName = 'idx_'+str(idx_loc)+'_vline_onBay_'+caseName
        self.plotMultiLocsOnBay(cellidxs, nSample, figureName)


    def plotMultiLocsOnBay(self, cellidxs, nSample, figureName):
        """

        :param cellidxs:
        :return:
        """
        plt.clf()
        idx_locs = cellidxs[:, 1].tolist()
        for isample in np.arange(nSample):
            asample, = plt.plot(self.C1Fields[isample, idx_locs], self.C2Fields[isample, idx_locs], '-o', markersize=5,alpha=transparency)
            plt.plot(self.C1Fields[isample, idx_locs[0]], self.C2Fields[isample, idx_locs[0]], '*', markersize=5, color='k')
            plt.plot(self.C1Fields[isample, idx_locs[-1]], self.C2Fields[isample, idx_locs[-1]], '^', markersize=5, color='blue')
        am, = plt.plot(self.C1Field_base[idx_locs], self.C2Field_base[idx_locs], '-o',markersize=5, c='red')
        #plt.legend([am, asample],["Baseline", "Samples"],prop={'size':10},numpoints=1)
        self._plotBaryTriangle()
        fname = figurefolder + figureName +'.pdf'
        plt.savefig(fname)


    def plotScalarStat_1Loc(self, scalarFields, scalarName, caseName, 
                            xRange='None', compareGaussian='True', shadedFlag='False'):
        """

        :param onePtCoord:
        :return:
        """
        plt.clf()
        idx_loc = int(self.cellidx[1])
        meanField, varField, stdField = pF.statEva_scalarField(scalarFields)
        #pdb.set_trace()
        mp.rc('xtick', labelsize=15)
        mp.rc('ytick', labelsize=15)
        sample_1Loc = scalarFields[:, idx_loc]
        # do statistics
        pdf = ss.kde.gaussian_kde(sample_1Loc)
        #pdb.set_trace()
        scalarMax = np.max(sample_1Loc)
        scalarMin = np.min(sample_1Loc)
        if xRange == 'None':
            scalarMax = np.max(sample_1Loc)+0.1*np.max(sample_1Loc)
            scalarMin = np.min(sample_1Loc)-0.2*abs(np.max(sample_1Loc))
        n_bins = 1000
        span = (scalarMax - scalarMin) / float(n_bins)
        bins = np.linspace(scalarMin, scalarMax, n_bins)
        cdf = np.cumsum(pdf(bins)*span)
        pdf_analytical = ss.norm.pdf(bins, loc=meanField[idx_loc], scale=stdField[idx_loc])
        cdf_analytical = ss.norm.cdf(bins, loc=meanField[idx_loc], scale=stdField[idx_loc])
        #pdb.set_trace()
        plt.figure(1)
        plt.clf()
        sp, = plt.plot(bins, pdf(bins), 'red', lw=3.0)
        at = plt.axvline(x=meanField[idx_loc], lw=3.0, color='black', ls='--')
        xbegin = meanField[idx_loc]-2*stdField[idx_loc]
        xend = meanField[idx_loc]+2*stdField[idx_loc]
        xx = np.linspace(xbegin, xend, 500)
        density_be = pdf(xx)
        density_an_be = ss.norm.pdf(xx, loc=meanField[idx_loc], scale=stdField[idx_loc])
        if shadedFlag == 'True':
            self._fillPDFInside(xx, density_be, 'darkred')
        if compareGaussian == 'True':
            spa, = plt.plot(bins, pdf_analytical, 'blue', ls='--', lw=3.0, dashes=(9, 2, 2, 2))
            if shadedFlag == 'True':
                self._fillPDFInside(xx, density_an_be, 'darkblue', 'black', 0.3, '\\\\')
            if legendFlag == 'True':
                if scalarName == 'c1':
                    c1_base = np.loadtxt('./resultData/RComponent_samples/c1_base')
                    aBase = plt.axvline(x=c1_base[idx_loc], lw=3.0, color='darkred', ls='--', dashes=(1, 2, 2))
                    plt.legend([sp, spa, at, aBase],[r'$C_1$', r'Gaussian', r'Sample Mean', r'Baseline'],prop={'size':15},numpoints=1)
                elif scalarName == 'c2':
                    c2_base = np.loadtxt('./resultData/RComponent_samples/c2_base')
                    aBase = plt.axvline(x=c2_base[idx_loc], lw=3.0, color='darkred', ls='--', dashes=(1, 2, 2))             
                    plt.legend([sp, spa, at, aBase],[r'$C_2$', r'Gaussian', r'Sample Mean', r'Baseline'],prop={'size':15},numpoints=1)               
                elif scalarName == 'c3':
                    c3_base = np.loadtxt('./resultData/RComponent_samples/c3_base')
                    aBase = plt.axvline(x=c3_base[idx_loc], lw=3.0, color='darkred', ls='--', dashes=(1, 2, 2))            
                    plt.legend([sp, spa, at, aBase],[r'$C_3$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)                          
                elif scalarName == 'Xi':
                    plt.legend([sp, spa, at],[r'$\Delta\xi$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'Eta':
                    plt.legend([sp, spa, at],[r'$\Delta\eta$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'K':
                    plt.legend([sp, spa, at],[r'$\Delta$ln$k$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'VA':
                    plt.legend([sp, spa, at],[r'$\Delta\varphi_1$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'VB':
                    plt.legend([sp, spa, at],[r'$\Delta\varphi_2$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'VC':
                    plt.legend([sp, spa, at],[r'$\Delta\varphi_3$', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)
                elif scalarName == 'TKE':
                    k_bar = np.loadtxt('./resultData/RComponent_samples/k_bar')
                    aBase = plt.axvline(x=k_bar[idx_loc], lw=3.0, color='darkred', ls='--', dashes=(1, 2, 2))            
                    plt.legend([sp, spa, at, aBase],[r'$k/U_b^2$', r'Gaussian', r'Sample Mean', r'Baseline'],prop={'size':15},numpoints=1)                
        #aL = plt.axvline(x=meanField[idx_loc]-2*stdField[idx_loc], lw=3.0, color='black', ls='--')
        #aR = plt.axvline(x=meanField[idx_loc]+2*stdField[idx_loc], lw=3.0, color='black', ls='--')
            if scalarName == 'Deltac1':
                plt.legend([sp, spa, at],[r'Sample', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)   
            elif scalarName == 'Deltac2':
                plt.legend([sp, spa, at],[r'Sample', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1)   
            elif scalarName == 'Deltac3':
                pass
                #plt.legend([sp, spa, at],[r'Sample', r'Gaussian', r'Sample Mean'],prop={'size':15},numpoints=1, loc='upper left')   
        if scalarName == 'c1':
            plt.xlabel(r'$C_1$', fontsize=20)
        elif scalarName == 'c2':
            plt.xlabel(r'$C_2$', fontsize=20)
        elif scalarName == 'c3':
            plt.xlabel(r'$C_3$', fontsize=20)
        elif scalarName == 'Xi':
            plt.xlabel(r'$\Delta\xi$', fontsize=20)
        elif scalarName == 'Eta':
            plt.xlabel(r'$\Delta\eta$', fontsize=20)
        elif scalarName == 'K':
            plt.xlabel(r'$\Delta$ln$k$', fontsize=20)
        elif scalarName == 'VA':
            plt.xlabel(r'$\Delta\varphi_1$ (degree)', fontsize=20)
        elif scalarName == 'VB':
            plt.xlabel(r'$\Delta\varphi_2$ (degree)', fontsize=20)
        elif scalarName == 'VC':
            plt.xlabel(r'$\Delta\varphi_3$ (degree)', fontsize=20)
        elif scalarName == 'TKE':
            plt.xlabel(r'$k/U_b^2$', fontsize=20)
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))        
        elif scalarName == 'Deltac1':
            plt.xlabel(r'$\Delta C_1$', fontsize=20)
        elif scalarName == 'Deltac2':
            plt.xlabel(r'$\Delta C_2$', fontsize=20)
        elif scalarName == 'Deltac3':
            plt.xlabel(r'$\Delta C_3$', fontsize=20)                        
        plt.ylabel(r'PDF', fontsize=20)
        #plt.legend([sp, at, aL],["PDF", "Mean", "2*std"],prop={'size':10},numpoints=1)
        #plt.title(scalarName + ' at (x, y, z) =  '+ str(onePtCoord))
        frame = plt.gca()
        frame.set_ylim([0.0, plt.ylim()[1]])
        if xRange != 'None':
            frame.set_xlim([xRange[0], xRange[1]])
            frame.xaxis.set_ticks(np.arange(xRange[0], xRange[1], 30))
            frame.set_ylim([0.0, plt.ylim()[1]])
            #pdb.set_trace()
        fname = figurefolder + scalarName+'_pdf_'+str(idx_loc)+'_'+ caseName +'.pdf'
        plt.savefig(fname)
        #pdb.set_trace()
        plt.figure(2)
        sc, = plt.plot(bins, cdf, 'red', lw=3.0)
        sca, = plt.plot(bins, cdf_analytical,  'blue', ls='--', lw=3.0, dashes=(9, 2, 2, 2))
        if scalarName == 'c1':
            plt.xlabel(r'$C_1$', fontsize=20)
            plt.legend([sc, sca],[r'$C_1$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'c2':
            plt.xlabel(r'$C_2$', fontsize=20)
            plt.legend([sc, sca],[r'$C_2$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'c3':
            plt.xlabel(r'$C_3$', fontsize=20)
            plt.legend([sc, sca],[r'$C_3$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'Xi':
            plt.xlabel(r'$\Delta\xi$', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta\xi$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'Eta':
            plt.xlabel(r'$\Delta\eta$', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta\eta$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'K':
            plt.xlabel(r'$\Delta$ln$k$', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta$ln$k$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'VA':
            plt.xlabel(r'$\Delta\varphi_1$ (degree)', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta\varphi_1$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'VB':
            plt.xlabel(r'$\Delta\varphi_2$ (degree)', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta\varphi_2$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'VC':
            plt.xlabel(r'$\Delta\varphi_3$ (degree)', fontsize=20)
            plt.legend([sc, sca],[r'$\Delta\varphi_3$', r'Gaussian'],prop={'size':15}, loc=4)
        elif scalarName == 'TKE':
            plt.xlabel(r'$k/U_b^2$', fontsize=20)
            plt.legend([sc, sca],[r'$k/U_b^2$', r'Gaussian'],prop={'size':15}, loc=4)
        plt.ylabel(r'CDF', fontsize=20)        
        #plt.title(scalarName + ' at (x, y, z) =  '+ str(onePtCoord))
        fname = figurefolder + scalarName+'_cdf_'+str(idx_loc)+'_'+ caseName +'.pdf'
        plt.savefig(fname)

    def CToDeltaC(self, c1Fields, c2Fields, c3Fields, baseFolder):
        """
        convert C1, C2 and C3 to DeltaC1, DeltaC2 and DeltaC3
        """
        c1_base = np.loadtxt(baseFolder+'c1_base')
        c2_base = np.loadtxt(baseFolder+'c2_base')
        c3_base = np.loadtxt(baseFolder+'c3_base')
        c1_baseFields = np.tile(c1_base, (self.nSample, 1))
        c2_baseFields = np.tile(c2_base, (self.nSample, 1))
        c3_baseFields = np.tile(c3_base, (self.nSample, 1))
        
        deltac1Fields = c1Fields - c1_baseFields
        deltac2Fields = c2Fields - c2_baseFields
        deltac3Fields = c3Fields - c3_baseFields
        
        return deltac1Fields, deltac2Fields, deltac3Fields

    def log2kTologek(self, deltaLog2k):
        """
        convert log2k to logek
        """
        kRatiokbase = 2**deltaLog2k
        deltaLogek = np.log(kRatiokbase)
        return deltaLogek
        
    def plotTauAllLoc1SampleOnBay(self, idx_sample):
        """
        plot ensemble of one location of Tau fields on Baycentric triangle
        (Normally not used)
        :return:
        None
        """
        plt.clf()
        pylab.rc('font', family='Times New Roman')
        plt.plot(self.C1Fields[idx_sample, :], self.C2Fields[idx_sample, :], 'o', color='#1b9e76',alpha=transparency)
        #plt.scatter(self.C1Field_base, self.C2Field_base, s=50, c='red')
        self._plotBaryTriangle()
        plt.show()

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

    def _fillPDFInside(self, x, density, facecolor, color='black', transparency=0.3, hatch='//'):
        """
        fill shaded region in PDF
        """

        plt.fill_between(x, density, facecolor=facecolor, lw=0, color=color, alpha=transparency, hatch=hatch)
        plt.fill_between(x, density, facecolor='None', lw=0, color=color, hatch=hatch)

# Main function
if __name__ == '__main__':
    mainInputFile = './mainInput.in'
    plotInputFile = './plotInfo.in'
    paramDict = readInputData(plotInputFile)
    caseName = paramDict['caseName']
    ## parse plot control for Bay plot
    plotAllLocs = ast.literal_eval(paramDict['plotAllLocs'])
    plotAllTauComponents = ast.literal_eval(paramDict['plotAllTauComponents'])
    plotVlineBay = ast.literal_eval(paramDict['plotVlineBay'])
    plotHlineBay = ast.literal_eval(paramDict['plotHlineBay'])
    ns_VH = int(paramDict['ns_VH'])
    onePtCoord = extractListFromDict(paramDict, 'onePtCoord')
    onePtCoord = np.array([[float(pn) for pn in onePtCoord]])
    pltPost = postPlot(mainInputFile, onePtCoord)
    idx_sample = int(paramDict['idx_sample'])
    compareGaussianFlag = paramDict['compareGaussian']
    legendFlag = paramDict['legendFlag']
    shadedFlag = paramDict['shadedFlag']
    # Scattering plot of Tau at one location on Barycentric triangle
    pltPost.plotTau1LocSamplesOnBay(caseName, legendFlag)
    # Contour plot of Tau at one location on Barycentric triangle
    pltPost.plotTau1LocContourOnBay(caseName)
    
    if plotVlineBay:
        pltPost.plotVlineonBay(ns_VH, caseName)

    if plotHlineBay:
        pltPost.plotHlineonBay(ns_VH, caseName)

    if plotAllTauComponents:
        deltaFolder = './resultData/deltaRComponent_samples/'
        fileNameCommon = 'delta'
        deltaRComponentPool = ['Xi', 'Eta', 'K', 'VA', 'VB', 'VC']
        for dR_c in deltaRComponentPool:
            print 'processing the deltaR component: ' + dR_c
            #pdb.set_trace()
            scalarFields_temp = np.loadtxt(deltaFolder+fileNameCommon+dR_c+'_s')
            if dR_c == 'K':
               scalarFields_temp = pltPost.log2kTologek(scalarFields_temp)         
            if (dR_c == 'VA' or dR_c == 'VB' or dR_c == 'VC'):
                xRange = [-90, 91]
                scalarFields_temp = 180.0*scalarFields_temp/np.pi
                pltPost.plotScalarStat_1Loc(scalarFields_temp, dR_c, caseName, xRange, compareGaussianFlag, shadedFlag)
            else:
                pltPost.plotScalarStat_1Loc(scalarFields_temp, dR_c, caseName, 'None', compareGaussianFlag, shadedFlag)

        RComFolder = './resultData/RComponent_samples/'
        c1Fields = np.loadtxt(RComFolder+'c1_s')
        c2Fields = np.loadtxt(RComFolder+'c2_s')
        c3Fields = np.loadtxt(RComFolder+'c3_s')
        deltac1Fields, deltac2Fields, deltac3Fields = pltPost.CToDeltaC(c1Fields, c2Fields, c3Fields, RComFolder)
        print 'processing the deltaR component: c1'
        pltPost.plotScalarStat_1Loc(c1Fields, 'c1', caseName, 'None', compareGaussianFlag, shadedFlag)
        print 'processing the deltaR component: c2'
        pltPost.plotScalarStat_1Loc(c2Fields, 'c2', caseName, 'None', compareGaussianFlag, shadedFlag)
        print 'processing the deltaR component: c3'
        pltPost.plotScalarStat_1Loc(c3Fields, 'c3', caseName, 'None', compareGaussianFlag, shadedFlag)
        print 'processing the deltaR component: Deltac1'
        pltPost.plotScalarStat_1Loc(deltac1Fields, 'Deltac1', caseName, 'None', compareGaussianFlag, shadedFlag)
        print 'processing the deltaR component: Deltac2'
        pltPost.plotScalarStat_1Loc(deltac2Fields, 'Deltac2', caseName, 'None', compareGaussianFlag, shadedFlag)
        print 'processing the deltaR component: Deltac3'
        pltPost.plotScalarStat_1Loc(deltac3Fields, 'Deltac3', caseName, 'None', compareGaussianFlag, shadedFlag)                
        
        
        
        TKEFields = np.loadtxt(RComFolder+'k_s')
        print 'processing TKE component'
        pltPost.plotScalarStat_1Loc(TKEFields, 'TKE', caseName, 'None', compareGaussianFlag, shadedFlag)

    if plotAllLocs:
        pltPost.plotTauAllLoc1SampleOnBay(idx_sample)




