#! /usr/bin/env python
import sys, os, os.path
import matplotlib as mp
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pdb
import shutil
from scipy.interpolate import Rbf, griddata, CloughTocher2DInterpolator
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from pylab import *
from utilities import readInputData



if not os.path.exists('figures/'):
    os.mkdir('figures/')

paramDict = readInputData('plotInfo.in')
timeDir = paramDict['timeDir']

paramDict_main = readInputData('MainInput.in')
Ns = int(paramDict_main['Ns'])
samplePlotList = range(1,Ns+1)

xLocVec = ['0p05', '0p10', '0p15', '0p20', '0p25', '0p30', '0p35', '0p40', '0p45', '0p50', '0p55', '0p60', '0p65', '0p70', '0p75', '0p80', '0p90', '0p95', '1p00', '1p25', '1p5', '1p75', '2p00']





starti = 0
c_T = 4.2538
smallValue = 1e-5
Uref = 27
nSampleLes = 100
dataDir = './LESdata'



for xLoc in xLocVec:
    if xLoc in ['0p80', '0p90', '0p95', '1p00', '1p25', '1p5', '1p75', '2p00']:
        yLocVec = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
    
    elif xLoc in ['0p05', '0p45', '0p50', '0p55', '0p60', '0p65', '0p70', '0p75']:
        yLocVec = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
    
    elif xLoc in ['0p10', '0p15', '0p20', '0p25', '0p30', '0p35', '0p40']:
        yLocVec = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

    if xLoc in ['0p05', '0p10', '0p15']:
        yRange = [0.10, 0.48, 0, 0.3]
        zRange = [0.05, 0.48, 0, 0.3]
        scaleVy = 0.1
        scaleVz = 1.0
    elif xLoc in ['0p25', '0p30']:
        yRange = [0.10, 0.48, 0, 0.3]
        zRange = [0.05, 0.48, 0, 0.3]
        scaleVy = 0.2
        scaleVz = 1.0
    elif xLoc in ['0p35', '0p40']:
        yRange = [0.10, 0.48, 0, 0.3]
        zRange = [0.05, 0.48, 0, 0.3]
        scaleVy = 0.3
        scaleVz = 1.0
    elif xLoc in ['0p45', '0p50']:
        yRange = [0.05, 0.48, 0, 0.3]
        zRange = [0.05, 0.48, 0, 0.3]
        scaleVy = 0.3
        scaleVz = 1.0
    elif xLoc in ['0p55', '0p60', '0p65', '0p70', '0p75', '0p80', '0p90', '0p95', '1p00', '1p25', '1p5', '1p75', '2p00']:
        yRange = [-0.05, 0.48, 0, 0.3]
        zRange = [0.0, 0.48, 0, 0.3]
        scaleVy = 0.5
        scaleVz = 2.0                        
    
    dataFileName = 'LES_VREMAN_Fine_XC' + xLoc + '_cellCenter'
    dataFileDir = os.path.join(dataDir, dataFileName)
    dataM = np.loadtxt(dataFileDir, delimiter=',')
    nVar = dataM.shape[1]

    RANSDir = 'wingBody_base-observation/postProcessing/sets/' + timeDir + '.000000'

    print "plot x = ", xLoc

    if nVar == 13:
        yIdx = 1;
        zIdx = 2;
        uyIdx = 5;
        uzIdx = 6;
    elif nVar == 12:
        yIdx = 0;
        zIdx = 1;
        uyIdx = 4;
        uzIdx = 5;

    yMatrix = dataM[:, yIdx]
    zMatrix = dataM[:, zIdx]
    uyMatrix = dataM[:, uyIdx]
    uzMatrix = dataM[:, uzIdx]

    yField = dataM[:, yIdx]

    plt.figure()
    ax1=plt.subplot(111)
    i = starti
    for yLoc in yLocVec:
        xSum = 0
        xMean = 0                
        for case in os.listdir('.'):
            if fnmatch.fnmatch(case, 'wingBody_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(case, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = case + '/postProcessing/sets/' + timeDir + '.000000/' + \
                                'Plane' + xLoc + '_' + str(i) + '_U.xy'
                    M = np.loadtxt(filename, comments = '%')        
                    zField_sample = M[:, 2]
                    uyField_sample = M[:, 4]
                    uzField_sample = M[:, 5]
                    
                    x = scaleVy * uyField_sample + yLocVec[i]
                    y = zField_sample                    
                    xSum = xSum + x                    
                    pS, = plt.plot(x, y, color='#1b9e76', alpha=1, lw = 2, mfc = 'none')
        
        xMean = xSum / len(samplePlotList)
        pMean, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))                    
        
        targetXLoc = yField[(np.abs(yField - yLoc)).argmin()]    
        idx_selec_yLoc = np.where((yField>=targetXLoc-smallValue) & (yField<=targetXLoc+smallValue))
        dataM_selec = dataM[idx_selec_yLoc, :]
        dataM_selec = dataM_selec[0]
        zField_selec_unsort = dataM_selec[:, zIdx]
        idxSortz = np.argsort(zField_selec_unsort)
        dataM_selec_sorted = dataM_selec[idxSortz, :]
        
        yField_LES = dataM_selec_sorted[:, yIdx]
        zField_LES = dataM_selec_sorted[:, zIdx]
        uyField_LES = dataM_selec_sorted[:, uyIdx]
        uzField_LES = dataM_selec_sorted[:, uzIdx]
        
        RANSFileName = 'Plane' + xLoc + '_' + str(i) + '_U.xy'
        RANSFileDir = os.path.join(RANSDir, RANSFileName)
        URANS = np.loadtxt(RANSFileDir)
        zField_RANS = URANS[:, 2]
        uyField_RANS = URANS[:, 4]
        uzField_RANS = URANS[:, 5]

        x = scaleVy * uyField_LES / Uref + yLocVec[i]
        y = zField_LES

        pT, = plt.plot(x, y, color='black', alpha=1, lw = 2, mfc = 'none')

        x = scaleVy * uyField_RANS + yLocVec[i]
        y = zField_RANS
        pB, = plt.plot(x, y, color='red', alpha=1, lw = 2, mfc = 'none')
        
        
        
        xzero = np.zeros(x.shape) + yLocVec[i]
        p_zero, = plt.plot(xzero, y, '--', color='black', alpha=1.0, lw = 1.0, mfc = 'none')    

        i = i + 1

    plt.title('Uy, x/C = ' + xLoc)
    plt.ylabel("$z/C$")
    plt.xlabel(r'$y/C;\quad$ $'+ str(scaleVy) + 'U_y/U_{ref}$ $+y/C$')
    plt.legend([pB, pT, pS, pMean], ["Baseline", "LES", "samples", "Mean"], prop={'size':10}, numpoints=1,
               bbox_to_anchor=(0.3, 1.01), loc=3, ncol=3)            

    plt.axis(yRange)

    plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    mp.rcParams.update({'font.size':12})

    figureName = './figures/Uy_xLoc_' + xLoc + '_time' + timeDir + '.pdf'
    plt.savefig(figureName, bbox_inches='tight')

    plt.figure()
    ax1=plt.subplot(111)
    i = starti
    for yLoc in yLocVec:
        xSum = 0
        xMean = 0                
        for case in os.listdir('.'):
            if fnmatch.fnmatch(case, 'wingBody_base-tmp*'):
                plotFlag = False
                for sample in samplePlotList:
                    if fnmatch.fnmatch(case, '*_' + str(sample) + '.0*'):
                        plotFlag = True
                if plotFlag ==True:
                    filename = case + '/postProcessing/sets/' + timeDir + '.000000/' + \
                                'Plane' + xLoc + '_' + str(i) + '_U.xy'
                    M = np.loadtxt(filename, comments = '%')        
                    zField_sample = M[:, 2]
                    uyField_sample = M[:, 4]
                    uzField_sample = M[:, 5]
                    
                    x = scaleVz * uzField_sample + yLocVec[i]
                    y = zField_sample                    
                    xSum = xSum + x                    
                    pS, = plt.plot(x, y, color='#1b9e76', alpha=1, lw = 2, mfc = 'none')
        
        xMean = xSum / len(samplePlotList)
        pMean, = plt.plot(xMean, y, '-',color ='blue',lw=2,dashes=(9,2,2,2))                    
        
        targetXLoc = yField[(np.abs(yField - yLoc)).argmin()]    
        idx_selec_yLoc = np.where((yField>=targetXLoc-smallValue) & (yField<=targetXLoc+smallValue))
        dataM_selec = dataM[idx_selec_yLoc, :]
        dataM_selec = dataM_selec[0]
        zField_selec_unsort = dataM_selec[:, zIdx]
        idxSortz = np.argsort(zField_selec_unsort)
        dataM_selec_sorted = dataM_selec[idxSortz, :]
        
        yField_LES = dataM_selec_sorted[:, yIdx]
        zField_LES = dataM_selec_sorted[:, zIdx]
        uyField_LES = dataM_selec_sorted[:, uyIdx]
        uzField_LES = dataM_selec_sorted[:, uzIdx]
        
        RANSFileName = 'Plane' + xLoc + '_' + str(i) + '_U.xy'
        RANSFileDir = os.path.join(RANSDir, RANSFileName)
        URANS = np.loadtxt(RANSFileDir)
        zField_RANS = URANS[:, 2]
        uyField_RANS = URANS[:, 4]
        uzField_RANS = URANS[:, 5]

        x = scaleVz * uzField_LES / Uref + yLocVec[i]
        y = zField_LES

        pT, = plt.plot(x, y, color='black', alpha=1, lw = 2, mfc = 'none')

        x = scaleVz * uzField_RANS + yLocVec[i]
        y = zField_RANS
        pB, = plt.plot(x, y, color='red', alpha=1, lw = 2, mfc = 'none')
        
        
        
        xzero = np.zeros(x.shape) + yLocVec[i]
        p_zero, = plt.plot(xzero, y, '--', color='black', alpha=1.0, lw = 1.0, mfc = 'none')    

        i = i + 1

    plt.title('Uz, x/C = ' + xLoc)
    plt.ylabel("$z/C$")
    plt.xlabel(r'$y/C;\quad$ $'+ str(scaleVz) + 'U_z/U_{ref}$ $+y/C$')
    plt.legend([pB, pT, pS, pMean], ["Baseline", "LES", "samples", "Mean"], prop={'size':10}, numpoints=1,
               bbox_to_anchor=(0.3, 1.01), loc=3, ncol=3)            

    plt.axis(zRange)
    #plt.axis([0.05, 0.48, 0, 0.3])

    plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    mp.rcParams.update({'font.size':12})

    figureName = './figures/Uz_xLoc_' + xLoc + '_time' + timeDir + '.pdf'
    plt.savefig(figureName, bbox_inches='tight')    
