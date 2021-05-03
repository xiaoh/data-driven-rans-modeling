#!/usr/bin/env python

# description        :visualize regression results

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Mar.27, 2016
# revision           :Oct.21, 2016
########################################################################################################################

# SYSTEM MODULE
import numpy as np
from sklearn.cluster import KMeans
import matplotlib as mp
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pdb
# LOCAL MODULE
import postFuns as pF

def getMeshIdxOf_xVec(caseMesh, xVec, y0, z0):
    """
    To obtain mesh idx of vertial profiles at xVec = [x0, x1, ..., xn]
    """
    meshIdx_dict = {} # define a dictionary for store mesh idx given x loc
    for x0 in xVec:   
        hCoord, vCoord = pF.getVerandHorCoord(x0, y0, z0, caseMesh)
        idx_v = pF.coord2cellIdx(vCoord, caseMesh)    
        idx_locs = idx_v[:, 1].tolist()        
        meshIdx_dict[str(x0)] = idx_locs
    
    return meshIdx_dict        
    
def extractPredictionFullField(caseMesh, deltaTau, xVec, y0, z0):
    """
    To extract profiles on give xVec = [x0, x1, ..., xn]
    
    Arg:
        deltaTau [nCell by nDisc], order = [deltaXi, deltaEta, deltaLog2K, deltaVA, deltaVB, deltaVC]
    """
    
    profiles_dict = {} # define a dictionary for store profiles at given x loc
    
    for x0 in xVec:   
        hCoord, vCoord = pF.getVerandHorCoord(x0, y0, z0, caseMesh)
        idx_v = pF.coord2cellIdx(vCoord, caseMesh)    
        idx_locs = idx_v[:, 1].tolist()        
        deltaTau_i = deltaTau[idx_locs, :]        
        profiles_dict[str(x0)] = deltaTau_i
    
    return profiles_dict

def getProfiles_MeshIdx(caseMesh, deltaTau, xVec, y0, z0):
    """
    To extract profiles on give xVec = [x0, x1, ..., xn]
    
    Arg:
        deltaTau [nCell by nDisc], order = [deltaXi, deltaEta, deltaLog2K, deltaVA, deltaVB, deltaVC]
    """
    
    meshIdx_dict = {} # define a dictionary for store profiles at given x loc
    
    for x0 in xVec:   
        hCoord, vCoord = pF.getVerandHorCoord(x0, y0, z0, caseMesh)
        idx_v = pF.coord2cellIdx(vCoord, caseMesh)    
        idx_locs = idx_v[:, 1].tolist()              
        meshIdx_dict[str(x0)] = idx_locs
    
    return meshIdx_dict


    
if __name__ == "__main__":    
    caseName = 'pehill'
    caseMesh = '../flowClustering/templateMesh/pehill/pehill10595'
    nFeatures = 10
    nClusters = 1
    Re_prediction = '10595'
    deltaTauCom = ['Xi', 'Eta', 'K', 'VA', 'VB', 'VC']
    xVec = np.arange(10) + 1    
    xVec = xVec.astype(str)
    xVec_openFoam = [0.0, 0.5, 1., 2., 3., 4., 5., 6., 7., 8.]
    xVec_arbituray = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
    
    scaleV = 0.6
    extractPredictionTo8lines()    
    plotDiscrepancyCompare(scaleV)
    
    #extractPredictionFullField()
    #plotDiscrepancyCompare_full(scaleV)

    

