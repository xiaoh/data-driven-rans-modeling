#!/usr/bin/env python
# description        :Generate Tau discrepancy
# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Oct.13, 2016
# revision           :Oct.19, 2016
########################################################################################################################

# SYSTEM MODULE
import numpy as np
import matplotlib as mp
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
import shutil
import subprocess
import pdb
# LOCAL MODULE
import utilities as uT
import ReynoldsStressRF as ReRF
import foamFileOperation as foamOp

deltaTau_headerNames_Euler = ['Index', 'ccx', 'ccy', 'ccz', 'deltaXb', 'deltaYb', 'deltaLog2K', 'deltaVA', 'deltaVB', 
                              'deltaVC']
tauPhysical_Euler = ['Xb', 'Yb', 'K', 'VA', 'VB', 'VC']
                              
deltaTau_headerNames_Quaternion = ['Index', 'ccx', 'ccy', 'ccz', 'deltaXb', 'deltaYb', 'deltaLog2K', 'theta', 'vx', 
                                   'vy', 'vz', 'indexList']

class getTrainingDataMain():
    """
    main function to obtain training data
    Arg:
        caseType:                   flow type
        Re:                         Reynold number
        nMesh:                      number of mesh cells
        deltaOrienParamerization:   method of parameterizing deltaV
                                    
                                    (Quaternion, Euler, rotationMatrix) 
    """
    def __init__(self, caseType, Re, nMesh, dataSparseness, dataSource, currentRANSDir, turbulenceModel, deltaOrienParamerization='Quaternion'): 
        """
        initialization
        
        """
        self.timeDir = '0.000000'    
        self.caseType = caseType
        self.Re = Re
        self.nMesh = nMesh
        self.deltaOrienParamerization = deltaOrienParamerization
        self.dataSparseness = dataSparseness
        self.dataSource = dataSource
        srcDir = os.path.normpath(os.path.join(os.path.abspath(__file__), '../../..'))
        self.baseDir = os.path.join(srcDir, 'database')
        self.currentRANSDir = currentRANSDir
        self.baseRANSDir = os.path.join(srcDir, 'RANSTemp', caseType, 'cell'+nMesh, 'cases_'+turbulenceModel)

    def getDiscrepancy(self, reflectCoord=False):
        """
        To get the discrepancy between DNS Tau and RANS Tau
        """
        caseDir = os.path.join(self.baseDir, self.caseType, self.dataSparseness)
        tauCurrentDir = os.path.join(self.currentRANSDir, self.caseType+'Re'+self.Re, '0.000000', 'Tau')
        tauRANSDir = os.path.join(self.currentRANSDir, self.caseType+'Re'+self.Re, '0.000000', 'TauRANS')

        if self.dataSparseness == 'fullField':
            tauDataDir = os.path.join(caseDir, 'Tau_'+self.caseType+'-Re'+self.Re+'-'+self.dataSource+'-cell'+self.nMesh)
            
            Tau_base = foamOp.readTurbStressFromFile(tauRANSDir)
            Tau_dns = foamOp.readTurbStressFromFile(tauDataDir)               
            
            if reflectCoord:
                Tau_base[:, 1] = -Tau_base[:, 1]
                Tau_base[:, 4] = -Tau_base[:, 4]
                Tau_dns[:, 1] = -Tau_dns[:, 1]
                Tau_dns[:, 4] = -Tau_dns[:, 4]
            
            nCell = int(self.nMesh)
            caseZeroDir = os.path.join(self.currentRANSDir, self.caseType+'Re'+self.Re, '0.000000')

            # obtain mesh information
            cellIdx = np.array([np.arange(nCell)]).T
            coord = foamOp.readTurbCoordinateFromFile(caseZeroDir)
            # generate delta tau
            mapTau = ReRF.ReynoldsStressRF(caseZeroDir, 'TauRANS', nCell, 1, correctInitTau='False')
            deltaXb, deltaYb, deltaLog2K = mapTau.getDeltaXYbaryLog2K(Tau_base, Tau_dns)
            # if you want deltaXi, deltaEta            
            # deltaXi, deltaEta, deltaLog2K = mapTau.getDeltaXiEtaLog2K(Tau_base, Tau_dns) 
            
            if self.deltaOrienParamerization == 'Euler': 
                deltaVA, deltaVB, deltaVC = mapTau.getDeltaThetaVABC(Tau_base, Tau_dns)
                deltaV = np.hstack([deltaVA, deltaVB, deltaVC])            
                # for verify: test if tau_DNS_reconstruct == Tau_dns?
                # tauDNS_reconstruct = mapTau.perturbTauInBary(deltaXb, deltaYb, deltaLog2K, deltaV)
                deltaTau = np.hstack([cellIdx, coord, deltaXb, deltaYb, deltaLog2K, deltaV])
                headerNames = deltaTau_headerNames_Euler
                    
            elif self.deltaOrienParamerization == 'Quaternion':
                deltaVecList, indexList = mapTau.deltaVec(Tau_base, Tau_dns)
                theta, vx, vy, vz = mapTau.getQuaternion(Tau_base, Tau_dns, indexList)
                # for verify: test if tau_DNS_reconstruct == Tau_dns?
                # tauDNS_reconstruct = mapTau.reconstructTauFromQuatInBary(Tau_base,deltaXb,deltaYb,deltaLog2K,vx,vy,vz,theta,indexList)
                deltaTau = np.hstack([cellIdx, coord, deltaXb, deltaYb, deltaLog2K, theta, vx, vy, vz, indexList])
                headerNames = deltaTau_headerNames_Quaternion
                                                      
            # polish header
            headStr = ''
            for headerName in headerNames:
                headStr = headStr + headerName + '      '
            headStr = headStr + '\n'    
            
            currentTrainDataDir = './trainData'
            if not os.path.exists(currentTrainDataDir):
                os.makedirs(currentTrainDataDir) 
            
            # Save learning
            deltaTauFileName = 'deltaTau_' + self.caseType + self.Re + '_' + self.deltaOrienParamerization
            dataPath = os.path.join(currentTrainDataDir, deltaTauFileName)
            np.savetxt(dataPath, deltaTau, header = headStr) 
            
            # save as openFOAM files
            self._writeTruthFieldToOpenFOAM(deltaTau, Tau_dns, headerNames)
            
        elif self.dataSparseness == 'coarseData':
            #TODO: How to unify the coarse data with full field data needs more thinking
            pass
        
        return headerNames
            

    def _writeTruthFieldToOpenFOAM(self, deltaTau, Tau_DNS, headerNames):
        """
        write true delta tau components and true tau to OpenFOAM field files and sample them for further visualization
        """
        print "writting all training data as OpenFOAM field files in 0.000000 folder\n" 

        currentRANSCaseDir = os.path.join(self.currentRANSDir, self.caseType+'Re'+self.Re)
        baseRANSCaseDir = os.path.join(self.baseRANSDir, 'Re'+self.Re)
        
        i = 4                            
        for component in headerNames[4:]:
            scalarField = deltaTau[:, i] # last column is the field
            i = i + 1            
            for ReCase in [currentRANSCaseDir, baseRANSCaseDir]:            
                componentPath = os.path.join(ReCase, self.timeDir, component+'_truth')
                templateFilePath = os.path.join(ReCase,self.timeDir, 'scalarFieldTemplate')
                if os.path.exists(componentPath):
                    os.remove(componentPath)
                shutil.copy(templateFilePath, componentPath)
                uT.replace(componentPath, '<Feature>', component)
                foamOp.writeScalarToFile(scalarField, componentPath)
                                        
                # sample the truth
                shutil.copy(os.path.join(ReCase,'system/sampleDict.template'),os.path.join(ReCase,'system/sampleDict'))
                uT.replace(os.path.join(ReCase,'system/sampleDict'), '<field>', component+'_truth')
                subprocess.check_call('sample >> log.sampleTruth', shell=True, cwd=ReCase)
            

         
         
         
         

            
def getDiscrepancy(caseType, Re):
    """
    To get the discrepancy between DNS Tau and RANS Tau
    """
    caseDir = baseFolder + caseType + '/Re' + Re
    Tau_base = foamOp.readTurbStressFromFile(caseDir+'/0/Tau')
    Tau_dns = foamOp.readTurbStressFromFile(caseDir+'/0/TauDNS')
    
    if reflectCoord:
        Tau_base[:, 1] = -Tau_base[:, 1]
        Tau_base[:, 4] = -Tau_base[:, 4]
        Tau_dns[:, 1] = -Tau_dns[:, 1]
        Tau_dns[:, 4] = -Tau_dns[:, 4]        
    
    nCell = Tau_base.shape[0]
    mapTau = ReRF.ReynoldsStressRF(caseDir + '/0/', 'Tau', nCell, 1, correctInitTau='False')
    
    deltaXi, deltaEta, deltaLog2K = mapTau.getDeltaXiEtaLog2K(Tau_base, Tau_dns)
    deltaVA, deltaVB, deltaVC = mapTau.getDeltaThetaVABC(Tau_base, Tau_dns)
    deltaV = np.hstack([deltaVA, deltaVB, deltaVC])
    
    cellIdx = np.array([np.arange(nCell)]).T
    coord = foamOp.readTurbCoordinateFromFile(caseDir+'/0/')
    
    deltaTau = np.hstack([cellIdx, coord, deltaXi, deltaEta, deltaLog2K, deltaV])
    
    np.savetxt('./processedData/deltaTau_Re'+Re+'_'+caseType, deltaTau, header='Index, x, y, z, deltaXi, deltaEta, deltaK, deltaVA, deltaVB, deltaVC')

def getDiscrepancy_modifyEuler(caseType, Re):
    """
    To get the discrepancy between DNS Tau and RANS Tau
    (The Euler angle is combined as deltaPhi_1 = (deltaVA + deltaVB), deltaPhi_2 = 
    (deltaVA - deltaVB), and deltaPhi_3 = deltaVC  
    """
    caseDir = baseFolder + caseType + '/Re' + Re
    Tau_base = foamOp.readTurbStressFromFile(caseDir+'/0/Tau')
    Tau_dns = foamOp.readTurbStressFromFile(caseDir+'/0/TauDNS')

    if reflectCoord:
        Tau_base[:, 1] = -Tau_base[:, 1]
        Tau_base[:, 4] = -Tau_base[:, 4]
        Tau_dns[:, 1] = -Tau_dns[:, 1]
        Tau_dns[:, 4] = -Tau_dns[:, 4]    
    
    nCell = Tau_base.shape[0]
    mapTau = ReRF.ReynoldsStressRF(caseDir + '/0/', 'Tau', nCell, 1, correctInitTau='False')
    
    deltaXi, deltaEta, deltaLog2K = mapTau.getDeltaXiEtaLog2K(Tau_base, Tau_dns)
    deltaVA, deltaVB, deltaVC = mapTau.getDeltaThetaVABC(Tau_base, Tau_dns)
    deltaPhi1 = deltaVA + deltaVB
    deltaPhi2 = deltaVA - deltaVB
    deltaPhi3 = deltaVC
    
    deltaPhi = np.hstack([deltaPhi1, deltaPhi2, deltaPhi3])
    
    cellIdx = np.array([np.arange(nCell)]).T
    coord = foamOp.readTurbCoordinateFromFile(caseDir+'/0/')
    
    deltaTau = np.hstack([cellIdx, coord, deltaXi, deltaEta, deltaLog2K, deltaPhi])
    
    np.savetxt('./processedData/deltaTau_Re'+Re+'_'+caseType+'_modifyEuler', deltaTau, 
                header='Index, x, y, z, deltaXi, deltaEta, deltaK, deltaPhi1, deltaPhi2, deltaPhi3')


def getDiscrepancy_rotation(caseType, Re):
    """
    To get the discrepancy (xi, eta, k, E) between DNS Tau and RANS Tau
    """
    caseDir = baseFolder + caseType + '/Re' + Re
    Tau_base = foamOp.readTurbStressFromFile(caseDir+'/0/Tau')
    Tau_dns = foamOp.readTurbStressFromFile(caseDir+'/0/TauDNS')
    nCell = Tau_base.shape[0]
    mapTau = ReRF.ReynoldsStressRF(caseDir + '/0/', 'Tau', nCell, 1, correctInitTau='False')
    deltaXi, deltaEta, deltaLog2K = mapTau.getDeltaXiEtaLog2K(Tau_base, Tau_dns)
    deltaVecList, indexList = mapTau.deltaVec(Tau_base, Tau_dns)
    E = mapTau.getDeltaCosine(Tau_base, Tau_dns, indexList)
    
    deltaE1F = np.array([E[:, 0]]).T
    deltaE2F = np.array([E[:, 1]]).T
    deltaE3F = np.array([E[:, 2]]).T
    deltaE4F = np.array([E[:, 3]]).T
    deltaE5F = np.array([E[:, 4]]).T
    deltaE6F = np.array([E[:, 5]]).T
    deltaE7F = np.array([E[:, 6]]).T
    deltaE8F = np.array([E[:, 7]]).T
    deltaE9F = np.array([E[:, 8]]).T
    
    foamOp.writeScalarToFile(deltaE1F, caseDir+'/0/deltaE1')
    foamOp.writeScalarToFile(deltaE2F, caseDir+'/0/deltaE2')
    foamOp.writeScalarToFile(deltaE3F, caseDir+'/0/deltaE3')
    foamOp.writeScalarToFile(deltaE4F, caseDir+'/0/deltaE4')
    foamOp.writeScalarToFile(deltaE5F, caseDir+'/0/deltaE5')
    foamOp.writeScalarToFile(deltaE6F, caseDir+'/0/deltaE6')
    foamOp.writeScalarToFile(deltaE7F, caseDir+'/0/deltaE7')
    foamOp.writeScalarToFile(deltaE7F, caseDir+'/0/deltaE8')
    foamOp.writeScalarToFile(deltaE9F, caseDir+'/0/deltaE9')
    cellIdx = np.array([np.arange(nCell)]).T
    coord = foamOp.readTurbCoordinateFromFile(caseDir+'/0/')

    deltaTau = np.hstack([cellIdx, coord, deltaXi, deltaEta, deltaLog2K, E])
    np.savetxt('./processedData/deltaTau_Re'+Re+'_'+caseType+'_rotation', deltaTau, header='Index, x, y, z, deltaXi, deltaEta, deltaLogK, deltaE')
    
if __name__ == "__main__":    
    baseFolder = '../../markerGen/templateCases/RANS/'
    caseType = 'duct'
    Re = '2200'
    modifyEulerAngle = True
    reflectCoord = True
    ####################################################################################################################    
    print "We are generating DNS discrepancy data for training: case = ", caseType, 'with Re = ', Re
    if modifyEulerAngle:
        print "The Euler angle is modified as phi_1 = alpha + gamma, and phi_2 = alpha - gamma, and phi_3 = beta"
        getDiscrepancy_modifyEuler(caseType, Re)
    else:
        print "The Euler angle is raw, no combination is imposed"
        getDiscrepancy(caseType, Re)
    
    #getDiscrepancy_rotation(caseType, Re)
    

