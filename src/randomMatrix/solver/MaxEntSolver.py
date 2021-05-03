#!/usr/bin/env python

# description        :Main drive of model-form uncertainty propagation by using Random matrix
#                     and maximum entropy theory

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Dec.03, 2015
# revision           :Dec.05, 2015

####################################################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.sparse as sp
import scipy.stats as ss
from scipy.stats import norm
import numpy.linalg as LA
# system, file operation
import pdb
import os
import os.path as ospt
import ast
import time
# plotting
import matplotlib.pyplot as plt  # for plotting
# Import local modules
from interfaceRANS import interfaceRANS
from RMFields import RMF
from ReynoldsStressRF import ReynoldsStressRF
from utilities import readInputData, extractListFromDict
import foamFileOperation as foamOp
try:
    import pp
    hasPP = True
    print "Parallel Python imported successfully"
except ImportError, e:
    print e, 'Warning: no pp module, I only can use one core to do the calculation'
    hasPP = False
    pass

class MaxEntSolver:
    """
    To generated a group of Reynolds stress samples based on RMF,
    and propagate this Reynolds stress ensemble to velocity.

    Arg:
        Nsample         : Number of required fields

    """

    def __init__(self, mainInputFile):
        ## parse parameters in dictionary
        paramDict = readInputData(mainInputFile)
        self.caseName = paramDict['caseName']
        self.caseNameRun = self.caseName + '-Run'
        self.caseSolver = paramDict['caseSolver']
        self.nSample = int(paramDict['nSample'])
        self.nSample_propagate = int(paramDict['nSample_propagate'])
        self.nCell_cfd = int(paramDict['nCell_cfd'])
        self.nCell_kl = int(paramDict['nCell_kl'])
        self.nCell_pce = int(paramDict['nCell_pce'])
        self.writeInterval = float(paramDict['writeInterval'])
        self.propagateFlag = ast.literal_eval(paramDict['propagateFlag'])
        self.kernelType = paramDict['kernelType']
        nKL = int(paramDict['nKL'])
        self.truncateTol = float(paramDict['truncateTol'])
        lenXYZ = extractListFromDict(paramDict, 'lenXYZ')
        self.lenXYZ = np.array([[float(pn) for pn in lenXYZ]])
        self.pOrder = int(paramDict['pOrder'])
        self.resultDir = './resultData/'
        ## specify path to baseline case, with 2 other sets of mesh (mesh for kl and pce)
        self.cfdDir = ospt.join(os.getcwd(), self.caseName)
        self.klDir = self.cfdDir + '-kl'
        self.pceDir = self.cfdDir + '-pce'
        ## construct coord
        xCoord = foamOp.readTurbCoordinateFromFile(self.klDir+'/0/')
        ## generate arguments dictionary for covariance
        Arg_covGen = self._prepareArg_covGen()
        ## generate arguments dictionary for KL expansion
        Arg_calKLModes = self._prepareArg_calKLModes(nKL, Arg_covGen['weightField'])
        ## initialize class of interfaceRANS
        self.intf = interfaceRANS(self.nSample, self.nCell_cfd)
        hostCellIdx_pce = np.loadtxt(self.pceDir+'/constant/indexHost.txt')
        #pdb.set_trace()
        stationaryFlag = ast.literal_eval(paramDict['stationaryFlag'])
        deltaStationaryValue = float(paramDict['deltaStationaryValue'])
        deltaFieldDir = paramDict['deltaFieldDir']
        deltaField = self._prepareDeltaField(stationaryFlag, deltaStationaryValue, deltaFieldDir)
        ## initialize class of random matrix
        self.randomMatrix = RMF(self.nCell_kl, self.nCell_pce, Arg_covGen, Arg_calKLModes, deltaField, xCoord, hostCellIdx_pce)
        self.d = 3
        ## open space for ndarray
        # ensemble of random matrix field R(x), and corresponding ensemble of Tau field
        self.RFields = np.zeros([self.nSample, self.nCell_cfd, self.d, self.d])
        self.TauSample = np.zeros([self.nSample, self.nCell_cfd, 6])
        self.deltaTauSample = np.zeros([self.nSample, self.nCell_cfd, 6])
        self.GSample = np.zeros([self.nSample, self.nCell_kl, 6])
        # mean of random matrix field Rbar, and corresponding Taubar field
        # read openFOAM Tau file to get Rbar (as mean Reynolds stress field)
        #tauFile = '/0/Tau'
        tauFile = '0/Tau'      ## JL: Fixed the path issue
        self.Rbar, self.Taubar = self.intf.getRbar(self.cfdDir, tauFile)

        # enanble PP only if user specified in mainInput file
        # the PP model has been imported successfully
        self.enablePP = ast.literal_eval(paramDict['enablePP']) and hasPP
        # initialize parallel python
        if(self.enablePP):
            print "Parallel Python will be used with ",
            self.job_server = pp.Server()
            if (int(paramDict.get('cpuToUse', 0)) > 0):
                n = int(paramDict['cpuToUse'])
                print n, " cores"
                self.job_server.set_ncpus(n)
            else:
                print "all available cores"
                self.job_server.set_ncpus() # use all CPUs available
            self.jobs = []

    def sampleTau(self):
        """
        draw nSample realizations (samples) of the Reynolds stress fields [TauFields]

        :arg:
            None
        :return:
            None
        """
        # perform Cholesky factorization Rbar = Lbar.T dot Lbar
        LbarField = self._choleskyRMFields(self.Rbar, self.nCell_cfd)
        # Generate nSample identity random matrix fields, (LFields)
        GFields = self.randomMatrix.sample(self.nSample, self.pOrder)
        for isamp in np.arange(self.nSample):
            GFlatSample_i = self.intf.rm2tau(GFields[isamp, :, :, :], self.nCell_kl)
            self.GSample[isamp, :, :] = GFlatSample_i
            np.savetxt(self.resultDir+'G_samples/G_s'+str(isamp), GFlatSample_i)
        #pdb.set_trace()
        #TODO, mapping mesh (from KL to CFD needed to be considered)
        # For loop to construct Reynolds stress samples (1 : nSamples)
        for isamp in np.arange(self.nSample):
            for icell in np.arange(self.nCell_cfd):
                #pdb.set_trace()
                L_temp = LbarField[icell, :, :]
                G_temp = GFields[isamp, icell, :, :]
                LTG_temp = np.dot(L_temp.T, G_temp)
                LTGL_temp = np.dot(LTG_temp, L_temp)
                self.RFields[isamp, icell, :, :] = LTGL_temp
        # convert rm fields to tau fields (OpenFOAM format)
        for isamp in np.arange(self.nSample):
            # convert R to Tau fields
            TauSample_i = self.intf.rm2tau(self.RFields[isamp, :, :, :], self.nCell_cfd)
            self.TauSample[isamp, :, :] = TauSample_i
            self.deltaTauSample[isamp, :, :] = TauSample_i - self.Taubar
            ## output data
            np.savetxt(self.resultDir+'R_samples/R_s'+str(isamp), TauSample_i)

    def writeKLModesToOF(self):
        self.randomMatrix.GaussianRF.kl.writeKLModesToOF(self.klDir)

    def mappingTau2Phy(self):
        """
        map Tau to physical variable:
        :return:
        """
        # Tau component field
        self.XC1Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.XC2Fields = np.zeros([self.nSample, self.nCell_cfd])
        
        self.c1Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.c2Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.c3Fields = np.zeros([self.nSample, self.nCell_cfd])
                
        self.XiFields = np.zeros([self.nSample, self.nCell_cfd])
        self.EtaFields = np.zeros([self.nSample, self.nCell_cfd])
        self.KFields = np.zeros([self.nSample, self.nCell_cfd])
        self.VAFields = np.zeros([self.nSample, self.nCell_cfd])
        self.VBFields = np.zeros([self.nSample, self.nCell_cfd])
        self.VCFields = np.zeros([self.nSample, self.nCell_cfd])
                
        # delta component field
        self.deltaXiFields = np.zeros([self.nSample, self.nCell_cfd])
        self.deltaEtaFields = np.zeros([self.nSample, self.nCell_cfd])
        self.deltaLog2KFields = np.zeros([self.nSample, self.nCell_cfd])
        self.deltaVAFields = np.zeros([self.nSample, self.nCell_cfd])
        self.deltaVBFields = np.zeros([self.nSample, self.nCell_cfd])
        self.deltaVCFields = np.zeros([self.nSample, self.nCell_cfd])
        
        # initialize Tau2Phy mapping class
        # Note: we shut down the capping here, since the plotting will plot the exactly where tau is
        self.mapTau = ReynoldsStressRF('None', self.Taubar, self.nCell_cfd, 1, 'False')
                
        # Baseline component
        k,V1,V2,V3,C,NP = self.mapTau._tau2PhysParams(self.Taubar)
        X = self.mapTau._C2X(C)
        RS = self.mapTau._phys2Natural(X)        
        VA_base, VB_base, VC_base = self.mapTau.getThetaVABC(self.Taubar)
        
        self.c1Field_base = C[:, 0]
        self.c2Field_base = C[:, 1]
        self.c3Field_base = C[:, 2]
        self.XC1Field_base = X[:, 0]
        self.XC2Field_base = X[:, 1]        
        self.XiField_base = RS[:, 0]
        self.EtaField_base = RS[:, 1]
        self.kField_base = k[:, 0]        
        self.VAField_base = VA_base[:, 0]
        self.VBField_base = VB_base[:, 0]
        self.VCField_base = VC_base[:, 0]
        
        for isamp in np.arange(self.nSample):
            
            # samples of perturbed Ta
            TauNew = self.TauSample[isamp, :, :]
            k,V1,V2,V3,C,NP = self.mapTau._tau2PhysParams(TauNew) # collapse time = 1.005s (3000 cells)
            X = self.mapTau._C2X(C) # collapse time = 0.02s (3000 cells)
            RS = self.mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
            VA, VB, VC = self.mapTau.getThetaVABC(TauNew) # collapse time = 1.005s (3000 cells)
            
            self.c1Fields[isamp, :] = C[:, 0]
            self.c2Fields[isamp, :] = C[:, 1]
            self.c3Fields[isamp, :] = C[:, 2]                        
            self.XC1Fields[isamp, :] = X[:, 0]
            self.XC2Fields[isamp, :] = X[:, 1]
            self.XiFields[isamp, :] = RS[:, 0]
            self.EtaFields[isamp, :] = RS[:, 1]
            self.KFields[isamp, :] = k[:, 0]
            self.VAFields[isamp, :] = VA[:, 0]
            self.VBFields[isamp, :] = VB[:, 0]
            self.VCFields[isamp, :] = VC[:, 0]                                   
            
            # get delta components
            deltaXi, deltaEta, deltaLog2K = self.mapTau.getDeltaXiEtaLog2K(self.Taubar, TauNew)
            deltaVA, deltaVB, deltaVC = self._Theta2deltaTheta(VA, VB, VC, VA_base, VB_base, VC_base)
            
            self.deltaXiFields[isamp, :] = deltaXi[:, 0]
            self.deltaEtaFields[isamp, :] = deltaEta[:, 0]
            self.deltaLog2KFields[isamp, :] = deltaLog2K[:, 0]
            self.deltaVAFields[isamp, :] = deltaVA[:, 0]
            self.deltaVBFields[isamp, :] = deltaVB[:, 0]
            self.deltaVCFields[isamp, :] = deltaVC[:, 0]
                    
        # output
        np.savetxt(self.resultDir+'RComponent_samples/c1_s', self.c1Fields)
        np.savetxt(self.resultDir+'RComponent_samples/c2_s', self.c2Fields)
        np.savetxt(self.resultDir+'RComponent_samples/c3_s', self.c3Fields)        
        np.savetxt(self.resultDir+'RComponent_samples/XC1_s', self.XC1Fields)
        np.savetxt(self.resultDir+'RComponent_samples/XC2_s', self.XC2Fields)        
        np.savetxt(self.resultDir+'RComponent_samples/Xi_s', self.XiFields)
        np.savetxt(self.resultDir+'RComponent_samples/Eta_s', self.EtaFields)
        np.savetxt(self.resultDir+'RComponent_samples/TKE_s', self.KFields)
        np.savetxt(self.resultDir+'RComponent_samples/VA_s', self.VAFields)
        np.savetxt(self.resultDir+'RComponent_samples/VB_s', self.VBFields)
        np.savetxt(self.resultDir+'RComponent_samples/VC_s', self.VCFields)

        np.savetxt(self.resultDir+'RComponent_samples/c1_base', self.c1Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/c2_base', self.c2Field_base)        
        np.savetxt(self.resultDir+'RComponent_samples/c3_base', self.c3Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/XC1_base', self.XC1Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/XC2_base', self.XC2Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/Xi_base', self.XiField_base)
        np.savetxt(self.resultDir+'RComponent_samples/Eta_base', self.EtaField_base)
        np.savetxt(self.resultDir+'RComponent_samples/TKE_base', self.kField_base)        
        np.savetxt(self.resultDir+'RComponent_samples/VA_base', self.VAField_base)
        np.savetxt(self.resultDir+'RComponent_samples/VB_base', self.VBField_base)
        np.savetxt(self.resultDir+'RComponent_samples/VC_base', self.VCField_base)
                        
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaXi_s', self.deltaXiFields)
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaEta_s', self.deltaEtaFields)
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaK_s', self.deltaLog2KFields)
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaVA_s', self.deltaVAFields)
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaVB_s', self.deltaVBFields)
        np.savetxt(self.resultDir+'deltaRComponent_samples/deltaVC_s', self.deltaVCFields)
        



    def _Theta2deltaTheta(self, VA, VB, VC, VA_bar, VB_bar, VC_bar):
        """

        :return:
        """
        deltaVA = VA - VA_bar
        deltaVB = VB - VB_bar
        deltaVC = VC - VC_bar
        #pdb.set_trace()
        # modify deltaV
        deltaVA = (deltaVA > np.pi/2)*-np.pi + (deltaVA > np.pi/2)*deltaVA + (deltaVA < np.pi/2)*deltaVA
        deltaVA = (deltaVA < -np.pi/2)*np.pi + (deltaVA < -np.pi/2)*deltaVA + (deltaVA > -np.pi/2)*deltaVA

        deltaVB = (deltaVB > np.pi/2)*-np.pi + (deltaVB > np.pi/2)*deltaVB + (deltaVB < np.pi/2)*deltaVB
        deltaVB = (deltaVB < -np.pi/2)*np.pi + (deltaVB < -np.pi/2)*deltaVB + (deltaVB > -np.pi/2)*deltaVB

        deltaVC = (deltaVC > np.pi/2)*-np.pi + (deltaVC > np.pi/2)*deltaVC + (deltaVC < np.pi/2)*deltaVC
        deltaVC = (deltaVC < -np.pi/2)*np.pi + (deltaVC < -np.pi/2)*deltaVC + (deltaVC > -np.pi/2)*deltaVC

        return deltaVA, deltaVB, deltaVC

    def getBarySamples(self):
        """
        get C1 and C2 ensembles for Baycentric triangle
        :return:
        """
        self.XC1Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.XC2Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.c1Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.c2Fields = np.zeros([self.nSample, self.nCell_cfd])
        self.c3Fields = np.zeros([self.nSample, self.nCell_cfd])

        for isamp in np.arange(self.nSample):
            TauNew = self.TauSample[isamp, :, :]
            k,V1,V2,V3,C,NP = self.mapTau._tau2PhysParams(TauNew)
            X = self.mapTau._C2X(C)
            
            self.c1Fields[isamp, :] = C[:, 0]
            self.c2Fields[isamp, :] = C[:, 1]
            self.c3Fields[isamp, :] = C[:, 2]
            
            
            self.XC1Fields[isamp, :] = X[:, 0]
            self.XC2Fields[isamp, :] = X[:, 1]
        
        # Get base C1 and C2 field
        k,V1,V2,V3,C,NP = self.mapTau._tau2PhysParams(self.Taubar)
        X = self.mapTau._C2X(C)
        
        self.c1Field_base = C[:, 0]
        self.c2Field_base = C[:, 1]
        self.c3Field_base = C[:, 2]

        self.XC1Field_base = X[:, 0]
        self.XC2Field_base = X[:, 1]                

        np.savetxt(self.resultDir+'RComponent_samples/XC1_s', self.XC1Fields)
        np.savetxt(self.resultDir+'RComponent_samples/XC2_s', self.XC2Fields)
        np.savetxt(self.resultDir+'RComponent_samples/XC1_base', self.XC1Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/XC2_base', self.XC2Field_base)

        np.savetxt(self.resultDir+'RComponent_samples/c1_s', self.c1Fields)
        np.savetxt(self.resultDir+'RComponent_samples/c2_s', self.c2Fields)
        np.savetxt(self.resultDir+'RComponent_samples/c3_s', self.c3Fields)                
        np.savetxt(self.resultDir+'RComponent_samples/c1_base', self.c1Field_base)
        np.savetxt(self.resultDir+'RComponent_samples/c2_base', self.c2Field_base)        
        np.savetxt(self.resultDir+'RComponent_samples/c3_base', self.c3Field_base)
         
    def plotTau1LocSamplesOnBay(self, idx_p):
        """
        plot ensemble of one location of Tau fields on Baycentric triangle
        :return:
        """
        plt.clf()
        #pdb.set_trace()
        plt.plot(self.XC1Fields[:, idx_p],self.XC2Fields[:, idx_p],'ko')
        plt.scatter(self.XC1Field_base[idx_p], self.XC2Field_base[idx_p], s=50, c='red')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.show()


    def propagateTau(self):
        """
        propagate nSample tau via OpenFOAM to velocity

        :return:
            None
        """

        # prepare cases folder
        TauSample_thin = self._thiningSamples(self.nSample_propagate)
        self.intf.genFolders(self.nSample_propagate, self.caseName, self.caseNameRun, TauSample_thin, self.writeInterval)
        caseCount = np.linspace(1, self.nSample_propagate, self.nSample_propagate)
        isample = 0
        for case in caseCount:
            tmpCaseName = self.caseName + "-tmp_" + str(case)
            print "#", case, "/", self.nSample_propagate, " solving PDE equations to propagate to velocity"
            self._writeAllFieldToOpenFOAM(tmpCaseName, isample)
            if (self.enablePP):
                self.jobs.append(
                    # invoke the dynamic core
                    self.job_server.submit(foamOp.callFoam, (tmpCaseName, self.caseSolver, False))
                )
            else:
                foamOp.callFoam(tmpCaseName, self.caseSolver, False)
            isample = isample + 1

        print "Waiting for all OpenFOAM runs to finish ..."
        # Barrier to make sure all cases are finished before moviing on
        if(self.enablePP):
            part_sum1 = [job() for job in self.jobs]
            self.jobs = []

        part_sum1 = 0

        # finished all ensemble propagation. print stats
        if(self.enablePP):
            self.job_server.print_stats()
        # propagate baseline
        foamOp.callFoam(self.caseNameRun, self.caseSolver, False)

    def _writeAllFieldToOpenFOAM(self, tmpCaseName, isample):
        """

        :return:
        """

        foamOp.writeScalarToFile(self.deltaXiFields[isample, :], tmpCaseName+'/0/deltaXi')
        foamOp.writeScalarToFile(self.deltaEtaFields[isample, :], tmpCaseName+'/0/deltaEta')
        foamOp.writeScalarToFile(self.deltaLog2KFields[isample, :], tmpCaseName+'/0/deltaK')
        foamOp.writeScalarToFile(self.deltaVAFields[isample, :], tmpCaseName+'/0/deltaVA')
        foamOp.writeScalarToFile(self.deltaVBFields[isample, :], tmpCaseName+'/0/deltaVB')
        foamOp.writeScalarToFile(self.deltaVCFields[isample, :], tmpCaseName+'/0/deltaVC')
        foamOp.writeScalarToFile(self.XC1Fields[isample, :], tmpCaseName+'/0/C1Field')
        foamOp.writeScalarToFile(self.XC1Fields[isample, :], tmpCaseName+'/0/C2Field')
        foamOp.writeTurbStressToFile(self.deltaTauSample[isample, :], tmpCaseName+'/0/deltaTau')

    def _thiningSamples(self, nSample_thin):
        """

        :param nSample_thin:
        :return:
        """

        # TODO: the thining is naive, more advanced way to be used
        TauSample_thin = self.TauSample[0:nSample_thin, :]
        np.savetxt(self.resultDir+'RComponent_samples/C1_s_thin', self.XC1Fields[0:nSample_thin])
        np.savetxt(self.resultDir+'RComponent_samples/C2_s_thin', self.XC2Fields[0:nSample_thin])
        np.savetxt(self.resultDir+'RComponent_samples/C1_base_thin', self.XC1Field_base[0:nSample_thin])
        np.savetxt(self.resultDir+'RComponent_samples/C2_base_thin', self.XC2Field_base[0:nSample_thin])
        return TauSample_thin

    def _choleskyRMFields(self, RMField, nCell):
        """

        :param RMFields:
        :return:
        """
        #pdb.set_trace()
        # a 3 by 3 matrix for numerical stability
        # TODO this needs to be changed
        smallVal = 1e-16
        M_eps = smallVal * np.identity(3)
        LField = np.zeros(RMField.shape)
        for icell in np.arange(nCell):
            #print icell
            RM_temp = RMField[icell, :, :] + M_eps
            LT = LA.cholesky(RM_temp)
            LField[icell, :, :] = LT.T
        return LField


    def _prepareArg_covGen(self):
        """
        prepare arguments for covariance generation (Arg_covGen), call in _init_ function

        :return:
            Arg_covGen
        """
        ## generate length scale field in x, y and z directions, respectively
        lenXField, lenYField, lenZField = self._genLenScale()
        ## read cellVolumn Field as the weight Field for KL expansion
        cellVField = foamOp.readTurbCellVolumeFromFile(self.klDir+'/0/')
        weightField = cellVField[:, 0]
        sigmaField = foamOp.readScalarFromFile(self.klDir+'/0/Sigma')
        #TODO: truncateTol can be specified by user
        #truncateTol = -np.log(1e-10)
        # TruncateTol can be specified by user
        truncateTol = -np.log(self.truncateTol) 
        ## assembly arguments dictionary for generation of covariance
        if self.kernelType == 'SqExp': 
            Arg_covGen = {
                    'sigmaField': sigmaField,
                    'kernelType': 'SqExp',
                    'lenXField': lenXField,
                    'lenYField': lenYField,
                    'lenZField': lenZField,
                    'weightField':weightField,
                    'truncateTol': truncateTol
                 }
        elif self.kernelType == 'givenStructure':
            CovStruct = np.loadtxt('./Corr.dat')
            Arg_covGen = {
                    'sigmaField': sigmaField,
                    'kernelType': 'givenStructure',
                    'CovStruct': CovStruct,
                    'weightField':weightField,
                    'truncateTol': truncateTol
                 }
        return Arg_covGen

    def _prepareArg_calKLModes(self, nKL, weightField):
        """
        prepare arguments for KL expansion (Arg_calKLModes), call in _init_ function

        :return:
            Arg_calKLModes
        """
        Arg_calKLModes = {
                    'nKL': nKL,
                    'weightField':weightField
                }
        return Arg_calKLModes

    def _prepareDeltaField(self, stationaryFlag, deltaStationary, deltaFieldDir='None'):
        """

        :return:
        """
        if (stationaryFlag or deltaFieldDir == 'None'):
            "constant delta field is used, value = ", deltaStationary
            deltaField = np.ones([self.nCell_pce, 1]) * deltaStationary
        else:

            deltaField = foamOp.readScalarFromFile(deltaFieldDir)
            deltaField = np.array([deltaField]).T
        #pdb.set_trace()
        return deltaField

    def _genLenScale(self):
        """
        generate length scale field for x, y and z

        :arg:
        :return:
            lenXField:
            lenYField:
            lenZField:
        """
        # TODO: this function is an interface to specify non-stationary length scale field
        #pdb.set_trace()
        lenXField = self.lenXYZ[0, 0] * np.ones(self.nCell_kl)
        lenYField = self.lenXYZ[0, 1] * np.ones(self.nCell_kl)
        lenZField = self.lenXYZ[0, 2] * np.ones(self.nCell_kl)

        return lenXField, lenYField, lenZField
