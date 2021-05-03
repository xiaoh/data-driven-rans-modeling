# Copyright: Jianxun Wang (vtwjx@vt.edu)
# Feb.26, 2016

"""
    
    This module is an Dynamic model interface to specific forward model:
    OpenFoam Reynolds stress Tau Solver
    
    It  consist of 3 founctions:
    1 generateEnsemble: generate ensemble
    2 forcastToTime: evolve ensemble to next time using forward model (VegetationFOAM)
    3 getBackgroundVars: Get observations and Kalman Gain Matrix
    
"""

# system import
import sys
import os
import os.path as ospt
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
import scipy.sparse as sp
import pdb
import ast
import time


try:
    import pp
    hasPP = True
    print "Parallel Python imported successfully"
except ImportError, e:
    print e, 'I will not use parallel python in the simulation!'
    hasPP = False
    pass

# local import
from utilities import readInputData, replace, extractListFromDict
from sigmaFieldOperations import computeSigmaField
import deltaTauRandomField as ranF


iDebug = 1
debugFolderName = './debugData/'


class FoamVegetationSolver():
    """
    A particular dynamic forward-model interface: vegetation solver
    The state variable include: [U uw] = [u_1,v_1,w_1,...,u_n,v_n,w_n,uw_1,uw_2,...uw_n]
    The parameters need to be augmented: coefficients [omegaV] for Beta(x)
    """
    
    # static variables
    NTensor = 6                         # dimension of physical tensor
    NVec = 3                            # dimension of physical vector
    NScalar = 1                         # dimension of physical scalar

    def __init__(self, Ns, DAInterval, Tend,  ModelInput, OutputControl='NotExist'):
        """
        Initialization
        """
        ## Extract forwar Model Input parameters
        paramDict = readInputData(ModelInput)
        self.caseName = paramDict['caseName']
        self.caseNameObservation = self.caseName + "-observation"
        # The number of cells 
        self.Ncell = int(paramDict['Ncell'])        
        # The number of cells which are observed
        self.NcellObs = int(paramDict['NcellObs'])
        # hyperparameter for KL expansion
        # --- length scales
        lenConstFlag = ast.literal_eval(paramDict['lenConstFlag'])
        lenVec = extractListFromDict(paramDict, 'lenVec')
        self.lenVec = np.array([[float(pn) for pn in lenVec]])
        # --- variance field (for Beta(x, y, z))
        sigmaConstFlag = ast.literal_eval(paramDict['sigmaConstFlag'])
        self.sigma = float(paramDict['sigmaVec'])
        # --- number of KL modes                 
        self.nModes = int(paramDict['Nmode'])
        # --- kernel type for KL expansion
        self.kernelType = paramDict['kernelType']        
        # enable paralization
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
        
        ## Start KL expansion for random fields                       
        # Read cell center coordinates
        self.baseCaseDir = ospt.join(os.getcwd(), self.caseName)
        meshCoord3D = np.loadtxt(self.baseCaseDir + '/grid.dat') 
        volCoord = np.loadtxt(self.baseCaseDir + '/cellVolume.dat')
        # construct variance field sigma(x) for Xi, Eta, K, Vs
        if sigmaConstFlag:
            print "\nConstant sigma fields are used" 
            print "sigma = ", self.sigma
            self.sigmaField = np.ones(self.Ncell) * self.sigma
        else:
            print "\nNonstationary sigma fields are used" 
            rbfLengthScale = float(paramDict['rbfLengthScale'])
            rbfKernel = paramDict['rbfKernel']       
            scatteredSigmaFile = ospt.join(self.caseName, 'scatSigma.dat')        
            if ospt.exists(scatteredSigmaFile):            
                self.sigma = np.absolute(computeSigmaField(scatteredSigmaFile, meshCoord3D, \
                                       rbfKernel, rbfLengthScale))
            else:
                print "WARNING: I cannot find file: scatSigma.dat, constant sigma is used"
                self.sigmaField = np.ones(self.Ncell) * self.sigma
        
        # construct length scale field        
        if lenConstFlag: 
            print "\nstationary length scale fields are used"             
            lenXField = np.ones(self.Ncell) * self.lenVec[0, 0]
            lenYField = np.ones(self.Ncell) * self.lenVec[0, 1]
            lenZField = np.ones(self.Ncell) * self.lenVec[0, 2]        
        else:
            print "\nNonstationary length scale fields are used"
            scatteredLxFile = ospt.join(self.caseName, 'scatLx.dat')
            scatteredLyFile = ospt.join(self.caseName, 'scatLy.dat')
            scatteredLzFile = ospt.join(self.caseName, 'scatLz.dat')         
        
            if ospt.exists(scatteredLxFile) and ospt.exists(scatteredLyFile) \
               and ospt.exists(scatteredLzFile):         
                lenXField = np.absolute(computeSigmaField(scatteredLxFile, meshCoord3D, \
                                               rbfKernel, rbfLengthScale))
                lenYField = np.absolute(computeSigmaField(scatteredLyFile, meshCoord3D, \
                                               rbfKernel, rbfLengthScale))
                lenYField = np.absolute(computeSigmaField(scatteredLyFile, meshCoord3D, \
                                               rbfKernel, rbfLengthScale))
            else:
                print "WARNING: scatLx.dat or scatLy.dat or scatLz.dat cannot be found!"
                lenXField = np.ones(self.Ncell) * self.lenVec[0, 0]
                lenYField = np.ones(self.Ncell) * self.lenVec[0, 1]
                lenZField = np.ones(self.Ncell) * self.lenVec[0, 2]                    

        # prepare the argument for generate covariance matrix for all components
        truncateTol = -np.log(1e-10)         
        Arg_covGen = {
                            'sigmaField': self.sigmaField,
                            'lenXField': lenXField,
                            'lenYField': lenYField,
                            'lenZField': lenZField,
                            'weightField':volCoord,
                            'truncateTol': truncateTol                            
                      }
        Arg_calModes = {
                            'nKL': self.nModes,
                            'weightField':volCoord
                        }
        
        # initialize the KL expansion class                                                    
        self.rf = ranF.randomField(meshCoord3D, Arg_covGen, Arg_calModes, 'beta')          
        # determine if perform KL expansion                                      
        klCalculate =  ast.literal_eval(paramDict['klCalculate'])        
        modeExist = os.path.exists('randomData_beta/KLModes.dat')        
        if (klCalculate):
            print "\nStart KL expansions ......"
            self.KLModes = self.rf.KLExpansion()
        elif (not modeExist):               
            print "\nWARNING: Cannot find KL files, performing KL expansion calculation"        
            self.KLModes = self.rf.KLExpansion()
        else:
            print "\nLoading KL expansion data, won't perform KL calculation"        
            self.KLModes = np.loadtxt('randomData_beta/KLModes.dat')
       

        # clear output data figures
        os.system('rm -rf ' + debugFolderName + '*')
        os.system('mkdir ' + debugFolderName +'init')
        os.system('mkdir ' + debugFolderName + 'init-propagate')        
        # Ns omega vector
        self.omegaEns = np.zeros([self.nModes, Ns]) # need to be part of parameters
        # Ns beta fields
        self.betaFieldEns = np.zeros([Ns, self.Ncell, FoamVegetationSolver.NScalar])        
        tic_in = time.time()
        for i in range(Ns):                 
            omega = self.rf.uncorrUniRandom(self.nModes)
            self.omegaEns[:, i] = omega[:,0]                                                      
            self.betaField = self.rf.reconstructField(omega, self.KLModes)       
            self.betaFieldEns[i, :, :] = self.betaField
            
            np.savetxt(debugFolderName+'init/'+'beta_init_s-'+str(i), self.betaField) 
        # unknown parameter ensemble
        self.XC = self.omegaEns
        pdb.set_trace()    
        self.Ns = Ns   # number of samples in the ensemble
        # number of variables in each cell (not including parameter)
        #TODO: check every cell if u, v, w, uw?
        self.Nvariable = FoamVegetationSolver.NVec + FoamVegetationSolver.NScalar                        
        self.Npara = self.nModes
        self.Nstate = (self.Ncell * self.Nvariable) + self.Npara # size of state variables
        #TODO: make sure observation in each cell is u, v, w, uw?
        self.NstateObs = self.NcellObs * self.Nvariable # size of observed state variables        

        # Ensemble state vectors stacked in matrices
        # Matrix (NVec*Ncell, Ns): state U, order: (u, v, w, uw)
        self.XU = np.zeros([self.Nvariable * self.Ncell, Ns])
        # Matrix (Npara, Ns): state parameter  gf

        # Ensemble states project to observation space
        # ensemble matrix in observation space (Ns, NVec*NcellObs), (u, v, w, uw)
        self.HXU = np.zeros([self.Nvariable * self.NcellObs, Ns])

        # Observation vectors stacked in matrices
        # state matirx of observation: part U (u,v,w)
        self.ObsU = np.zeros([self.Nvariable * self.NcellObs, Ns])

        # Observation error covariance (could be given or default generate)
        if self.pseudoObs == 0:
            #This is for real experiment Data            
            self.rmu = float(paramDict['ObserveMean'])
            ObsSigmaFixedVec = extractListFromDict(paramDict, 'ObsSigmaFixedVec')                         
            self.ObsSigmaFixedVec = [float(pn) for pn in ObsSigmaFixedVec]

            ObsRelCoeffVec = extractListFromDict(paramDict, 'ObsRelCoeffVec')                
            self.ObsRelCoeffVec = [float(pn) for pn in ObsRelCoeffVec]
            
            self.ObsRmaCoeff = float(paramDict['ObsRmaCoeff'])
            self.ObsErrCoeff = float(paramDict['ObsErrCoeff'])                
            
            rsigmaVecSq = [pn**2 for pn in self.ObsSigmaFixedVec]
            self.Robs = sp.diags(np.array(rsigmaVecSq * self.NcellObs),0)

        self.DAInterval = DAInterval                # DA step interval
    
    def generateEnsemble(self, DAInterval):
        """ Generate OF case folders and X, HX

        Args:
        DAInterval: DA step interval

        Returns:
        X: ensemble matrix of whole states
        HX: ensemble matrix of whole states in observation space
        """
        
    def forecastToTime(self, X, nextEndTime):
        """ Call the dynamic core to evolve the ensemble states to time 

        Args:
        X: ensemble matrix of whole state at current time
        nextEndTime: target time of this forecast

        Returns:
        X: ensemble matrix of whole states at next DA time
        HX: projection of X (above) in observation space
        """

    def getBackgroundVars(self, HX, XP, nextEndTime):

        """ Function is to generate observation and get kalman Gain Matrix

        Args:
        HX: ensemble matrix of whole state in observation space
        P: covariance matrix of ensemble
        nextEndTime: next DA interval

        Returns:
        Obs: state matrix of observation
        KalmanGainMatrix
        """
                
