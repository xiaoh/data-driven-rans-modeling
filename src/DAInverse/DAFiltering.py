# copyright: Jianxun Wang (vtwjx@vt.edu)
# Mar. 31, 2015
# Modified by Carlos Michelen-Strofer 2017
# 
# This module contains different filtering methods, each implemented as its own 
# child class of the DAFilter class.
# Currently, only the ensemble Kalman filter (EnKF) is implemented.
# Requirements for child classes: 
#   * Initialize with the use of an input file (class dependent).
#   * Report basic instance information with the '__str__' method.
#   * Have a 'solve()' method requiring no extra inputs.
#   * Have a 'report()' method requiring no extra inputs.
#   * Inherit the 'plot()' method from the parent class (DAFilter).
#
# TODO (CM):
#   * Move most of the content of _correctForecast() to solve()
#   * Calculate Kalman Gain here, not in the dynamic model (dyn model provides X, HX, R, Obs)
#   * Move all the misfit calculations/reporting and debugging out of solve.

# system modules
import numpy as np
import os
import sys
import ast
import importlib
from numpy import linalg as LA

# local modules
from utilities import readInputData
import DAInverse.dynModels as dynModels


# parent class
class DAFilter:

    def __init__(self, InputFile):
        ''' Parse input file and assign values to class attributes.
        '''
        self.dynModel = None

    def __str__(self):
        """ Print basic summary information of the filtering technique.
        """
        s = 'An empty data assimilation filtering technique ...'
        return s

    def solve(self):
        ''' Implementation of the filtering technique.
        '''
        pass

    def report(self):
        ''' Report summary information at each step.
        '''
        pass

    def plot(self):
        """ Call the dynamic model to plot mean, variance, error etc.
        """
        self.dynModel.plot()

# specific filtering techniques (child classes)
class EnKF(DAFilter):
    """ Ensemble Kalman Filter.
        Class to implement functionalities of Ensemble Kalman Filter
    """
    
    def __init__(self, InputFile):
        """ Parse input file and initialize parameters """
        ## Parse InputFile
        paramDict = readInputData(InputFile)        
        
        # EnKF inputs
        self.Tend = float(paramDict["Tend"]) # total run time     
        self.DAInterval = float(paramDict['DAInterval']) # DA time step interval                 
        self.Ns = int(paramDict['Ns']) # number of sample
        self.convergenceResi = float(paramDict['convergenceResi'])
        self.reachmaxiteration = ast.literal_eval(paramDict['reachmaxiteration'])
        
        # EnKF private variables
        self._sensitivityOnly = ast.literal_eval(paramDict['sensitivityOnly'])
        self._iDebug = paramDict['iDebug']
        self._debugFolderName = paramDict['debugFolderName']
        if not os.path.exists(self._debugFolderName):
            os.makedirs(self._debugFolderName)  

        # Create instance of dynamic model class 
        self.forwardModel = paramDict['forwardModel']
        self.forwarModelInput = paramDict['forwardModelInput']
        dynModel = getattr(importlib.import_module('dynModels.' + self.forwardModel), self.forwardModel)
        self.dynModel = dynModel(self.Ns, self.DAInterval, self.Tend, self.forwarModelInput)
        
        ## Initialize time
        self.time = 0 # current time
        self.T = np.arange(0, self.Tend, self.DAInterval)
        
        ## Initialize states 
        self.X = np.zeros([self.dynModel.Nstate, self.Ns]) # ensemble matrix (Ns, Nstate)
        self.HX = np.zeros([self.dynModel.NstateObs, self.Ns]) # ensemble matrix projected to observed space (Ns, NstateSample)
        self.Obs = np.zeros([self.dynModel.NstateObs, self.Ns]) # observation matrix (Ns, NstateSample)
        
        ## Initialize misfit       
        self.misfitX_L1 = []; 
        self.misfitX_L2 = [];
        self.misfitX_inf = [];
        self.obsSigma = [];
        self.obsSigmaInf = [];
        self.sigmaHX = [];                 
    
    def __str__(self):
        """ Print basic summary information of the ensemble, e.g.,
            model name, number of samples, DA interval
        """
        s = 'Ensemble Kalman Filter' + \
            '\n   Number of samples: {}'.format(self.Ns) + \
            '\n   Run time:          {}'.format(self.Tend) +  \
            '\n   DA interval:       {}'.format(self.DAInterval) + \
            '\n   Forward model:     {}'.format(self.forwardModel)
        return s

    def solve(self):
        """ Solves the parameter estimation problem.
        """
        # Generate initial state Ensemble
        (self.X, self.HX) = self.dynModel.generateEnsemble()
        if(self._sensitivityOnly):
            print "Sensitivity study completed."
            sys.exit(0)
        
        # Main DA loop
        ii = 0
        for t in self.T:
            nextEndTime = 2 * self.DAInterval + t
            DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
            print "#######################################################################"
            print "\n Data Assimilation step = ", DAstep, "\n"   
            
            # propagate the state ensemble to next DA time
            self.X, self.HX = self.dynModel.forecastToTime(self.X, nextEndTime)
            if (self._iDebug):
                np.savetxt(self._debugFolderName+'X_'+ str(DAstep) + '.txt', self.X)
                np.savetxt(self._debugFolderName+'HX_'+ str(DAstep) + '.txt', self.HX)
            
            # correct the propagated results     
            self.X, sigmaHXNorm, misfitMax = self._correctForecasts(self.X, self.HX, nextEndTime)
            
            #Check iteration convergence and report misfits                        
            Robs = self.dynModel.get_Robs()
            sigmaInf = np.max(np.diag(np.sqrt(Robs)));         
            sigma = ( LA.norm(np.sqrt(np.diag(Robs))) / self.dynModel.NstateObs)
            print "Std of ensemble = ", sigmaHXNorm
            print "Std of observation error", sigma
            self.obsSigma.append(sigma)
            self.obsSigmaInf.append(sigmaInf)
            if  misfitMax > sigmaInf:
                print "infinit misfit(", misfitMax, ") is larger than Inf norm of observation error(", sigmaInf, "), so iteration continuing\n\n"
            else:
                print "infinit misfit of ensemble reaches Inf norm of observation error, considering converge\n\n"
                self.misfitX_L1 = np.array(self.misfitX_L1)
                self.misfitX_L2 = np.array(self.misfitX_L2)
                self.misfitX_inf = np.array(self.misfitX_inf)
                self.obsSigma = np.array(self.obsSigma)
            if ii > 0:
                iterativeResidual = (abs(self.misfitX_L2[ii] - self.misfitX_L2[ii-1])) /abs(self.misfitX_L2[0])    
                print "relative Iterative residual = ", iterativeResidual
                print "relative Convergence criterion = ", self.convergenceResi 
                if iterativeResidual  < self.convergenceResi:    
                    if self.reachmaxiteration:
                        pass
                    else:
                        print "Iteration converges, finished \n\n"
                        break
                np.savetxt('./misfit_L1.txt', self.misfitX_L1)
                np.savetxt('./misfit_L2.txt', self.misfitX_L2)
                np.savetxt('./misfit_inf.txt', self.misfitX_inf)
                np.savetxt('./obsSigma.txt', self.obsSigma)
                np.savetxt('./obsSigmaInf.txt', self.obsSigmaInf)
            ii = ii + 1
        
        # Save misfits
        self.misfitX_L1 = np.array(self.misfitX_L1)
        self.misfitX_L2 = np.array(self.misfitX_L2)
        self.misfitX_inf = np.array(self.misfitX_inf)
        self.obsSigma = np.array(self.obsSigma)
        self.sigmaHX = np.array(self.sigmaHX)

        np.savetxt('./misfit_L1.txt', self.misfitX_L1)
        np.savetxt('./misfit_L2.txt', self.misfitX_L2)
        np.savetxt('./misfit_inf.txt', self.misfitX_inf)
        np.savetxt('./obsSigma.txt', self.obsSigma)
        np.savetxt('./sigmaHX.txt', self.sigmaHX)
        # os.system('plotIterationConvergence.py')            
    
    def report(self):
        """ Report summary information at each step
        """
        raise NotImplementedError

    def clean(self):
        """ Call the dynamic model to do any necessary cleanup before exiting.
        """
        self.dynModel.clean()
    
    ############################ Priviate function ################################    
    def _correctForecasts(self, X, HX, nextEndTime):
        """ Filtering step: Correct the propagated ensemble X 
            via EnKF filtering procedure
            
            Arg:
            X: state ensemble matrix (Nstate by Ns)
            HX: state ensemble matrix projected to observed space
            nextEndTime: next DA time spot
            
            Return:
            X: state ensemble matrix
        """
        DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
        # get XMean and tile XMeanVec to a full matrix (Nsate by Ns)
        XMeanVec = np.mean(X, axis=1) # vector mean of X
        XMeanVec = np.array([XMeanVec])
        XMeanVec = XMeanVec.T
        XMean = np.tile(XMeanVec, (1, self.Ns))
        
        # get coveriance matrix P
        XP = X - XMean
        
        # calculate Kalman Gain matrix and get observation
        (Obs, KalmanGainMatrix) = self.dynModel.getBackgroundVars(HX, XP, nextEndTime)        
        if (self._iDebug):
            np.savetxt(self._debugFolderName+'Obs_'+str(DAstep), Obs)    

        # analyze step
        dX = np.dot(KalmanGainMatrix, Obs - HX)
        X = X + dX
        if (self._iDebug):
            np.savetxt(self._debugFolderName+'kalmanGain_'+str(DAstep), KalmanGainMatrix  )
            np.savetxt(self._debugFolderName+'updateX_'+str(DAstep), X) 
        
        # calculate misfits
        Nnorm = self.dynModel.NstateObs * self.Ns
        diff = abs(Obs - HX)
        misFit_L1 = np.sum(diff) / Nnorm
        misFit_L2 = np.sqrt(np.sum(diff**2)) / Nnorm
        misFit_Inf = LA.norm(diff, np.inf)                
        misFitMax = np.max([misFit_L1, misFit_L2, misFit_Inf])
        sigmaHX = np.std(HX, axis=1)
        sigmaHXNorm = LA.norm(sigmaHX) / self.dynModel.NstateObs
        
        self.misfitX_L1.append(misFit_L1)
        self.misfitX_L2.append(misFit_L2)
        self.misfitX_inf.append(misFit_Inf)
        self.sigmaHX.append(sigmaHXNorm)
        
        print "\nThe misfit between predicted QoI with observed QoI is:"
        print "L1 norm = ", misFit_L1,", \nL2 norm = ", misFit_L2, "\nInf norm = ", misFit_Inf, "\n\n"
    
            
     
        return X, sigmaHXNorm, misFitMax


# specific filtering techniques (child classes)
class EnRML(DAFilter):
    """ Ensemble Randomized Maximum Likelihood.
        Class to implement functionalities of EnRML Method
    """

    def __init__(self, InputFile):
        ## Parse InputFile
        paramDict = readInputData(InputFile)        
        
        # EnRML inputs
        self.Tend = float(paramDict["Tend"]) # total run time     
        self.DAInterval = float(paramDict['DAInterval']) # DA time step interval                 
        self.Ns = int(paramDict['Ns']) # number of sample
        self.beta = float(paramDict["beta"])
        self.constbeta = ast.literal_eval(paramDict['constbeta'])
        self.convergenceResi = float(paramDict['convergenceResi'])
        self.reachmaxiteration = ast.literal_eval(paramDict['reachmaxiteration'])
        # EnRML private variables
        self._sensitivityOnly = ast.literal_eval(paramDict['sensitivityOnly'])
        self._iDebug = paramDict['iDebug']
        self._debugFolderName = paramDict['debugFolderName']
        if not os.path.exists(self._debugFolderName):
            os.makedirs(self._debugFolderName) 
        # Create instance of dynamic model class 
        self.forwardModel = paramDict['forwardModel']
        self.forwarModelInput = paramDict['forwardModelInput']
        dynModel = getattr(importlib.import_module('dynModels.' + self.forwardModel), self.forwardModel)
        self.dynModel = dynModel(self.Ns, self.DAInterval, self.Tend, self.forwarModelInput)
    
        ## Initialize time
        self.time = 0 # current time
        self.T = np.arange(0, self.Tend, self.DAInterval)
        
        ## Initialize states 
        self.Xpr = np.zeros([self.dynModel.Nstate, self.Ns]) # prior ensemble matrix (Ns, Nstate)
        self.Xf = np.zeros([self.dynModel.Nstate, self.Ns]) # forcast ensemble matrix (Ns, Nstate)
        self.Xa = np.zeros([self.dynModel.Nstate, self.Ns]) # analysis ensemble matrix (Ns, Nstate)
        self.HX = np.zeros([self.dynModel.NstateObs, self.Ns]) # ensemble matrix projected to observed space (Ns, NstateSample)
        self.Obs = np.zeros([self.dynModel.NstateObs, self.Ns]) # observation matrix (Ns, NstateSample)
        
        ## Initialize misfit       
        self.misfitX_L1 = []; 
        self.misfitX_L2 = [];
        self.misfitX_inf = [];
        self.obsSigma = [];
        self.obsSigmaInf = [];
        self.sigmaHX = [];
        self.steplength = [];

    def __str__(self):
        """ Print basic summary information of the ensemble, e.g.,
            model name, number of samples, DA interval
        """
        s = 'Ensemble Randomized Maximum Likelihood Method' + \
            '\n   Number of samples: {}'.format(self.Ns) + \
            '\n   Run time:          {}'.format(self.Tend) +  \
            '\n   DA interval:       {}'.format(self.DAInterval) + \
            '\n   Forward model:     {}'.format(self.forwardModel)
        return s

    def solve(self):
        """ Solves the parameter estimation problem.
        """
        # Generate initial state Ensemble
        (self.Xa, self.HX) = self.dynModel.generateEnsemble()
        if(self._sensitivityOnly):
            print "Sensitivity study completed."
            sys.exit(0)
        ii = 0
        self.Xpr = self.Xa.copy()
        Obs = self.dynModel._Observe(self.DAInterval)
        Robs = self.dynModel.get_Robs()
        for t in self.T:
            nextEndTime = 2 * self.DAInterval + t
            DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
            print "#######################################################################"
            print "\n Data Assimilation step = ", DAstep, "\n"   
            
            # propagate the state ensemble to next DA time
            self.Xf, self.HX = self.dynModel.forecastToTime(self.Xa.copy(), nextEndTime)
            if (self._iDebug):
                np.savetxt(self._debugFolderName+'X_'+ str(DAstep) + '.txt', self.Xf)
                np.savetxt(self._debugFolderName+'HX_'+ str(DAstep) + '.txt', self.HX)
            
            # correct the propagated results
            self.Xa, sigmaHXNorm, misfitMax = self._correctForecasts(self.Xpr, self.Xf, self.HX,Obs, Robs, nextEndTime)
            #Check iteration convergence and report misfits (Todo:move to report function)                      
            sigmaInf = np.max(np.diag(np.sqrt(Robs)));
            sigma = ( LA.norm(np.sqrt(np.diag(Robs))) / self.dynModel.NstateObs)
            print "Std of ensemble = ", sigmaHXNorm
            print "Std of observation error", sigma
            self.obsSigma.append(sigma)
            self.obsSigmaInf.append(sigmaInf)
            if  misfitMax > sigmaInf:
                print "infinit misfit(", misfitMax, ") is larger than Inf norm of observation error(", sigmaInf, "), so iteration continuing\n\n"
            else:
                print "infinit misfit of ensemble reaches Inf norm of observation error, considering converge\n\n"
                self.misfitX_L1 = np.array(self.misfitX_L1)
                self.misfitX_L2 = np.array(self.misfitX_L2)
                self.misfitX_inf = np.array(self.misfitX_inf)
                self.obsSigma = np.array(self.obsSigma)
            if ii > 0:
                iterativeResidual = (abs(self.misfitX_L2[ii] - self.misfitX_L2[ii-1])) /abs(self.misfitX_L2[0])    
                print "relative Iterative residual = ", iterativeResidual
                print "relative Convergence criterion = ", self.convergenceResi 
        
                if iterativeResidual  < self.convergenceResi:
                    if self.reachmaxiteration:
                        pass
                    else:
                        print "Iteration converges, finished \n\n"
                        break
                
                np.savetxt('./misfit_L1.txt', self.misfitX_L1)
                np.savetxt('./misfit_L2.txt', self.misfitX_L2)
                np.savetxt('./misfit_inf.txt', self.misfitX_inf)
                np.savetxt('./obsSigma.txt', self.obsSigma)
                np.savetxt('./obsSigmaInf.txt', self.obsSigmaInf)
                if self.constbeta:
                    print "Use constant step length beta = " ,self.beta
                    pass
                else:
                    print "Use self-adjusting step length \n"
                    if self.misfitX_L2[ii] < self.misfitX_L2[ii-1]:
                        self.beta = 1.2*self.beta
                        print "Iteration Converging, increase step length beta = ", self.beta
                    else:
                        self.beta = 0.8*self.beta
                        self.Xa = self.Xf.copy()
                        print "Iteration Diverging, reduce step length beta = ", self.beta
            self.steplength.append(self.beta)
            ii = ii + 1
        # Save misfits
        self.misfitX_L1 = np.array(self.misfitX_L1)
        self.misfitX_L2 = np.array(self.misfitX_L2)
        self.misfitX_inf = np.array(self.misfitX_inf)
        self.obsSigma = np.array(self.obsSigma)
        self.sigmaHX = np.array(self.sigmaHX)
        self.steplength = np.array(self.steplength)

        np.savetxt('./misfit_L1.txt', self.misfitX_L1)
        np.savetxt('./misfit_L2.txt', self.misfitX_L2)
        np.savetxt('./misfit_inf.txt', self.misfitX_inf)
        np.savetxt('./obsSigma.txt', self.obsSigma)
        np.savetxt('./sigmaHX.txt', self.sigmaHX)
        np.savetxt('./steplength.txt',self.steplength)

    def report(self):
        """ Report summary information at each step
        """
        raise NotImplementedError

    def clean(self):
        """ Call the dynamic model to do any necessary cleanup before exiting.
        """
        self.dynModel.clean()

    ############################ Priviate function ################################    
    def _correctForecasts(self, Xpr, X, HX, Obs,Robs, nextEndTime):

        """ Filtering step: Correct the propagated ensemble X 
            via EnRML filtering procedure
            
            Arg:
            Xpr: Prior state ensemble matrix (Nstate by Ns)
            X: state ensemble matrix (Nstate by Ns)
            HX: state ensemble matrix projected to observed space
            Obs: Observation
            nextEndTime: next DA time spot
            
            Return:
            X: state ensemble matrix
        """
        DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
        # get XMean and tile XMeanVec to a full matrix (Nsate by Ns)
        XMeanVec = np.mean(X, axis=1) # vector mean of X
        XMeanVec = np.array([XMeanVec])
        XMeanVec = XMeanVec.T
        XMean = np.tile(XMeanVec, (1, self.Ns))
        # get prior XMean and tile XMeanVec to a full matrix (Nsate by Ns)
        XprMeanVec = np.mean(Xpr, axis=1) # vector mean of X
        XprMeanVec = np.array([XprMeanVec])
        XprMeanVec = XprMeanVec.T
        XprMean = np.tile(XprMeanVec, (1, self.Ns))
        # get coveriance matrix
        XP0 = Xpr - XprMean
        XP =  X - XMean
        
        (GNGainMatrix, penalty) = self.dynModel.getGaussNewtonVars(HX, X,Xpr,XP, XP0, Obs, Robs, nextEndTime)

        # analyze step
        # Robs = self.dynModel.get_Robs()
        Cdd = np.diag(np.diag(np.sqrt(Robs)))
        #dX = - np.dot(Cdd.dot(LA.inv(Hess)),(X-Xpr)) - GNGainMatrix.dot(HX-Obs)
        dX = - GNGainMatrix.dot(HX-Obs) + penalty
        X = self.beta*Xpr + (1.0-self.beta)* X+ self.beta*dX
        if (self._iDebug):
            np.savetxt(self._debugFolderName+'dX_'+str(DAstep), dX)
            np.savetxt(self._debugFolderName+'updateX_'+str(DAstep), X) 
            np.savetxt(self._debugFolderName+'Obs_'+str(DAstep), Obs)  
        
        # calculate misfits
        Nnorm = self.dynModel.NstateObs * self.Ns
        diff = abs(Obs - HX)
        misFit_L1 = np.sum(diff) / Nnorm
        misFit_L2 = np.sqrt(np.sum(diff**2)) / Nnorm
        misFit_Inf = LA.norm(diff, np.inf)                
        misFitMax = np.max([misFit_L1, misFit_L2, misFit_Inf])
        sigmaHX = np.std(HX, axis=1)
        sigmaHXNorm = LA.norm(sigmaHX) / self.dynModel.NstateObs
        
        self.misfitX_L1.append(misFit_L1)
        self.misfitX_L2.append(misFit_L2)
        self.misfitX_inf.append(misFit_Inf)
        self.sigmaHX.append(sigmaHXNorm)
        
        print "\nThe misfit between predicted QoI with observed QoI is:"
        print "L1 norm = ", misFit_L1,", \nL2 norm = ", misFit_L2, "\nInf norm = ", misFit_Inf, "\n\n"
    
        return X, sigmaHXNorm, misFitMax


# specific filtering techniques (child classes)
class EnKFMDA(DAFilter):
    """ Ensemble kalman filter-Multi Data Assimilation.
        Class to implement functionalities of EnKF Multiple-Data-Assimilation
    """
    
    def __init__(self, InputFile):
        """ Parse input file and initialize parameters """
        ## Parse InputFile
        paramDict = readInputData(InputFile)        
        
        # EnKF-MDA inputs
        self.Tend = float(paramDict["Tend"]) # total run time     
        self.DAInterval = float(paramDict['DAInterval']) # DA time step interval                 
        self.Ns = int(paramDict['Ns']) # number of sample
        self.convergenceResi = float(paramDict['convergenceResi'])
        self.reachmaxiteration = ast.literal_eval(paramDict['reachmaxiteration'])
        
        # EnKF-MDA private variables
        self._sensitivityOnly = ast.literal_eval(paramDict['sensitivityOnly'])
        self._iDebug = paramDict['iDebug']
        self._debugFolderName = paramDict['debugFolderName']
        if not os.path.exists(self._debugFolderName):
            os.makedirs(self._debugFolderName)  

        # Create instance of dynamic model class 
        self.forwardModel = paramDict['forwardModel']
        self.forwarModelInput = paramDict['forwardModelInput']
        dynModel = getattr(importlib.import_module('dynModels.' + self.forwardModel), self.forwardModel)
        self.dynModel = dynModel(self.Ns, self.DAInterval, self.Tend, self.forwarModelInput)
        
        ## Initialize time
        self.time = 0 # current time
        self.T = np.arange(0, self.Tend, self.DAInterval)
        
        ## Initialize states 
        self.X = np.zeros([self.dynModel.Nstate, self.Ns]) # ensemble matrix (Ns, Nstate)
        self.HX = np.zeros([self.dynModel.NstateObs, self.Ns]) # ensemble matrix projected to observed space (Ns, NstateSample)
        self.Obs = np.zeros([self.dynModel.NstateObs, self.Ns]) # observation matrix (Ns, NstateSample)
        self.PertObs = np.zeros([self.dynModel.NstateObs, self.Ns]) # observation perturbation matrix (Ns, NstateSample)
        
        ## Initialize misfit       
        self.misfitX_L1 = []; 
        self.misfitX_L2 = [];
        self.misfitX_inf = [];
        self.obsSigma = [];
        self.obsSigmaInf = [];
        self.sigmaHX = [];                 
    
    def __str__(self):
        """ Print basic summary information of the ensemble, e.g.,
            model name, number of samples, DA interval
        """
        s = 'Ensemble Kalman Filter' + \
            '\n   Number of samples: {}'.format(self.Ns) + \
            '\n   Run time:          {}'.format(self.Tend) +  \
            '\n   DA interval:       {}'.format(self.DAInterval) + \
            '\n   Forward model:     {}'.format(self.forwardModel)
        return s

    def solve(self):
        """ Solves the parameter estimation problem.
        """
        # Generate initial state Ensemble
        (self.X, self.HX) = self.dynModel.generateEnsemble()
        if(self._sensitivityOnly):
            print "Sensitivity study completed."
            sys.exit(0)
        # Main DA loop
        ii = 0
        for t in self.T:
            nextEndTime = 2 * self.DAInterval + t
            DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
            print "#######################################################################"
            print "\n Data Assimilation step = ", DAstep, "\n"   
            
            # propagate the state ensemble to next DA time
            self.X, self.HX = self.dynModel.forecastToTime(self.X, nextEndTime)
            if (self._iDebug):
                np.savetxt(self._debugFolderName+'X_'+ str(DAstep) + '.txt', self.X)
                np.savetxt(self._debugFolderName+'HX_'+ str(DAstep) + '.txt', self.HX)

            # correct the propagated results     
            self.X, sigmaHXNorm, misfitMax = self._correctForecasts(self.X, self.HX, nextEndTime)
            
            #Check iteration convergence and report misfits                        
            Robs = self.dynModel.get_Robs()
            sigmaInf = np.max(np.diag(np.sqrt(Robs)));
            sigma = ( LA.norm(np.sqrt(np.diag(Robs))) / self.dynModel.NstateObs)
            print "Std of ensemble = ", sigmaHXNorm
            print "Std of observation error", sigma
            self.obsSigma.append(sigma)
            self.obsSigmaInf.append(sigmaInf)
            if  misfitMax > sigmaInf:
                print "infinit misfit(", misfitMax, ") is larger than Inf norm of observation error(", sigmaInf, "), so iteration continuing\n\n"
            else:
                print "infinit misfit of ensemble reaches Inf norm of observation error, considering converge\n\n"
                self.misfitX_L1 = np.array(self.misfitX_L1)
                self.misfitX_L2 = np.array(self.misfitX_L2)
                self.misfitX_inf = np.array(self.misfitX_inf)
                self.obsSigma = np.array(self.obsSigma)
            if ii > 0:
                iterativeResidual = (abs(self.misfitX_L2[ii] - self.misfitX_L2[ii-1])) /abs(self.misfitX_L2[0])    
                print "relative Iterative residual = ", iterativeResidual
                print "relative Convergence criterion = ", self.convergenceResi 
                if iterativeResidual  < self.convergenceResi:    
                    if self.reachmaxiteration:
                        pass
                    else:
                        print "Iteration converges, finished \n\n"
                        #break

                np.savetxt('./misfit_L1.txt', self.misfitX_L1)
                np.savetxt('./misfit_L2.txt', self.misfitX_L2)
                np.savetxt('./misfit_inf.txt', self.misfitX_inf)
                np.savetxt('./obsSigma.txt', self.obsSigma)
                np.savetxt('./obsSigmaInf.txt', self.obsSigmaInf)

            ii = ii + 1
        
        # Save misfits
        self.misfitX_L1 = np.array(self.misfitX_L1)
        self.misfitX_L2 = np.array(self.misfitX_L2)
        self.misfitX_inf = np.array(self.misfitX_inf)
        self.obsSigma = np.array(self.obsSigma)
        self.sigmaHX = np.array(self.sigmaHX)

        np.savetxt('./misfit_L1.txt', self.misfitX_L1)
        np.savetxt('./misfit_L2.txt', self.misfitX_L2)
        np.savetxt('./misfit_inf.txt', self.misfitX_inf)
        np.savetxt('./obsSigma.txt', self.obsSigma)
        np.savetxt('./sigmaHX.txt', self.sigmaHX)
        # os.system('plotIterationConvergence.py')            
    
    def report(self):
        """ Report summary information at each step
        """
        raise NotImplementedError

    def clean(self):
        """ Call the dynamic model to do any necessary cleanup before exiting.
        """
        self.dynModel.clean()
    
    ############################ Priviate function ################################    
    def _correctForecasts(self, X, HX, nextEndTime):
        """ Filtering step: Correct the propagated ensemble X 
            via EnKF-MDA filtering procedure
            
            Arg:
            X: state ensemble matrix (Nstate by Ns)
            HX: state ensemble matrix projected to observed space
            nextEndTime: next DA time spot
            
            Return:
            X: state ensemble matrix
        """
        Nmda = self.Tend/self.DAInterval
        DAstep = (nextEndTime - self.DAInterval) / self.DAInterval
        # get XMean and tile XMeanVec to a full matrix (Nsate by Ns)
        XMeanVec = np.mean(X, axis=1) # vector mean of X
        XMeanVec = np.array([XMeanVec])
        XMeanVec = XMeanVec.T
        XMean = np.tile(XMeanVec, (1, self.Ns))
        
        # get coveriance matrix P
        XP = X - XMean
        
        # calculate Kalman Gain matrix and get observation
        (Obs,pertObs, KalmanGainMatrix) = self.dynModel.getBackgroundVarsMDA(Nmda, HX, XP, nextEndTime)

        ObsInf=Obs+np.sqrt(Nmda)*pertObs           
        # analyze step
        dX = np.dot(KalmanGainMatrix, ObsInf - HX)
        X = X + dX
        if (self._iDebug):
            np.savetxt(self._debugFolderName+'kalmanGain_'+str(DAstep), KalmanGainMatrix  )
            np.savetxt(self._debugFolderName+'updateX_'+str(DAstep), X) 
            np.savetxt(self._debugFolderName+'pertObs_'+str(DAstep), pertObs)
            np.savetxt(self._debugFolderName+'Obs_'+str(DAstep), Obs)

        # calculate misfits
        Nnorm = self.dynModel.NstateObs * self.Ns
        diff = abs(Obs - HX)
        misFit_L1 = np.sum(diff) / Nnorm
        misFit_L2 = np.sqrt(np.sum(diff**2)) / Nnorm
        misFit_Inf = LA.norm(diff, np.inf)                
        misFitMax = np.max([misFit_L1, misFit_L2, misFit_Inf])
        sigmaHX = np.std(HX, axis=1)
        sigmaHXNorm = LA.norm(sigmaHX) / self.dynModel.NstateObs
        
        self.misfitX_L1.append(misFit_L1)
        self.misfitX_L2.append(misFit_L2)
        self.misfitX_inf.append(misFit_Inf)
        self.sigmaHX.append(sigmaHXNorm)
        
        print "\nThe misfit between predicted QoI with observed QoI is:"
        print "L1 norm = ", misFit_L1,", \nL2 norm = ", misFit_L2, "\nInf norm = ", misFit_Inf, "\n\n"

        return X, sigmaHXNorm, misFitMax
