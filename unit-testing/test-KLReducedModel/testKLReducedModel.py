# Unit testing for KLReducedModel.py
####################################################
# import standard modules
import unittest
import pdb
import numpy as np
import os
import sys

# import own modules
sys.path.append("../../src/python/")
import KLReducedModel as kl

global unitTest
unitTest = True
   
class TestKLReducedModel(unittest.TestCase):
    """ Unit test the testKLReducedModel function
    """    
    
    def test_generateMesh(self):
        
        self.coor1D = np.array([[0, 0.1, 0.4, 0.8, 1.0]]) 
        self.coor1D = self.coor1D.T     
        
        self.coor2D = np.array([[0, 0.1, 0.4, 0.8, 1.0],
                           [0, 2.0, 4.0, 8.0, 10.0]])
        self.coor2D = self.coor2D.T 
                
        self.coor3D = np.array([[0, 0.1, 0.4, 0.8, 1.0],
                           [0, 2.0, 4.0, 8.0, 10],
                           [1, 1.1, 1.4, 1.8, 2.0]])
        self.coor3D = self.coor3D.T                            
        
        self.kernelType = 'SqExp'
        self.hyperPara = np.array([5, 0.05]) # sigma, len
        self.m = 5    
        
        self.KL1D = kl.KLReducedModel(self.coor1D, self.kernelType,
                                      self.hyperPara, self.m)        
        self.KL2D = kl.KLReducedModel(self.coor2D, self.kernelType, 
                                      self.hyperPara, self.m)    
        self.KL3D = kl.KLReducedModel(self.coor3D, self.kernelType, 
                                      self.hyperPara, self.m)
        
        # generate meshgrid files
        self.KL1D._generateMesh()
        self.KL2D._generateMesh()
        self.KL3D._generateMesh()
     
    def test_calKLModes1D(self):
        # parameters specification
        m = 64
        kernelType = "SqExp"
        hyperPara = np.array([5, 0.05])
        
        # generate 1D meshing
        coor1D = np.linspace(0, 1, 128) 
        coor1D = np.reshape(coor1D, (128, 1))       
        
        # Initial KL class and call the solver
        KL1D = kl.KLReducedModel(coor1D, kernelType, hyperPara, m)
        cov, eig, KLModes = KL1D.calKLModes()
        
        # Plotting
        KL1D.plotCov(cov)
        os.system('mv *.eps klExpansionData1D')
    
    def test_calKLModes2D(self):
        # parameters specification
        m = 64
        kernelType = "SqExp"                        
        hyperPara = np.array([5, 0.5])
        
        # generate 2D meshing and cell areas
        coor2D = np.loadtxt('./dataForTesting/cellCenter.dat')
        area2D = np.loadtxt('./dataForTesting/cellArea.dat')
        
        # Initial KL class and call the solver
        KL2D = kl.KLReducedModel(coor2D, kernelType, hyperPara, m, area2D, area2D)
        cov, eig, KLModes = KL2D.calKLModes()
        
        # Plotting
        KL2D.plotCov(cov)
        os.system('mv *.eps klExpansionData2D')

    def test_calKLModes3D(self):
        # parameters specification
        m = 64
        kernelType = "SqExp"                        
        hyperPara = np.array([5, 0.5, 0.5, 0.5])
        
        # generate 3D meshing and cell areas
        coor3D = np.loadtxt('./dataForTesting/cellCenter3D.dat')
        area3D = np.loadtxt('./dataForTesting/cellArea.dat')
        
        # Initial KL class and call the solver
        KL3D = kl.KLReducedModel(coor3D, kernelType, hyperPara, m, area3D, area3D)
        cov, eig, KLModes = KL3D.calKLModes()
        
        # Plotting
        KL3D.plotCov(cov)
        os.system('mv *.eps klExpansionData3D')      

if __name__ == '__main__':
    unittest.main()
