# Unit testing plotting for KLReducedModel.py
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
np.random.seed(10)
class TestKLReducedModelNonRun(unittest.TestCase):
    """ Unit test the plot function in testKLReducedModel
    """ 
    def test_recomputeFieldReduced(self):
        # parameters specification
        m = 64
        kernelType = "SqExp"
        hyperPara = np.array([5, 0.5])                

        # generate 2D meshing and cell areas
        coor2D = np.loadtxt('./dataForTesting/cellCenter.dat')
        area2D = np.loadtxt('./dataForTesting/cellArea.dat')
        
        # initialize KLReducedModel class
        KL2D = kl.KLReducedModel(coor2D, kernelType, hyperPara, m)        
        
        fileDir = './klExpansionData2D/KLmodes.dat'
        recField, coord = KL2D.recomputeFieldReduced(0, fileDir)
        np.savetxt('reconstructFieldReduced.dat', recField)  
        KL2D.plotField(recField, coord, 20)
        os.system('mv field.eps fieldReduced.eps')
        
    def test_recomputeFieldFull(self):
        m = 64
        kernelType = "SqExp"
        hyperPara = np.array([5, 0.5])                

        # generate 2D meshing and cell areas
        coor2D = np.loadtxt('./dataForTesting/cellCenter.dat')
        area2D = np.loadtxt('./dataForTesting/cellArea.dat')
        
        # initialize KLReducedModel class
        KL2D = kl.KLReducedModel(coor2D, kernelType, hyperPara, m)        
        
        fileDir = './klExpansionData2D/cov.dat'
        recField = KL2D.recomputeFieldFull(0, fileDir)
        np.savetxt('reconstructFieldFull.dat', recField)  
        KL2D.plotField(recField, coor2D, 20)
        
        
    def test_plotKLModes1D(self):
        # parameters specification
        m = 64
        kernelType = "SqExp"
        hyperPara = np.array([5, 0.05])
        
        # generate 1D meshing
        coor1D = np.linspace(0, 1, 128) 
        coor1D = np.reshape(coor1D, (128, 1))       
        
        # Initial KL class and call the solver
        KL1D = kl.KLReducedModel(coor1D, kernelType, hyperPara, m)
        
        fileDir = "./klExpansionData1D/KLmodes.dat"
        nkl = 4
        KL1D.plotKLModes1D(nkl, fileDir)
            

if __name__ == '__main__':
    unittest.main()        
