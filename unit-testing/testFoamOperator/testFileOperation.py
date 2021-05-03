# Unit testing for File operator
####################################################
# import standard modules
import unittest
import pdb
import numpy as np
import sys

# import own modules
sys.path.append("../../src/python/")
import foamFileOperation as foamOp
global unitTest 
unitTest = True;

tauFile = "./testCase/Tau"
tau = foamOp.readTurbStressFromFile(tauFile)
   
class TestReadTurbStressFromFile(unittest.TestCase):
    """ Unit test the readTurbStressFromFile function
    """
    def test_extracTensor(self):
        tauBench = "(0.000134029 6.52603e-06 0 0.000222434 0 0.000171712)\n\
(0.000158926 -2.48196e-05 0 0.00020589 0 0.000180039)\n\
(0.000210832 -1.80677e-05 0 0.000188057 0 0.000201275)\n\
(0.00026382 1.8599e-05 0 0.000162855 0 0.000216637)\n"
        
        tauFile = "./testCase/Tau"
        resMid = foamOp.extractTensor(tauFile)
        
        self.assertEqual(tauBench, resMid.group())
    
    def test_extracScalar(self):
        xBench = "(\n0.0025\n0.0075\n0.0125\n0.0175\n)"
        scalarFile = "./testCase/ccx"
        resMid = foamOp.extractScalar(scalarFile)
        
        self.assertEqual(xBench, resMid.group())        
        
    
            
    def test_readTurbStressFromFile(self):
        tauBench = np.array([[  1.34029000e-04,   6.52603000e-06,   0.00000000e+00,
                                2.22434000e-04,   0.00000000e+00,   1.71712000e-04],
                             [  1.58926000e-04,  -2.48196000e-05,   0.00000000e+00,
                                2.05890000e-04,   0.00000000e+00,   1.80039000e-04],
                             [  2.10832000e-04,  -1.80677000e-05,   0.00000000e+00,
                                1.88057000e-04,   0.00000000e+00,   2.01275000e-04],
                             [  2.63820000e-04,   1.85990000e-05,   0.00000000e+00,
                                1.62855000e-04,   0.00000000e+00,   2.16637000e-04]])
        #pdb.set_trace()
        tauFile = "./testCase/Tau"
        tau = foamOp.readTurbStressFromFile(tauFile)
        
        self.assertTrue(np.allclose(tauBench, tau, rtol=1e-06, atol=1e-08))
    
    def test_readTurbCoordinateFromFile(self):
        coordinateBench = np.array([[ 0.0025,  0.0575,  0.005 ],
                                    [ 0.0075,  0.0875,  0.005 ],
                                    [ 0.0125,  0.0475,  0.005 ],
                                    [ 0.0175,  0.0375,  0.005 ]])
        fileDir = "./testCase/"
        coordinate = foamOp.readTurbCoordinateFromFile(fileDir)
        
        self.assertTrue(np.allclose(coordinateBench, coordinate, rtol=1e-06, atol=1e-08))
        
    def test_extractTensorPattern(self):
        tauFile = "./testCase/Tau"
        (resStart, resEnd) = foamOp.extractTensorPattern(tauFile)
    
    def test_writeTurbStressToFile(self):
        tauUpdate = np.array([[ 1, 1, 1, 2, 2, 3],
                              [ 3, 2, 1, 1, 2, 3],
                              [ 5, 5, 5, 3, 4, 3],
                              [ 4, 4, 4, 0, 3, 3]])  
        tauFile = "./testCase/Tau"
        foamOp.writeTurbStressToFile(tauUpdate, tauFile)  

            

                

if __name__ == '__main__':
    unittest.main()