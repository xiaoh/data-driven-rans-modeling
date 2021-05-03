#!/usr/bin/python
# Unit testing for KLReducedModel.py
####################################################
# import standard modules
import unittest
import pdb
import numpy as np
import sys
import matplotlib.pylab as plt


# import own modules
sys.path.append("../../src/python/")
import ReynoldsStressRF as RS
from foamFileOperation import *

global unitTest
unitTest = True

class TestReynoldsStressRF(unittest.TestCase):
    """ Unit test the ReynoldsStressRF function
    """

    def setUp(self):
        self.R0 = RS.ReynoldsStressRF(os.getcwd(),"/Tau.org",3,1);

    def testC2X(self):
        X = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(X[0,0],0.50425592)
        self.assertAlmostEqual(X[0,1],0.00065027)
        self.assertAlmostEqual(X[1,0],0.53050669)
        self.assertAlmostEqual(X[1,1],0.00124285)
        self.assertAlmostEqual(X[2,0],0.5594895)
        self.assertAlmostEqual(X[2,1],0.00209133)

    def testX2C(self):
        X = self.R0._C2X(self.R0.C)
        C = self.R0._X2C(X)

        self.assertAlmostEqual(self.R0.C[0,0],C[0,0])
        self.assertAlmostEqual(self.R0.C[0,1],C[0,1])
        self.assertAlmostEqual(self.R0.C[1,0],C[1,0])
        self.assertAlmostEqual(self.R0.C[1,1],C[1,1])
        self.assertAlmostEqual(self.R0.C[2,0],C[2,0])
        self.assertAlmostEqual(self.R0.C[2,1],C[2,1])

    def testTau2PhysParam(self):
        tau = self.R0.tau
        self.R0._tau2PhysParams(tau)
        self.assertAlmostEqual(self.R0.C[0,0],self.R0.COrg[0,0])
        self.assertAlmostEqual(self.R0.C[0,1],self.R0.COrg[0,1])

        self.assertAlmostEqual(self.R0.k[0],self.R0.kOrg[0])
        self.assertAlmostEqual(self.R0.k[1],self.R0.kOrg[1])

        self.assertAlmostEqual(self.R0.V1[0,0],self.R0.V1Org[0,0])
        self.assertAlmostEqual(self.R0.V1[0,1],self.R0.V1Org[0,1])
        self.assertAlmostEqual(self.R0.V2[0,0],self.R0.V2Org[0,0])
        self.assertAlmostEqual(self.R0.V2[0,1],self.R0.V2Org[0,1])
        self.assertAlmostEqual(self.R0.V3[0,0],self.R0.V3Org[0,0])
        self.assertAlmostEqual(self.R0.V3[0,1],self.R0.V3Org[0,1])

    def testC2Tau(self):
        Cs = self.R0.C
        tau = self.R0._C2Tau(Cs)
        self.assertAlmostEqual(self.R0.tau[0,0],tau[0,0])
        self.assertAlmostEqual(self.R0.tau[1,1],tau[1,1])
        self.assertAlmostEqual(self.R0.tau[2,2],tau[2,2])
        self.assertAlmostEqual(self.R0.tau[0,3],tau[0,3])
        self.assertAlmostEqual(self.R0.tau[1,4],tau[1,4])
        self.assertAlmostEqual(self.R0.tau[2,5],tau[2,5])

    def testPhys2Natural(self):
        X = self.R0._C2X(self.R0.C)
        RSs = self.R0._phys2Natural(X)

        self.assertAlmostEqual(RSs[0,0],0.00851824)
        self.assertAlmostEqual(RSs[0,1],-0.99849828)
        self.assertAlmostEqual(RSs[1,0],0.06110107)
        self.assertAlmostEqual(RSs[1,1],-0.99712975)
        self.assertAlmostEqual(RSs[2,0],0.11926701)
        self.assertAlmostEqual(RSs[2,1],-0.99517027)

        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(RSs[:,0],RSs[:,1],'k+')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([-2,2])
        plt.ylim([-2,2])
        plt.savefig("a.pdf")

    def testNatural2PhysP0(self):
        X = self.R0._C2X(self.R0.C)
        RSs = self.R0._phys2Natural(X)
        XNew = self.R0._natural2Phys(RSs)
        self.assertAlmostEqual(X[0,0],XNew[0,0])
        self.assertAlmostEqual(X[1,1],XNew[1,1])
        self.assertAlmostEqual(X[2,0],XNew[2,0])

    def testNatural2PhysP01(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 0.1
            deltaEta[i,0] = 0.1

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(XNew[0,0],0.55150542)
        self.assertAlmostEqual(XNew[1,1],0.04454412)
        self.assertAlmostEqual(XNew[2,0],0.60388708)

        self.assertAlmostEqual(tau[0,0],8.47676111e-02)
        self.assertAlmostEqual(tau[1,1],-1.19118022e-03)
        self.assertAlmostEqual(tau[2,2],1.50102413e-04)

        plt.clf()
        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(XNew[:,0],XNew[:,1],'k+')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.savefig("b.pdf")

    def testNatural2PhysP10(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 1.0
            deltaEta[i,0] = 1.0

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(XNew[0,0],0.74962457)
        self.assertAlmostEqual(XNew[1,1],0.43425556)
        self.assertAlmostEqual(XNew[2,0],0.74879257)

        plt.clf()
        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(XNew[:,0],XNew[:,1],'k+')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.savefig("c.pdf")

    def testNatural2PhysP_10(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = -1.0
            deltaEta[i,0] = -1.0

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(XNew[0,0],0.00425912)
        self.assertAlmostEqual(XNew[1,1],0.0)
        self.assertAlmostEqual(XNew[2,0],0.0596335)

        plt.clf()
        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(XNew[:,0],XNew[:,1],'k+')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([-0.1,1.1])
        plt.ylim([-0.1,1])
        plt.savefig("d.pdf")

    def testNatural2PhysP10_10(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 1.0
            deltaEta[i,0] = -1.0

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(XNew[0,0],1.0)
        self.assertAlmostEqual(XNew[1,1],0.0)
        self.assertAlmostEqual(XNew[2,0],1.0)

        plt.clf()
        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(XNew[:,0],XNew[:,1],'k+')
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([-0.1,1.1])
        plt.ylim([-0.1,1])
        plt.savefig("e.pdf")

    def testNatural2PhysP0_10(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 0.0
            deltaEta[i,0] = -1.0

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(XNew[0,1],0.0)
        self.assertAlmostEqual(XNew[1,1],0.0)
        self.assertAlmostEqual(XNew[2,1],0.0)

        plt.clf()
        plt.plot(X[:,0],X[:,1],'ko')
        plt.plot(XNew[:,0],XNew[:,1],'k+',ms = 20)
        plt.plot([0,1,0.5,0.5,0],[0,0,3**0.5/2.0,3**0.5/2.0,0],'b-')
        plt.xlim([-0.1,1.1])
        plt.ylim([-0.1,1])
        plt.savefig("f.pdf")

    def testPerturbK(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        deltaK = np.zeros([self.R0.nMesh,1])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 0.0
            deltaEta[i,0] = -1.0
            deltaK[i,0] = 1.0

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta,deltaK)
        XNew = self.R0._C2X(self.R0.C)

        self.assertAlmostEqual(tau[0,0],1.68978647e-01)

    def testAngles_01(self):
        deltaXi = np.zeros([self.R0.nMesh,1])
        deltaEta = np.zeros([self.R0.nMesh,1])
        deltaK = np.zeros([self.R0.nMesh,1])
        deltaThetaV = np.zeros([self.R0.nMesh,3])
        for i in range(self.R0.nMesh):
            deltaXi[i,0] = 0
            deltaEta[i,0] = 0
            deltaK[i,0] = 0
            deltaThetaV[i,0] = -0.1
            deltaThetaV[i,1] = -0.1
            deltaThetaV[i,2] = -0.1

        tauOld = self.R0.tau
        TAOld, TBOld, TCOld = self.R0.getThetaVABC(tauOld)

        X = self.R0._C2X(self.R0.C)
        tau = self.R0.perturbTau(deltaXi,deltaEta,deltaK,deltaThetaV)
        TA, TB, TC = self.R0.getThetaVABC(tau)

        self.assertAlmostEqual(TAOld[0],1.56605528)
        self.assertAlmostEqual(TBOld[1],1.57024287)
        self.assertAlmostEqual(TCOld[2],1.35408205e-05)
        self.assertAlmostEqual(TA[0],1.46605528)
        self.assertAlmostEqual(TB[1],1.47024287)
        self.assertAlmostEqual(TC[2],3.04160619)


def main():
    suite = unittest.TestLoader().loadTestsFromTestCase(TestReynoldsStressRF)
    unittest.TextTestRunner(verbosity = 2).run(suite)

if __name__ == '__main__':
    main()
