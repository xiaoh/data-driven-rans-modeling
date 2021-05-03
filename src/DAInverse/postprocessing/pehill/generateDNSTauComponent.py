#!/aoe/bin/python27

# description        :Generate DNS Tau Component

# author             :Jianxun Wang (vtwjx@vt.edu)
# copyright          :Heng Xiao's Group
# date               :Feb.24, 2016
# revision           :Feb.24, 2016
###############################################################################

## Import system modules
# sci computing
import numpy as np
import scipy.stats as ss
# system, file operation
import pdb
import sys
import ast
# plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib as mp
import pylab
from scipy import interpolate
# Import local modules
from utilities import readInputData, extractListFromDict
from ReynoldsStressRF import ReynoldsStressRF


def readDNS():
    """
    get DNS component (c1, c2, c3, XC1, XC2, Xi, Eta, k, VA, VB, VC)
    """
    xArray = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
    for x in xArray:
        tau_vline = np.loadtxt('sets/X/ransTauFine' + x)
        y_vline = tau_vline[:, 0]
        tau_dns_vline = tau_vline[:, 1:7]
                
        tau_rans_vline = np.loadtxt('sets/X/line_x' + x +'_Tau.xy')
        y_vline_rans = tau_rans_vline[:, 0]       
        tau_rans_vline = tau_rans_vline[:, 1:]
                
        mapTau = ReynoldsStressRF('None', tau_rans_vline, 400, 1, 'True')        
        k,V1,V2,V3,C,NP = mapTau._tau2PhysParams(tau_dns_vline)
        X = mapTau._C2X(C)
        RS = mapTau._phys2Natural(X) # collapse time = 0.02s (3000 cells)
        VA, VB, VC = mapTau.getThetaVABC(tau_dns_vline) # collapse time = 1.005s (3000 cells)
            
        np.savetxt(dns_folder + 'C_dns_x_' + x, C)
        np.savetxt(dns_folder + 'XC_dns_x_' + x, X)
        np.savetxt(dns_folder + 'XiEta_dns_x' + x, RS)
        np.savetxt(dns_folder + 'k_dns_x' + x, k)
        np.savetxt(dns_folder + 'VA_dns_x' + x, VA)
        np.savetxt(dns_folder + 'VB_dns_x' + x, VB)
        np.savetxt(dns_folder + 'VC_dns_x' + x, VC)
    
if __name__ == '__main__':
    dns_folder = './sets/X/dns_decompose/'
    readDNS()    
