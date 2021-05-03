#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from pylab import *
import pandas as pd
from scipy.stats import gaussian_kde
from utilities import readInputData
from scipy.interpolate import UnivariateSpline

# pdf plot
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
##fig = plt.figure()
#df = pd.DataFrame(scatSample[:,1],columns=['x'])
#df['x'].hist(ax=ax1, bins=10, color = 'grey')
#df['x'].plot(ax=ax2, kind='kde', linewidth=3, color='#d95f02')
#plt.xlabel('reattachment point')
##plt.xlim([0,6])

def getBound(x,CDF,confidence):
    for iter in range(0,len(CDF)-1):
        if (CDF[iter] - 1 + confidence) * (CDF[iter+1] -1 + confidence) <= 0:
            xLow = x[iter+1]
            print 'at ', xLow, ' the probability is ', CDF[iter+1] 
        elif (CDF[iter] - confidence) * (CDF[iter+1] - confidence) <= 0:
            xUp = x[iter]
            print 'at ', xUp, ' the probability is ', CDF[iter] 
    return [xLow,xUp]

def confRange(x,density, confidence=0.975):
    dx = x[1]-x[0]
    density /= (dx*density).sum()
    CDF = np.cumsum(density*dx)
    #plt.plot(x,CDF)
    #plt.show()
    xBound = getBound(x,CDF,confidence)
    return xBound

def intervalPlot(x,density,cl,alpha,ha):
    xb = confRange(x,density)
    min_index = np.where(x==xb[0])[0][0]
    max_index = np.where(x==xb[1])[0][0]
    plt.fill_between(x[min_index:max_index],0,density[min_index:max_index],lw=0,facecolor=cl,color='black',hatch=ha,alpha=alpha)
    plt.fill_between(x[min_index:max_index],0,density[min_index:max_index],facecolor='None',lw=0,hatch=ha)
    plt.plot((x[min_index],x[min_index]),(0,density[min_index]),color=cl,lw=2)
    plt.plot((x[max_index],x[max_index]),(0,density[max_index]),color=cl,lw=2)

xprior = np.loadtxt('reattch_3000')
xpost = np.loadtxt('reattch_300000')
fig = plt.figure()
densityPrior = gaussian_kde(xprior)
densityPost = gaussian_kde(xpost)
x1 = np.linspace(0,6,500)
x2 = np.linspace(2,6,500)
p1, = plt.plot(x1,densityPrior(x1),'--', color='darkred', dashes=(12,3),lw=3)
p2, = plt.plot(x2,densityPost(x2),color='blue', dashes = (9,2,2,2),lw=3)
p3, = plt.plot([4.615,4.615],[0,2],'k-',lw=3)
intervalPlot(x1,densityPrior(x1),'darkred',0.2,'\\\\')
intervalPlot(x2,densityPost(x2),'darkblue',0.5,'//')
#gca().add_patch(Rectangle((0.15,0.71),0.5,0.03,fc='darkred',alpha=0.2,lw=1,hatch='\\\\'))
#gca().add_patch(Rectangle((0.15,0.71),0.5,0.03,fc='None',lw=0,hatch='\\\\'))
#text(0.7,0.71,'95% probability(Prior)',fontsize=18)
#gca().add_patch(Rectangle((0.15,0.62),0.5,0.03,fc='darkblue',alpha=0.5,lw=1,hatch='//'))
#gca().add_patch(Rectangle((0.15,0.62),0.5,0.03,fc='None',lw=0,hatch='//'))
#text(0.7,0.57,'95% probability\n(Posterior)',fontsize=18)
plt.annotate("",xytext=(2.9,0.67),xy=(4,0.55),arrowprops=dict(arrowstyle='-|>'))
#text(1.6,0.65,'95% probability\n(posterior)',fontsize=16)
text(2.1,0.7,'95% probability\n(posterior)',fontsize=16)
plt.annotate("",xytext=(0.9,0.47),xy=(2,0.35),arrowprops=dict(arrowstyle='-|>'))
text(0.1,0.45,'95% probability\n(prior)',fontsize=16)
lg = plt.legend([p1,p2,p3],["Prior","Posterior","DNS (Breuer et al. 2009)"],loc = 0)
lg.draw_frame(False)
plt.ylim([0,2])
plt.xlabel(r'$x_{attach}/H$',fontsize=16)
plt.ylabel('Probability density',fontsize=16)
matplotlib.rcParams.update({'font.size':15})
plt.setp(gca().get_legend().get_texts(), fontsize = '16')
fig.savefig("pehill-pdf.pdf")
plt.show()
