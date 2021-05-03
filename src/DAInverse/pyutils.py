# Utility 
# Copyright: Jianxun Wang (vtwjx@vt.edu)
# Mar.27, 2015
import numpy as np

def readfile(filename):
    d0  = []
    nlines = 0
    for line in file(filename):
        line = line.rstrip('\n')
        line_list = [float(x) for x in line.split()]
        d0.append(line_list)
        nlines = nlines+1
    return d0,nlines

def column(matrix, i):
    return [row[i] for row in matrix]

def checkPtInside(xc,yc,x,y):
    nc=len(xc);
    nm=len(x);
    umask=np.array(np.zeros(nm))
    for i in range(nm):
        b=y[i]-x[i];
        nint = 0;
        for j in range(nc-1):
            if abs(xc[j]-xc[j+1])<1.0e-30:
                yi=xc[j]+b;
                if (xc[j]<=x[i]) & ((yi-yc[j])*(yi-yc[j+1])<=0.0):
                    nint = nint+1;
            else:
                a1=(yc[j+1]-yc[j]) / (xc[j+1]-xc[j]);
                b1=yc[j]-a1*xc[j];
                xi=(b1-b)/(1.0-a1);
                if (xi<=x[i]) & ((xi-xc[j])*(xi-xc[j+1])<=0.0):
                    nint = nint+1;
        if ( nint%2 == 0 ):
            umask[i]=1;
    return umask;

