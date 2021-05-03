#!/aoe/bin/python27
import sys, os, os.path
import matplotlib.pyplot as plt
import numpy as np
import math
import fnmatch
import hillShape
import pdb
from utilities import readInputData
import foamFileOperation as fp

def plotDomain():
    # Plot the simulation domain
    y=np.arange(0, 9, 0.01)
    yext = np.array([9, 9, 0, 0])

    h=hillShape.profile(y)
    hext = np.array([1, 3.036, 3.036, 1])
    y = np.append(y, yext)
    h = np.append(h, hext)

    plt.plot(y,h,'k-',lw=3)
    #plt.axis([-1, 10, -1, 3.5])

def plotPositions():
    Coord = fp.readVelocityFromFile('pehill/constant/obsLocations')
    Coord = Coord.reshape(-1,3)
    plt.plot(Coord[:,0],Coord[:,1], 'kx',markersize=8,markeredgewidth=1)

def main(iShow=False):

    plt.figure();
    ax1=plt.subplot(111)
    plotDomain()
    plotPositions()
    plt.axis([-0.5, 10, 0, 3.05])
    ax1.set_aspect(1.5)
    rect=plt.Rectangle((2, 0), 8, 1.2, facecolor="#aaaaaa", alpha=0.5, edgecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(8, 3.7)
    #lg = plt.legend(loc = 7)
    #lg.draw_frame(False)
    plt.savefig("Positions"".pdf")

    if iShow:
        plt.show()

if __name__ == "__main__":
   main(True)
