import numpy as np
from pyraf import iraf
import os
from astropy.io import fits
from sys import exit

def sextract(image, output):
   blah = True

def getPhotCal(image_g, image_i):
   blah = True

def makeBoxes(n=1000):
    xPos = 10900.0 * np.random.random_sample(size=n) + 50.0
    yPos = 10900.0 * np.random.random_sample(size=n) + 50.0
    xWidth = 80.0 * np.random.random_sample(size=n) + 30.0
    yWidth = 80.0 * np.random.random_sample(size=n) + 30.0
    
    with open('bgRegions.txt','w+') as f1:
        for i in range(n):
            print >> f1, "logical; box", xPos[i], yPos[i], xWidth[i], yWidth[i], 0

def bgCounts(image):
    if not os.path.exists("bgvals_"+image[4:-5]+".txt"):
        bg_file= open("bgvals_"+image[4:-5]+".txt", 'w+')

        b3,b4,b5,b6 = np.loadtxt('bgRegions.txt',usecols=(2,3,4,5),unpack=True)
        for i in range(len(b3)) :
            bx1 = b3[i] - (b5[i]/2.)
            bx2 = b3[i] + (b5[i]/2.)
            by1 = b4[i] - (b6[i]/2.)
            by2 = b4[i] + (b6[i]/2.)

            iraf.images.imstat(image[0:-5]+'['+repr(int(bx1))+':'+repr(int(bx2))+','+repr(int(by1))+':'+repr(int(by2))+']', fields="image,npix,mean,midpt,stddev,min,max", Stdout=bg_file)
        bg_file.close()

    bgmean, bgsig = np.loadtxt("bgvals_"+image[4:-5]+".txt",usecols=(3,4),unpack=True)

    bgSig = np.std(bgmean)
    bgMean = np.mean(bgmean)
    return bgMean, bgSig
