import numpy as np
from pyraf import iraf
import os
from astropy.io import fits
from sys import exit

def getImageSum(images, output):
    data = {}
    for image in images:
        hdulist = fits.open(image)
        scidata = hdulist[0].data
        data[image] = scidata
        
    check = data.values()[0] + data.values()[1] + data.values()[2] + data.values()[3]
    # print check[3795][1213]
    hdu = fits.PrimaryHDU(check)
    hdu.writeto(output)    
        
def getDivImage(images, output='check.fits'):
    gaps = {}
    for image in images:
        hdulist = fits.open(image)
        scidata = hdulist[0].data
        gaps[image] = (scidata != 0.0).astype(float)
        print (gaps[image])[3795][1213]
        hdulist.close()
    check = gaps.values()[0] + gaps.values()[1] + gaps.values()[2] + gaps.values()[3]
    print check[3795][1213]
    hdu = fits.PrimaryHDU(check)
    hdu.writeto(output)

def bgCounts(image):
    if image.startswith('sh') or image.startswith('subsh'):
        if not os.path.exists("bgvals_"+image[0:-5]+".txt"):
            bg_file= open("bgvals_"+image[0:-5]+".txt", 'w+')

            b3,b4,b5,b6 = np.loadtxt('bgRegions.txt',usecols=(2,3,4,5),unpack=True)
            for i in range(len(b3)) :
                bx1 = b3[i] - (b5[i]/2.)
                bx2 = b3[i] + (b5[i]/2.)
                by1 = b4[i] - (b6[i]/2.)
                by2 = b4[i] + (b6[i]/2.)

                iraf.images.imstat(image[0:-5]+'['+repr(int(bx1))+':'+repr(int(bx2))+','+repr(int(by1))+':'+repr(int(by2))+']', fields="image,npix,mean,midpt,stddev,min,max", Stdout=bg_file)
            bg_file.close()

        bgmean, bgsig = np.loadtxt("bgvals_"+image[0:-5]+".txt",usecols=(3,4),unpack=True)
    else:
        if not os.path.exists("bgvals_"+image[4:-5]+".txt"):
            bg_file= open("bgvals_"+image[4:-5]+".txt", 'w+')

            b3,b4,b5,b6 = np.loadtxt('bgr_'+image[4:-5]+'.txt',usecols=(2,3,4,5),unpack=True)
            for i in range(len(b3)) :
                bx1 = b3[i] - (b5[i]/2.)
                bx2 = b3[i] + (b5[i]/2.)
                by1 = b4[i] - (b6[i]/2.)
                by2 = b4[i] + (b6[i]/2.)

                iraf.images.imstat(image[0:-5]+'['+repr(int(bx1))+':'+repr(int(bx2))+','+repr(int(by1))+':'+repr(int(by2))+']', fields="image,npix,mean,midpt,stddev,min,max", Stdout=bg_file)
            bg_file.close()

        bgmean, bgsig = np.loadtxt("bgvals_"+image[4:-5]+".txt",usecols=(3,4),unpack=True)

    bgSig = np.mean(bgsig)
    bgMean = np.mean(bgmean)
    return bgMean, bgSig

##############################################################################

def msc_getscale(images, sigmas, refimage=False, rootname='', prefix='sc', applyscale=False, verbose=False):
    # if you don't specify a reference image, use the first one in the list by default
    bgSigmaDict = sigmas
    if not refimage:
        ref = images[0]
    else :
        ref = refimage
        
    try : 
        expDict = {}
        for i,image in enumerate(images):
            print "Processing",image+"..."
            scimage = prefix + image
            
            # delete a previously scaled image if it exists
            if applyscale:
                if os.path.isfile(scimage):
                    iraf.imutil.imdelete(scimage, go_ahead='yes', verify='no')
                    
            # get the exposure time from the image header and add it to the list
            hdulist = fits.open(image)
            filterName = hdulist[0].header['filter']
            expDict[image] = hdulist[0].header['exptime']
            hdulist.close()
            
            coords = 'photRefStars'+filterName[0]+'.sort.txt'
            iraf.datapars.setParam('sigma',bgSigmaDict[image])
            iraf.apphot.phot('subsh'+image, coords=coords, interactive='no', output='subsh'+image[0:-5]+'.mag.1')
            
            with open('subsh'+image[0:-5]+'.txdump', 'w+') as tdFile:
                iraf.ptools.txdump(textfiles='subsh'+image[0:-5]+'.mag.1', fields="XCENTER, YCENTER, FLUX, MAG", expr='yes', Stdout=tdFile)
            os.remove('subsh'+image[0:-5]+'.mag.1')
    
        # print expDict[ref]    
        # read in the photometry
        # use a dictionary to keep track in a dynamic fashion
        magDict = {}
        fluxDict = {}
        xDict = {}
        yDict = {}
        for i,image in enumerate(images):
            xpos, ypos, fluxes, mags = np.loadtxt('subsh'+image[0:-5]+'.txdump',usecols=(0,1,2,3) ,dtype=str,unpack=True)
            magsFloat = []
            flux = []
            for j,mag in enumerate(mags):
                if mag == 'INDEF':
                    magsFloat.append(99.99)
                    flux.append(float(fluxes[j]))
                else:
                    magsFloat.append(float(mag))
                    flux.append(float(fluxes[j]))
                
            key = image
            value = np.array(magsFloat)
            
            magDict[key] = value
            fluxDict[key] = np.array(flux)
            xDict[key] = xpos
            yDict[key] = ypos
                
        scale = []
        std = []
        scaleDict = {}
        stdDict = {}
        sigThreshold = 0.005
        for j,image in enumerate(images):
            print 'Calculating scaling factor for',image
            n = 1
            expRatio = (float(expDict[ref])/float(expDict[image]))
            magTempA = magDict[image]
            magTempRef = magDict[ref]
            xTempA = xDict[image]
            yTempA = yDict[image]
            magA = magTempA[np.where(magTempA<99.99)]
            magRef = magTempRef[np.where(magTempA<99.99)]
            xA = xTempA[np.where(magTempA<99.99)]
            yA = yTempA[np.where(magTempA<99.99)]
            
            rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
            if verbose:
                for i,r in enumerate(rat):
                    print image, i, r
                    
            sigTest = np.std(rat)
            print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
            if sigTest <= sigThreshold:
                scale.append(np.mean(rat))
                scaleDict[image] = np.mean(rat)
                std.append(np.std(rat))
                with open(image+'.scalepos','w+') as posfile:
                    for m in range(len(xA)):
                        print >> posfile, xA[m], yA[m]
            else:
                while sigTest > sigThreshold:
                    magTempA = magA
                    magTempRef = magRef
                    xTempA = xA
                    yTempA = yA
                    magA = magTempA[np.where(abs(rat-np.median(rat))<sigTest)]
                    xA = xTempA[np.where(abs(rat-np.median(rat))<sigTest)]
                    yA = yTempA[np.where(abs(rat-np.median(rat))<sigTest)]
                    magRef = magTempRef[np.where(abs(rat-np.median(rat))<sigTest)]
                    rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
                    if verbose:
                        for i,r in enumerate(rat):
                            print image, i, r
                    sigTest = np.std(rat)
                    
                    n = n + 1
                    if n > 10:
                        print "Iteration did not converge to sigma <", repr(sigThreshold),"for", image
                        print "Quitting..."
                        exit()
                    print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
                with open(image+'.scalepos','w+') as posfile:
                    for m in range(len(xA)):
                        print >> posfile, xA[m], yA[m]
                scale.append(np.mean(rat))
                scaleDict[image] = np.mean(rat)
                std.append(np.std(rat))
        if applyscale:
            for image in magDict:
                iraf.imutil.imarith('subsh'+image,'/',scaleDict[image],prefix+'subsh'+image, verbose='yes')
                
        return scale, std
        
    except IOError:
        print 'Make sure you have an image and a list of ref star coordinates!'
        print 'Quitting...'
        exit()