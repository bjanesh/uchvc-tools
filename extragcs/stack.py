from pyraf import iraf
from astropy.io import fits
import numpy as np
import os
import glob
from msc_helpers import msc_getscale, bgCounts, getDivImage

###########
# first set up a few things based on the directory we're in and what files we will need
objectName = os.getcwd()[-5:].lower() # the directory we are in
mscList = 'msc_'+objectName+'.list'   # a list of the single ext. fits images


###########
# check for / add header keywords for setairmass
# ST <- LSTHDR, EPOCH <- 2000
# !!!! the mosaic 1.1 images already have the ST keyword, so just set EPOCH to 2000
# iraf.imutil.hedit(images='mscList', fields='EPOCH', value=2000.0, add='yes', verify='no')

###########
# run setairmass

if not os.path.isfile('setairmass.done'):
    iraf.astutil.setairmass.setParam('images', "msc*fits")          # Input images
    iraf.astutil.setairmass.setParam('intype', "beginning")    # Input keyword time stamp
    iraf.astutil.setairmass.setParam('outtype', "effective")    # Output airmass time stamp\n
    iraf.astutil.setairmass.setParam('ra', "ra")           # Right acsension keyword (hours)
    iraf.astutil.setairmass.setParam('dec', "dec")          # Declination keyword (degrees)
    iraf.astutil.setairmass.setParam('equinox', "radeceq")        # Equinox keyword (years)
    iraf.astutil.setairmass.setParam('st', "st")           # Local siderial time keyword (hours)
    iraf.astutil.setairmass.setParam('ut', "time-obs")     # Universal time keyword (hours)
    iraf.astutil.setairmass.setParam('date', "date-obs")     # Observation date keyword
    iraf.astutil.setairmass.setParam('exposure', "exptime")      # Exposure time keyword (seconds)
    iraf.astutil.setairmass.setParam('airmass', "airmass")      # Airmass keyword (output)
    iraf.astutil.setairmass.setParam('utmiddle', "utmiddle")     # Mid-observation UT keyword (output)
    iraf.astutil.setairmass.setParam('scale', 750.)           # The atmospheric scale height\n
    iraf.astutil.setairmass.setParam('show', 'yes')            # Print the airmasses and mid-UT?
    iraf.astutil.setairmass.setParam('update', 'yes')            # Update the image header?
    iraf.astutil.setairmass.setParam('override', 'yes')            # Override previous assignments?
    iraf.astutil.setairmass()
    with open('setairmass.done','w+') as f1:
        print >> f1, True
else:
    print 'setairmass already done'

###########
# align the images to each other (imtranspose/imrotate/rotate?)
# might need to use getshfts.cl here -- see kathy notes for more info
# finally use imalign to do the aligning

# first make @-lists for the files in each filter
# the reference image should be first

if not os.path.isfile('shift.refimage'):
    refImageSh = raw_input('Enter the shifting reference image name: ')
    with open('shift.refimage', 'w+') as f1:
        print >> f1, refImageSh
else:
    refImageSh = str(np.loadtxt('shift.refimage', usecols=(0,), dtype=str, unpack=True))
    # refImageSh = refImageShin[0]
    print "reference image: ", refImageSh

filesSh = glob.glob('msc_*.fits')
filesSh.pop(filesSh.index(refImageSh))
filesSh.append(refImageSh)
filesSh.reverse()

if not os.path.isfile('mscFilesSh.list'):   
    with open('mscFilesSh.list','w+') as f1:
        for i,mfile in enumerate(filesSh):
            print >> f1, mfile  

# import the getshfts task as a pyraf task
iraf.task(getshfts = "home$scripts/getshfts.cl")

# then run getshfts
if not os.path.isfile('mscFilesSh.shft'):
    iraf.getshfts(images='@mscFilesSh.list', rootname='mscFilesSh')

if not os.path.isfile('sh'+refImageSh):    
    iraf.immatch.imalign.setParam('input','@mscFilesSh.list')
    iraf.immatch.imalign.setParam('reference',refImageSh)
    iraf.immatch.imalign.setParam('coords','mscFilesSh.reg')
    iraf.immatch.imalign.setParam('output','sh//@mscFilesSh.list')
    iraf.immatch.imalign.setParam('shifts','mscFilesSh.shft')
    iraf.immatch.imalign.setParam('niterate',10)
    iraf.immatch.imalign.setParam('interp_type','poly5')
    iraf.immatch.imalign.setParam('boundary_type','constant')
    iraf.immatch.imalign.setParam('constant',0.0)
    iraf.immatch.imalign.setParam('verbose','yes')
    iraf.immatch.imalign.setParam('trimimages','yes')

    iraf.immatch.imalign()

###########
# print out a table of the images for notes
# for each image in the single ext list, read the header and print out relevant info
filterName = []
expTime = []
airmass = []
bgMean = []
bgSigma = []

bgMeanDict = {}
bgSigmaDict = {}

for i,file_ in enumerate(glob.glob('msc*.fits')):
    hdulist = fits.open(file_)
    filterName.append(hdulist[0].header['filter'])
    expTime.append(hdulist[0].header['exptime'])
    airmass.append(hdulist[0].header['airmass'])
    hdulist.close()
    bgMeanx, bgSigmax = bgCounts('sh'+file_)
    bgMean.append(bgMeanx)
    bgSigma.append(bgSigmax)
    bgSigmaDict[file_] = bgSigmax
    bgMeanDict[file_] = bgMeanx
    
###########
# figure out what image should be the reference image for each filter set
# you want it to be the the one with the most flux
# with photometric conditions, this is the image with the lowest airmass
# you might need to do imexam to check this
# !!!! use the background in the images. pick the one with the lowest background counts
# so we need to measure the background in the images

filterName, airmass, expTime, mscFile, bgMean, bgSigma = zip(*sorted(zip(filterName, airmass, expTime, glob.glob('msc*.fits'), bgMean, bgSigma)))

bgMeanB = [bg for i,bg in enumerate(bgMean) if (filterName[i][0]=='B')]
bgMeanV = [bg for i,bg in enumerate(bgMean) if (filterName[i][0]=='V')]
bgMeanR = [bg for i,bg in enumerate(bgMean) if (filterName[i][0]=='R')]

bgSigmaB = [bg for i,bg in enumerate(bgSigma) if (filterName[i][0]=='B')]
bgSigmaV = [bg for i,bg in enumerate(bgSigma) if (filterName[i][0]=='V')]
bgSigmaR = [bg for i,bg in enumerate(bgSigma) if (filterName[i][0]=='R')]

refImages = [np.argmin(bgMeanB),np.argmin(bgMeanR)+len(bgMeanB),np.argmin(bgMeanV)+len(bgMeanB)+len(bgMeanR)]
rIm = ['']*len(filterName)
rIm[refImages[0]],rIm[refImages[1]],rIm[refImages[2]] = '<--','<--','<--'

refImageB = mscFile[refImages[0]]
refImageR = mscFile[refImages[1]]
refImageV = mscFile[refImages[2]]

filesB = [mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='B')]
filesV = [mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='V')]
filesR = [mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='R')]

filesB.pop(np.argmin(bgMeanB))
filesV.pop(np.argmin(bgMeanV))
filesR.pop(np.argmin(bgMeanR))

filesB.append(refImageB)
filesV.append(refImageV)
filesR.append(refImageR)

filesB.reverse()
filesV.reverse()
filesR.reverse()

with open('mscFilesB.list','w+') as f1:
    for i,mfile in enumerate(filesB):
        print >> f1, mfile

with open('mscFilesV.list','w+') as f1:
    for i,mfile in enumerate(filesV):
        print >> f1, mfile

with open('mscFilesR.list','w+') as f1:
    for i,mfile in enumerate(filesR):
        print >> f1, mfile

print '##############################################################'
print '# Image       Airmass   Filter  Exp  "Mean" BG     BG sig  ref'
for i,x in enumerate(filterName):
    print '{0:12s} {1:9.6f}  {2:1s}     {3:5.0f}  {4:9.3f}  {5:9.3f}  {6:3s}'.format(mscFile[i][0:-5], airmass[i], filterName[i][0], expTime[i], bgMean[i], bgSigma[i], rIm[i])

with open('imageData.txt', 'w+') as f1:
    print >> f1, '##############################################################'
    print >> f1, '# Image       Airmass   Filter  Exp  "Mean" BG     BG sig  ref'
    for i,x in enumerate(filterName):
        print >> f1, '{0:12s} {1:9.6f}  {2:1s}     {3:5.0f}  {4:9.3f}  {5:9.3f}  {6:3s}'.format(mscFile[i][0:-5], airmass[i], filterName[i][0], expTime[i], bgMean[i], bgSigma[i], rIm[i])
    
###########
# now do the actual stacking

###########
# subtract background counts
# measure background counts using empty regions
# use imarith to actually subtract the counts
# mean median value should be ZERO

for i,file_ in enumerate(glob.glob('msc*.fits')):
    # bgMeanshx, bgSigmashx = bgCounts(file_)
    # bgMeansh.append(bgMeanshx)
    # bgSigmash.append(bgSigmashx)
    if not os.path.isfile('subsh'+file_):
        iraf.imutil.imarith('sh'+file_,'-',bgMeanDict[file_],'subsh'+file_, verbose='yes')
        bgMeantest, bgSigmatest = bgCounts('subsh'+file_)
        print bgMeantest, bgSigmatest
       
# filesBsh = ['sh'+mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='B')]
# getDivImage(filesBsh, output='checkB.fits')
# filesVsh = ['sh'+mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='V')]
# getDivImage(filesVsh, output='checkV.fits')
# filesRsh = ['sh'+mfile for i,mfile in enumerate(mscFile) if (filterName[i][0]=='R')] 
# getDivImage(filesRsh, output='checkR.fits')

###########
# scale the images to a common flux scale
# need to know phot parameters (typical fwhm, saturation, pick an aperture/annulus size)
# use getscale script to find the scale factors
# make sure sigma values < 0.01

# first pick 10-20 reference stars using imexam in each filter (use the reference image)
if not os.path.isfile('photRefStarsB.txt'):
    iraf.tv.display('sh'+refImageB, frame=1)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('rplot',20.0)
    iraf.tv.imexam('sh'+refImageB, frame=1, logfile='photRefStarsB.txt', keeplog='yes')

if not os.path.isfile('photRefStarsV.txt'):
    iraf.tv.display('sh'+refImageV, frame=1)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('rplot',20.0)
    iraf.tv.imexam('sh'+refImageV, frame=1, logfile='photRefStarsV.txt', keeplog='yes')
    
if not os.path.isfile('photRefStarsR.txt'):
    iraf.tv.display('sh'+refImageR, frame=1)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('rplot',20.0)
    iraf.tv.imexam('sh'+refImageR, frame=1, logfile='photRefStarsR.txt', keeplog='yes')

sortFiles = glob.glob('photRefStars*.sort.txt')
if len(sortFiles) != 3:
    print 'Edit IMEXAM output files to remove saturated stars and extended sources...'
    raw_input("Press Enter when finished:")

fwhmB = np.loadtxt('photRefStarsB.sort.txt', usecols=(3,), unpack=True)
fwhmV = np.loadtxt('photRefStarsV.sort.txt', usecols=(3,), unpack=True)
fwhmR = np.loadtxt('photRefStarsR.sort.txt', usecols=(3,), unpack=True)

fwhmB_AVG = np.mean(fwhmB)
fwhmV_AVG = np.mean(fwhmV)
fwhmR_AVG = np.mean(fwhmR)
# 
# iraf.tv.display('sh'+refImageB, frame=1)
# iraf.unlearn(iraf.tv.tvmark)
# iraf.tv.tvmark.setParam('label',"no")
# iraf.tv.tvmark.setParam('pointsize',7)
# iraf.tv.tvmark.setParam('mark',"circle")
# mark_radii = str( repr(int(4.0*fwhmB_AVG)-2)+','+repr(int(4.0*fwhmB_AVG)-1)+','+repr(int(4.0*fwhmB_AVG))+','+repr(int(4.0*fwhmB_AVG)+1)+','+repr(int(4.0*fwhmB_AVG)+2))
# iraf.tv.tvmark(frame=1, coords='photRefStarsB.txt', radii=mark_radii, color=205)
# mark_radii2 = str( repr(int(5.5*fwhmB_AVG)-2)+','+repr(int(5.5*fwhmB_AVG)-1))
# mark_radii3 = str( repr(int(5.5*fwhmB_AVG)+11)+','+repr(int(5.5*fwhmB_AVG)+12))
# iraf.tv.tvmark(frame=1, coords='photRefStarsB.txt', radii=mark_radii2, color=207)
# iraf.tv.tvmark(frame=1, coords='photRefStarsB.txt', radii=mark_radii3, color=207)

# iraf.task(msc_getscale = "home$scripts/msc_getscale.cl")

iraf.unlearn(iraf.apphot.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)

iraf.apphot.phot.setParam('interactive','no')
iraf.apphot.phot.setParam('verify','no')

iraf.datapars.setParam('gain',"gain")
iraf.datapars.setParam('ccdread',"rdnoise")
iraf.datapars.setParam('exposure',"exptime")
iraf.datapars.setParam('airmass',"airmass")
iraf.datapars.setParam('filter',"filter")
iraf.datapars.setParam('obstime',"time-obs")
iraf.datapars.setParam('datamax',218000.)
iraf.photpars.setParam('zmag',0.)
iraf.centerpars.setParam('cbox',9.)
iraf.centerpars.setParam('maxshift',3.)
iraf.fitskypars.setParam('salgorithm',"median")
iraf.fitskypars.setParam('dannulus',10.)

iraf.datapars.setParam('fwhmpsf',fwhmB_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmB_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmB_AVG)
scaleB, stdB = msc_getscale(filesB, bgSigmaDict, refimage=refImageB, rootname='', prefix='sc', applyscale=False, verbose=False)

iraf.datapars.setParam('fwhmpsf',fwhmV_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmV_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmV_AVG)
scaleV, stdV = msc_getscale(filesV, bgSigmaDict, refimage=refImageV, rootname='', prefix='sc', applyscale=False, verbose=False)

iraf.datapars.setParam('fwhmpsf',fwhmR_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmR_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmR_AVG)
scaleR, stdR = msc_getscale(filesR, bgSigmaDict, refimage=refImageR, rootname='', prefix='sc', applyscale=False, verbose=False)

print '### FIRST PASS SCALING FACTORS'
for i, file_ in enumerate(filesB):
    print file_, bgSigmaB[i], scaleB[i], stdB[i]  

for i, file_ in enumerate(filesV):
    print file_, bgSigmaV[i], scaleV[i], stdV[i]
    
for i, file_ in enumerate(filesR):
    print file_, bgSigmaR[i], scaleR[i], stdR[i]

# if there are scale factors bigger than 1, we have the wrong reference image
# so we need to resort images by scaling factor and try again   
if any(filter(lambda x: x > 1.0, scaleB)) :
    scaleB, filesB = zip(*sorted(zip(scaleB, filesB), reverse=True))
    print scaleB, filesB
    refImageB = filesB[0]
    iraf.datapars.setParam('fwhmpsf',fwhmB_AVG)
    iraf.photpars.setParam('apertures',4.0*fwhmB_AVG)
    iraf.fitskypars.setParam('annulus',5.5*fwhmB_AVG)
    scaleB, stdB = msc_getscale(filesB, bgSigmaDict, refimage=refImageB, rootname='', prefix='sc', applyscale=False, verbose=False)

if any(filter(lambda x: x > 1.0, scaleV)) :
    scaleV, filesV = zip(*sorted(zip(scaleV, filesV), reverse=True))
    refImageV = filesV[0]
    iraf.datapars.setParam('fwhmpsf',fwhmV_AVG)
    iraf.photpars.setParam('apertures',4.0*fwhmV_AVG)
    iraf.fitskypars.setParam('annulus',5.5*fwhmV_AVG)
    scaleV, stdV = msc_getscale(filesV, bgSigmaDict, refimage=refImageV, rootname='', prefix='sc', applyscale=False, verbose=False)

if any(filter(lambda x: x > 1.0, scaleR)) :
    scaleR, filesR = zip(*sorted(zip(scaleR, filesR), reverse=True))
    refImageR = filesR[0]
    iraf.datapars.setParam('fwhmpsf',fwhmR_AVG)
    iraf.photpars.setParam('apertures',4.0*fwhmR_AVG)
    iraf.fitskypars.setParam('annulus',5.5*fwhmR_AVG)
    scaleR, stdR = msc_getscale(filesR, bgSigmaDict, refimage=refImageR, rootname='', prefix='sc', applyscale=False, verbose=False)

print '### SECOND PASS SCALING FACTORS, these were applied'


iraf.datapars.setParam('fwhmpsf',fwhmB_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmB_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmB_AVG)
scaleB, stdB = msc_getscale(filesB, bgSigmaDict, refimage=refImageB, rootname='', prefix='sc', applyscale=True, verbose=False)

iraf.datapars.setParam('fwhmpsf',fwhmV_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmV_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmV_AVG)
scaleV, stdV = msc_getscale(filesV, bgSigmaDict, refimage=refImageV, rootname='', prefix='sc', applyscale=True, verbose=False)

iraf.datapars.setParam('fwhmpsf',fwhmR_AVG)
iraf.photpars.setParam('apertures',4.0*fwhmR_AVG)
iraf.fitskypars.setParam('annulus',5.5*fwhmR_AVG)
scaleR, stdR = msc_getscale(filesR, bgSigmaDict, refimage=refImageR, rootname='', prefix='sc', applyscale=True, verbose=False)

with open('filesB.scl','w+') as f1:
    for i, file_ in enumerate(filesB):
        print >> f1, file_, bgSigmaB[i], scaleB[i], stdB[i]  

with open('filesV.scl','w+') as f1:
    for i, file_ in enumerate(filesV):
        print >> f1, file_, bgSigmaV[i], scaleV[i], stdV[i]

with open('filesR.scl','w+') as f1:
    for i, file_ in enumerate(filesR):
        print >> f1, file_, bgSigmaR[i], scaleR[i], stdR[i]

###########
# finally combine images
# imcombine with average combining, ccdclip
# add back in the background of the REFERENCE IMAGE with imarith

# write to a file what the mean background and bg sigma are
with open(objectName+'.bginfo', 'w+') as f1:
    print >> f1, 'B: ', bgMeanDict[refImageB], bgSigmaDict[refImageB]
    print >> f1, 'V: ', bgMeanDict[refImageV], bgSigmaDict[refImageV]
    print >> f1, 'R: ', bgMeanDict[refImageR], bgSigmaDict[refImageR]


iraf.immatch.imcombine.setParam('reject','ccdclip')
iraf.immatch.imcombine.setParam('gain','gain')
iraf.immatch.imcombine.setParam('rdnoise','rdnoise')
iraf.immatch.imcombine.setParam('combine','average')
if not os.path.isfile(objectName+'_B.fits'):
    iraf.immatch.imcombine('scsubsh//@mscFilesB.list', objectName+'_B.1.fits')
    iraf.imutil.imarith(objectName+'_B.1.fits','+',bgMeanDict[refImageB],objectName+'_B.fits', verbose='yes')
    # iraf.imutil.imarith('checkB.fits','*',bgMeanDict[refImageB],'checkB.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_B.1.fits','+','checkB.2.fits',objectName+'_B.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_B.2.fits','/','checkB.fits',objectName+'_B.fits', verbose='yes')
if not os.path.isfile(objectName+'_V.fits'):
    iraf.immatch.imcombine('scsubsh//@mscFilesV.list', objectName+'_V.1.fits')
    iraf.imutil.imarith(objectName+'_V.1.fits','+',bgMeanDict[refImageV],objectName+'_V.fits', verbose='yes')
    # iraf.imutil.imarith('checkV.fits','*',bgMeanDict[refImageV],'checkV.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_V.1.fits','+','checkV.2.fits',objectName+'_V.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_V.2.fits','/','checkV.fits',objectName+'_V.fits', verbose='yes')
if not os.path.isfile(objectName+'_R.fits'):
    iraf.immatch.imcombine('scsubsh//@mscFilesR.list', objectName+'_R.1.fits')
    iraf.imutil.imarith(objectName+'_R.1.fits','+',bgMeanDict[refImageB],objectName+'_R.fits', verbose='yes')
    # iraf.imutil.imarith('checkR.fits','*',bgMeanDict[refImageR],'checkR.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_R.1.fits','+','checkR.2.fits',objectName+'_R.2.fits', verbose='yes')
    # iraf.imutil.imarith(objectName+'_R.2.fits','/','checkR.fits',objectName+'_R.fits', verbose='yes')

