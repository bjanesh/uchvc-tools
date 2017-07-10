#!/usr/bin/env python

import os
import numpy as np
from pyraf import iraf
import glob

#######################
objname = 'AGC238626'
filter_ = 'i'   ####### change me!
max_mag = -4.0
step = 0.1
nartstars = 100
fwhm = 5.671
sigma = 7.445809
threshold = 3.0
#######################
img = objname+'_'+filter_+'_crop.fits'
psf_img = objname+'_'+filter_+'.fits.psf.1.fits'

# make the cropped image
iraf.imcopy(objname+'_'+filter_+'.fits[4500:6500,4500:6500]', img)

iraf.daophot(_doprint=0)
iraf.ptools(_doprint=0)
iraf.tables(_doprint=0)

iraf.unlearn(iraf.ptools.pselect,iraf.tables.tcreate,iraf.tables.tmatch,iraf.tables.tinfo,iraf.tables.tdump)

iraf.datapars.datamax=55000.
iraf.daopars.psfrad=int(round(2.0*fwhm))
iraf.datapars.sigma=sigma
iraf.findpars.threshold=threshold

if os.path.isfile('mask.reg'):
    m3,m4,m5,m6 = np.loadtxt('mask.reg',usecols=(2,3,4,5),unpack=True)

inst_mags = np.arange(max_mag, 1.0, step)

cTable = open('ctable_'+filter_+'.out','w+')
print >> cTable, "# Inst mag bin    % Complete"
print >> cTable, "#   "

for mag in inst_mags:
    art_img = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'.fits'
    dao_coo = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'.coo.1'
    phot_out = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'.mag.1'
    art_coo = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'.fits.art'
    art_tab = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'_art.tab'
    phot_tab = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'_phot.tab'
    match_tab = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'_tmatch.tab'
    match_file = objname+'_'+filter_+'_crop_add.'+'{:5.3f}'.format(mag)+'_tmatch.tdump'
    next_mag = mag + step
    
    print 'adding stars'
    iraf.daophot.addstar(img, "", psf_img, art_img, next_mag, mag, nartstars, verbose='no', verify='no', mode='h')
    print 'detecting sources'
    iraf.daophot.daofind(art_img, dao_coo, verbose='no', verify='no', mode='h')
    print 'photing sources'
    iraf.phot(art_img, dao_coo, phot_out, "", verbose='no', verify='no', mode='h')
    iraf.images.imdel(art_img, "yes", verify='no', mode='h')

# g image pselects
# for file_ in glob.glob('*add*.mag.1'):
    if not os.path.isfile(phot_out+'a'):
        print 'pselect-ing', phot_out
        if os.path.isfile('temp2') :    
            os.remove('temp2')
        iraf.ptools.pselect(infi=phot_out, outfi='temp1', expr="MAG != INDEF")
        if os.path.isfile('mask.reg'):
            for i in range(len(m3)) :
                mx1 = m3[i] - (m5[i]/2.)
                mx2 = m3[i] + (m5[i]/2.)
                my1 = m4[i] - (m6[i]/2.)
                my2 = m4[i] + (m6[i]/2.)
                iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
                os.rename('temp2', 'temp1')
        os.rename('temp1', phot_out+'a')
        if os.path.isfile('temp1') :    
            os.remove('temp1')
        
# for file_ in glob.glob('*.art'):

    if not os.path.isfile(art_coo+'.1a'):
        print 'pselect-ing', art_coo
        if os.path.isfile('temp2') :    
            os.remove('temp2')
        iraf.ptools.pselect(infi=art_coo, outfi='temp1', expr="MAG != INDEF")
        if os.path.isfile('mask.reg'):
            for i in range(len(m3)) :
                mx1 = m3[i] - (m5[i]/2.)
                mx2 = m3[i] + (m5[i]/2.)
                my1 = m4[i] - (m6[i]/2.)
                my2 = m4[i] + (m6[i]/2.)
                iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
                os.rename('temp2', 'temp1')
        os.rename('temp1', art_coo+'.1a')
        if os.path.isfile('temp1') :    
            os.remove('temp1')

# for file_ in glob.glob('*.pselect'):
    if not os.path.isfile(art_tab):
        print 'creating artificial star table from', art_coo+'.1a'
        iraf.tables.tcreate.setParam("nlines",0)
        iraf.tables.tcreate(art_tab, "/Users/wjanesh/projects/uchvc-tools/art_table.description", art_coo+'.1a')
    
# for file_ in glob.glob('*.mag.1a'):
    print 'reformatting', phot_out
    f1 = open(phot_out, 'r')
    f2 = open('temp1', 'w')
    for line in f1:
        if line[0] == 'A' :
            line = line[0:22]+' '+line[23:72]+' '+line[73:-1]+'\n'
        f2.write(line.replace('\\', ''))
    f1.close()
    f2.close()
    
    os.remove(phot_out)
    os.rename('temp1', phot_out)
    if not os.path.isfile(phot_tab):
        print 'creating phot table from', phot_out
        iraf.tables.tcreate.setParam("nlines",5)
        iraf.tables.tcreate(phot_tab, "/Users/wjanesh/projects/uchvc-tools/phot_table.description", phot_out)
    
# artFiles = glob.glob('*_art.tab')
# photFiles = glob.glob('*_phot.tab')
# namedFiles = glob.glob('*add*.mag.1a')
    iraf.tables.tmatch.setParam('incol1',"ID, XCEN, YCEN, MAG")
    iraf.tables.tmatch.setParam('incol2',"ID, XINIT, YINIT, MAG, MERR, IFILTER")
# 
# for i in range(len(artFiles)):
#     if not os.path.isfile(namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch.tab"):
    print 'matching',art_tab, '&', phot_tab
    iraf.tables.tmatch(art_tab, phot_tab, match_tab, "xcen,ycen", "xinit,yinit", 6.)
    if not os.path.isfile(match_file):
        iraf.tables.tdump.setParam('datafile',match_file)
        iraf.tables.tdump(match_tab)
# 
    iraf.tables.tinfo.setParam('ttout',"no")
# 
# for i, file_ in enumerate(glob.glob('*tmatch.tab')):
    iraf.tables.tinfo(art_tab)
    artk = iraf.tables.tinfo.getParam('nrows')

    iraf.tables.tinfo(match_tab)
    matchk = iraf.tables.tinfo.getParam('nrows')
    pct = float(matchk)/float(artk)
    print mag, pct
    print >> cTable, " ",mag,"        ",pct
cTable.close()