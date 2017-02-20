#!/usr/bin/env python

import os
import numpy as np
from pyraf import iraf
import glob
#######################
filter_ = 'odi_g'   ######## change me!
#######################


iraf.ptools(_doprint=0)
iraf.tables(_doprint=0)

iraf.unlearn(iraf.ptools.pselect,iraf.tables.tcreate,iraf.tables.tmatch,iraf.tables.tinfo,iraf.tables.tdump)
m3,m4,m5,m6 = np.loadtxt('mask.reg',usecols=(2,3,4,5),unpack=True)

for 
    addstar("AGC249525_"//(s1)//"_sh", "", "default", "AGC249525_"//(s1)//"_sh_add."//(s3), y, z, nartstars)
    daofind("AGC249525_"//(s1)//"_sh_add."//(s3)//".fits", "default")
    phot("AGC249525_"//(s1)//"_sh_add."//(s3), "AGC249525_"//(s1)//"_sh_add."//(s3)//".fits.coo.1", "default", "")
    imdel("AGC249525_"//(s1)//"_sh_add."//(s3), yes)

# g image pselects
# for file_ in glob.glob('*add*.mag.1'):
    if not os.path.isfile(file_+'a'):
        print 'pselect-ing', file_
        if os.path.isfile('temp2') :    
            os.remove('temp2')
        iraf.ptools.pselect(infi=file_, outfi='temp1', expr="MAG != INDEF")
        for i in range(len(m3)) :
            mx1 = m3[i] - (m5[i]/2.)
            mx2 = m3[i] + (m5[i]/2.)
            my1 = m4[i] - (m6[i]/2.)
            my2 = m4[i] + (m6[i]/2.)
            iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
            os.rename('temp2', 'temp1')
        os.rename('temp1', file_+'a')
        if os.path.isfile('temp1') :    
            os.remove('temp1')
        
# for file_ in glob.glob('*.art'):
    if not os.path.isfile(file_[0:-4]+'.pselect'):
        print 'pselect-ing', file_
        if os.path.isfile('temp2') :    
            os.remove('temp2')
        iraf.ptools.pselect(infi=file_, outfi='temp1', expr="MAG != INDEF")
        for i in range(len(m3)) :
            mx1 = m3[i] - (m5[i]/2.)
            mx2 = m3[i] + (m5[i]/2.)
            my1 = m4[i] - (m6[i]/2.)
            my2 = m4[i] + (m6[i]/2.)
            iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
            os.rename('temp2', 'temp1')
        os.rename('temp1', file_[0:-4]+'.pselect')
        if os.path.isfile('temp1') :    
            os.remove('temp1')

# for file_ in glob.glob('*.pselect'):
    if not os.path.isfile("art_star_list."+filter_+"."+file_[-12:-8]+".tab"):
        print 'creating artificial star table from', file_
        iraf.tables.tcreate.setParam("nlines",0)
        iraf.tables.tcreate("art_star_list."+filter_+"."+file_[-12:-8]+".tab", "art_table.description", file_)
    
# for file_ in glob.glob('*.mag.1a'):
    print 'reformatting', file_
    f1 = open(file_, 'r')
    f2 = open('temp1', 'w')
    for line in f1:
        if line[0] == 'A' :
            line = line[0:22]+' '+line[23:72]+' '+line[73:-1]+'\n'
        f2.write(line.replace('\\', ''))
    f1.close()
    f2.close()
    
    os.remove(file_)
    os.rename('temp1', file_)
    if not os.path.isfile("phot_output."+filter_+"."+file_[-11:-7]+".tab"):
        print 'creating phot table from', file_
        iraf.tables.tcreate.setParam("nlines",5)
        iraf.tables.tcreate("phot_output."+filter_+"."+file_[-11:-7]+".tab", "phot_table.description", file_)
    
artFiles = glob.glob('art*.tab')
photFiles = glob.glob('phot*.tab')
namedFiles = glob.glob('*add*.mag.1a')
iraf.tables.tmatch.setParam('incol1',"ID, XCEN, YCEN, MAG")
iraf.tables.tmatch.setParam('incol2',"ID, XINIT, YINIT, MAG, MERR, IFILTER")

for i in range(len(artFiles)):
    if not os.path.isfile(namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch.tab"):
        print 'matching',artFiles[i], '&', photFiles[i]
        iraf.tables.tmatch(artFiles[i], photFiles[i], namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch", "xcen,ycen", "xinit,yinit", 6.)
    if not os.path.isfile(namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch.tdump"):
        iraf.tables.tdump.setParam('datafile',namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch.tdump")
        iraf.tables.tdump(namedFiles[i][0:9]+'.'+filter_+"."+artFiles[i][-8:-4]+"_tmatch")

iraf.tables.tinfo.setParam('ttout',"no")

cTable = open('ctable.out','w+')
print >> cTable, "# Inst mag bin    % Complete"
print >> cTable, "#   "
for i, file_ in enumerate(glob.glob('*tmatch.tab')):
    iraf.tables.tinfo(artFiles[i])
    artk = iraf.tables.tinfo.getParam('nrows')

    iraf.tables.tinfo(file_)
    matchk = iraf.tables.tinfo.getParam('nrows')
    pct = float(matchk)/float(artk)
    print >> cTable, " ",file_[12:16],"        ",pct
cTable.close()