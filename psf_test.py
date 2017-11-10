#! /usr/local/bin/python
import os, sys, shutil
import numpy as np
from subprocess import call
from astropy.io import fits
from pyraf import iraf


def main(argv):
	home_root = os.environ['HOME']
	iraf.images(_doprint=0)
	iraf.tv(_doprint=0)
	iraf.ptools(_doprint=0)
	iraf.noao(_doprint=0)
	iraf.digiphot(_doprint=0)
	iraf.photcal(_doprint=0)
	iraf.apphot(_doprint=0)  
	iraf.imutil(_doprint=0)

	for file_ in os.listdir("./"):
		if file_.endswith("_i.fits"):
			fits_file_i = file_
			title_string = file_[0:9]
		if file_.endswith("_g.fits"):
			fits_file_g = file_
				
	
	fits_h_i = fits.open(fits_file_i)
	fits_h_g = fits.open(fits_file_g)
	
	fwhm_i = fits_h_i[0].header['FWHMPSF']
	fwhm_g = fits_h_g[0].header['FWHMPSF']
	print 'Image FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(fwhm_g,fwhm_i)
	
	# if not os.path.isdir('psf'): # make a psf subfolder if it doesn't exist
	# 	os.mkdir('psf')
	
	# # and then copy in the data files if they don't already exist as files
	# if not os.path.isfile('psf/'+fits_file_i):
	# 	iraf.images.imcopy(fits_file_i,'psf/'+fits_file_i,verbose="yes")
	# 
	# if not os.path.isfile('psf/'+fits_file_g) :
	# 	iraf.images.imcopy(fits_file_g,'psf/'+fits_file_g,verbose="yes")
	# 
	# # also copy in the apcor star lists
	# if not (os.path.isfile('psf/apcor_stars_i.txt') or os.path.islink('psf/apcor_stars_i.txt')) :
	# 	shutil.copyfile('apcor_stars_i.txt','psf/apcor_stars_i.txt')	
	# 
	# if not (os.path.isfile('psf/apcor_stars_g.txt') or os.path.islink('psf/apcor_stars_g.txt')) :
	# 	shutil.copyfile('apcor_stars_g.txt','psf/apcor_stars_g.txt')	
	# 	
	# # and change to the psf folder
	# os.chdir('psf')

	iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
	
	iraf.apphot.phot.setParam('interactive',"no")
	iraf.apphot.phot.setParam('verify',"no")
	iraf.datapars.setParam('datamax',"INDEF")
	iraf.datapars.setParam('gain',"gain")
	iraf.datapars.setParam('ccdread',"rdnoise")
	iraf.datapars.setParam('exposure',"exptime")
	iraf.datapars.setParam('airmass',"airmass")
	iraf.datapars.setParam('filter',"filter")
	iraf.datapars.setParam('obstime',"time-obs")
	iraf.datapars.setParam('sigma',"INDEF")
	iraf.photpars.setParam('zmag',0.)
	iraf.centerpars.setParam('cbox',9.)
	iraf.centerpars.setParam('maxshift',3.)
	iraf.fitskypars.setParam('salgorithm',"median")
	iraf.fitskypars.setParam('dannulus',10.)
	
	iraf.datapars.setParam('fwhmpsf',fwhm_g)
	iraf.photpars.setParam('apertures',iraf.nint(4.*fwhm_g))
	iraf.fitskypars.setParam('annulus',(6.*fwhm_g))
	
	iraf.apphot.phot(image=fits_file_g, coords='apcor_stars_g.txt', output="g_psfstars.mag.1")
	
	iraf.datapars.setParam('fwhmpsf',fwhm_i)
	iraf.photpars.setParam('apertures',iraf.nint(4.*fwhm_i))
	iraf.fitskypars.setParam('annulus',(6.*fwhm_i))
	
	iraf.apphot.phot(image=fits_file_i, coords='apcor_stars_i.txt', output="i_psfstars.mag.1")
	
	with open("apcor_stars_i.txt") as foo:
		lines = len(foo.readlines())
	
	iraf.daophot(_doprint=0)
	iraf.pstselect.setParam('image',fits_file_i)
	iraf.pstselect.setParam('photfile',"i_psfstars.mag.1")
	iraf.pstselect.setParam('pstfile',"default")
	iraf.pstselect.setParam('maxnpsf',lines)
	iraf.pstselect.setParam('plottype',"mesh")
	
	iraf.daophot.pstselect(verify='no', mode='h')
	
	with open("apcor_stars_g.txt") as foo:
		lines = len(foo.readlines())
	
	iraf.pstselect.setParam('image',fits_file_g)
	iraf.pstselect.setParam('photfile',"g_psfstars.mag.1")
	iraf.pstselect.setParam('pstfile',"default")
	iraf.pstselect.setParam('maxnpsf',lines)
	iraf.pstselect.setParam('plottype',"mesh")
	
	iraf.daophot.pstselect()
	
	iraf.datapars.setParam('datamax',"INDEF")
	iraf.datapars.setParam('gain',"gain")
	iraf.datapars.setParam('ccdread',"rdnoise")
	iraf.datapars.setParam('exposure',"exptime")
	iraf.datapars.setParam('airmass',"airmass")
	iraf.datapars.setParam('filter',"filter")
	iraf.datapars.setParam('obstime',"time-obs")
	iraf.datapars.setParam('sigma',"INDEF")
	iraf.daopars.setParam('psfrad',iraf.nint(3.*fwhm_i))
	iraf.daopars.setParam('fitrad',fwhm_i)
	
	iraf.psf.setParam('image',fits_file_i)
	iraf.psf.setParam('photfile',"i_psfstars.mag.1")
	iraf.psf.setParam('pstfile',fits_file_i+".pst.1")
	iraf.psf.setParam('interactive', 'no')
	iraf.daophot.psf()
	
	iraf.psf.setParam('image',fits_file_g)
	iraf.psf.setParam('photfile',"g_psfstars.mag.1")
	iraf.psf.setParam('pstfile',fits_file_g+".pst.1")
	iraf.psf.setParam('interactive', 'no')
	iraf.daophot.psf()
	
if __name__ == "__main__":
	main(sys.argv[1:])	
