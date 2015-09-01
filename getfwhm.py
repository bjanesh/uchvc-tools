from pyraf import iraf
from pyraf.irafpar import makeIrafPar, IrafParList
from stsci.tools.irafglobals import *
from pyraf.pyrafglobals import *

def getfwhm(images='@imh.lis', coordlist='ref_stars', outfile='getfwhm.log',radius=4.0,buffer=7.0,width=5.0,rplot=15.0,center='no',verbose='no',imagelist=None,mode='al',DOLLARnargs=0,taskObj=None):
	Vars = IrafParList('getfwhm')
	Vars.addParam(makeIrafPar(images, datatype='string', name='images',mode='a',prompt='Input image(s)'))
	Vars.addParam(makeIrafPar(coordlist, datatype='string', name='coordlist',mode='a',prompt='List of object positions'))
	Vars.addParam(makeIrafPar(outfile, datatype='string', name='outfile',mode='a',prompt='Output file'))
	Vars.addParam(makeIrafPar(radius, datatype='real', name='radius', mode='h',prompt='Object radius'))
	Vars.addParam(makeIrafPar(buffer, datatype='real', name='buffer', mode='h',prompt='Background buffer width'))
	Vars.addParam(makeIrafPar(width, datatype='real', name='width', mode='h',prompt='Background width'))
	Vars.addParam(makeIrafPar(rplot, datatype='real', name='rplot', mode='h',prompt='Plotting radius'))
	Vars.addParam(makeIrafPar(center, datatype='bool', name='center', mode='h',prompt='Center object in aperture?'))
	Vars.addParam(makeIrafPar(verbose, datatype='bool', name='verbose',mode='h',prompt='Verbose output?'))
	Vars.addParam(makeIrafPar(imagelist, datatype='struct', name='imagelist',list_flag=1,mode='h'))
	Vars.addParam(makeIrafPar(mode, datatype='string', name='mode', mode='h'))
	Vars.addParam(makeIrafPar(DOLLARnargs, datatype='int', name='$nargs',mode='h'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='len', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='image', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='imagefile',mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='imlist', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='coords', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='outputfile',mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='rad', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='buff', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='wid', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='rpl', mode='u'))

	Vars.imlist = Vars.images
	Vars.imagefile = iraf.mktemp('tmp$getfwhm')
	iraf.sections(Vars.imlist, option = 'fullname', Stdout=Vars.imagefile)
	Vars.imagelist = Vars.imagefile
	Vars.coords = Vars.coordlist
	Vars.outputfile = Vars.outfile
	Vars.rad = Vars.radius
	Vars.buff = Vars.buffer
	Vars.wid = Vars.width
	Vars.rpl = Vars.rplot
	iraf.rimexam.radius = Vars.rad
	iraf.rimexam.buffer = Vars.buff
	iraf.rimexam.width = Vars.wid
	iraf.rimexam.rplot = Vars.rpl
	if (Vars.center == yes):
		iraf.rimexam.center = yes
	else:
		iraf.rimexam.center = no
	iraf.rimexam.fittype = 'gaussian'
	iraf.rimexam.iterati = 1
	iraf.clPrint('# RIMEXAM Parameter Settings:  Radius=', Vars.rad,', Buffer=',Vars.buff,', Width=',Vars.wid,StdoutAppend=Vars.outputfile)
	while (iraf.fscan(locals(), 'Vars.imagelist', 'Vars.image') != EOF):
		Vars.len = iraf.strlen(Vars.image)
		if (iraf.substr(Vars.image, Vars.len - 3, Vars.len) == '.imh'):
			Vars.image = iraf.substr(Vars.image, 1, Vars.len - 4)
		if (Vars.verbose):
			iraf.clPrint('...processing ', Vars.image)
		iraf.clPrint(Vars.image)
		iraf.imexamine(Vars.image, logfile = Vars.outputfile, keeplog = yes,defkey = 'a',imagecur = Vars.coords,wcs = 'world',use_display = no)