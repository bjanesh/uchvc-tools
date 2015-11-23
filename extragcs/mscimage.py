from pyraf import iraf
from pyraf.irafpar import makeIrafPar, IrafParList
from stsci.tools.irafglobals import *
from pyraf.pyrafglobals import *

def mscimage(input=None, output=None, format='image', pixmask=no,verbose=')_.verbose',wcssource='image',reference='',ra=INDEF,dec=INDEF,scale=INDEF,rotation=INDEF,blank=0.0,interpolant='poly5',minterpolant='linear',boundary='reflect',constant=0.0,fluxconserve=no,ntrim=8,nxblock=INDEF,nyblock=INDEF,interactive=no,nx=10,ny=20,fitgeometry='general',xxorder=4,xyorder=4,xxterms='half',yxorder=4,yyorder=4,yxterms='half',fd_in='',fd_ext='',fd_coord='',mode='ql',DOLLARnargs=0,taskObj=None):

	Vars = IrafParList('mscimage')
	Vars.addParam(makeIrafPar(input, datatype='string', name='input', mode='a',prompt='List of input mosaic exposures'))
	Vars.addParam(makeIrafPar(output, datatype='string', name='output',mode='a',prompt='List of output images'))
	Vars.addParam(makeIrafPar(format, datatype='string', name='format',enum=['image', 'mef'],mode='h',prompt='Output format (image|mef)'))
	Vars.addParam(makeIrafPar(pixmask, datatype='bool', name='pixmask',mode='h',prompt='Create pixel mask?'))
	Vars.addParam(makeIrafPar(verbose, datatype='bool', name='verbose',mode='h',prompt='Verbose output?\n\n# Output WCS parameters'))
	Vars.addParam(makeIrafPar(wcssource, datatype='string', name='wcssource',enum=['image', 'parameters', 'match'],mode='h',prompt='Output WCS source (image|parameters|match)'))
	Vars.addParam(makeIrafPar(reference, datatype='file', name='reference',mode='h',prompt='Reference image'))
	Vars.addParam(makeIrafPar(ra, datatype='real', name='ra', max=24.0,min=0.0,mode='h',prompt='RA of tangent point (hours)'))
	Vars.addParam(makeIrafPar(dec, datatype='real', name='dec', max=90.0,min=-90.0,mode='h',prompt='DEC of tangent point (degrees)'))
	Vars.addParam(makeIrafPar(scale, datatype='real', name='scale', mode='h',prompt='Scale (arcsec/pixel)'))
	Vars.addParam(makeIrafPar(rotation, datatype='real', name='rotation',max=360.0,min=-360.0,mode='h',prompt='Rotation of DEC from N to E (degrees)\n\n# Resampling parmeters'))
	Vars.addParam(makeIrafPar(blank, datatype='real', name='blank', mode='h',prompt='Blank value'))
	Vars.addParam(makeIrafPar(interpolant, datatype='string',name='interpolant',mode='h',prompt='Interpolant for data'))
	Vars.addParam(makeIrafPar(minterpolant, datatype='string',name='minterpolant',mode='h',prompt='Interpolant for mask'))
	Vars.addParam(makeIrafPar(boundary, datatype='string', name='boundary',enum=['nearest', 'constant', 'reflect', 'wrap'],mode='h',prompt='Boundary extension'))
	Vars.addParam(makeIrafPar(constant, datatype='real', name='constant',mode='h',prompt='Constant boundary extension value'))
	Vars.addParam(makeIrafPar(fluxconserve, datatype='bool',name='fluxconserve',mode='h',prompt='Preserve flux per unit area?'))
	Vars.addParam(makeIrafPar(ntrim, datatype='int', name='ntrim', min=0,mode='h',prompt='Edge trim in each extension'))
	Vars.addParam(makeIrafPar(nxblock, datatype='int', name='nxblock',mode='h',prompt='X dimension of working block size in pixels'))
	Vars.addParam(makeIrafPar(nyblock, datatype='int', name='nyblock',mode='h',prompt='Y dimension of working block size in pixels\n\n# Geometric mapping parameters'))
	Vars.addParam(makeIrafPar(interactive, datatype='bool', name='interactive',mode='h',prompt='Fit mapping interactively?'))
	Vars.addParam(makeIrafPar(nx, datatype='int', name='nx', mode='h',prompt='Number of x grid points'))
	Vars.addParam(makeIrafPar(ny, datatype='int', name='ny', mode='h',prompt='Number of y grid points'))
	Vars.addParam(makeIrafPar(fitgeometry, datatype='string',name='fitgeometry',enum=['shift', 'xyscale', 'rotate', 'rscale', 'rxyscale', 'general'],mode='h',prompt='Fitting geometry'))
	Vars.addParam(makeIrafPar(xxorder, datatype='int', name='xxorder', min=2,mode='h',prompt='Order of x fit in x'))
	Vars.addParam(makeIrafPar(xyorder, datatype='int', name='xyorder', min=2,mode='h',prompt='Order of x fit in y'))
	Vars.addParam(makeIrafPar(xxterms, datatype='string', name='xxterms',mode='h',prompt='X fit cross terms type'))
	Vars.addParam(makeIrafPar(yxorder, datatype='int', name='yxorder', min=2,mode='h',prompt='Order of y fit in x'))
	Vars.addParam(makeIrafPar(yyorder, datatype='int', name='yyorder', min=2,mode='h',prompt='Order of y fit in y'))
	Vars.addParam(makeIrafPar(yxterms, datatype='string', name='yxterms',mode='h',prompt='Y fit cross terms type\n\n'))
	Vars.addParam(makeIrafPar(fd_in, datatype='struct', name='fd_in',list_flag=1,mode='h',prompt=''))
	Vars.addParam(makeIrafPar(fd_ext, datatype='struct', name='fd_ext',list_flag=1,mode='h',prompt=''))
	Vars.addParam(makeIrafPar(fd_coord, datatype='struct', name='fd_coord',list_flag=1,mode='h',prompt=''))
	Vars.addParam(makeIrafPar(mode, datatype='string', name='mode', mode='h',prompt=''))
	Vars.addParam(makeIrafPar(DOLLARnargs, datatype='int', name='$nargs',mode='h'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='in', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='out', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='ref', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='pl', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='image', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='trimsec', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='outsec', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='plsec', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='inlists', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='extlist', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='pllist', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='coord', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='db', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='wcsref', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='outtemp', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='file', name='pltemp', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nc', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nl', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='ncref', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nlref', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='cmin', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='cmax', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='lmin', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='lmax', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nimage', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nimages', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nxblk', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='int', name='nyblk', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='x', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='y', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='rval', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='xmin', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='xmax', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='ymin', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='ymax', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='crpix1', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='real', name='crpix2', mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='extname',mode='u'))
	Vars.addParam(makeIrafPar(None, datatype='string', name='str', mode='u'))

	iraf.cache('mscextensions', 'mscgmask')
	Vars.inlists = iraf.mktemp('tmp$iraf')
	Vars.extlist = iraf.mktemp('tmp$iraf')
	Vars.pllist = iraf.mktemp('tmp$iraf')
	Vars.coord = iraf.mktemp('tmp$iraf')
	Vars.db = iraf.mktemp('tmp$iraf')
	Vars.outtemp = iraf.mktemp('tmp')
	Vars.wcsref = iraf.mktemp('tmp')
	Vars.pltemp = iraf.mktemp('tmp')
	iraf.joinlists(Vars.input, Vars.output, output = Vars.inlists, delim = ' ',short=yes,type = 'image')
	Vars.fd_in = Vars.inlists
	while (iraf.fscan(locals(), 'Vars.fd_in', 'Vars.PYin', 'Vars.out') != EOF):
		if (iraf.imaccess(Vars.out)):
			iraf.printf('Warning: Image already exists (%s)\n', Vars.out)
			continue
		if (Vars.pixmask):
			Vars.pl = Vars.out
			Vars.nc = iraf.strlen(Vars.pl)
			if (Vars.nc > 5 and iraf.substr(Vars.pl, Vars.nc - 4, Vars.nc) == '.fits'):
				Vars.pl = iraf.substr(Vars.pl, 1, Vars.nc - 5)
			elif (Vars.nc > 4 and iraf.substr(Vars.out, Vars.nc - 3, Vars.nc) == '.imh'):
				Vars.pl = iraf.substr(Vars.pl, 1, Vars.nc - 4)
			Vars.pl = Vars.pl + '_bpm'
			if (Vars.format == 'image' and iraf.imaccess(Vars.pl)):
				iraf.printf('Warning: Mask already exists (%s)\n', Vars.pl)
				continue
		else:
			Vars.pl = ''
		iraf.mscextensions(Vars.PYin, output = 'file', index = '0-',extname = '',extver = '',lindex = no,lname = yes,lver = no,ikparams = '',Stdout=Vars.extlist)
		Vars.nimages = int(iraf.mscextensions.nimages)
		Vars.nimage = 0
		if (Vars.nimages < 1):
			iraf.printf("WARNING: No input image data found in `%s'.\
",Vars.PYin)
			iraf.delete(Vars.extlist, verify = no)
			continue
		if (not iraf.imaccess(Vars.wcsref)):
			Vars.ref = Vars.reference
			if (Vars.wcssource == 'match'):
				Vars.wcsref = Vars.ref
			else:
				iraf.mscwtemplate('@' + Vars.extlist, Vars.wcsref,wcssource = Vars.wcssource,reference = Vars.ref,ra = Vars.ra,dec = Vars.dec,scale = Vars.scale,rotation = Vars.rotation,projection = '',verbose = Vars.verbose)
		Vars.fd_ext = Vars.extlist
		while (iraf.fscan(locals(), 'Vars.fd_ext', 'Vars.image') != EOF):
			Vars.nimage = Vars.nimage + 1
			if (Vars.nimages > 1):
				Pipe1 = iraf.hselect(Vars.image, 'extname', yes, Stdout=1)
				iraf.scan(locals(), 'Vars.extname', Stdin=Pipe1)
				del Pipe1
				if (iraf.nscan() == 0):
					Vars.extname = 'im' + str(Vars.nimage)
				Pipe1 = iraf.printf('%s[%s,append]\n', Vars.outtemp,Vars.extname,Stdout=1)
				iraf.scan(locals(), 'Vars.outsec', Stdin=Pipe1)
				del Pipe1
				Pipe1 = iraf.printf('%s%s\n', Vars.pl, Vars.extname, Stdout=1)
				iraf.scan(locals(), 'Vars.plsec', Stdin=Pipe1)
				del Pipe1
			else:
				Vars.extname = ''
				Vars.outsec = Vars.outtemp
				Vars.plsec = Vars.pl
			if (Vars.pixmask and iraf.imaccess(Vars.plsec)):
				iraf.delete(Vars.coord, verify = no)
				iraf.delete(Vars.db, verify = no)
				iraf.printf('Warning: Mask already exists (%s)\n', Vars.plsec)
				continue
			if (Vars.verbose):
				iraf.printf('Resampling %s ...\n', Vars.image)
			Pipe1 = iraf.hselect(Vars.image, 'naxis1,naxis2', yes, Stdout=1)
			iraf.scan(locals(), 'Vars.nc', 'Vars.nl', Stdin=Pipe1)
			del Pipe1
			Vars.cmin = 1 + Vars.ntrim
			Vars.cmax = Vars.nc - Vars.ntrim
			Vars.lmin = 1 + Vars.ntrim
			Vars.lmax = Vars.nl - Vars.ntrim
			Pipe1 = iraf.printf('[%d:%d,%d:%d]\n', Vars.cmin, Vars.cmax,Vars.lmin,Vars.lmax,Stdout=1)
			iraf.scan(locals(), 'Vars.trimsec', Stdin=Pipe1)
			del Pipe1
			if (Vars.wcssource == 'match'):
				Pipe1 = iraf.hselect(Vars.ref, 'naxis1,naxis2', yes, Stdout=1)
				iraf.scan(locals(), 'Vars.ncref', 'Vars.nlref', Stdin=Pipe1)
				del Pipe1
				Vars.xmin = (Vars.ncref - 1.) / (Vars.nx - 1.)
				Vars.ymin = (Vars.nlref - 1.) / (Vars.ny - 1.)
				Vars.ymax = 1
				while (Vars.ymax <= Vars.nlref + 1):
					Vars.xmax = 1
					while (Vars.xmax <= Vars.ncref + 1):
						iraf.clPrint(Vars.xmax, Vars.ymax, Vars.xmax,Vars.ymax,StdoutAppend=Vars.coord)
						Vars.xmax = Vars.xmax + Vars.xmin
					Vars.ymax = Vars.ymax + Vars.ymin
				iraf.mscctran(Vars.coord, Vars.db, Vars.ref, 'logical','world',columns = '3 4',units = '',formats = '%.4H %.3h',min_sigdigit = 10,verbose = no)
				iraf.delete(Vars.coord, verify=no)
				iraf.wcsctran(Vars.db, Vars.coord, Vars.image + Vars.trimsec,inwcs = 'world',outwcs = 'logical',columns = '3 4',units = 'hours native',formats = '',min_sigdigit = 10,verbose = no)
				iraf.delete(Vars.db, verify=no)
			else:
				Vars.nc = Vars.cmax - Vars.cmin + 1
				Vars.nl = Vars.lmax - Vars.lmin + 1
				Vars.xmin = (Vars.nc - 1.) / (Vars.nx - 1.)
				Vars.ymin = (Vars.nl - 1.) / (Vars.ny - 1.)
				Vars.ymax = 1
				while (Vars.ymax <= Vars.nl + 1):
					Vars.xmax = 1
					while (Vars.xmax <= Vars.nc + 1):
						iraf.clPrint(Vars.xmax, Vars.ymax, Vars.xmax,Vars.ymax,StdoutAppend=Vars.coord)
						Vars.xmax = Vars.xmax + Vars.xmin
					Vars.ymax = Vars.ymax + Vars.ymin
				iraf.mscctran(Vars.coord, Vars.db, Vars.image + Vars.trimsec,'logical','world',columns = '1 2',units = '',formats = '%.4H %.3h',min_sigdigit = 10,verbose = no)
				iraf.delete(Vars.coord, verify=no)
				iraf.wcsctran(Vars.db, Vars.coord, Vars.wcsref,inwcs = 'world',outwcs = 'logical',columns = '1 2',units = 'hours native',formats = '',min_sigdigit = 10,verbose = no)
				iraf.delete(Vars.db, verify=no)
			Vars.xmax = 0.
			Vars.xmin = 1.
			Vars.ymax = 0.
			Vars.ymin = 1.
			Vars.fd_coord = Vars.coord
			while (iraf.fscan(locals(), 'Vars.fd_coord', 'Vars.x', 'Vars.y') != EOF):
				if (iraf.nscan() < 2):
					continue
				if (Vars.xmax < Vars.xmin):
					Vars.xmin = Vars.x
					Vars.xmax = Vars.x
					Vars.ymin = Vars.y
					Vars.ymax = Vars.y
				else:
					Vars.xmin = float(iraf.minimum(Vars.x, Vars.xmin))
					Vars.xmax = float(iraf.maximum(Vars.x, Vars.xmax))
					Vars.ymin = float(iraf.minimum(Vars.y, Vars.ymin))
					Vars.ymax = float(iraf.maximum(Vars.y, Vars.ymax))
			Vars.fd_coord = ''
			if (Vars.xmax <= Vars.xmin or Vars.ymax <= Vars.ymin):
				iraf.error(1, 'No overlap for matching reference')
			Vars.cmin = int(iraf.nint(Vars.xmin - 1.5))
			Vars.cmax = int(iraf.nint(Vars.xmax + 1.5))
			Vars.lmin = int(iraf.nint(Vars.ymin - 1.5))
			Vars.lmax = int(iraf.nint(Vars.ymax + 1.5))
			iraf.geomap(Vars.coord, Vars.db, Vars.cmin, Vars.cmax, Vars.lmin,Vars.lmax,transforms = '',results = '',fitgeometry = Vars.fitgeometry,function = 'chebyshev',xxorder = Vars.xxorder,xyorder = Vars.xyorder,xxterms = Vars.xxterms,yxorder = Vars.yxorder,yyorder = Vars.yyorder,yxterms = Vars.yxterms,reject = INDEF,calctype = 'double',verbose = no,interactive = Vars.interactive,graphics = 'stdgraph',cursor = '')
			if (Vars.wcssource == 'match'):
				Vars.cmin = 1
				Vars.lmin = 1
				Vars.cmax = Vars.ncref
				Vars.lmax = Vars.nlref
			if (Vars.nxblock == INDEF):
				Vars.nxblk = Vars.cmax - Vars.cmin + 3
			else:
				Vars.nxblk = Vars.nxblock
			if (Vars.nyblock == INDEF):
				Vars.nyblk = Vars.lmax - Vars.lmin + 3
			else:
				Vars.nyblk = Vars.nyblock
			iraf.geotran(Vars.image + Vars.trimsec, Vars.outsec, Vars.db,Vars.coord,geometry = 'geometric',xin = INDEF,yin = INDEF,xshift = INDEF,yshift = INDEF,xout = INDEF,yout = INDEF,xmag = INDEF,ymag = INDEF,xrotation = INDEF,yrotation = INDEF,xmin = Vars.cmin,xmax = Vars.cmax,ymin = Vars.lmin,ymax = Vars.lmax,xsample = 10.,ysample = 10.,xscale = 1.,yscale = 1.,ncols = INDEF,nlines = INDEF,interpolant = Vars.interpolant,boundary = 'constant',constant = Vars.constant,fluxconserve = Vars.fluxconserve,nxblock = Vars.nxblk,nyblock = Vars.nyblk,verbose = no)
			iraf.wcscopy(Vars.outsec, Vars.wcsref, verbose=no)
			Vars.xmin = 0.
			Vars.ymin = 0.
			Pipe1 = iraf.hselect(Vars.outsec, 'crpix1,crpix2', yes, Stdout=1)
			iraf.scan(locals(), 'Vars.xmin', 'Vars.ymin', Stdin=Pipe1)
			del Pipe1
			Vars.xmin = Vars.xmin - Vars.cmin + 1
			Vars.ymin = Vars.ymin - Vars.lmin + 1
			if (Vars.nimage == 1):
				Vars.crpix1 = Vars.xmin
				Vars.crpix2 = Vars.ymin
			else:
				Vars.crpix1 = float(iraf.maximum(Vars.crpix1, Vars.xmin))
				Vars.crpix2 = float(iraf.maximum(Vars.crpix2, Vars.ymin))
			iraf.hedit(Vars.outsec, 'crpix1', Vars.xmin, add=yes, verify=no,show=no,update=yes)
			iraf.hedit(Vars.outsec, 'crpix2', Vars.ymin, add=yes, verify=no,show=no,update=yes)
			if (Vars.pixmask):
				Pipe1 = iraf.printf('%s%s\n', Vars.pl, Vars.extname, Stdout=1)
				iraf.scan(locals(), 'Vars.plsec', Stdin=Pipe1)
				del Pipe1
				iraf.mscgmask(Vars.image + Vars.trimsec, Vars.pltemp + '.pl','BPM',mval = 10000)
				iraf.geotran(Vars.pltemp, Vars.plsec + '.fits', Vars.db,Vars.coord,geometry = 'geometric',xin = INDEF,yin = INDEF,xshift = INDEF,yshift = INDEF,xout = INDEF,yout = INDEF,xmag = INDEF,ymag = INDEF,xrotation = INDEF,yrotation = INDEF,xmin = Vars.cmin,xmax = Vars.cmax,ymin = Vars.lmin,ymax = Vars.lmax,xsample = 10.,ysample = 10.,interpolant = Vars.minterpolant,boundary = 'constant',constant = 20000.,fluxconserve = no,nxblock = Vars.nxblk,nyblock = Vars.nyblk,verbose = no)
				iraf.imdelete(Vars.pltemp, verify=no)
				iraf.mscpmask(Vars.plsec + '.fits', Vars.plsec + '.pl')
				iraf.imdelete(Vars.plsec + '.fits', verify=no)
				iraf.hedit(Vars.outsec, 'BPM', Vars.plsec + '.pl', add=yes,show=no,verify=no,update=yes)
				iraf.wcscopy(Vars.plsec, Vars.outsec, verbose=no)
				iraf.clPrint(Vars.plsec, StdoutAppend=Vars.pllist)
			else:
				iraf.hedit(Vars.outsec, 'BPM', PYdel=yes, add=no, addonly=no,show=no,verify=no,update=yes)
			iraf.delete(Vars.coord, verify = no)
			iraf.delete(Vars.db, verify = no)
		Vars.fd_ext = ''
		iraf.delete(Vars.extlist, verify = no)
		if (Vars.nimages > 1 and Vars.format == 'image'):
			if (Vars.verbose):
				iraf.printf('Creating image %s ...\n', Vars.out)
			iraf.mscextensions(Vars.outtemp, output = 'file', index = '',extname = '',extver = '',lindex = no,lname = yes,lver = no,ikparams = '',Stdout=Vars.extlist)
			if (Vars.pixmask):
				iraf.combine('@' + Vars.pllist, Vars.pltemp + '.pl',headers = '',bpmasks = Vars.pl,rejmasks = '',nrejmasks = '',expmasks = '',sigmas = '',imcmb = '',ccdtype = '',amps = no,subsets = no,delete = no,combine = 'average',reject = 'none',project = no,outtype = 'real',outlimits = '',offsets = 'wcs',masktype = 'none',maskvalue = '0',blank = 0.,scale = 'none',zero = 'none',weight = 'none',statsec = '',lthreshold = INDEF,hthreshold = 0.99,nlow = 1,nhigh = 1,nkeep = 1,mclip = yes,lsigma = 3.,hsigma = 3.,rdnoise = '0.',gain = '1.',snoise = '0.',sigscale = 0.1,pclip =  - 0.5,grow = 0.,Stdout='dev$null')
				iraf.imdelete(Vars.pltemp, verify=no)
				iraf.combine('@' + Vars.extlist, Vars.out, headers = '',bpmasks = '',rejmasks = '',nrejmasks = '',expmasks = '',sigmas = '',imcmb = '',ccdtype = '',amps = no,subsets = no,delete = no,combine = 'average',reject = 'none',project = no,outtype = 'real',outlimits = '',offsets = 'wcs',masktype = 'badvalue',maskvalue = '2',blank = 0.,scale = 'none',zero = 'none',weight = 'none',statsec = '',lthreshold = INDEF,hthreshold = INDEF,nlow = 1,nhigh = 1,nkeep = 1,mclip = yes,lsigma = 3.,hsigma = 3.,rdnoise = '0.',gain = '1.',snoise = '0.',sigscale = 0.1,pclip =  - 0.5,grow = 0.,Stdout='dev$null')
				iraf.hedit(Vars.out, 'BPM', Vars.pl, add=yes, verify=no,show=no,update=yes)
				iraf.hedit(Vars.pl, 'IMCMB???,PROCID??', add=no, addonly=no,PYdel=yes,update=yes,verify=no,show=no)
			else:
				iraf.combine('@' + Vars.extlist, Vars.out, headers = '',bpmasks = '',rejmasks = '',nrejmasks = '',expmasks = '',sigmas = '',imcmb = '',ccdtype = '',amps = no,subsets = no,delete = no,combine = 'average',reject = 'none',project = no,outtype = 'real',outlimits = '',offsets = 'wcs',masktype = 'none',maskvalue = '2',blank = 0.,scale = 'none',zero = 'none',weight = 'none',statsec = '',lthreshold = INDEF,hthreshold = INDEF,nlow = 1,nhigh = 1,nkeep = 1,mclip = yes,lsigma = 3.,hsigma = 3.,rdnoise = '0.',gain = '1.',snoise = '0.',sigscale = 0.1,pclip =  - 0.5,grow = 0.,Stdout='dev$null')
			Pipe2 = iraf.hselect('@' + Vars.extlist, 'gain', yes, Stdout=1)
			Pipe1 = iraf.average(data_value = 0., Stdin=Pipe2, Stdout=1)
			del Pipe2
			iraf.scan(locals(), 'Vars.rval', Stdin=Pipe1)
			del Pipe1
			iraf.hedit(Vars.out, 'gain', Vars.rval, add=yes, PYdel=no,update=yes,verify=no,show=no)
			Pipe2 = iraf.hselect('@' + Vars.extlist, 'rdnoise', yes, Stdout=1)
			Pipe1 = iraf.average(data_value = 0., Stdin=Pipe2, Stdout=1)
			del Pipe2
			iraf.scan(locals(), 'Vars.rval', Stdin=Pipe1)
			del Pipe1
			iraf.hedit(Vars.out, 'rdnoise', Vars.rval, add=yes, PYdel=no,update=yes,verify=no,show=no)
			iraf.hedit(Vars.out, 'IMCMB???,PROCID??', add=no, addonly=no,PYdel=yes,update=yes,verify=no,show=no)
			iraf.hedit(Vars.out,'NEXTEND,DETSEC,CCDSEC,AMPSEC,IMAGEID,DATASEC,TRIMSEC,BIASSEC',add=no,addonly=no,PYdel=yes,update=yes,verify=no,show=no)
			iraf.imdelete(Vars.outtemp, verify=no)
			if (iraf.access(Vars.pllist)):
				iraf.imdelete('@' + Vars.pllist, verify=no)
				iraf.delete(Vars.pllist, verify=no)
			iraf.delete(Vars.extlist, verify = no)
		elif (Vars.nimages > 1):
			iraf.imrename(Vars.outtemp, Vars.out, verbose=no)
			iraf.mscextensions(Vars.out, output = 'file', index = '',extname = '',extver = '',lindex = no,lname = yes,lver = no,ikparams = '',Stdout=Vars.extlist)
			Vars.fd_ext = Vars.extlist
			while (iraf.fscan(locals(), 'Vars.fd_ext', 'Vars.image') != EOF):
				Pipe1 = iraf.hselect(Vars.image, 'naxis1,naxis2,crpix1,crpix2',yes,Stdout=1)
				iraf.scan(locals(), 'Vars.nc', 'Vars.nl', 'Vars.xmin','Vars.ymin',Stdin=Pipe1)
				del Pipe1
				Vars.cmin = int(iraf.nint(Vars.crpix1 - Vars.xmin + 1))
				Vars.lmin = int(iraf.nint(Vars.crpix2 - Vars.ymin + 1))
				Vars.cmax = Vars.nc + Vars.cmin - 1
				Vars.lmax = Vars.nl + Vars.lmin - 1
				Pipe1 = iraf.printf('[%d:%d,%d:%d]\n', Vars.cmin, Vars.cmax,Vars.lmin,Vars.lmax,Stdout=1)
				iraf.scan(locals(), 'Vars.str', Stdin=Pipe1)
				del Pipe1
				iraf.hedit(Vars.image, 'DETSEC', Vars.str, add=yes, verify=no,show=no,update=yes)
				iraf.hedit(Vars.image, 'DTM1_1', 1., add=yes, verify=no,show=no,update=yes)
				iraf.hedit(Vars.image, 'DTM2_2', 1., add=yes, verify=no,show=no,update=yes)
				Vars.cmin = Vars.cmin - 1
				Vars.lmin = Vars.lmin - 1
				iraf.hedit(Vars.image, 'DTV1', Vars.cmin, add=yes, verify=no,show=no,update=yes)
				iraf.hedit(Vars.image, 'DTV2', Vars.lmin, add=yes, verify=no,show=no,update=yes)
				iraf.hedit(Vars.image,'CCDSUM,CCDSEC,AMPSEC,ATM1_1,ATM2_2,ATV1,ATV2',PYdel=yes,add=no,addonly=no,verify=no,show=no,update=yes)
			Vars.fd_ext = ''
			iraf.delete(Vars.extlist, verify=no)
		else:
			iraf.imrename(Vars.outsec, Vars.out, verbose=no)
		if (iraf.access(Vars.pllist)):
			iraf.delete(Vars.pllist, verify=no)
	Vars.fd_in = ''
	iraf.delete(Vars.inlists, verify = no)
	if (Vars.wcssource != 'match' and iraf.imaccess(Vars.wcsref)):
		iraf.imdelete(Vars.wcsref, verify=no)
		
mscimage()
