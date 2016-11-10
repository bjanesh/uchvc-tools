apphot
unlearn phot
unlearn datapars
unlearn photpars
unlearn centerpars
unlearn fitskypars
#
phot.interactive=no
phot.verify=no
# I know that all of the aperture-correction stars are not saturated,
#so I can leave datamax as INDEF
#datapars.datamax=INDEF
datapars.gain="gain"
datapars.ccdread="rdnoise"
datapars.exposure="exptime"
datapars.airmass="airmass"
datapars.filter="filter"
datapars.obstime="time-obs"
datapars.sigma=INDEF
photpars.zmag=0.
centerpars.cbox=9.
centerpars.maxshift=3.
fitskypars.salgorithm="median"
fitskypars.dannulus=10.
#
phot image="AGC198606_g_sh.fits" coords="g_apcor_stars.dat" datapars.fwhmpsf=6.6 photpars.apertures="6,9,12,15,18,21,24,27,30,33,36" fitskypars.annulus=39 output="g.apcor.mag.1"
#
phot image="AGC198606_i_sh.fits" coords="i_apcor_stars.dat" datapars.fwhmpsf=6.95 photpars.apertures="7,10,14,17,21,24,28,31,35,38,42" fitskypars.annulus=45 output="i.apcor.mag.1"





