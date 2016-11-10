apphot
unlearn phot
unlearn datapars
unlearn photpars
unlearn centerpars
unlearn fitskypars
#
phot.interactive=no
phot.verify=no
datapars.datamin=INDEF
datapars.datamax=50000.
datapars.gain="gain"
datapars.ccdread="rdnoise"
datapars.exposure="exptime"
datapars.airmass="airmass"
datapars.filter="filter"
datapars.obstime="time-obs"
datapars.sigma=INDEF
photpars.zmag=0.
centerpars.calgorithm="centroid"
centerpars.cbox=9.
centerpars.maxshift=3.
fitskypars.salgorithm="median"
fitskypars.dannulus=10.
#
# Use an aperture that is 1 x <fwhm>, because an aperture correction
# will be applied in the calc_calib_mags step
# Using a sky annulus that begins at 6 x <fwhm> should be fine
# g-band
phot image="AGC198606_g_sh.fits" coords="escut_g.pos" datapars.fwhmpsf=6.6 photpars.apertures="6" fitskypars.annulus=36 output="agc198606_sources_g.mag.1"
# i-band
phot image="AGC198606_i_sh.fits" coords="escut_i.pos" datapars.fwhmpsf=6.95 photpars.apertures="7" fitskypars.annulus=42 output="agc198606_sources_i.mag.1"
