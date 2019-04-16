import numpy as np
import os
from astropy.io import fits
import aplpy
from PIL import Image

targets = ["AGC174540", "AGC198511", "AGC198606", "HI0959+19", "HI1037+21", "HI1050+23", "AGC215417", "HI1151+20", "AGC226067", "AGC227987", "AGC229326", "AGC238626", "AGC238713", "AGC249000", "AGC249282", "AGC249320", "AGC249323", "AGC249525", "AGC258237", "AGC258242", "AGC258459", "AGC268069", "AGC268074"]

path = os.getcwd()

for agc in targets:
    folder = agc.lower()
    img_g = path+'/'+folder+'/'+agc+'_g.fits'
    img_i = path+'/'+folder+'/'+agc+'_i.fits'
      
    fits_g = fits.open(img_g)
    fits_i = fits.open(img_i)

    # first make the green channel from an average of the g and i images
    if not os.path.isfile(path+'/'+folder+'/'+agc+'_avg.fits'):
        newheader = fits_i[0].header
        newimage = (fits_g[0].data+fits_i[0].data)/2.0
        newhdu = fits.PrimaryHDU(data=newimage,header=newheader)
        newhdu.writeto(path+'/'+folder+'/'+agc+'_avg.fits') 
    # use aplpy to generate the actual png image  
    aplpy.make_rgb_image([img_i,path+'/'+folder+'/'+agc+'_avg.fits',img_g],path+'/'+folder+'/'+agc+'_rgb.png')
    if not os.path.isfile(path+'/'+folder+'/'+agc+'_rgb_thumb.png'):
        thumb_size = 1024,1024
        im = Image.open(path+'/'+folder+'/'+agc+'_rgb.png')
        im.thumbnail(thumb_size)
        im.save(path+'/'+folder+'/'+agc+'_rgb_thumb.png', "PNG")
# print os.listdir('.')         

# fig = mpl.figure(figsize=(10,10))

# f24 = aplpy.FITSFigure("agc268069/20130414T012307.1_AGC_268069_odi_g.7019.fits", hdu='OTA24.SCI', figure=fig, subplot=(3,3,1))
# f24.show_grayscale(vmin=145.0, vmax=285.0)
# f24.set_theme('publication')
# f24.frame.set_color('none')
# f24.axis_labels.hide()
# f24.tick_labels.hide()
# f24.ticks.hide()


