import numpy as np
import matplotlib.pylab as plt
from astropy.stats import sigma_clipped_stats
from astropy.io import fits as pyfits
from astropy.convolution import Gaussian2DKernel
from photutils.detection import detect_sources
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation
def bkg_boxes(frame,nboxes,length,sources=False):
	"""
	Function to calculate the sigma clipped statistics
	of a number of randomly generated boxes
	Variables:
	frame: fits image 
	nboxes: number of boxes to generate
	length: length of side of box in pixels
	sources: if sources = True, the sources in each box will be detected and masked
	       if sources = False, no masking is done 
	"""
	hdu = pyfits.open(frame,memmap=True)
	image = hdu[0].data
	
	#Get length of image in each axis
	naxis1 = float(hdu[0].header['NAXIS1'])
	naxis2 = float(hdu[0].header['NAXIS2'])
	
	#generate the centers of 1000 random boxes.
	#np.random.seed(1234)
	box_centers = np.random.random_integers(0,np.min([naxis1,naxis2]),size=(nboxes,2))
	
	#divide length by 2
	side = float(length)/2.0
	
	bg_stats = []
	centers = []
	for center in range(len(box_centers)):
		x1 = box_centers[center][0]-side
		x2 = box_centers[center][0]+side
		y1 = box_centers[center][1]-side
		y2 = box_centers[center][1]+side
	
		#Check to ensure that box is within image
		if (x1 > 0.0 and x2 < naxis1) and (y1 > 0.0 and y2 < naxis2):
		  centers.append(box_centers[center])
		  """
		  The centers that are within the image bounds are returned
		  in case you need to examine the regions used.
		  """      
		  box = image[x1:x2,y1:y2]
		  
		if (box >= 0).all() == True:
			"""
			Only boxes with non-negative values are kept.
			This should help deal with cell gaps
			The sigma and iter values might need some tuning.
			"""
			mean, median, std = sigma_clipped_stats(box, sigma=3.0, iters=10)
			if sources == False:
				bg_stats.append((mean, median, std))
			if sources == True:
				threshold = median + (std * 2.)
				segm_img = detect_sources(box, threshold, npixels=50)
				mask = segm_img.data.astype(np.bool)# turn segm_img into a mask
				selem = np.ones((15, 15))    # dilate using a 25x25 box
				mask2 = binary_dilation(mask, selem)
				mask_mean = np.mean(mask2)
				mean_mask, median_mask, std_mask = sigma_clipped_stats(box, sigma=3.0, mask=mask2)
				bg_stats.append((mean_mask, median_mask, std_mask))
	  
	bg_stats = np.reshape(np.array(bg_stats),(len(bg_stats),3))
	centers = np.reshape(np.array(centers),(len(centers),2))
	
	#Calculate median std of Background
	med = np.mean(bg_stats[:,1])
	#calculate standard deviation of the std values
	std = np.mean(bg_stats[:,2])
	
	#Locate the box that had the largest std
	#Array will be returned for plotting if wanted
	max_std = np.argmax(bg_stats[:,2])
	max_center = centers[max_std]
	max_box = image[max_center[0]-side:max_center[0]+side,max_center[1]-side:max_center[1]+side]
	#plt.imshow(max_box,origin='lower', cmap='Greys_r')
	#plt.show()
	
	return med,std,centers

# example of how to use function
# frame= 'gcp-f1tm.g.sh.fits'
# stats,med_std,std_std,centers,max_box = bkg_boxes(frame,500,200.0,sources=True)