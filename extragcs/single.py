from pyraf import iraf
from astropy.io import fits
import numpy as np
import os
import glob

os.putenv('iraf','/iraf/iraf')
iraf.stsdas(_doprint=1)

objectName = os.getcwd()[-5:].lower()

print objectName
# iraf.stsdas.analysis.gasp.xyeq()

os.putenv('iraf','/iraf.216/iraf')
iraf.stsdas(_doprint=1)

os.putenv('iraf','/iraf/iraf')
iraf.stsdas(_doprint=1)
