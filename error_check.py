# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:18:44 2021

@author: roub
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
from astropy.wcs import WCS
from matplotlib.pyplot import loglog
from photutils.centroids import centroid_com
from photutils.centroids import centroid_1dg, centroid_2dg
from photutils.aperture import CircularAperture, CircularAnnulus
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy import units as u
from sklearn.metrics import mean_squared_error as MSE

from photutils.centroids import centroid_1dg, centroid_2dg
from photutils.aperture import RectangularAperture
import matplotlib.pyplot as plt
data_arr = np.array([[1, 2], [3,4]])
error_arr = np.array([[1, 2], [3,4]])
aper = RectangularAperture([0.5,0.5], 1, 1)

plt.imshow(data_arr)
aper.plot()
plt.show()


phot_table = aperture_photometry(data_arr, aper, error=error_arr)
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
print(phot_table)

