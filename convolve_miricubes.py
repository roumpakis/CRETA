#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 13:15:59 2021

@author: tanio
"""

import sys, glob, os
sys.path.append('/Users/tanio/Sync/pywork/pys')
import pdb
import numpy as np
from numpy import unravel_index
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
import mylmfit2dfun
from lmfit import Parameters


def convolve_miricubes():
    
    # Default directories
    rootdir = '/Users/tanio/lowz/jwst/jwebbinars/stage3/'
    inidir = rootdir+''
    outdir = rootdir+'analysis/'

    # Temporary arrays to lable sub-band cubes
    channames = ['ch1','ch2','ch3','ch4']
    nchans = len(channames)
    bandnames = ['short','medium','long']
    nbands = len(bandnames)
    
    # Default data to be passed is:
    # > A single sub-band cube and the following info on the cube: 
    #   > The wavelength array
    #   > The spaxel size

    # This for loops will not exist since only one sub-band cube is given
    # These are just to double check evreything is done correctly
    for channame in channames:
        for bandname in bandnames:
            
            # Read subcube fits file
            fitsfile = fits.open(rootdir+'Level3_'+channame+'-'+bandname+'_s3d.fits')
            # Science frame
            fitsfile[1].data[fitsfile[1].data == 0] = np.nan
            # Error frame
            fitsfile[2].data[fitsfile[2].data == 0] = np.nan
            # Dimensions
            nwaveinds, nys, nxs = fitsfile[1].data.shape
            
            # Create 2D mesh grid
            x, y = np.meshgrid(np.arange(nxs), np.arange(nys))
            # Stack them (x,y,2)
            #xy = np.stack((x,y))
            
            # Initialize the list where all the PSF sigmas will be stored
            psf_sigma_eff = []
            for waveind in range(nwaveinds):
                
                # Read slice from cube
                frame = fitsfile[1].data[waveind]
                frame_err = fitsfile[2].data[waveind]
                # Extract subcube 10x10 centered at the center of the image. This probably can be improved!
                subframe = np.asarray(frame[int(np.rint(nys/2))-10:int(np.rint(nys/2))+10, int(np.rint(nxs/2))-10:int(np.rint(nxs/2))+10])
                # Select the indices of the pixel where the maximum flux is and add the offset from the extracted subcube
                peakind = tuple(map(sum, zip(np.where(subframe == np.nanmax(subframe)), [int(np.rint(nys/2))-10, int(np.rint(nxs/2))-10])))
                # Peak flux at the indices
                peak = frame[peakind]

                if waveind == 0:
                    # Set up parameters only at the first wavelength in each subcube
                    params = Parameters()
                    # Initial peak of gaussian set to spaxel with maximum value within the central 20x20 spaxels
                    params.add('par0', value = peak[0], min=0., max=np.inf)
                    # Initial position of the peak
                    params.add('par1', value = peakind[1][0]) # X coordinate pixel where the peak of the PSF is
                    params.add('par2', value= peakind[0][0]) # Y coordinate pixel where the peak of the PSF is
                    params.add('par3', value=1.5)
                    params.add('par4', value=1.)
                    params.add('par5', value=np.pi/4., min=0., max=2.*np.pi)
                else:
                    # Otherwise use the parameters obtained from the previous frame fit
                    params = gauss2d_fit.lmfitres.params
                
                # Remove nanas and flatten arrays x,y and frame
                valinds = np.where(~np.isnan(frame))
                xflat = x[valinds]
                yflat = y[valinds]
                xyflat = np.stack((xflat,yflat))

                # Set up fitting class
                gauss2d_fit = mylmfit2dfun.mylmfit2dfun(xyflat, frame[valinds], params, zerr=frame_err[valinds])
                # Run ML fit to get the best fit
                gauss2d_fit.lmfit(tilt=True)
                
                # Add effective sigma to the list
                current_channel_pixel_scale = 0.3
                psf_sigma_eff.append(np.sqrt(gauss2d_fit.pars[3] * gauss2d_fit.pars[4]) * current_channel_pixel_scale) # Units of arcsecond
                print(psf_sigma_eff[waveind])

            print('#############')
            # Write the sigma_eff list to a file

            # Do all sub-bands and channels and finish the function

            # Start the loop in the channels that are going to be extracted

            # Read the sigma_eff list form the file
            index_of_the_last_channel_last_wave = nwaveinds-1 # placeholder
            this_channel_pixel_scale = 0.3 # placeholder
            for waveind in range(nwaveinds):
                
                # Read sclice from cube
                frame = fitsfile[1].data[waveind]
                frame_err = fitsfile[2].data[waveind]
                
                # This will not be needed in the implementation (see below)
                valinds = np.where(~np.isnan(frame))
                xflat = x[valinds]
                yflat = y[valinds]
                xyflat = np.stack((xflat,yflat))
                
                # Calculate the kernel sigma used to convolve as the sqrt of the subtractions of the final and current sigmas squared
                kernel_sigma = np.sqrt(psf_sigma_eff[index_of_the_last_channel_last_wave]**2 - psf_sigma_eff[waveind]**2) / this_channel_pixel_scale
                print(kernel_sigma * this_channel_pixel_scale)
                # Generate the kernel
                gauss2d_kernel = Gaussian2DKernel(x_stddev=kernel_sigma)
                
                # Convolve the frame and the error frame
                convframe = convolve(frame, gauss2d_kernel, normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)
                convframe_err = convolve(frame_err, gauss2d_kernel, normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)
                
                # This fits again the convolved image. NOT NEEDED IN THE IMPLEMENTATION.
                # This is just to check that the PSF sigma of the convolved image is similar to the target PSF sigma of the final channel and wavelength
                conv_gauss2d_fit = mylmfit2dfun.mylmfit2dfun(xyflat, convframe[valinds], gauss2d_fit.lmfitres.params, zerr=convframe_err[valinds])
                conv_gauss2d_fit.lmfit(tilt=True)
                
                # Print the sigma of the current frame, the fitted sigma of the convolved frame, and the target sigma
                print(channame+bandname, psf_sigma_eff[waveind], np.sqrt(conv_gauss2d_fit.pars[3] * conv_gauss2d_fit.pars[4]) * this_channel_pixel_scale, psf_sigma_eff[index_of_the_last_channel_last_wave]) # all in arcsec
                
                pdb.set_trace()

                
if __name__ == "__main__":
    out = convolve_miricubes()
