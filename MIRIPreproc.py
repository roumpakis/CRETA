# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 13:06:50 2021

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
from photutils.aperture import RectangularAperture
from photutils.centroids import centroid_1dg, centroid_2dg
import os
import glob
current_path = os.path.abspath(os.getcwd())
from photutils.aperture import CircularAperture
class MIRIPreproc:
    def __init__(self):
        print()
        
        
        
    def getTargetR(self,base_l, base_r, target_l,pixel_scale,base_pixel_scale):
        return base_r * (target_l/base_l) * (base_pixel_scale/pixel_scale)


    #
    def getApertureSum(self,x,y,z,r, data,w,error, aprture_type, DQ):
       
        apers = []
        dq_slice = []
        for i in range(len(r)):
            #get the rectangular or circular aperture
            if aprture_type == 0:
                aper = CircularAperture([x,y], r[i])
            else:    
                aper = RectangularAperture([x,y], r[i], r[i])
            apers.append(aper)

        phot_table = aperture_photometry(data, apers, wcs=w,error = error)
        rs_photometry = []
        rs_error = []
        rs_area = []
        for i in range(len(r)):
            rs_photometry.append(phot_table["aperture_sum_"+str(i)][0])
            rs_error.append(phot_table["aperture_sum_err_"+str(i)][0])
            rs_area.append(apers[i])
        return [rs_photometry,rs_error,rs_area]
   
    # calculates the wavelength 
    #CDELT3 pixel scale of z axe
    def getSimulatedL(self,CRVAL3,CDELT3 ,i):
        return CRVAL3 + (i) * CDELT3


    #Load a subcube 
    def getFITSData(self,image_file):
        print('Load file: '+image_file)
        hdu_list = fits.open(image_file)
        # hdu_list.info()
        SCI = hdu_list['SCI']
        ERR = hdu_list['ERR']
        DQ = hdu_list['DQ'].data
        err_data = ERR.data
        CRPIX3 = hdu_list['SCI'].header['CRPIX3'] -  1
        CRVAL3 = hdu_list['SCI'].header['CRVAL3']
        CDELT1 = hdu_list['SCI'].header['CDELT1'] * 3600 #arcsec / pixel
        CDELT2 = hdu_list['SCI'].header['CDELT2'] * 3600 #arcsec/ pixel
        CDELT3 = hdu_list['SCI'].header['CDELT3'] #um / 'pixel'
        pixelScale = hdu_list['SCI'].header['CDELT1'] * 3600 # arcsec / pixel
        if hdu_list['PRIMARY'].header['INSTRUME'] == 'NIRSPEC':
                cube_name = hdu_list['PRIMARY'].header['GRATING'] 
        else:
                cube_name = 'ch_'+ hdu_list['PRIMARY'].header['CHANNEL'] + '_'+hdu_list['PRIMARY'].header['BAND']
        output_file_name =   hdu_list['PRIMARY'].header['OBS_ID']
        hdu_list.close()
        
        image_data = fits.getdata(image_file)
        headers = SCI.header
        headers_txt = repr(headers)
        headers_list= headers_txt.split("\n")
        
        primaryDict = {}
        for i in range(len(headers_list)):
            line = headers_list[i]
            if line.find('=') != -1:
                split1 = line.split("=")
                kkey = split1[0]
                kkey = "".join(kkey.split())
                
                split2 = split1[1].split('/')        
                vvalue = split2[0]
                vvalue = "".join(vvalue.split())
                primaryDict[kkey]=vvalue 
                
        res = []
        res.append(image_data)
        res.append(primaryDict)
        res.append(headers)
        res.append(CRPIX3)
        res.append(CRVAL3)
        res.append(CDELT3)
        res.append(pixelScale)
        res.append(err_data)
        res.append(CDELT1)
        res.append(CDELT2)  
        res.append(cube_name)
        res.append(output_file_name)
        res.append(DQ)
        return res
         
    def getChannelLs(self,subcube):
        all_li = []

        for i in range(subcube.image_before.shape[0]):
            l_i =  self.getSimulatedL(subcube.CRVAL3,subcube.CDELT3 ,i)
            all_li.append(l_i)

        return all_li
    
         #%%   
    def getPointSourceRs(self,subcube,base_r):
        all_ri = []

        for i in range(subcube.image_before.shape[0]):
            
            r_i = self.getRpixPointSource(base_r,subcube.pixel_scale, subcube.ls[i],subcube.base_l)
            all_ri.append(r_i)

        return all_ri
    
             #%%   
    def getExtendedSourceRs(self,subcube, base_r):
        all_ri = []

        for i in range(subcube.image_before.shape[0]):
            r_i = self.getRpixExtendedSource(base_r,subcube.pixel_scale)
            all_ri.append(r_i)

        return all_ri


#%%
    ##### Function for PSF sub-channel centering. Use a sub-image in order to avoid bad pixels#####
    ###############################################################################  
    # @subcube: The corresponding PSF sub-channel. (SubeCube)
    # @image: The data that we will use, with or without background subtraction based on user option. (np.array)
    ########### --> Return res   ############################
    # @res: A list with resulting centroids in pixels. (list)   
    ###############################################################################
    def PSFCentering(self,subcube,image):
        res = []
        for i in range(image.shape[0]):
            sliceIm = image[i,:,:]
            jj, kk= self.imageCentroid(sliceIm,subcube.bkg_y,subcube.bkg_x)
            res.append([jj,kk])
        return res   


    def PSFCenteringSky(self,subcube):
        res = []
        for i in range(len(subcube.ls)):
            [jj,kk] = subcube.xys[i]
            
            sky = subcube.wcs.pixel_to_world(jj,kk,subcube.ls[i])
            res.append(sky)
        return res 
#%% 
    ##### Function that creates a list centroids   #####
    ###############################################################################  
    # @subcube: Sub-channel for centroids calculation. (SubeCube)       
    # @image: The data that we will use, with or without background subtraction based on user option. (np.array)
    ########### --> Returns [x,y]   ############################
    # @y: Centroid Y coordinate. (int)   
    # @x: Centroid X coordinate. (int)  
    ###############################################################################
    def DataUsersCentering(self,subcube,image):
        res = []
        for i in range(image.shape[0]):
           # defaults to      
            c1 = SkyCoord(subcube.user_ra,subcube.user_dec, unit="deg")  
            x,y,z = subcube.wcs.world_to_pixel(c1, subcube.ls[i]*u.um)
            
            res.append([x,y])
        return res    
        # %% Aperture Photometry
    def AperturePhotometry(self,subcube,image):
        aperture_sum = []
        all_apers = []
        aper_elem = []
        all_area = []
        all_error = []
        for i in range(image.shape[0]):
            # print('BAND: ',subcube.name_band , ' PLANE : ', i)
            sliceIm = image[i,:,:]
            jj, kk= subcube.xys[i]

            [photometries, errors,areas] = self.getApertureSum(jj,kk,i,subcube.rs[i],sliceIm,subcube.wcs,subcube.error_data[i,:,:],subcube.aperture_type, subcube.DQ)
            
            aperture_sum.append(photometries)  
            all_area.append(areas)
            all_error.append(errors) 
        return [aperture_sum, all_area,all_error]
    


  
    
    #%% 
    ##### Calculate the image Centroid using an 11x11 box at the 'midle'of the image  #####
    ###############################################################################  
    # @image: The data that we will use, with or without background subtraction based on user option. (np.array)
    # @xx: Cordinate of X-axe, used as sub-image center. (float)
    # @yy: Cordinate of Y-axe, used as sub-image center. (float)    
    ########### --> Returns [column, rows]   ############################
    # @column: Centroid Y coordinate. (int)   
    # @row: Centroid X coordinate. (int)  
    ###############################################################################
    def imageCentroid(self, image,xx,yy):
           NY = int(image.shape[0]/2)
           NX = int (image.shape[1]/2)
           
            #  # 25% to 75% sub Image
           img = image.copy()
           start_X = int(NX/2)
           start_Y = int (NY/2)
           subImg = img[start_Y:start_Y+NY, start_X:start_X+NX]
           # print('subImg ',subImg.shape)           
           #the 11X11 image
           round_yy = int(yy)
           round_xx = int(xx)
           zoom2img = subImg[round_yy-5: round_yy+6, round_xx -5:round_xx  +6]
           #FLAG
           
           #if there are NaNs within the sub-image do not center
           if  len(np.where(zoom2img  == np.NaN)) != 0 : 
                 row = int(round_xx)
                 column= int(round_yy)
                 row = row+start_Y+ round_yy -5 
                 column = column+start_X+round_xx - 5
                 return [column,row]                 


           else:
               columns, rows= centroid_2dg(zoom2img) 
               # print(rows,columns)
               rows = rows+start_Y+ round_yy -5 
               columns = columns+start_X+round_xx - 5
               return [columns,rows]


    #%%   
    
    def totalImageFlux(self, image):

      res = []
      flux = 0
      for i in range(len(image)):
          flux = np.nansum(image[i])
          res.append(flux)
      return res      
  
  #%%% Calculate the PSF infinite aperture as the total of the elements 
  # @image:    the 3-D image cube of PSF
 ###
    def PSFInfFlux(self, image, delta_factor):

      res = []
      flux = 0
      for i in range(len(image)): 
          flux = np.nansum(image[i])
          res.append(flux)
      # res = list(np.array(res) * 10**6 * (delta_factor/206265**2))    
      return res      
  
   #%% Get the radius in pixels 
   # @r_ap: the user radius in arc sec
   # @ps: the channel's pixel scale
   # @l_i: the target lambda
   ##
  
    def getRpixPointSource(self,r_ap, ps , l_i, l_ap):
        return (np.array(r_ap)/ps) * (l_i/l_ap)
    
   #%% Get the radius in pixels 
   # @r_ap: the user radius in arc sec
   # @ps: the channel's pixel scale
   # @l_i: the target lambda
   ##
    
    def getRpixExtendedSource(self,r_ap, ps):
        return (np.array(r_ap)/ps)  
    
  #%% Infinite Correction Point Source 
    def PSFPointSourceCorrection(self, image):
        psf_inf = self.PSFInfFlux(image)
        
        
#%%
    def subtractUserBackground(self, subcube, r_in, r_out):
        res = []
        annulus = []
        annulus_aperture_list = []
        annulus_centroid = []
        res_rout = []
        anImg = subcube.image_before.copy()
        
        for z in range(len(anImg)):
           
            img = anImg[z,:,:]
            j,k= subcube.xys[z]
            annulus_aperture = CircularAnnulus([j,k], r_in=r_in[z], r_out=r_out[z])
            # if r_in[z] > subcube.rs[z]:
            #     print("=== WARNING ===  "+subcube.name_band+" [ lambda"+str(subcube.ls[z])+"] Annulus inner r("+str(r_in[z])+") is greater than aperture("+str(subcube.rs[z])+")")
            annulus_aperture_list.append(annulus_aperture)
            annulus_masks = annulus_aperture.to_mask(method='center')
            annulus_data = annulus_masks.multiply(img)
            
            ww = np.where(annulus_data != 0)
            annulus_data_1d = annulus_data[ww]
            mask2 = np.where(~np.isnan(annulus_data_1d))
            annulus_data_1d = annulus_data_1d[mask2] #exclude the NaN 
            mean, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            anImg[z,:,:] = anImg[z,:,:]  -  median_sigclip
            annulus.append(annulus_masks)
            annulus_centroid.append([j,k])
            res.append(median_sigclip)
            res_rout.append(r_out)

        return  [anImg,res,annulus,annulus_centroid,annulus_aperture_list,res_rout]             
        
    
#%%
    def addMaxValue(self,image, prece):
            for z in range(len(image)):
                maxV = np.nanmax(image[z])
                maxAdd = prece * maxV / 5
                image[z] = image[z]+maxAdd
            return image

       #%% Calculate the image Centroid using an 11x11 box at the 'midle'of the image 
   # @image: the 2D input image 
   # @xx:    coordinate of 11x11 of max flux
   # @yy:    coordinate of 11x11 of max flux
   ###
    def userCentroid(self, image,xx,yy):
           
         
           #the 11X11 image
           round_yy =  int(yy)
           round_xx = int(xx)
           
           zoom2img = image[round_yy - 5 : round_yy+6, round_xx-5 : round_xx+6]
           # print('ROUNDS ', round_xx,  '  ', round_yy)
           # if round_xx <= 0 or round_yy<=0  or round_yy>zoom2img.shape[1] or round_xx > zoom2img.shape[0]:
           #     print("EIMASTE STO IF , ",yy,'  ', xx )
           #     return[yy,xx]
           # elif   len(np.where(zoom2img  == np.NaN)) != 0 : 
           #       print("EIMASTE STO elIF , ",yy,'  ', xx ) 
           #       return[yy,xx]
           columns, rows= centroid_2dg(zoom2img) 
           # print('C: ', columns, ' R: ', rows)
           #fix the 11x11 coordinatesr
           rows = rows+round_yy-5
           columns = columns+round_xx-5

           
           return [columns,rows]

#%%
    ##### Function for centering at a specific wavelength that user defines, using 3 slices from data array.     #####
    ###############################################################################  
    # @images: List of all avaliable sub-channels. (list of SubeCube)
    # @l_c: Wavelength used for centering. (float)
    # @RA: The user defined RA before centering. (arcsec)
    # @dec: The user defined dec before centering. (arcsec)    
    ########### --> Returns sky   ############################
    # @sky: The sky coordinates after centering at wavelenth l_c. (SkyCoord)   
    ###############################################################################      
    def lambdaBasedCentering(self,images, l_c, RA,dec, dxdy = False, theImage = []):
           res_cubes = []
           ls_min = []
           ls_max = []
           found = False
           for i in range(len(images)):
               ls_min.append(images[i].ls[0])
               ls_max.append(images[i].ls[len(images[i].ls)-1])
                # print(images[i].name_band, "  ", ls_min, " ", ls_max, " l_c: ", l_c)
               if l_c >= ls_min[i] and l_c<= ls_max[i]:
                   res_cubes.append(i)     
                   found = True
                   
           c1 = SkyCoord(RA,dec, unit="deg")  # defaults to   
           the_image = images[res_cubes[0]] 
           # print('Centering with ', the_image.name_band)             
           x,y,z = the_image.wcs.world_to_pixel(c1, l_c*u.um)
           
           if dxdy:
               # print( the_image.image_before.shape)
               the_image = theImage
               z,y,x = the_image.image_before.shape
               x = x/2
               y = y/2
               sky_list = []
               res_cubes_all = []
               for i in range(len(the_image.image_before)):
                   plane = the_image.image_before[i,:,:]
                   jj, kk= self.userCentroid(plane,x,y)
                   # print('NEw Center in pixels: ', jj, '  ',kk)
                   sky = the_image.wcs.pixel_to_world(jj,kk,l_c)
                   sky_list.append(sky)
                   res_cubes_all.append( res_cubes[0] )
               return sky_list,res_cubes_all            
               # print('KANOUME TO DEFAULT RE MNIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaaa')
           else:    
               x = x.tolist()
               y = y.tolist()
               z = round(z.tolist())
           plane = the_image.image_before[z,:,:]
           if z-1>=0 and z+1<len(the_image.image_before):
               plane  = plane + the_image.image_before[z-1,:,:] + the_image.image_before[z+1,:,:]  #add one before and one plane after l_c
           elif     z-1<0:
                   plane  = plane + the_image.image_before[z+1,:,:] + the_image.image_before[z+2,:,:]  #add two next  plane after l_c
           elif     z+1>=len(the_image.image_before):
                   plane  = plane + the_image.image_before[z-1,:,:] + the_image.image_before[z-2,:,:]  #add two previous  plane before l_c

           # plt.imshow(plane)
           jj, kk= self.userCentroid(plane,x,y)
           # print('NEw Center in pixels: ', jj, '  ',kk)
           sky = the_image.wcs.pixel_to_world(jj,kk,l_c)
           return sky, res_cubes[0]

  #%% Get lambdas overlappinf correction 
    def getLambdasOverlappingCorrection(self, ch1_data,ch1_ls,ch2_data,ch2_ls,delta1,delta2):
        print('Stitching channels and bands')
        ch1_start = np.where(np.array(ch1_ls)>= ch2_ls[0] - delta2/2 )[0][0]
        ch2_stop =  np.where(np.array(ch2_ls)<= ch1_ls[len(ch1_ls)-1]+ delta1/2)[0]
        ch2_stop =ch2_stop[len(ch2_stop)-1]                                                                     
        ch1_overlapping_before = ch1_data[ch1_start:]
        ch1_over_ls= ch1_ls[ch1_start:]
        ch2_overlapping = ch2_data[:ch2_stop+1]
        # print(len(ch1_overlapping_before),len(ch2_overlapping))

        
        ch1_mean = np.mean(ch1_overlapping_before)
        ch2_mean = np.mean(ch2_overlapping)
        ratio = ch2_mean /ch1_mean
        # print(ratio,ch1_mean ,ch2_mean)
        ch1_fixed = []
        
        ch2_over_ls =  ch2_ls[:ch2_stop+1]
        for i in range(len(ch1_data)):
            ch1_fixed.append(ch1_data[i]*ratio)
        ch1_overlapping = ch1_fixed[ch1_start:]  


        # print(ch1_over_ls[0])
        # print(ch2_over_ls[0])
        # plt.plot(ch1_over_ls,ch1_overlapping_before, 'o' ,markersize=1, color='black',label='1st Subchannel Overlapping Part Before Scaling')
        # plt.plot(ch1_over_ls,ch1_overlapping,  'o' ,markersize=1, color='red', label='1st Subchannel Overlapping Part After Scaling')
        # plt.plot(ch2_over_ls,ch2_overlapping,  'o' ,markersize=1, color='green', label='2nd Subchannel Overlapping Part')
        # plt.xlabel('λ(μm)')
        # plt.ylabel('Flux')
        # plt.legend()
        # plt.show()        
        # # print('MSE = ',MSE(ch1_over_ls,ch2_over_ls))
        # plt.plot(ch1_ls,ch1_data, 'o' ,markersize= 0.5, color='black',label='Before Scaling')
        # plt.plot(ch1_ls,ch1_fixed,  'o' ,markersize=0.5, color='red', label='After Scaling')
        # plt.plot(ch2_ls,ch2_data,  'o' ,markersize=0.5, color='green', label='Following Sub-Channel')
        # plt.xlabel('λ(μm)')
        # plt.ylabel('Flux')
        # plt.legend()
        # plt.show()

        res_data = []
        res_ls = []
        # if len(ch1_overlapping) == len(ch2_overlapping):
        # print('fiiix')
        # res_data.extend(ch1_fixed[:ch1_start])
        # res_ls.extend(ch1_ls[:ch1_start])
        # for i in range(len(ch1_overlapping)):
        #         res_data.append((ch1_overlapping[i]+ch2_overlapping[i])/2)
        #         res_ls.append(ch1_over_ls[i])
        # res_ls.extend(ch2_ls[ch2_stop+1:])    
        # res_data.extend(ch2_data[ch2_stop+1:])
        # return [res_ls,res_data]

        # else:
            

        import pandas as pd
        res_all = []
        for i in range(len(ch1_fixed)):
                res_all.append([ch1_ls[i],ch1_fixed[i]])
                
        for i in range(len(ch2_data)):
                res_all.append([ch2_ls[i],ch2_data[i]])            
        df = pd.DataFrame(res_all, columns = ['ls', 'Flux']) #add here everithing
        df = df.sort_values(by=['ls'])

        return [list(df.ls),list(df.Flux),ratio]
       
    def fixSpectrumLambdas(self,all_data,all_ls,all_delta):
        print('Stitching band and channels....')
        ch1_data = all_data[0]
        ch1_ls = all_ls[0]
        tanio = []
      
        for i in range(len(all_data)-1):

            ch2_data = all_data[i+1]
            ch2_ls = all_ls[i+1]
            [ls,apers,ratio ] =self.getLambdasOverlappingCorrection(ch1_data,ch1_ls,ch2_data,ch2_ls,all_delta[i],all_delta[i+1])
            ch1_data = apers
            ch1_ls = ls
            tanio.append(ratio)
  
        return [ls,apers,tanio ]
    #%%LOAD 1D EXTRACTION
    def Load1DFile(self,filename):
        import math
        hdu_list = fits.open(filename)
        data = hdu_list['EXTRACT1D'].data
        lambdas = data[:]['WAVELENGTH']
        flux = data[:]['FLUX']
        error =  data[:]['ERROR']
        bright = data[:]['SURF_BRIGHT']
        background = data[:]['BACKGROUND']
        backgroundError =  data[:]['BERROR']
        area =  data[:]['NPIXELS']
        r = np.sqrt(area / (math.pi))
        # fits.close()
        return [lambdas,flux,error,bright,background,backgroundError,area,r]

    def LoadAll1D(self,files):
        data_all = []
        lambdas_all = []
        flux_all = []
        error_all = []
        bright_all = []
        background_all = []
        backgroundError_all = []
        area_all = []
        r_all = []
   
        for i in range(len(files)):
            [lambdas,flux,error,bright,background,backgroundError,area,r] = self.Load1DFile(files[i])
            lambdas_all.extend(lambdas)
            flux_all.extend(flux)
            error_all.extend(error)
            bright_all.extend(bright)
            background_all.extend(background)
            backgroundError_all.extend(backgroundError)
            area_all.extend(area)
            r_all.extend(r)
        
        return [lambdas_all,flux_all,error_all,bright_all,background_all,backgroundError_all,area_all,r_all]   
#%% 
    def getSubcubesAll(self,subcubes,background, aperture_correction):
        all_rs_arcsec = []
        all_ls = []
        all_apers = []
        all_xys = []
        all_area_pix = []
        all_bright = []
        all_error_spectrum =[]
        all_corrected_spectrum = []
        all_delta = []
        all_background = []
        all_names = []
        all_error_corrected = []
        all_unit_ratio = []
        all_r_in = []
        all_rs = []
        all_ps = []
        all_psc_flux = []
        all_psc_err= []
        for i in range(len(subcubes)):
            all_rs_arcsec.extend(subcubes[i].rs_arcsec)
            all_rs.extend(subcubes[i].rs)
            all_ls.extend(subcubes[i].ls)
            all_apers.extend(subcubes[i].apers)
            all_xys.extend(subcubes[i].xys)
            all_error_spectrum.extend(subcubes[i].error)
            all_corrected_spectrum.extend(subcubes[i].corrected_spectrum)
            
            
            all_delta.extend(subcubes[i].CDELT3L)
            all_names.extend(subcubes[i].nameL)
            all_unit_ratio.extend(subcubes[i].unit_ratio)
            
            ps_list = [subcubes[i].pixel_scale] * len(subcubes[i].rs)
            all_ps.extend(ps_list)
            
            if (background): 
                all_background.extend(subcubes[i].background_spectrum)
                all_r_in .extend(subcubes[i].bckg_rs)    
            if(aperture_correction):
                all_psc_flux.extend(subcubes[i].spectrum_PSF_corrected)
                all_psc_err.extend(subcubes[i].error_PSF_corrected)
        return [all_rs_arcsec,all_ls,all_apers,all_xys,all_area_pix,all_bright,all_error_spectrum,\
                all_corrected_spectrum,all_delta,all_names,all_unit_ratio,all_background,all_r_in,all_rs,all_ps,all_psc_flux,all_psc_err]    
    
    
    #%%    
    def getSubcubesAllAppended(self,subcubes,background):
        all_rs = []
        all_ls = []
        all_apers = []
        all_xys = []
        all_area_pix = []
        all_bright = []
        all_error_spectrum =[]
        all_corrected_spectrum = []
        all_delta = []
        all_background = []
        all_names = []
        
        for i in range(len(subcubes)):
            all_rs.append(subcubes[i].rs)
            all_ls.append(subcubes[i].ls)
            all_apers.append(subcubes[i].apers)
            all_xys.append(subcubes[i].xys)
            # all_area_pix.append(subcubes[i].area_pix)
            # all_bright.append(subcubes[i].bright)
            all_error_spectrum.append(subcubes[i].error)
            all_corrected_spectrum.append(subcubes[i].corrected_spectrum)
            
            all_delta.append(subcubes[i].CDELT3)
            all_names.extend(subcubes[i].name_band)
            if (background) : all_background.append(subcubes[i].background_spectrum)
            
        return [all_rs,all_ls,all_apers,all_xys,all_area_pix,all_bright,all_error_spectrum,all_corrected_spectrum,all_delta,all_names,all_background]     
    #%%
    def listMJSR2Jy(self,data, ratio):
        res = []
        for i in range(len(data)):
            res.append(data[i] ** 10** 6 * (ratio[i] / 206265**2))
            
        return res    
    

    #%%Grid in Arcseconds 
    def crateGridInArcSec(self, user_ra, user_dec, gridPoints_dist, gridPointsX, gridPointsY, cube, r, subchannels, pointSource, l_ap):
        NX = np.arange(0,gridPointsX)
        NY = np.arange(0,gridPointsY)
        
        # r_pix = gridPoints /2 
        gridPoints_pix = gridPoints_dist / cube.pixel_scale
        if r == -1:
            r = gridPoints_dist/2
            r_pix = ((gridPoints_pix/2)) 
        else:    
            r_pix = r / cube.pixel_scale
        c1 = SkyCoord(user_ra,user_dec, unit="deg")  # defaults to      
        user_x,user_y,user_z = cube.wcs.world_to_pixel(c1, cube.ls[0]*u.um)
        grids_xs = user_x +(NX - (gridPointsX-1)/2) * gridPoints_pix
        grids_ys = user_y +(NY - (gridPointsY-1)/2) * gridPoints_pix
        
        
        
        sky_list = []
        pixels_list = []
        coord_grid = []
        names = []
        sky_ra = []
        sky_dec = []
        for i in range(len(grids_xs)):
            for j in range(len(grids_ys)):
                sky = cube.wcs.pixel_to_world(grids_xs[i],grids_ys[j],0)
                coord_grid.append(sky)   
                sky_list.append(sky)
                pixels_list.append([i,j])
                names.append(str(i)+"_"+str(j))
                sky_ra.append(sky[0].ra)
                sky_dec.append(sky[0].dec)
                
        # for i in range(len(subchannels)):                
        #     self.plotGridSubchanel( user_ra, user_dec, gridPoints_dist, gridPointsX, gridPointsY, subchannels[i], r)
        params_path = current_path+"\\Params"
        # self.delteFilesatPath(params_path)
        # self.writeParamsFiles(coord_grid,r,l_ap,pointSource)        
        
        return sky_list,pixels_list,names, sky_ra, sky_dec
        
#%%   
    def plotGridSubchanel(self, user_ra, user_dec, gridPoints_dist, gridPointsX, gridPointsY, cube, r):
        NX = np.arange(0,gridPointsX)
        NY = np.arange(0,gridPointsY)
        from matplotlib.patches import Rectangle
        
        gridPoints_pix = gridPoints_dist / cube.pixel_scale
        if r == -1:
            r = gridPoints_dist/2
            r_pix = ((gridPoints_pix/2)) 
            # print("EXOUME grid_points: ", gridPoints_pix, " r: ", r_pix)
        else:    
            r_pix = r / cube.pixel_scale
            # print("EXOUME grid_points: ", gridPoints_pix, " xeirokinhto r: ", r_pix)
        c1 = SkyCoord(user_ra,user_dec, unit="deg")  # defaults to      
        user_x,user_y,user_z = cube.wcs.world_to_pixel(c1, cube.ls[0]*u.um)
        grids_xs = user_x +(NX - (gridPointsX-1)/2) * gridPoints_pix
        grids_ys = user_y +(NY - (gridPointsY-1)/2) * gridPoints_pix
        
        sky_list = []
        pixels_list = []
        coord_grid = []
        names = []
        for i in range(len(grids_xs)):
            for j in range(len(grids_ys)):
                sky = cube.wcs.pixel_to_world(grids_xs[i],grids_ys[j],0)
                coord_grid.append(sky)   
                sky_list.append(sky)
                pixels_list.append([ grids_xs[i],grids_ys[j] ])
                names.append(str(i)+"_"+str(j))
        

        img = cube.image_before[0,:,:]
        for i in range(1,len(cube.image_before)):
                    img = img + cube.image_before[i,:,:]
        
        plt.imshow(img, origin='lower')
        plt.plot(user_x,user_y, 'o',color="red", label="User Input Centroid")
        for i in range(len(pixels_list)):
            # xx = pixels_list[i][0] - (r_pix/2)
            xx = pixels_list[i][0] - r_pix
            yy = pixels_list[i][1] - r_pix
            # yy = pixels_list[i][1] - (r_pix/2)
            plt.gca().add_patch(Rectangle([xx,yy],2*r_pix,2*r_pix,linewidth=1,edgecolor='r',facecolor='none'))
            plt.plot(pixels_list[i][0],pixels_list[i][1], 'bo', markersize=3)
            plt.title(cube.name_band)
        plt.legend()    
        plt.show()          
         
        return sky_list,pixels_list,names            
            
        
        

    
    #%%
    def writeParamsFiles(self,sky_list,user_r_ap,lambda_ap, pointSource):
       print('Ok prepei na nai')
       print(repr(sky_list[0]))
       for i in range(len(sky_list)):
            f = open("Params//params_"+str(i)+".csv", "w")
            f.write('user_r_ap = '+str(user_r_ap)+"\n" )
            f.write('user_ra = '+str(sky_list[i][0].ra) +"\n" ) 
            f.write('user_dec = '+str(sky_list[i][0].dec) +"\n" )
            f.write('point_source = '+str(pointSource)+"\n" )
            f.write('lambda_ap = '+str(lambda_ap)+"\n" )
            f.write('aperture_correction = '+str(False)+"\n" ) 
            f.write('centering = '+str(False)+"\n" ) 
            f.write('lambda_cent = '+str(4.89049986650)+"\n" ) 
            f.write('background_sub  = '+str(False)+"\n" )
            f.write('r_ann_in = '+str(1.23)+"\n" )
            f.write('ann_width = '+str(1.23)+"\n" )
            f.write('PSF_path  = C:\\Users\\roub\\Desktop\\finale\\PSF\\'+"\n" )
            f.write('data_path   = C:\\Users\\roub\\Desktop\\finale\\Data\\'+"\n" )
            f.write('output_path   =C:\\Users\\roub\\Desktop\\finale\\Results\\'+"\n" )
            
    #%%
    def getApertureDQList(self,cube):
       # print('DQ List ')
       res = []
       for i in range(len(cube.DQ)):
               # for j in range(len(cube.area[i])):
                   aper = cube.area[i][0]
                   # print("DQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ AAAAAAAAAAAA A A A A ", cube.area[i])
                   # aperstats2 = ApertureStats(cube.DQ[i,:,:], aper)
                   mask = aper.to_mask()
                   dq_masked = mask.cutout(cube.DQ[i,:,:])
                   dqv = np.max(dq_masked)
                   # print(aperstats2.max) 
                   res.append(dqv)
       return res     
   #%%
   
   
   ##### Function that calculates the stitching ratio between all posible sub - channels ##### 
   #@realData: list of data SubCube elements 
   #@aperture_correction: boolean value 
    def stichingRatioCalculation(self, realData,aperture_correction, idx, grid):
        # 'G140H',
        cubesNames = [ 'G140H', 'ch_1_SHORT','ch_1_MEDIUM', 'ch_1_LONG' ,\
                      'ch_2_SHORT','ch_2_MEDIUM', 'ch_2_LONG' ,
                      'ch_3_SHORT','ch_3_MEDIUM', 'ch_3_LONG' ,
                      'ch_4_SHORT','ch_4_MEDIUM', 'ch_4_LONG' ]
        stichRatio = []
        channelExist = []
        dct = {}
        # print(len(realData))
        for i in range(len(realData)):
            print(realData[i].name_band)
            dct[realData[i].name_band] = realData[i]
        # print("to leksiko einai ", str(repr(dict)))
        for i in cubesNames:
            exists = False
            for cube in realData:
                if(cube.name_band == i ):
                    exists = True
            if(exists == True): #if sub-channel exists
               channelExist.append(1) #put 1
            else:
               channelExist.append(0) #else put 0
    
        cube_idx = 0
        allRatio = []
        print("KAI TO CHANNELS EXISTS PERIEXEIIIIIIIIIIIIIIIIIIII ", str(channelExist))
        for i in range(len(channelExist)-1):
            if(channelExist[i] == 1 and channelExist[i+1] == 1) and dct[cubesNames[i]].ls[len(dct[cubesNames[i]].ls)-1] > dct[cubesNames[i+1]].ls[0]:
                print('Calculate stitching ration', cubesNames[i], 'and ', cubesNames[i+1])
                # self.calculateStitchRatio(realData[cube_idx],realData[cube_idx+1])
                
                ratio = self.calculateStitchRatio(dct[cubesNames[i]],dct[cubesNames[i+1]],aperture_correction,idx, grid)
                cube_idx = cube_idx+1
                allRatio.append(ratio)
            else:
                print('Den mporoume na upologisoume sta ', cubesNames[i], 'and ', cubesNames[i+1])
                allRatio.append(np.NaN)
        # print(allRatio)
      
        
        return allRatio        
                
    def calculateStitchRatio(self, ch1,ch2, aperture_correction , idx, grid):

        if aperture_correction: 
                print('aderfia edw eimaste complettt')
                if grid:
                    ch1_data = np.array(ch1.spectrum_PSF_corrected)
                    ch2_data = np.array(ch2.spectrum_PSF_corrected)
                else:
                    ch1_data = np.array(ch1.spectrum_PSF_corrected)[idx,:]
                    ch2_data = np.array(ch2.spectrum_PSF_corrected)[idx, :]

        else:   
            if grid:
                ch1_data = np.array(ch1.corrected_spectrum)[idx,:]
                ch2_data = np.array(ch2.corrected_spectrum)[idx, :]
            else:    
                ch1_data = np.array(ch1.corrected_spectrum)[:,idx]
                ch2_data = np.array(ch2.corrected_spectrum)[:,idx]
            

        # plt.plot(ch1_data, label = "ch1")
        # plt.plot(ch2_data, label="ch2")
        # plt.legend()
        # plt.show()
        ch1_ls = ch1.ls
        ch2_ls = ch2.ls
        delta1 = ch1.CDELT3
        delta2 = ch2.CDELT3
        ch1_start = np.where(np.array(ch1_ls)>= ch2_ls[0] - delta2/2 )[0][0]
        ch2_stop =  np.where(np.array(ch2_ls)<= ch1_ls[len(ch1_ls)-1]+ delta1/2)[0]
        ch2_stop =ch2_stop[len(ch2_stop)-1]                                                                     
        ch1_overlapping_before = ch1_data[ch1_start:]
        ch1_over_ls= ch1_ls[ch1_start:]
        ch2_overlapping = ch2_data[:ch2_stop+1]

        ch1_mean = np.mean(ch1_overlapping_before)
        ch2_mean = np.mean(ch2_overlapping)
        ratio = ch2_mean /ch1_mean 
        
        return ratio
        
    #%%
    ##### Function for grid extraction, multiple points with square apertures #####
    ###############################################################################    
    # @ratio_list:  RA of grid central point in arcsec. (float)
    # @idx: dec of grid central point in arcsec. (float) 
    # @spectrum: Number of points in X pixel coordinates. (int)    
    ########### --> Return [df_res,realData_all,spec1ds]   ############################
    # @stitched_spectrum: The spectrum 1D element. (Spectrum1D)           
    ###############################################################################    
    def stichSpectrum(self, ratio_list, idx, spectrum):
        print("fix ... is ", ratio_list[idx])
        if np.isnan(ratio_list[idx]): #if there is no spectrum there
            return [np.NaN] * len(spectrum)
        else:
            print("den einai nan ", ratio_list[idx], ratio_list[idx:])
            if np.isnan(np.sum(ratio_list[idx+1:]) ):#if there is nan at some next point
                print('just one ratio ', ratio_list[idx], np.array(spectrum).shape)
                return np.array(spectrum )* ratio_list[idx]
            else:
                ratio = 1 
                print('multiple ratio')
                for i in range(idx, len(ratio_list)):
                    ratio = ratio* ratio_list[i]
                print("ratio is: ", ratio)    
                return np.array(spectrum) * ratio    
            
            
  #%% 
    def delteFilesatPath(self, path):
        print()         

        
        files = glob.glob(path+'/*')
        for f in files:
            os.remove(f)          
    
    
# #%%
# def gridExtraction(self, data_list, r_arcsec, background_sub ):
    
            
    
    
    def readGridParamsFile(self, path,filename):
        print('reading Grid Parameters')
    
    