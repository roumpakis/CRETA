# -*- coding: utf-8 -*-

"""

Created on Tue Jun 15 13:10:35 2021



@author: roub

"""





# import ipdb

import numpy as np

import matplotlib.pyplot as plt

from photutils.aperture import CircularAperture

from astropy.wcs import WCS

from MIRIPreproc import MIRIPreproc

from photutils.aperture import RectangularAperture

import os

import time

from astropy.coordinates import SkyCoord

from astropy import units as u

from photutils.aperture import SkyCircularAperture

from photutils.aperture import aperture_photometry

current_path = os.path.abspath(os.getcwd())

import sys, glob, os

sys.path.append('/Users/tanio/Sync/pywork/pys')

import pdb

import numpy as np

from numpy import unravel_index

from astropy.io import fits

from astropy.convolution import convolve, Gaussian2DKernel, convolve_fft

import mylmfit2dfun

from lmfit import Parameters

class SubCube:

    #%%

    def __init__(self,path, file_name,base_r,base_l, pixel_scale,base_pixel_scale, point_source,PSF, centering, back, r_in,width, aperture_type, convolve):

        

        self.preprocess = MIRIPreproc ()

        self.channel = self.preprocess.getFITSData(path+file_name)

        

        self.DQ_lista = []

        

        self.DQ = self.channel[12]

        

        self.pixel_scale = pixel_scale                  #Sub-cube's pixel scale

        self.base_pixel_scale = base_pixel_scale        #Pixel scale of sub-cube with the shortest wavelength

        self.image_before= self.channel[0].copy() 

        # self.zero_mask = self.image_before == 0         #Create the zero-mask 

        self.zero_mask = self.DQ != 0 



        self.image_before[self.zero_mask] = np.NaN      #Replace zero with NaN

        self.aperture_type = aperture_type

        

        self.primaryDict= self.channel[1]

        self.headers= self.channel[2]

        self.CRPIX3 = self.channel[3]

        self.CRVAL3 = self.channel[4]

        self.CDELT3 = self.channel[5]

        self.error_data = self.channel[7]               #The flux Error data 

        self.CDELT1_pix= self.channel[8]                    #First axis increment per pixel                 

        self.CDELT2_pix= self.channel[9]                    #Second axis increment per pixel                 

        self.CDELT1_arcsec = self.channel[8]                    #First axis increment per pixel                 

        self.CDELT2_arcsec = self.channel[9]                    #Second axis increment per pixel 

        self.CRVAL1 = self.channel[13]
        self.CRVAL2 = self.channel[14]
        

        self.name_band = self.channel[10]

        self.error_data[self.zero_mask] = np.NaN

        self.wcs = WCS( self.headers)

        self.total_flux_before = self.preprocess.totalImageFlux(self.image_before)

        self.base_r = base_r

        self.base_l = base_l

        self.bkg_x, self.bkg_y = self.calculateXYMax()      # Coordinates of pixel with max flux

        self.ls = self.preprocess.getChannelLs(self)        # Cube's wavelength list

        self.NZ,self.NY,self.NX = self.image_before.shape   # 

        self.nameL = []

        self.CDELT3L = []

        self.unit_ratio = []

        self.stiched_spectrum = []

        self.stiched_error = []

        for i in range(len(self.ls)):

            self.nameL.append(self.name_band)

            self.CDELT3L.append(self.CDELT3)

            self.unit_ratio.append(self.CDELT1_pix * self.CDELT2_pix)

        if point_source:

            self.rs =  self.preprocess.getPointSourceRs(self,base_r) 

            

        else:

            self.rs =  self.preprocess.getExtendedSourceRs(self,base_r)

   

        self.rs_arcsec = []

        self.rs_arcsec = np.array(self.rs) * self.pixel_scale 

       

        ## CHECK IF SIGMA EFF FILES EXISTS BEFORE CALLING THE FUNCTION !!!!!!!!!!!!!!

        seff_filename ="Sigma_Eff/seff_"+self.name_band+".csv"

        if convolve:

            print('>>>>>> Call the convolution Function')

            if PSF:

                if   os.path.isfile(seff_filename):

                                print('Just load Sigma Effective')

                                self.psf_sigma_eff =self.readSigmaEffective(seff_filename) #read PSF centroids from file

                else:

                                print('Calculating Sigma Effictive!!!!')

                                self.convolve_subband(PSF)

        



#%% Den xrhsimopoieitai

    def doAreaCalculations(self):

        self.area_pix = []

        for i in range(len(self.apers)):

           self.area_pix.append(np.pi*self.rs[i]**2)  

#%%

    ##### Function for  Flux unit correction. Convert flux spectrum to My  #####

    ###############################################################################  

          #Change sub-channel's flux and error unit

    ###############################################################################  

    def doFluxUnitCorrection(self):

        

        delta_factor = self.CDELT1_pix*self.CDELT2_pix

        self.delta_factor_list = []

        self.corrected_spectrum = (np.array(self.apers)* 10**6 * (delta_factor/206265**2)) # Convert flux spectrum to My  

        self.delta_factor_list.append(delta_factor)

        self.delta_factor_list * len(self.rs)

        self.error= np.array(self.error) * 10**6 * (delta_factor/206265**2) #Convert error spectrum to My

        

        if self.name_band == 'G140H'   :

            self.error = np.array(self.error )/ 206265**2



#%%

    def doCenters(self,user_ra,user_dec,PSF):

            self.user_ra = user_ra

            self.user_dec = user_dec

            if PSF:

                self.xys = self.preprocess.PSFCentering(self, self.image_before)

                # self.xys_sky = self.preprocess.PSFCenteringSky(self, self.image_before)

            else:

                # TDS modified to NOT center, just simply get the RA, dec and x, y

                #sky = self.preprocess.lambdaBasedCentering([self], self.ls[0],self.user_ra,self.user_dec)

                #self.user_ra = sky[0][0].ra

                #self.user_dec = sky[0][0].dec

                self.user_ra = user_ra

                self.user_dec = user_dec

                self.xys = self.preprocess.DataUsersCentering(self, self.image_before)

                

#%%                

    def doBackgroundSubtraction(self,point_source,r_in,width):

            print('>> '+self.name_band+':Background subtraction')

            

            if point_source:

                # self.bckg_rs =  self.preprocess.getPointSourceRs(self,r_in)

                # self.width_list =  self.preprocess.getPointSourceRs(self,width)

                self.bckg_rs =  self.preprocess.getExtendedSourceRs(self,r_in)

                self.width_list =  self.preprocess.getExtendedSourceRs(self,width)                

                

            else:

                self.bckg_rs =  self.preprocess.getExtendedSourceRs(self,r_in)

                self.width_list =  self.preprocess.getExtendedSourceRs(self,width)

            self.subtractBackground( self.bckg_rs, width)

            

#%% 

    def doPhotometry(self,PSF,back):

        time_photometry = time.time()

        print('>> '+self.name_band+':Perfoming photometry')     

        if back:

            [self.apers,self.area,self.error] = self.preprocess.AperturePhotometry(self, self.image_after)

        else:            

            [self.apers,self.area,self.error] = self.preprocess.AperturePhotometry(self, self.image_before )

            

            

        if self.name_band == 'G140H'   :

            print("NIRSpec Data")

            self.apers = np.array(self.apers )/ 206265**2

        if PSF:

            

            self.PSF_correction = []

            for i in range(len(self.apers[1])):

                self.PSF_correction.append(np.array(np.array(self.PSF_inf_flux) / np.array(self.apers)[:,i]).T)       
                
            # for i in range(len(self.apers[0])):

            #         self.PSF_correction.append(np.divide(np.array(self.PSF_inf_flux) , np.array(self.apers)[:,i]))

        print("+++++ Photometry Time: %s seconds +++++" % (time.time() - time_photometry))        

        #%%

    def subtractBackground(self, r_in, width ):   

        time_background= time.time()

        r_out = []



        for i in range(len(r_in)):

               width_pix = self.width_list[i]

             # width_pix = width / self.pixel_scale

               r_out.append(r_in[i]+width_pix)

        self.image_after, self.med_sigma, self.annulus, self.annulus_centroid, self.annulus_aperture, self.routs = self.preprocess.subtractUserBackground(self, r_in, r_out)

        self.image_after[self.zero_mask] = np.NaN

        

        # create the background flux spectrum

        self.background_spectrum = []

        for i in range(len(self.med_sigma)):

            # self.background_spectrum.append(self.med_sigma[i] * self.annulus_aperture[i].area)

            self.background_spectrum.append(self.med_sigma[i])

        self.inf_prece = []

        self.inf_apers = []

        self.const_apers = []

        self.inf_rs = []



        print("+++++ Background Subtraction Time: %s seconds +++++" % (time.time() - time_background))

        

        

    #%%

    def ComputeSourceInfAperture(self, r):

        print ('Compute INF ')

        r_start = r

        prece = []

        dr = 1

        max_v_l = []

        total_image_flux= self.preprocess.totalImageFlux(self.image_before)

        all_bright = []

        inf_ratio = []

        self.inf_rs=[]

        for z in range(len(self.image_before)):      

            # print('for')

            inf_aperture = []

            rr = []

            data = self.image_before[z,:,:]

            bright = []

                

            x,y = self.preprocess.imageCentroid(data,self.bkg_x,self.bkg_y) 

            x =self.xys[z][0]

            y =self.xys[z][1]

            # print(x,y,r)

            aperture_list = []

            while r < data.shape[0] and r<data.shape[1]:

                 aper = CircularAperture([x,y] , r)

                 # print(data.shape,aper.area)

                 aperture_list.append(aper)

                 # na ftiaksw ta lista me ta aperture kai na kalesw mia fora to photometry

                 a = self.preprocess.getApertureSum(x,y,z,r, data,self.wcs)



                 inf_aperture.append(a["aperture_sum"][0])

                 bright.append(a["aperture_sum"][0] / aper.area) 

                 rr.append(r)

                 r =r+dr

                 # print(z,x,y,r)

    

            r =    r_start  



            max_r =  np.where(np.isnan(inf_aperture))[0][0] - 1

            max_v = inf_aperture[max_r]

            max_v_l.append(max_v)

            inf_ratio.append(max_v_l[z]/total_image_flux[z])

            # prece.append(self.apers[z] / max_v)

            # inf_aperture = []  

            self.inf_rs.append(rr[max_r])

            all_bright.append(bright)



        self.inf_apers = max_v_l

        self.inf_ratio = inf_ratio

        self.inf_perce= prece

        return 

# #%%        

#     def ComputeExtendedSourceAperture(self,image,r):

#         print ('Compute R')



#         self.const_apers = \

#           self.preprocess.extendedSourcePhotometry(self,image,r)

# #%%        

#     def ComputePointSourceAperture(self,image):

#         self.apers = \

#           self.preprocess.pointSourcePhotometry(self,image)

          

          #%%        

    def PSFCorrection(self,ratio, PSF_ls):

        self.spectrum_PSF_corrected = []

        self.error_PSF_corrected = []

        # print('========= PSF CORERECETT --> ',ratio,len(ratio))



        for i in range(len(self.apers[0])):

            # PSF IS NOT FUNCTIONING CORRECTLY BECAUSE THERE IS NO OBSERVED PSF YET

            # ipdb.set_trace()
            plt.title('Corrected Spectrum')
            plt.plot(self.corrected_spectrum)
            plt.show()
            
            
            plt.title('RATIO')
            plt.plot(ratio[0])
            plt.show()
            print(np.array(ratio).shape)
            
            plt.plot(self.ls, label = 'Data ls')
            plt.plot(PSF_ls, label = 'PSF ls')
            plt.legend()
            plt.show()
            self.spectrum_PSF_corrected.append((np.array(self.corrected_spectrum)[:,i]* (np.array(ratio)[i,:])).T)

            self.error_PSF_corrected.append((np.array(self.error)[:,i]*(np.array(ratio)[i,:])).T)
          
            # self.ls --> refers to real data wavelength list


        if self.name_band == 'G140H':

            print("NIRSpec Data")

            self.spectrum_PSF_corrected = np.array(self.spectrum_PSF_corrected )/ 206265**2 

            self.error_PSF_corrected = np.array(self.error_PSF_corrected )/ 206265**2 

            

            

            

            

            

            

            

            

            

#%%

    def gridExtractionPhotometry(self, sky_ra, sky_dec, r_arcsec, data, plot):

        r_pix = r_arcsec /  self.pixel_scale

        # print("r pix is ", r_pix)

        positions = []

        for i in range(len(sky_ra)):

            sky = SkyCoord(sky_ra[i], sky_dec[i], unit = "deg")

            x,y,z  = self.wcs.world_to_pixel(sky,self.ls[0] * u.um )

            positions.append((x,y))

        aperture =RectangularAperture(positions, r_pix,r_pix)  

        res = []

        DQ_phot = []

        for i in range(len(data)):

            phot_table = aperture_photometry(data[i,:,:], aperture, wcs=self.wcs, error = self.error_data[i,:,:])

            DQ_phot_table = aperture_photometry(self.DQ[i,:,:], aperture, wcs=self.wcs)

            # if plot:

                # plt.plot(phot_table["aperture_sum"])    

                # plt.show()

            res.append(phot_table)

            DQ_phot.append(DQ_phot_table["aperture_sum"])

        aps = [aperture]*len(data)  

        self.DQ_lista = DQ_phot

        # print("META TO PHOTOMETRY KOPELA KI EXOUME ", np.array(DQ_phot).shape)

        return res , aps , DQ_phot 

#%%

    def sliceGridExtractionPhotometry(self, sky_ra, sky_dec, r_arcsec, data, plot, idx):

        r_pix = r_arcsec /  self.pixel_scale

        # print("r pix is ", r_pix)

        positions = []

        for i in range(len(sky_ra)):

            sky = SkyCoord(sky_ra[i], sky_dec[i], unit = "deg")

            x,y,z  = self.wcs.world_to_pixel(sky,self.ls[0] * u.um )

            positions.append((x,y))

        aperture =RectangularAperture(positions, r_pix,r_pix)  

        res = []

        DQ_phot = []

        

        phot_table = aperture_photometry(data[idx,:,:], aperture, wcs=self.wcs, error = self.error_data[idx,:,:])

        DQ_phot_table = aperture_photometry(self.DQ[idx,:,:], aperture, wcs=self.wcs)

            # if plot:

                # plt.plot(phot_table["aperture_sum"])    

                # plt.show()

        res.append(phot_table)

        DQ_phot.append(DQ_phot_table["aperture_sum"])

        aps = [aperture]*len(data)  

        self.DQ_lista.append(DQ_phot)

        # print("META TO PHOTOMETRY KOPELA KI EXOUME ", np.array(DQ_phot).shape)

        return res , aps , DQ_phot 

    

#%%

    def centroid2sky(self):

            print('pixels to sky')

            self.sky_xys = []

            for i in range(len(self.xys)):

                self.sky_xys.append(self.wcs.pixel_to_world(self.xys[i][0], self.xys[i][1], self.ls[i]))

      #%%   Plotting Functions

    def plotAllApers(self):

        plt.loglog(self.ls,self.apers, label = 'Point Source')

        plt.loglog(self.ls,self.const_apers,label = 'const r')

        plt.loglog(self.ls,self.inf_apers, label = 'inf aperture')

        plt.xlabel('λ(μm)')

        plt.ylabel('Aperture Flux')

        plt.legend()

        plt.show()

        

    def plotApers(self, crackL, cracksV, maxV, color):

        plt.loglog(crackL,cracksV,"+r")

        plt.loglog(self.ls,self.apers, label = 'Point Source', color = color)

        plt.xlabel('λ(μm)')

        plt.ylabel('Aperture Flux')

        plt.legend()

        plt.show()        

        



    def calculateXYMax(self):

           image = self.image_before[0,:,:]

           NX = int(image.shape[0]/2)

           NY = int (image.shape[1]/2)

           # print(NX,NY)



            # 25% to 75% sub Inage

           img = image

           start_Y = int(NY/2)

           start_X = int (NX/2)

           subImg = img[start_Y:start_Y+NY, start_X:start_X+NX]



           # max flux point

           xys = np.where(subImg == np.nanmax(subImg))

           yy = xys[0][0]

           xx = xys[1][0]



           return [yy,xx]

           # return centroid_com(subImg)



    def plotAllAperPLT(self):

        # plt.plot(self.ls,self.apers_before, label = 'Point Source Pefore Bck')

        plt.plot(self.ls,self.apers, label = 'Extented r')

        plt.plot(self.ls,self.const_apers,label = 'const r')

        plt.plot(self.ls,self.inf_apers, label = 'inf aperture')

        plt.plot(self.ls,self.diffs, label = 'Diffs')

        plt.xlabel('λ(μm)')

        plt.ylabel('Aperture Flux')

        plt.legend()

        plt.show()       

        

        

        

    def plotCircles(self, background):

        

        colors = ['white', 'red', 'orange','magenta','lime','cyan','yellow','crimson','goldenrod','darkviolet','deepskyblue','chocolate','yellowgreen']

        img = self.image_before[0,:,:]

        for i in range(1,len(self.image_before)):

                    img = img + self.image_before[i,:,:]



        img = self.image_before[0,:,:]            

        plt.figure()

        plt.imshow(img) 

        apers = []

        patches = []



        path = current_path+"/Results/"

        

        for i in range(len(self.rs[0])):

            aper = CircularAperture(self.xys[0], self.rs[0][i]) 

            apers.append(aper)

            ap_patches = aper.plot(color=colors[i], lw=2,label='Photometry aperture_'+str(i))

            patches.append(ap_patches[0])

            

           

        if background:    

            ann_patches = self.annulus_aperture[0].plot(color='red', lw=2,label='Background annulus')

            patches.append(ann_patches[0])

        handles = ap_patches

        # plt.ylim([0,self.NY])

        #plt.ion() 

        plt.legend(handles=handles)

        plt.title(self.name_band)

        plt.savefig(path+self.name_band+".png")

        #plt.show()

        plt.close()

        



        

        

    def plotInformationMask(self,i):

            plt.imshow(self.image_after[i,:,:] != 0)

            aper = CircularAperture([self.xys[i]], self.rs[i]) 

            aper_inf = CircularAperture([self.xys[i]], self.rs[i])

            ap_patches = aper.plot(color='white', lw=2,label='Photometry aperture')

            ap_patches_inf = aper_inf.plot(color='green', lw=2,label='Inf aperture')

            ann_patches = self.annulus_aperture[i].plot(color='red', lw=2,label='Background annulus')

            handles = (ap_patches[0], ann_patches[0],ap_patches_inf[0])

            plt.ylim([0,self.NY])

            plt.legend(handles=handles)

            plt.show()



    def pltLR(self):

        plt.plot(self.ls,self.rs)

        plt.ylabel('r(pix)')

        plt.xlabel('λ(μm)')        

        plt.show()           

        

    def pltLXY(self):

        plt.plot(self.ls,self.xys)

        plt.ylabel('centroid(x,y)')

        plt.xlabel('λ(μm)')

        plt.show()        

        

        

    def PSFFlux(self):

        plt.plot(self.ls,self.PSF_inf_flux, label = 'PSF Inf Flux')

        plt.plot(self.ls,self.apers, label = 'PSF Aperture Flux')

        plt.ylabel('Flux')

        plt.xlabel('λ(μm)')

        plt.legend()

        plt.show()          

        

        plt.plot(self.ls,self.PSF_correction)

        plt.ylabel('PSF_inf / PSF_aper (Flux)')

        plt.xlabel('λ(μm)')

        plt.legend()

        plt.show()  

    def PSFCorrectedFlux(self,PSFCube):

        for i in range(len(PSFCube.PSF_inf_flux)):

            PSFCube.PSF_inf_flux[i] = PSFCube.PSF_inf_flux[i] * 5 

        plt.plot(self.ls,self.apers, label = 'Before PSF Correction')

        plt.plot(self.ls,self.apers_corrected, label = 'After PSF Correction')

        plt.plot(self.ls,PSFCube.PSF_inf_flux,label = 'PSF INF Flux')

        plt.ylabel('Flux')

        plt.xlabel('λ(μm)')

        plt.legend()

        plt.show()      

        

        

    def plotGrid(self, pixels,r,x,y, cubes):

        from matplotlib.patches import Rectangle

      

        

        # for i in range(1,len(self.image_before)):

        #     for j in range(self.image_before.shape[1]):

        #                     for k in range(self.image_before.shape[2]):

        #                             img[j,k] = img[j,k] + self.image_before[i,j,k]

            

        img = cubes.image_before[0,:,:]

        plt.ion() 

        for i in range(1,len(cubes.image_before)):

                    img = img + cubes.image_before[i,:,:]

        

        plt.imshow(img)

        plt.plot(x,y, 'o',color="red", label="User Input Centroid")

        for i in range(len(pixels)):

            xx = pixels[i][0] - (r/2)

            yy = pixels[i][1] - (r/2)

            plt.gca().add_patch(Rectangle([xx,yy],r,r,linewidth=1,edgecolor='r',facecolor='none'))

            plt.plot(pixels[i][0],pixels[i][1], 'bo')

            plt.title(cubes.name_band)

        plt.legend()    

        plt.show() 





        

    def plotCenters(self):

        from matplotlib.patches import Rectangle

        plt.ion() 

        plt.imshow(self.image_before[0,:,:])

        plt.plot(self.xys[0][0],self.xys[0][1], 'br')

        plt.show()    

        last = len(self.xys)-1

        plt.ion() 

        plt.imshow(self.image_before[last,:,:])

        plt.plot(self.xys[last][0],self.xys[last][1], 'bo')

        plt.show()         

        

#%%         

    def plotApertureInfo(self, background):

        from matplotlib.patches import Rectangle

 

        img = self.image_before[0,:,:]

        for i in range(1,len(self.image_before)):

                    img = img + self.image_before[i,:,:]



        plt.imshow(img)

        plt.plot(self.xys[0],self.xys[1], 'o',color="red", label="User Input Centroid")       

#%%

    def convolve_subband(self, PSF):

             import numpy as np

             # Dimensions

             nwaveinds, nys, nxs = self.image_before.shape

             # Create 2D mesh grid

             x, y = np.meshgrid(np.arange(nxs), np.arange(nys))

             # Initialize the list where all the PSF sigmas will be stored

             self.psf_sigma_eff = []       

             # Read slice from cube

             for waveind in range(nwaveinds):

                 frame = self.image_before[waveind]

                 frame_err =  self.error_data[waveind]

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

                 current_channel_pixel_scale = self.pixel_scale

                 self.psf_sigma_eff.append(np.sqrt(gauss2d_fit.pars[3] * gauss2d_fit.pars[4]) * current_channel_pixel_scale) # Units of arcsecond

                #### FOR NOW IT IS HARDCODED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                 # self.psf_sigma_eff.append( self.ls[waveind] / 5 * 0.33 ) # Units of arcsecond

                 

                 

                 

             import numpy as np

             from sklearn.linear_model import LinearRegression

             model = LinearRegression()

             model_x = np.array(self.ls).reshape(-1,1)

             # print(model_x.shape)

             # if self.name_band == 'ch_4_LONG' and PSF:

             #     print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa')

             #     model_y =  self.psf_sigma_eff[: int(0.5 * len(self.psf_sigma_eff))]

             #     model_x = np.array(self.ls[: int(0.5 * len(self.ls))]).reshape(-1,1)

             #     model_x2 = np.array(self.ls[int(0.5 * len(self.ls))  :]).reshape(-1,1)

             # else:

             model_y =  self.psf_sigma_eff

             model.fit(model_x, model_y)

             model_pred = model.predict(model_x)

             # all_model_pred = []

             # all_model_pred.extend(model_pred)

             # if self.name_band == 'ch_4_LONG' and PSF:

             #     model_pred2 =  model.predict(model_x2)

             #     all_model_pred.extend(model_pred2)

             #     print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa', len(all_model_pred), len(self.ls))

             plt.ion()    

             plt.plot(self.ls, self.psf_sigma_eff, label = 'Sigma Effective')  

             plt.plot(self.ls, model_pred, label = 'Fitted') 

             plt.title(self.name_band)

             plt.legend()

             plt.show() 

             print(np.array(self.psf_sigma_eff).shape, np.array(model_pred).shape)

             self.psf_sigma_eff = list(model_pred)

             #Regression

             # a,b

             # psf_sigma_eff = (a+ b+ls) * current_pixel_scale 

             if PSF:

                 self.writeSigmaEffective()







#%%

    def writeSigmaEffective(self):

            print('Write Sigma Effective to File')

            f = open("Sigma_Eff/seff_"+self.name_band+".csv", "w")

            for j in range(len(self.ls)):



                line = str(self.psf_sigma_eff) +'\n'

                f.write(line)

            f.close()   

    def readSigmaEffective(self,file):

             # print(file)

             res = []       

             f = open(file, "r")

             for line in f:





                 line = line.split('[')[1]

                 line = line.split(']')[0]

                 lines = line.split(',')



                 for i in lines:



                     # print(i)

                     res.append(float(i))

             f.close()    

             return res 



#%%

    def fixConvolved(self, reference_psf_sigma_eff,psf_sigma_eff_list):

        print('>> Fix Data With conv')            

        this_channel_pixel_scale = self.pixel_scale # placeholder

        nwaveinds, nys, nxs = self.image_before.shape

        plt.ion() 

        plt.imshow(self.image_before[0])

        plt.title(self.name_band+' \nData Before Conv')

        plt.show()

       

        for waveind in range(len(self.ls)):

                

                # Read sclice from cube

                frame = self.image_before[waveind,:,:]

                frame_err = self.error_data[waveind,:,:]



                # Create 2D mesh grid

                x, y = np.meshgrid(np.arange(nxs), np.arange(nys))

                # This will not be needed in the implementation (see below)

                valinds = np.where(~np.isnan(frame))

                xflat = x[valinds]

                yflat = y[valinds]

                xyflat = np.stack((xflat,yflat))

                

                # Calculate the kernel sigma used to convolve as the sqrt of the subtractions of the final and current sigmas squared

                kernel_sigma = np.sqrt(reference_psf_sigma_eff **2 -psf_sigma_eff_list[waveind]**2) / self.pixel_scale

                                                                    # USE PSF sigma effective !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                # print(kernel_sigma * this_channel_pixel_scale)

                # Generate the kernel

                gauss2d_kernel = Gaussian2DKernel(x_stddev=kernel_sigma)

                

                # Convolve the frame and the error frame

                self.image_before[waveind, :,:] = convolve(frame, gauss2d_kernel, normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)

                self.error_data[waveind,:,:] = convolve(frame_err, gauss2d_kernel, normalize_kernel=True, nan_treatment='interpolate', preserve_nan=True)

                

                # This fits again the convolved image. NOT NEEDED IN THE IMPLEMENTATION.

                # This is just to check that the PSF sigma of the convolved image is similar to the target PSF sigma of the final channel and wavelength

                # conv_gauss2d_fit = mylmfit2dfun.mylmfit2dfun(xyflat, convframe[valinds], gauss2d_fit.lmfitres.params, zerr=convframe_err[valinds])

                # conv_gauss2d_fit.lmfit(tilt=True)

                

                # Print the sigma of the current frame, the fitted sigma of the convolved frame, and the target sigma

                # print(channame+bandname, psf_sigma_eff[waveind], np.sqrt(conv_gauss2d_fit.pars[3] * conv_gauss2d_fit.pars[4]) * this_channel_pixel_scale, psf_sigma_eff[index_of_the_last_channel_last_wave]) # all in arcsec

                

                # pdb.set_trace()

        plt.ion()         

        plt.imshow(self.image_before[0])

        plt.title(self.name_band+' \nData After Conv')

        plt.show()

        plt.ion() 

        plt.imshow(self.error_data[0])

        plt.title(self.name_band+ ' \nError After Conv')

        plt.show()       



    def PSFGridPhotometry(self):

        print('>> Fix Data With conv')

        delta_factor = self.CDELT1_pix*self.CDELT2_pix

        apers = []

        photometries = []

        print('AAAAAAAAAAAAAAAAAAAAAAAAAAAA')

        for i in range(len(self.rs[0])):

            # print(self.xys[i][0], self.xys[i][1],  self.rs[0][0], np.array( self.rs)[0,0] )
            # x_pix,y_pix,z_pix  = self.wcs.world_to_pixel(sky,self.ls[0] * u.um )
            # EDW EINAI TO PROBLHMA TOU GRID
            aper = RectangularAperture([self.xys[i][0],self.xys[i][1]], self.rs[0][i], self.rs[0][i])

            apers.append(aper)

            phot_table = aperture_photometry(self.image_before[i,:,:], aper, wcs=self.wcs,error = self.error_data[i,:,:])

            photometries.append(phot_table["aperture_sum"])

        print(phot_table)

        print('exw mazepsei sunolika ', len(apers)) 

        return photometries

       



        

