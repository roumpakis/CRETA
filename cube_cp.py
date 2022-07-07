
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 17:02:01 2021

@author: roub
"""
import numpy as np
import matplotlib.pyplot as plt
from MIRIPreproc import MIRIPreproc
from userAPI import userAPI
import pandas as pd
import time
from astropy import units as u
import os
from specutils import Spectrum1D
import astropy
from astropy.nddata import StdDevUncertainty
import datetime
from astropy.wcs import WCS
current_path = os.path.abspath(os.getcwd())
start_time = time.time()
from astropy.table import Table
preprocess = MIRIPreproc ()
user = userAPI()
class cube_cp:


    def __init__(self):
        print('User API Created')


#%%
    ##### Function for single point extraction                                #####
    ###############################################################################  
        
    # @aperture_type: Aperture type: 0 for Circular, 1 for Rectangular. (int)
    # @convolve: Fix resolution option. (Boolean)
    # @parameters_file: Use the parameters file or the command execution option. (boolean)
    # @user_ra: Center RA in degrees. (float)
    # @user_dec: Center Dec in degrees. (float)
    # @user_r_ap: user defined radius in arcsec. (float)
    # @point_source: Point or extended source extraction option. (boolean)
    # @lambda_ap: Wavelength that aperture is defined, only for point source (float)
    # @apperture_currection: Apperture correction option (boolean)
    # @centering: Center user input with a 11x11 box (boolean)
    # @lambda_cent: Wavelength of centering (float)
    # @background: Background subtraction option (boolean)
    # @r_ann_in: Inner annulus radius (float)
    # @width: Width of annulus (float)    
    ########### --> Return [df_res,data,,pec1d]   ############################
    # @df_res: A dataframe with extraction information. (pandas.DataFrame)
    # @data: A list of data elements. (list of SubeCube)
    # @meta: Dictionary with metadata information. (dict)
    # @pec1d: The spectrum 1D element. (Spectrum1D)        
    ###############################################################################   
        
    def singlePointExtraction(self, aperture_type=1, convolve=False, parameters_file = True, user_ra=0, user_dec=0,\
                         user_r_ap=[0.25], point_source=True, lambda_ap = 5, apperture_correction=False,centering=False,\
                        lambda_cent=5,background=False, r_ann_in=1.23, ann_width=0.2):
        path = current_path+"\\Results\\"
        filename = current_path+'\\params.txt'
        import time 
        [df_res,data,meta] = self.CRETA(filename,"",aperture_type,[],0,convolve,  \
                            user_ra, user_dec,user_r_ap, point_source, lambda_ap, apperture_correction,centering, \
                                lambda_cent,background, r_ann_in, ann_width,parameters_file)
        
        pec1d = self.create1DSpectrum(df_res,meta)
        
        ts = time.time()
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d %H:%M:%S")
        now_str = str(now)
        now_str = now_str.replace(':', '-')
        now_str = now_str.replace(' ', '_')
        outname = path+"JWST_"+str(now_str)+"singe_extraction_spec1d.fits"


        namesList = (df_res['Band_name'].values.tolist())  
        for i in range(len(namesList)):
            namesList[i] = str(namesList[i])
        plt.plot(namesList)
        plt.show()
        t = self.customFITSWriter([df_res], False, outname,[pec1d], False,  namesList)  
        
        self.plotStoreApertures(data, background)
        print('>> Single Point extraction execution time: ' + str(time.time() - ts ))
        return [df_res,data, meta,pec1d]
#%%
    
    ##### CRETA is the main function where the spectrum extraction is implemented #####
    ###############################################################################    
    # @params_file_name:  The filename where the parameters set is stored. (string)
    # @fid: The file identity, a name for output file separation. (string) 
    # @aperture_type: Aperture type: 0 for Circular, 1 for Rectangular. (int)  
    # @convolve: Fix resolution option. (Boolean)
    # @parameters_file: Use the parameters file or the command execution option. (boolean)
    # @user_ra: Center RA in degrees. (float)
    # @user_dec: Center Dec in degrees. (float)
    # @user_r_ap: user defined radius in arcsec. (float)
    # @point_source: Point or extended source extraction option. (boolean)
    # @lambda_ap: Wavelength that aperture is defined, only for point source (float)
    # @apperture_currection: Apperture correction option (boolean)
    # @centering: Center user input with a 11x11 box (boolean)
    # @lambda_cent: Wavelength of centering (float)
    # @background: Background subtraction option (boolean)
    # @r_ann_in: Inner annulus radius (float)
    # @width: Width of annulus (float)    
    # @grid_distance: distance between 2 points on grid in arcseconds. (float) 
    ########### --> Return [df,realData_all,meta_dict]############################
    # @df: A dataframe with extraction information. (pandas.DataFrame)
    # @realData_all: A list of data elements. (list of SubeCube)
    # @meta_dict: Dictionary with metadata information. (dict)
    ###############################################################################
    def CRETA(self, params_file_name, fid,aperture_type,names, grid_extraction, convolve,\
             user_ra, user_dec,\
                         user_r_ap, point_source, lambda_ap, aperture_correction,centering,\
                        lambda_cent,background, r_ann_in, ann_width,parameters_file = True ):
        preprocess = MIRIPreproc ()
        
        #%%        print ('ok')
        ###
        # Step 1: Load Paramaters from params.txt file
        ###
        pixel_scale = [] #pixel scale to arc sec
        base_l_list = [] #the first plane lambda, R
        base_r_list = []
        
        time_parameters_loading = time.time()
        REAL_DATA = False
        PSF = True
        from pathlib import Path
        user = userAPI()
        # print(params_file_name)
        
        
        
        if parameters_file:
            params = user.loadUserParams( params_file_name) #Load User Parameters
            aper_rs = params[0].split(",")
            # print(repr(aper_rs))
            user_rs_arcsec = []
            for i in range(len(aper_rs)):
                user_rs_arcsec.append(float(aper_rs[i]))   
            if 'm' in params[1] and 'm' in params[2]:
                   from astropy.coordinates import SkyCoord
                   Stringc = SkyCoord(params[1], params[2], frame='icrs')
                   user_ra = float(repr(Stringc.ra).split(" ")[1])
                   user_dec = float(repr(Stringc.dec).split(" ")[1])
            else:    
                user_ra = float(params[1])
                user_dec = float(params[2])
            
            params[3] = params[3].replace(" ","")
            params[5] = params[5].replace(" ","")
            params[6] = params[6].replace(" ","")
            params[8] = params[8].replace(" ","")
            point_source = params[3] == 'True'
            lambda_ap = float(params[4])
            aperture_correction = params[5] == 'True'
            centering = params[6] == 'True'
            l_c =  float(params[7])
            background =  params[8] == 'True'
            r_in = float(params[9])
            width = float(params[10])
            
        else:
            user_rs_arcsec = user_r_ap
            l_c = lambda_cent
            r_in = r_ann_in
            width = ann_width
            params = []
            params.append(str(user_r_ap))
            params.append(str(user_ra))
            params.append(str(user_dec))
            params.append(str(point_source))
            params.append(str(lambda_ap))
            params.append(str(aperture_correction))
            params.append(str(centering))
            params.append(str(lambda_cent))
            params.append(str(background))
            params.append(str(r_ann_in))
            params.append(str(ann_width))
            
            
            
        PSF_path = current_path+"\\PSF\\"
        data_path = current_path+"\\Data\\"
        output_path = current_path+"\\Results\\"
              
        #%% 
        ###
        # Step 2: Print user parameters
        ###       
        
        print(' PSF:', PSF_path)
        print(' Data:', data_path)
        #Print user parameters
        print('########################################')
        print('     Load User Parameters ')
        print('########################################')
        print('r: '+str(user_rs_arcsec)+'(arcsec)')
        print('[RA,δ]: ['+str(user_ra)+','+str(user_dec)+'](degrees)')
        print('Point Source: '+str(point_source))
        print('Aperture Correction: '+str(aperture_correction)+'(PSF Correction)')
        print('Centering: '+str(centering))
        print('Centering lambda: '+str(l_c)+'(μm)')
        print('Background Subtraction: '+str(background))
        if background:
            print('Background Inner Radious, Annulus Width: '+str(r_in)+','+str(width)+'(arcsec,arcsec)')
        print('PSF sub-cubes Path: '+PSF_path)
        print('Data sub-cubes Path: '+data_path)
        print('########################################')
        print("!!!!! Loading User's Parameters': %s seconds !!!!!" % (time.time() - time_parameters_loading))    
        
        
        #%% 
        ###
        # Step 3: Create the metadata Dictionary that we will use it for the Spectrum1D output file
        ###
        if point_source:
            aper_type = "point source"
        else:
            aper_type = "extended source"
        from astropy.coordinates import SkyCoord    
        c = SkyCoord(ra=user_ra*u.degree, dec=user_dec*u.degree, frame='icrs')    
        meta_dict = {"centroid":c.to_string('hmsdms'),"exrtaction_type":aper_type , \
                     "aperture_correction":aperture_correction, "background_subtraction":background, \
                     "Background Inner Radious":r_ann_in, "Annulus Width":ann_width,   
                     "Centering":centering, "Centering lambda":l_c                         
                      }
        #%% Load PSF: PSF_all is a list with all PSF sub-cubes sorted by wavelength
  #%%Load Real Data
        print('\nLoading Data..')    
        files = os.listdir(current_path+"\\Data\\")
        for i in files: #exclude hidden files from mac
            if i.startswith('.'):
                files.remove(i)      
        [realData_all, pixel_scale, base_l_list,base_r_list,output_file_name, read_only] = user.getSubCubes(data_path,user_rs_arcsec,point_source,  REAL_DATA, centering, background,r_in,width,lambda_ap,aperture_type, [], files, False)        
        print('\nLoading PSF..')
        timePSF_loading = time.time()
        [PSF_all, pixel_scale, base_l_list,base_r_list,output_file_name, zzz] = user.getSubCubes(PSF_path,user_rs_arcsec,point_source,  PSF, centering, background,r_in,width,lambda_ap,aperture_type, read_only, files, convolve)
        print("!!!!! Loading PSF Cubes': %s seconds !!!!!" % (time.time() - timePSF_loading))              
        for i in range(len(realData_all)):
            if convolve:
                      realData_all[i].fixConvolved(PSF_all[-1].psf_sigma_eff[-1],PSF_all[i].psf_sigma_eff) 
                    

        #%% Centering
        time_centering = time.time()
        if centering:
            new_sky =  preprocess.lambdaBasedCentering(PSF_all, l_c, user_ra, user_dec)
            print('New Sky is: ', new_sky)
            user_ra = new_sky[0][0].ra
            user_dec = new_sky[0][0].dec
        else:
            from astropy.coordinates import SkyCoord
            c = SkyCoord(ra=user_ra*u.degree, dec=user_dec*u.degree, frame='icrs')    
            user_ra = c.ra
            user_dec = c.dec
            
        print("!!!!! Centering exec time': %s seconds !!!!!" % (time.time() - time_centering))  

        #%% Create Centroids  of Apertures PSF/ Real Data
        if aperture_correction:
            for i in range(len(PSF_all)):

                
                filename = current_path+"\\Centroids\\xys_"+PSF_all[i].name_band+".csv"
                PSF_inf_filename = current_path+"\\PSF_INF\\inf_"+filename+".csv"
                
                # print(filename)
                if   os.path.isfile(filename):
                    print('Just load Centers')
                    PSF_all[i].xys = user.readCubeCentroids(filename) #read PSF centroids from file
                else:
                    
                    PSF_all[i].doCenters(user_ra,user_dec,PSF) #centering PSF cube 
                    user.writeCubeCentroids(PSF_all[i],i)  #PSF centroids in file

                #INF FLUX
                if   os.path.isfile(PSF_inf_filename):
                    print('Just load PSF Inf Flux')
                    PSF_all[i].PSF_inf_flux = user.readPSFInfFlux(PSF_inf_filename) #read PSF centroids from file
                else:
                    user.writePSFInfFlux(PSF_all)
                     
            
        for i in range(len(realData_all)):
            realData_all[i].doCenters(user_ra,user_dec,REAL_DATA)  
            
        #%% PSF Photometry
        time_PSF_photometry_all = time.time()    
        if aperture_correction:    ## APERTURE CORRECTIO !!!! ###
        
            for i in range(len(PSF_all)):
                
                if background:
                    PSF_all[i].doBackgroundSubtraction(point_source,r_in,width)  ## Background Subtraction if needed
                PSF_all[i].doPhotometry(PSF,background) 
                PSF_all[i].doFluxUnitCorrection()
        print("!!!!! PSF Photometry': %s seconds !!!!!" % (time.time() - time_PSF_photometry_all))      
        
        #%%Real Data Photometry
        time_data_photometry_all = time.time() 
        for i in range(len(realData_all)):

            if background:
                realData_all[i].doBackgroundSubtraction(point_source,r_in,width)   ## Background Subtraction if needed
        # Photometry
            realData_all[i].doPhotometry(REAL_DATA,background)
            # realData_all[i].doAreaCalculations() #calculate the area 
            realData_all[i].doFluxUnitCorrection() #change the photometry unit to MJ/sr
        print("!!!!! Real Data Photometry  time': %s seconds !!!!!" , (time.time() - time_data_photometry_all))    
        
        #%% PSF CORRECTION
        if aperture_correction:
            for i in range(len(PSF_all)):
                realData_all[i].PSFCorrection(PSF_all[i].PSF_correction)
        
        
        #%% Create all all_lists 
        #    
        time_create_list_all  = time.time() 
        [all_rs_arcsec,all_ls,all_apers,all_xys,all_area_pix,all_bright,all_error_spectrum,all_corrected_spectrum,all_delta,\
         all_names,all_unit_ratio,all_background,all_r_in,all_rs,all_ps, all_psc_flux, all_psc_err] =\
            preprocess.getSubcubesAll(realData_all,background, aperture_correction)

        print("!!!!! Create list all': %s seconds !!!!!" , (time.time() - time_create_list_all))
        
        #%%aperture_correction 
        time_create_list_all  = time.time() 
        if aperture_correction:
            PSF_ratio = []
            spectrum_PSF_corrected = []
            error_PSF_corrected = []
            for i in range(len(PSF_all)):
                PSF_ratio.append([])
                spectrum_PSF_corrected.append([])
                error_PSF_corrected.append([]) 
                for j in range(len(PSF_all[i].rs[0])):
                    PSF_ratio[i].append(np.array(PSF_all[i].PSF_correction)[j,:])
                    spectrum_PSF_corrected[i].append(np.array(realData_all[i].spectrum_PSF_corrected)[j,:])
                    error_PSF_corrected[i].append(np.array(realData_all[i].error_PSF_corrected)[j,:])
        
            
        print("!!!!! Aperture Correction': %s seconds !!!!!" % (time.time() - time_create_list_all))        
                
        #%%
        print("!!!!! Stitching Process!!!!!")
        cubesNames = [ 'G140H', 'ch_1_SHORT','ch_1_MEDIUM', 'ch_1_LONG' ,\
                      'ch_2_SHORT','ch_2_MEDIUM', 'ch_2_LONG' ,
                      'ch_3_SHORT','ch_3_MEDIUM', 'ch_3_LONG' ,
                      'ch_4_SHORT','ch_4_MEDIUM', 'ch_4_LONG' ]
            
            
        dct = {}
        for i in range(len(realData_all)):
            dct[realData_all[i].name_band] = realData_all[i]
            
        all_ratio_list = []
        print(realData_all[0].rs[0])
        for qq in range(len(realData_all[0].rs[0])):
            print("QQQQQ = ",qq)
            rl = realData_all[0].preprocess.stichingRatioCalculation(realData_all,aperture_correction,qq, False)
            all_ratio_list.append(rl)
        data_idx = 0  
       
        for j in range(len(PSF_all[0].rs[0])):          #for every radius
             print("... for radius ", np.array(PSF_all[0].rs)[0,j])
             for i in range(len(cubesNames)-1):         # for every band name that would exist
                        print('i == ', i , "j === ",j)
                        if cubesNames[i] in dct:        # if the datacube is avaliable
                            data = dct[cubesNames[i]]

                            if cubesNames[i+1] in dct:  # if we can calculate the stittching ratio
                                
                                ratio = np.array(all_ratio_list)[j,i]
                                         
                                if aperture_correction: #if PSC, use stich corrected spectrum
                                    beforeStich = np.array(data.spectrum_PSF_corrected)[j,:]
                                    meta = preprocess.stichSpectrum( list(np.array(all_ratio_list)[j,:]), i, beforeStich) # stitch aperture
                                    data.stiched_spectrum.append(meta)
                                    
                                    beforeStich_error = np.array(data.error_PSF_corrected)[j,:]
                                    stitched_error= preprocess.stichSpectrum( list(np.array(all_ratio_list)[j,:]), i, beforeStich_error) #stitch aperture
                                    data.stiched_error.append( stitched_error )#stitched spectrum  
                                    
                                    data_idx = data_idx+1

                                else: #use 
                                    beforeStich = np.array(data.corrected_spectrum)[:,j]
                                    stitched_flux = preprocess.stichSpectrum( list(np.array(all_ratio_list)[j,:]), i, beforeStich) #stitch aperture
                                    data.stiched_spectrum.append( stitched_flux )#stitched spectrum
                                    beforeStich_error = np.array(data.error)[:,j]
                                    stitched_error= preprocess.stichSpectrum( list(np.array(all_ratio_list)[j,:]), i, beforeStich_error) #stitch aperture
                                    data.stiched_error.append( stitched_error )#stitched spectrum                                    


                                    data_idx = data_idx+1
                                    
                            else: #if cube does not exists

                                        data.stiched_spectrum.append([np.NaN] * len(data.apers))
                                        data.stiched_error.append([np.NaN] * len(data.apers))
                        
             all_stittched_spectrum = []
             all_sttitched_error = []
             final_apers = []
             final_ls = []
             for i in range(len(realData_all)-1):

                 
                    final_apers.extend(np.array(realData_all[i].apers)[j,:]) 
                    final_ls.extend(np.array(realData_all[i].ls)) 
                    # print(realData_all[i].name_band , "  exei stitched ", np.array(realData_all[i].stiched_spectrum)[j,0])
                    all_stittched_spectrum.extend(np.array(realData_all[i].stiched_spectrum)[j,:])
                    all_sttitched_error.extend(np.array(realData_all[i].stiched_error)[j,:]) #if aperture correction error user corrected error


            #Check if r_ap is out the FoV, Photometry contains NaNs
             if     np.isnan(np.sum(final_apers)):
                    print('ERROR: r_ap is larger than FoV')


             spectrum_PSF_corrected_all = [] 
             error_PSF_corrected_all = [] 
             PSF_ratio_all = []      
            
            #PSF CORRECTION
             if aperture_correction:
                for i in range(len(spectrum_PSF_corrected)):
                      spectrum_PSF_corrected_all.extend(np.array(spectrum_PSF_corrected[i])[j,:])  
                      error_PSF_corrected_all.extend(np.array(error_PSF_corrected[i])[j,:])  
                      PSF_ratio_all.extend(np.array(PSF_ratio[i])[j,:])
                      

            
            
  #%%          
             time_stich =     time.time()
             res_all = []
             res_all.append(all_ls)
             res_all.append(all_names)
             res_all.append(np.array(all_corrected_spectrum)[:,j])
             res_all.append(np.array(all_error_spectrum)[:,j])
             res_all.append(np.array(all_rs_arcsec)[:,j])
            
             if background:
                    res_all.append(all_background)
               
             if aperture_correction:
                            
                            res_all.append(spectrum_PSF_corrected_all)
                            res_all.append((error_PSF_corrected_all))
                            res_all.append((PSF_ratio_all))
                
             if len(np.array(final_apers).shape)!=1:
                    res_all.append(np.array(all_stittched_spectrum)[j,:])
                    res_all.append(np.array(all_sttitched_error)[j,:])
             else:
                    res_all.append(all_stittched_spectrum)
                    res_all.append(np.array(all_sttitched_error))
                    
                
             # print("ERROR SHAPE: ",res_all)
             all_DQ_list = []
             for i in range(len(realData_all)):
                        cube = realData_all[i]
                        temp = cube.preprocess.getApertureDQList(cube)

                        all_DQ_list.extend(temp)
             res_all.append(all_DQ_list) 

             print("!!!!!!!!! STICHING: %s seconds ---" , (time.time() - time_stich))            
        #%%Create DF

             time_writting_output =     time.time()
            
             column_names = ['Wave', 'Band_name','Flux_ap','Flux_err_ap','R_ap']
            
             if background:
                column_names.append('Background')
             if aperture_correction:
                column_names.append('Flux_ap_PSC')
                column_names.append('Flux_Err_ap_PCS')
                column_names.append('PSC')
                
             column_names.append('Flux_ap_st')    
             column_names.append('Flux_err_ap_st')
             column_names.append('DQ')

             # print(background,aperture_correction,len(res_all))
                
             df = pd.DataFrame( res_all)
            
             path = current_path+"\\Results\\"
             df = df.T
             df.columns = column_names
             df = df.sort_values(by=['Wave'])  
             #%%
             #CHANGE DF dType
             df['Wave']= df['Wave'].astype(float)
             df['Band_name']= df['Band_name'].astype(str)
             df['Flux_ap']= df['Flux_ap'].astype(float)
             df['Flux_err_ap']= df['Flux_err_ap'].astype(float)
             df['R_ap']= df['R_ap'].astype(float)
             if aperture_correction:
                 df['Flux_ap_PSC']= df['Flux_ap_PSC'].astype(float)
                 df['Flux_Err_ap_PCS']= df['Flux_Err_ap_PCS'].astype(float)
                 df['PSC']= df['PSC'].astype(float)                 
             df['Flux_ap_st']= df['Flux_ap_st'].astype(float)
             df['Flux_err_ap_st']= df['Flux_err_ap_st'].astype(float)
             df['DQ']= df['DQ'].astype(float)

            #%%
             plt.loglog(df['Wave'],df['Flux_ap'],label = 'Flux Before PSC')
             
             if aperture_correction:
                plt.loglog(df['Wave'],df['Flux_ap_PSC'],label = 'Flux After PSC')
               
             plt.loglog(df['Wave'],df['Flux_ap_st'],label = 'Flux Stiched')
             
             plt.xlabel("Wavelength (μm)")
             plt.ylabel("Flux (Jy)")
             plt.legend()

             plt.loglog(df['Wave'],df['Flux_err_ap'], '--' ,markersize=1,label = 'Flux Error')
             if aperture_correction:
                   plt.loglog(df['Wave'],df['Flux_Err_ap_PCS'], '--' ,markersize=1,label = 'Flux Error PSC')
             plt.loglog(df['Wave'],df['Flux_err_ap_st'], '--' ,markersize=1,label = 'Error Stiched')
             
             plt.xlabel("Wavelength (μm)")
             plt.ylabel("Flux (Jy)")
             plt.legend()
             plt.savefig(path+"flux_and_err_spectrum_"+str(j)+".png")
             plt.show()
             
             aperture_lamda_issue = -1
             if background:
                
                if  len(np.where(np.array(all_rs)[:,j] > np.array(all_r_in))[0]) != 0 : 
                    index_with_issue = np.where(np.array(all_rs)[:,j] > np.array(all_r_in))[0][0]
                    aperture_lamda_issue = all_ls[index_with_issue]
                    
                    
             #create output file name based on timestamp       
             now = datetime.datetime.now()
             now = now.strftime("%Y-%m-%d %H:%M:%S")
             now_str = str(now)
             now_str = now_str.replace(':', '-')
             now_str = now_str.replace(' ', '_')    
             
             
             user.writeResultsFile("JWST_"+str(now_str)+'_'+str(fid)+"_"+str(user_rs_arcsec[j])+'.csv',params,df,all_ratio_list,output_path,user_ra,user_dec,aperture_lamda_issue,grid_extraction, 0, 0, 0, PSF_path, data_path )   
            

             print("---Writting Output : %s seconds ---" % (time.time() - time_writting_output))
        return [df,realData_all,meta_dict]
        print("---Execution Time: %s seconds ---" % (time.time() - start_time))
        
        
#%%
    ##### Function for Spectrum1D output file creation                          ###
    ###############################################################################    
    # @df_res: A dataframe with extraction information. (pandas.DataFrame)
    # @meta: Dictionary with metadata information. (dict)
    ########### --> Return [df_res,data, meta,pec1d]   ############################
    # @pec1d: The spectrum 1D element. (Spectrum1D)        
    ###############################################################################    
        
    def create1DSpectrum(self,df_res,meta):
            fluxes = []
            errors = []
            df = df_res
            fluxes.append(df['Flux_ap'].values * u.Jy)
            errors.append(df['Flux_err_ap'].values * u.Jy)
            DQ = df['DQ']
            if 'Flux_Err_ap_PCS' in df.columns:
                fluxes.append(df['Flux_ap_PSC'].values * u.Jy)
                errors.append(df['Flux_Err_ap_PCS'].values * u.Jy)
                
            
            fluxes.append(df['Flux_ap_st'].values * u.Jy)
            errors.append(df['Flux_err_ap_st'].values * u.Jy)
            fluxes.append(DQ)
            errors.append(DQ)
                
            wave = df_res['Wave'].values * u.um    
            wave_all = []
            for i in range(len(fluxes)):
                wave_all.append(wave)
            q = astropy.units.Quantity(np.array(fluxes), unit=u.Jy) 
            # q = astropy.units.Quantity(np.array(fluxes[1]), unit=u.Jy) 
            wave_all = np.array(wave_all).T * u.um
            unc = StdDevUncertainty(np.array(errors))
            pec1d = Spectrum1D(spectral_axis=wave, flux=q ,uncertainty = unc, meta=meta)   
            return pec1d
        
#%%       
        
# #%% Create a spectrum 1D with wave x[ band_name, pixel_scale,radius, data_quality ]
#     def createInfoSpectrum(self, df_res, dq, ps):
#          wave = df_res['Wave'].values * u.um
#          fluxes = []
#          channels = df_res['Band Name'].values 
#          rs = df_res['R_ap'].values 
#          fluxes.append(channels)
#          fluxes.append(rs)
#          fluxes.append(ps)
#          print(fluxes)
         
         
    def plotStoreApertures(self,data, background):
           
            for i in range(len(data)):
                data[i].plotCircles(background)
                
#%%
    def customFITSWriter(self, df_res, PSC, filename, spec1ds, aperture_correction, band_names):
       waves = []
       Flux_ap = []
       Flux_err_ap = []
       Flux_ap_st = []
       Flux_err_ap_st = []
       DQ = []
       Flux_ap_PSC = []
       Flux_Err_ap_PCS = []
       plt.plot(band_names)
       plt.show()
       names = ["Wave" ,'Band_Name', "Flux_ap", "Flux_err_ap", "Flux_ap_st", "Flux_err_ap_st", "DQ" ]
       if aperture_correction:
           names.append('Flux_ap_PSC')
           names.append('Flux_Err_ap_PCS')
       for i in range(len(df_res)):
           waves.append(df_res[i]["Wave"])
           # Band_name.append(df_res[i]["Band_name"])
           Flux_ap.append(df_res[i]["Flux_ap"])
           Flux_err_ap.append(df_res[i]["Flux_err_ap"])
           # R_ap.append(df_res[i]["R_ap"])
           
           if aperture_correction:
                Flux_ap_PSC.append(df_res[i]["Flux_ap_PSC"])
                Flux_Err_ap_PCS.append(df_res[i]["Flux_Err_ap_PCS"])
      
           Flux_ap_st.append(df_res[i]["Flux_ap_st"])
           Flux_err_ap_st.append(df_res[i]["Flux_err_ap_st"])
           DQ.append(df_res[i]["DQ"])
       band_names_fixed = []
       for i in range(len(waves)):
           band_names_fixed.append(band_names)           
       # all_data = [waves, Band_name, Flux_ap, Flux_err_ap, R_ap ]
       all_data = [waves ,band_names_fixed, Flux_ap, Flux_err_ap, Flux_ap_st, Flux_err_ap, DQ]
       if aperture_correction:
           all_data.append(Flux_ap_PSC)
           all_data.append(Flux_err_ap_st)
       # if PSC:
       #     all_data.append(Flux_ap_PSC)
       #     all_data.append(Flux_Err_ap_PCS)
       # all_data.append(Flux_ap_st)
       # all_data.append(Flux_err_ap_st)
       # all_data.append(DQ)           
       # tab = Table(all_data, names=df_res[0].columns)
       dct = {}
       for i in range(len(spec1ds)):
           print(str(spec1ds[i].meta))
           meta_str = str(spec1ds[i].meta)
           meta_str = meta_str.replace("{", "")
           meta_str = meta_str.replace("}", "")
           dct[str(i)] = meta_str
           

       print(len(band_names), len(waves[0]), len(band_names), len(all_data) )    
       plt.plot(band_names)
       plt.title('ta onomata')
       plt.show()
       print( type(df_res), type(PSC), type(filename), type(spec1ds), type(aperture_correction), type(band_names))
       print(dct)
       tab = Table(all_data, names = names, meta = dct)
       # import pdb
       # pdb.set_trace()
       tab.write(filename, format="fits")  
       return tab
           
           
       
#%%
    ##### Function that reads a grid extraction fits file.           ###
    ###############################################################################    
    # @filename: The fits filename . (string)
    ########### --> Return res_spec1d  ############################
    # @res_spec1d: A list of extracted spectra. (list of Spectrum1D)        
    ###############################################################################   

    def customFITSReader(filename): 
           from astropy.io import fits 
           hdu_list = fits.open(filename)
           res_spec1d = []
           for i in range(len(hdu_list[1].data)):
               table = hdu_list[1].data
               wave = table["Wave"] * u.um 
               Flux_ap = table["Flux_ap"] * u.Jy
               Flux_err_ap = table["Flux_err_ap"] * u.Jy
               Flux_ap_st = table["Flux_ap_st"] * u.Jy
               Flux_err_ap_st = table["Flux_err_ap_st"] * u.Jy
               DQ = table["DQ"]
               metad =  hdu_list[1].header[str(i)]
               dict_list = metad.split(",")
               dct={}
               for j in range( len(dict_list)):
                        line = dict_list[j]
                        key = line.split(":")[0]
                        value = line.split(":")[1]
                        dct[key] = value
                        # print(line, "   ")
                       
    
               fluxes = [Flux_ap[i], Flux_ap_st[i], DQ[i]]
               errors = [Flux_err_ap[i],Flux_err_ap_st[i]]
               errors.append(len(DQ[i]) * [0])
               q = astropy.units.Quantity(np.array(fluxes), unit=u.Jy) 
               unc = StdDevUncertainty(np.array(errors))
               pec1d = Spectrum1D(spectral_axis=wave[i].T, flux=q ,uncertainty = unc, meta = dct) 
               res_spec1d.append(pec1d)
           return res_spec1d
          
         #%%
    ##### Function that create grid extraction with default set of parameters.           ###
    ###############################################################################    
    # @path: Path to data files. (string)
    ########### --> Return cube_data  ############################
    # @res_spec1d: A list of data sub-channels. (list of SubCube)        
    ###############################################################################     
    def gridExtraction(self, point_source = False,lambda_ap=0,  centering = False, lambda_cnt = 0, parameters_file = False, plots=False, first_subband = 'G140H', last_subband = 'ch_4_LONG', x_steps = -1, y_steps = -1, r = -1, distance = -1, user_ra = 0, user_dec = 0,  user_centroid=False, aperture_correction = False, convolve=False):
        
        #%% create the read only list 
         cubesNames = [ 'G140H', 'ch_1_SHORT','ch_1_MEDIUM', 'ch_1_LONG' ,\
                      'ch_2_SHORT','ch_2_MEDIUM', 'ch_2_LONG' ,
                      'ch_3_SHORT','ch_3_MEDIUM', 'ch_3_LONG' ,
                      'ch_4_SHORT','ch_4_MEDIUM', 'ch_4_LONG' ]
         
             
         path = current_path+'\\Data\\'    
         if parameters_file:
                 all_grid_params = userAPI.loadUserParams(userAPI,current_path+'\\grid_params.txt')    
                 first_subband = all_grid_params[0]
                 first_subband = first_subband.replace(" ",'')
                 last_subband = all_grid_params[1]
                 last_subband = last_subband.replace(" ",'')
                 x_steps = int(all_grid_params[2])
                 y_steps = int(all_grid_params[3])
                 r = float(all_grid_params[4])
                 distance = float(all_grid_params[5])
                 user_ra = float(all_grid_params[6])
                 user_dec = float(all_grid_params[7])
                 user_centroid = all_grid_params[8] == 'True'
       
         
         
         
         if first_subband in cubesNames and last_subband in cubesNames:
             first_subband_idx = cubesNames.index(first_subband)
             last_subband_idx = cubesNames.index(last_subband)
             read_only = cubesNames[first_subband_idx : last_subband_idx+1]
         else:
             print("WARNING:Some names do not match,double check input")
             #we have to stop the execution!!!!!!!
             read_only = cubesNames
             
         all_DFs = []
         all_dcts = []
         spec1ds = []
         import time 
         ts = time.time()
         files = os.listdir(path)
         for i in files: #exclude hidden files from mac
            if i.startswith('.'):
                files.remove(i)  
                
         print('MPL MPL MPL ',read_only)       
         cubes = []
         cubes_lambda = []
         dct_cube_file = {}
         dct_cube_lambda = {}
         new_files = []
         for i in range(len(files)):
                cube =preprocess.getFITSData(path+files[i])
                if cube[10]  in read_only: #if sub ban 
                    cubes.append(cube)
                    cubes_lambda.append(cube[4])
                    # print("ANTISOTITOOT" , files[i], cube[10])
                    new_files.append(files[i])
                    # dct_cube_file[cube[10]] = files[i]
                    # dct_cube_lambda[cube[10]] = cubes_lambda[i]
         files = new_files   
         #%% Exclude the sub bands that are not on the 
         ii = 0
         print("BEFORE LEN: ", len(cubes)) 
         # cubes_new = []
         # files_new = []
         # cubes_lambda_new = []
         # for i in range(len(cubes)):
            
         # #     if i[10] not in read_only: #if sub band is not in read only list, remove it
         # #         print("Removiiiing : ", i[10], "  ", files[ii])
         # #         cubes.remove(i)
         # #         files.remove(files[ii])
         # #         cubes_lambda.remove(cubes_lambda[ii])
         # #     ii = ii+1

         # # print("AFTER LEN: ",len(cubes))   
        
         # # files = []
         # # cubes_lambda = []
         # # read_only_cp = read_only.copy()
         # # for ix in range(len(read_only)):
         # #    if read_only[ix] in dct_cube_file.keys(): 
         # #        files.append(dct_cube_file[read_only[ix]])   
         # #        cubes_lambda.append(dct_cube_lambda[read_only[ix]]) 
         # #    else:
         # #        read_only_cp.remove(read_only[ix])
         
         # # read_only = read_only_cp
             
         #     if cubes[i][10]  in read_only: #if sub ban 
         #         cubes_new.append(cubes[i])
         #         files_new.append(files[i])
         #         cubes_lambda_new.append(cubes_lambda[i])
             
         # print("READ ONLUUUUU ", read_only, files) 
         # cubes = cubes_new
         # cubes_lambda = cubes_lambda_new
         # files = files_new
         
         #%%   
         [cubes_sorted ,asx]= user.sortCubesByLambda(cubes,cubes_lambda,files)   
         for i in cubes:
             print('a onomata ', i[10])
         for i in files:
             print("AOUAOAUA",i)
         cube_data = cubes_sorted[-1]
         print('Kai xrhsimopoioume ws cube ton: ', cube_data[10])
         image = cube_data[0].copy()
         
         CDELT1 = cube_data[8]
         CDELT2 = cube_data[9]
         pixel_scale = np.sqrt(CDELT1*CDELT2)
         DQ = cube_data[12]
         CRVAL3 = cube_data[4]
         CDELT3 = cube_data[5]
         
         nan_mask = DQ!=0
         image[nan_mask] = np.NaN 
         wcs = WCS(cube_data[2])
         ls = []
         for i in range(image.shape[0]):
            l_i =  preprocess.getSimulatedL(CRVAL3,CDELT3 ,i)
            ls.append(l_i)
         l_ap = ls[-1]   #use as l_ap the last lambda of longer wavelength 
         
         
         if r < 0:  
             #correlated errors?
             r = pixel_scale/2
             
         if distance < 0:
            distance = 2*r 
         
            
         print(image.shape)
         [NZ, NY, NX ] = image.shape
         
         #IF the user does not define values for grid points in X or Y coordinate
         if x_steps < 0 or y_steps<0:
             
             if distance<0:
                 x_steps = NX
                 y_steps = NY
             else:
                 
                 x_steps = (NX*pixel_scale)/distance
                 y_steps = (NY*pixel_scale)/distance
                 # print('super mpiaaax ', x_steps, y_steps)
                 r = distance/2
         x_pix = float((NX-1)/2 )
         y_pix = float((NY-1)/2)
         print('x_pisx: ', x_pix, 'y_pix: ', y_pix)

         print('>> Grid Extraction Parameters Set <<<<')
         print('First sub-band: ', first_subband)
         print('Last subband: ', last_subband)
         print('X Steps: ', x_steps)
         print('Y Steps: ', y_steps)
         print('Distance: ', distance)
         # print(NX, "oupla upla ",NY)

#%% Load PSF   
         # point_source = True #SSSSSSSSSSSSSSSSSSSNKLFENLF$H)$H)$G$NG$B$U($U(JJ$(GJOP$NG$BG$$GK)))
         if point_source or convolve:
            print('POINT SOURCEE')         
            PSF_path = current_path+"\\PSF\\"
            PSF_files = os.listdir(PSF_path)
            for i in PSF_files : #exclude hidden files from mac
                if i.startswith('.'):
                    PSF_files .remove(i)    
            [PSF_all, pxs, bsl,base_r_list,ofn, zzz] = user.getSubCubes(PSF_path,r,point_source,  True, centering, False,0,0,lambda_ap,1, read_only, PSF_files, convolve)         
            
        #%%Load Real Data
         print('\nLoading Data..')  

         [realData_all, pixel_scale, base_l_list,base_r_list,output_file_name, read_only_data] = user.getSubCubes(path,r,point_source,  False, False, False,0.0,0.0,l_ap,1, read_only, files, convolve)
         for i in range(len(realData_all)):
              realData_all[i].rs = [realData_all[i].rs] 
              if convolve:
                  realData_all[i].fixConvolved(PSF_all[-1].psf_sigma_eff[-1],PSF_all[i].psf_sigma_eff)
              print(PSF_all[i].name_band)    
         print("Ta real data einaiiiii :", len(realData_all), " kai points: ", x_steps*y_steps)
         

         
        #%% Centering Process
         if user_centroid:
              print('user centroid')
              if centering:
                print('centering')
                new_sky =  preprocess.lambdaBasedCentering(realData_all,lambda_cnt, user_ra, user_dec) #center by labda, using the 11x11 box
                ra = new_sky[0][0].ra
                dec = new_sky[0][0].dec  
                print('Center coords after centering: ra = ', ra, ' dec = ',dec)
              else:   
                 ra = user_ra
                 dec = user_dec
                 print("User Defined RA,DEC without Centering:  ", ra, " dec ", dec)
         else: 
             sky = wcs.pixel_to_world(x_pix, y_pix, ls[-1]*u.um)
             ra = sky[0].ra
             dec = sky[0].dec
           
             
         #%%Create Grid points    
         sky_list,pixels_indices,names, sky_ra, sky_dec = realData_all[-1].preprocess.crateGridInArcSec(ra,dec,distance, x_steps, y_steps, realData_all[-1], r, realData_all, False, l_ap)
         all_photometries = []
         all_aps = []
         # for each avaliable sub-bands perform grid Photometry
         for i in range(len(realData_all)):
             subchanel_photometry, aps, DQ_list = realData_all[i].gridExtractionPhotometry( sky_ra, sky_dec, r, realData_all[i].image_before, plots)
             all_photometries.append(subchanel_photometry)
             all_aps.append(aps)
         ts2 = time.time()
         print("!!!!! Grid Extraction until Photometry': %s seconds !!!!!" % (ts2 -ts))
         # self.multipleParamsFiles(ra, dec, distance, x_steps,y_steps, path,1, r, False, l_ap) 
         
         for i in range(len(realData_all)):                
            preprocess.plotGridSubchanel( ra, dec, distance, x_steps, y_steps, realData_all[i], r)   
            
            
#%% PQD
          
         if aperture_correction:   
             for i in range(len(PSF_all)):
                  PSF_all[i].rs = [PSF_all[i].rs]    
             if aperture_correction:
                for i in range(len(PSF_all)):                
                    filename = current_path+"\\Centroids\\xys_"+PSF_all[i].name_band+".csv"
                    PSF_inf_filename = current_path+"\\PSF_INF\\inf_"+PSF_all[i].name_band+".csv"
                   
                    print(filename)
                    if   os.path.isfile(filename):
                        print('Just load Centers')
                        PSF_all[i].xys = user.readCubeCentroids(filename) #read PSF centroids from file
                    else:
                        
                        PSF_all[i].doCenters(user_ra,user_dec,True) #centering PSF cube 
                        user.writeCubeCentroids(PSF_all[i],i)  #PSF centroids in file
                        
                        
                      
    
                    #INF FLUX
                    if   os.path.isfile(PSF_inf_filename):
                        print('Just load PSF Inf Flux')
                        PSF_all[i].PSF_inf_flux = user.readPSFInfFlux(PSF_inf_filename) #read PSF centroids from file
                    else:
                        user.writePSFInfFlux(PSF_all)
    
             time_PSF_photometry_all = time.time()   
             #Center grid to PSF
             
             for i in range(len(PSF_all)):
                 PSF_sky_cnt_filename =  current_path+"\\Centroid_Sky\\sky_"+PSF_all[i].name_band+".csv"
                    
                 if   os.path.isfile(PSF_sky_cnt_filename):
                        print('Just load SKYY')
                        PSF_all[i].xys_sky = user.readCentroidSky(PSF_sky_cnt_filename) #read PSF centroids from file
                 else:
                        print("WRITE PSF SKY")
                        user.writeCentroidSky(PSF_all)  #PSF centroids in file  
             
             if aperture_correction:    ## APERTURE CORRECTIO !!!! ###
                PSF_correction_ratio = []
                subband_correction_ratio = []
                # calculate 1 correction value and apply for each spaxel
    
                for i in range(len(PSF_all)):
                    # PSF_cnt_sky_list =  preprocess.PSFCenteringSky(PSF_all[i])
                    
                    # for j in range(len(PSF_cnt_sky_list)): #all subband slices
                                # PSF_cnt_sky = PSF_cnt_sky_list[j]
                                
                                # PSF_sky_PSF_list,pixels_indices,PSF_names, PSF_sky_ra, PSF_sky_dec = preprocess.crateGridInArcSec(PSF_cnt_sky[0].ra,PSF_cnt_sky[0].dec,distance, x_steps, y_steps, PSF_all[i], r, PSF_all, point_source, l_ap)
                                # PSF_photometry, PSF_aps, PSF_DQ_list = PSF_all[i].sliceGridExtractionPhotometry(PSF_sky_ra, PSF_sky_dec, r, PSF_all[i].image_before, plots, j)       
                                PCR = []
                                band_photos = PSF_all[i].PSFGridPhotometry()
                                # print('Band Photozzz ',np.array(band_photos).shape)
                                for j in range(len(band_photos)):
                                    PCR.append(PSF_all[i].PSF_inf_flux[j] / band_photos[j])
                                # print(PCR)
                            
                                PSF_correction_ratio.append(PCR)
                    # preprocess.plotGridSubchanel(PSF_cnt_sky[0].ra,PSF_cnt_sky[0].dec, distance, x_steps, y_steps, PSF_all[i], r)                            
                # PSF_correction_ratio.append(np.array(subband_correction_ratio))
                # subband_correction_ratio = []
             PSC = np.array(PSF_correction_ratio) 
             # print(PSC[0][0])
             print('A!@#$%',np.array(PSC[0]).shape)
             for i in range(len(PSF_all)):
                 print(PSF_all[i].name_band, realData_all[i].name_band)

                 #%% 
        ###
        # Step 3: Create the metadata Dictionary that we will use it for the Spectrum1D output file
        ###

         aper_type = "extended_source"
         from astropy.coordinates import SkyCoord  
         # print(ra, dec)
         for az in range(len(sky_list)):
              extraction_c= SkyCoord(ra=sky_ra[az], dec=sky_dec[az])    

              
              meta_dict = { 'extraction_RA':sky_ra[az],'exrtaction_type':aper_type , 'aperture_correction':aperture_correction, 'background_subtraction':False, 
                           'Background Inner Radious':0, 'Annulus Width':0,  
                           'Centering':False, 'Centering lambda':l_ap, 'x_steps' : x_steps,   'y_steps' : y_steps,\
                            'distance':distance, 'r' : r , 'grid_x': pixels_indices[az][0], 'grid_y': pixels_indices[az][1], 'CDELT1':distance, 
                            'CDELT1':distance, 'extraction_DEC':sky_dec[az], 'first_wave':cubes[0][4]
                            }              
              all_dcts.append(meta_dict)    
             
             
          #%%   
         
            
         # flux_cube =  np.array  [len(all_photometries[i]),NY, NX]
         # flux_cube = np.NaN
         for grid_point_idx in range(len(all_photometries[0][0])):
             all_apers = []
             all_error = []
             PSC_flux = []
             PSC_err = []
             PSC_ratio = []
             for i in range(len(all_photometries)): #for each sub-channel
                 cube_apers = []
                 cube_error = []
                 cube_xys = []
                 realData_all[i].spectrum_PSF_corrected = []
                 realData_all[i].error_PSF_corrected = []
                 # sub-channel i
                 for j in range(len(all_photometries[i])):
                     #j is wavelength
                                 #[sub-channel][wavelength][grid_point]
                     photo = all_photometries[i][j][grid_point_idx]["aperture_sum"]
                     err = all_photometries[i][j][grid_point_idx]["aperture_sum_err"]
                     xx = all_photometries[i][j][grid_point_idx]["xcenter"]
                     yy = all_photometries[i][j][grid_point_idx]["ycenter"]
                     all_apers.append(photo)
                     all_error.append(err)
                     
                     cube_apers.append(photo)
                     cube_error.append(err)
                     cube_xys.append([xx,yy])
                     # cube_area.append()

                 realData_all[i].apers =  [cube_apers]
                 realData_all[i].error  = [cube_error]
                 realData_all[i].xys = [cube_xys]
                 realData_all[i].area = all_aps[i]
                 realData_all[i].doFluxUnitCorrection()


                 print( np.array(realData_all[i].corrected_spectrum).shape, (len(all_photometries[i])))
                 print(realData_all[i].name_band)
                 if aperture_correction:
                     for j in range(len(PSC[i])):
                         # print(PSF_all[j].name_band, realData_all[j].name_band)
                         realData_all[i].spectrum_PSF_corrected.append(PSC[i][j] * np.array(realData_all[i].corrected_spectrum[0,j]))
                         realData_all[i].error_PSF_corrected.append(PSC[i][j] * np.array(realData_all[i].error[0,j]))


             background = False
             # aperture_correction = False
             time_create_list_all  = time.time() 
             [all_rs_arcsec,all_ls,all_apers,all_xys,all_area_pix,all_bright,all_error_spectrum,all_corrected_spectrum,all_delta,\
              all_names,all_unit_ratio,all_background,all_r_in,all_rs,all_ps, all_psc_flux, all_psc_err] =preprocess.getSubcubesAll(realData_all,background, aperture_correction)



        
        #%%
             
             print("!!!!! Stitching Process!!!!!")
             #create a dictionary that contains the data cubes
             dct = {}
             for i in range(len(realData_all)):
                 dct[realData_all[i].name_band] = realData_all[i]
                 
             #calculates the stitching ratio between all avaliable sub-bands                
             all_ratio_list = []
             rl = realData_all[0].preprocess.stichingRatioCalculation(realData_all,aperture_correction,0, True)
             all_ratio_list = rl
             data_idx = 0  
  
    
             for i in range(len(cubesNames)-1):         # for every band name that  exists
                            print('i == ', i , "j === ",j)
                            if cubesNames[i] in dct:        # if the datacube is avaliable
                                data = dct[cubesNames[i]]
                                    
                                if cubesNames[i+1] in dct  :  # if we can calculate the stittching ratio
                                    ratio = np.array(all_ratio_list)

                                    #if aperture correction option is on, we have to stitch the psc corrected spectrum
                                    if aperture_correction:
                                        beforeStich = np.array(data.spectrum_PSF_corrected)
                                        beforeStich_error = np.array(data.error_PSF_corrected)
                                        
                                    #otherwise we have to stitch the original photometry spectrum    
                                    else:    
                                        beforeStich_error = np.array(data.error)[0]    
                                        beforeStich = np.array(data.corrected_spectrum)[0]
                                        
                                    #perform the stitching and assign it back to the cubes                                        
                                    stitched_flux = preprocess.stichSpectrum( list(np.array(all_ratio_list)), i, beforeStich) #stitch aperture
                                    data.stiched_spectrum = stitched_flux #stitched spectrum
                                    stitched_error= preprocess.stichSpectrum( list(np.array(all_ratio_list)), i, beforeStich_error) #stitch aperture
                                    data.stiched_error = stitched_error 
                                
                                    data_idx = data_idx+1

                                        
                                else: #if cube does not exists

                                            data.stiched_spectrum = []
                                            data.stiched_error = []
                                            data.stiched_spectrum.extend([np.NaN] * len(data.corrected_spectrum[0,:]))
                                            data.stiched_error.extend([np.NaN] * len(data.corrected_spectrum[0,:]))
                                            
                            
             all_stittched_spectrum = []
             all_sttitched_error = []
             final_apers = []
             final_ls = []
             for i in range(len(realData_all)-1):
    
                        final_apers.extend(np.array(realData_all[i].apers)[0,:]) 
                        final_ls.extend(np.array(realData_all[i].ls)[:]) 
                        all_stittched_spectrum.extend(realData_all[i].stiched_spectrum)
                        all_sttitched_error.extend((realData_all[i].stiched_error))#if aperture correction error user corrected error

    
             #attach the speuctum of last sub-channel @stitched spectrum 47
             final_apers.extend(np.array(realData_all[-1].apers)[0,:])
             final_ls.extend(np.array(realData_all[-1].ls)[:])                         
             if aperture_correction:
                 print('Ta teleutaia nabwww')
                 all_stittched_spectrum.extend( np.array(realData_all[-1].spectrum_PSF_corrected))  
                 all_sttitched_error.extend( np.array(realData_all[-1].error_PSF_corrected))

             else:                       
                 all_stittched_spectrum.extend(np.array(realData_all[-1].corrected_spectrum)[0])  
                 all_sttitched_error.extend(np.array(realData_all[-1].error)[0])                       
             final_error= []   
             PSF_ratio_all = []    
             plt.plot(all_stittched_spectrum)
             plt.plot(all_sttitched_error)
             plt.show()

                
        #%%  create lists that contains  the extracted results, for all data cubes
             all_psc_flux = []
             all_psc_err = []
             all_corrected_spectrum = []
             all_error_spectrum = []
             for qqq in range(len(realData_all)):
                all_psc_flux.extend(realData_all[qqq].spectrum_PSF_corrected)
                all_psc_err.extend(realData_all[qqq].error_PSF_corrected)
                all_corrected_spectrum.extend(realData_all[qqq].corrected_spectrum[0])
                all_error_spectrum.extend(realData_all[qqq].error[0])                
                 

             time_stich =     time.time()
             res_all = []
             res_all.append(all_ls)
             res_all.append(all_names)
             res_all.append((all_corrected_spectrum))
             res_all.append((all_error_spectrum))
             res_all.append((all_rs_arcsec))

             res_all.append(all_stittched_spectrum)
             res_all.append(all_sttitched_error)
        
             if aperture_correction:
                 res_all.append(np.array(all_psc_flux))
                 res_all.append(np.array(all_psc_err))
                 res_all.append(np.array(PSC_ratio))

             all_DQ_list = []
             for i in range(len(realData_all)):
                            all_DQ_list.extend(list(np.array(realData_all[i].DQ_lista) // 513)) 

             res_all.append(np.array(all_DQ_list)[:,grid_point_idx])
    
             print("!!!!!!!!! STICHING: %s seconds ---" , (time.time() - time_stich))            
            #%%Create a data frame that contains the extracted information based on the above lists
    
             time_writting_output =     time.time()
                
             column_names = ['Wave', 'Band_name','Flux_ap','Flux_err_ap','R_ap']
             column_names.append('Flux_ap_st')    
             column_names.append('Flux_err_ap_st')
             
             if aperture_correction:
                 column_names.append('Flux_ap_PSC')
                 column_names.append('Flux_Err_ap_PCS')
                 column_names.append('PSC')
                 
             # print(background,aperture_correction,len(res_all))
             column_names.append('DQ')       
             df = pd.DataFrame( res_all)
                
             path = current_path+"\\Results\\"
             df = df.T
             df.columns = column_names
             df = df.sort_values(by=['Wave']) 
             df = df.fillna(value=np.nan)
             

                 #CHANGE DF data Type
             df['Wave']= df['Wave'].astype(float)
             df['Band_name']= df['Band_name'].astype(str)
             df['Flux_ap']= df['Flux_ap'].astype(float)
             df['Flux_err_ap']= df['Flux_err_ap'].astype(float)
             df['R_ap']= df['R_ap'].astype(float)
             if aperture_correction:
                             df['Flux_ap_PSC']= df['Flux_ap_PSC'].astype(float)
                             df['Flux_Err_ap_PCS']= df['Flux_Err_ap_PCS'].astype(float)
                             df['PSC']= df['PSC'].astype(float)                 
             df['Flux_ap_st']= df['Flux_ap_st'].astype(float)
             df['Flux_err_ap_st']= df['Flux_err_ap_st'].astype(float)
             df['DQ']= df['DQ'].astype(float)
    
                # %% Plot the resulting spectra
             plt.loglog(df['Wave'],df['Flux_ap'],label = 'Flux Before PSC')  
             if aperture_correction:
                         plt.loglog(df['Wave'],np.array(all_psc_flux),label = 'Flux After PSC')
                   
             plt.loglog(df['Wave'],df['Flux_ap_st'],label = 'Flux Stiched')
                 
             plt.xlabel("Wavelength (μm)")
             plt.ylabel("Flux (Jy)")
             plt.legend()
             plt.show()
    
             plt.loglog(df['Wave'], df['Flux_err_ap'], '--' ,markersize=1,label = 'Flux Error')
             if aperture_correction:
                            plt.loglog(df['Wave'],df['Flux_Err_ap_PCS'], '--' ,markersize=1,label = 'Flux Error PSC')
             plt.loglog(df['Wave'],df['Flux_err_ap_st'], '--' ,markersize=1,label = 'Error Stiched')
                 
             plt.xlabel("Wavelength (μm)")
             plt.ylabel("Flux (Jy)")
             plt.legend()
             plt.show()
             
 #%%                
             aperture_lamda_issue = -1
             if background:
                    
                    if  len(np.where(np.array(all_rs)[:,j] > np.array(all_r_in))[0]) != 0 : 
                        index_with_issue = np.where(np.array(all_rs)[:,j] > np.array(all_r_in))[0][0]
                        aperture_lamda_issue = all_ls[index_with_issue]
                        
                 #create output file name based on timestamp       
             now = datetime.datetime.now()
             now = now.strftime("%Y-%m-%d %H:%M:%S")
             now_str = str(now)
             now_str = now_str.replace(':', '-')
             now_str = now_str.replace(' ', '_') 
             all_DFs.append(df)
             pec1d = self.create1DSpectrum(df,all_dcts[grid_point_idx])
             spec1ds.append(pec1d)    
                 
    
             print("---Writting Output : %s seconds ---" % (time.time() - time_writting_output))        
         ts = time.time()
         now = datetime.datetime.now()
         now = now.strftime("%Y-%m-%d %H:%M:%S")
         now_str = str(now)
         now_str = now_str.replace(':', '-')
         now_str = now_str.replace(' ', '_')
         results_path = current_path+"\\Results\\"
         output_file_name = results_path+"JWST_"+str(now_str)+"_"+str(i)+"_Grid_spec1d.fits"
                # pec1d.write() # write the fits file
         t = self.customFITSWriter(all_DFs, False, output_file_name,spec1ds, aperture_correction, all_names)  
         # self.customFITSReader(output_file_name) 
         
         return all_DFs, res_all, realData_all,all_dcts
             
         
            
            
            
            
            
            
            
            
            
            
             
             
             