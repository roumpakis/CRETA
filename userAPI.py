"""
Created on Mon Aug 30 12:07:56 2021

@author: roub
"""
from astropy.coordinates import SkyCoord
import glob
import os
from astropy.io import fits
from MIRIPreproc import MIRIPreproc
from SubCube import SubCube
import pandas as pd
import numpy as np  
preprocess = MIRIPreproc ()
current_path = os.path.abspath(os.getcwd())
class userAPI:
     #%%
    def __init__(self):
        print('User API Created')
     #%%   
     # get PSF/REAL DATA cubes from API   
    def getSubCubes(self,path,user_r_arcsec,point_source,  PSF, centering, background,r_in,width,lambda_ap,aperture_type, read_only, files, convolve):
        # print("Eimais to subbbbbbbbbbbbbbbbbbbbbbbb ", PSF,"    ", read_only)
        # files = os.listdir(path)
        # for i in files: #exclude hidden files from mac
        #     if i.startswith('.'):
        #         files.remove(i)
        #     if not i in read_only:
        #         print("!O!O!O!O!O!O!O Den diabazw apo ta PSF to ", i , " giati den uparxei sta data")
        #         files.remove(i)
                
        # print("OPA TOUMPA :", files)        
        PSFCubes = []
        cubes_lambda = []
        print(read_only)
        new_files = []
        for i in range(len(files)):
                
                
                cube =preprocess.getFITSData(path+files[i])
                if cube[10] in read_only and PSF:
                    print('Loading and .......... ', cube[10])
                    PSFCubes.append(cube)
                    cubes_lambda.append(cube[4])
                    new_files.append(files[i])
                else:
                    PSFCubes.append(cube)
                    cubes_lambda.append(cube[4])
                    new_files.append(files[i])
 
        [PSFs , files ]= self.sortCubesByLambda(PSFCubes,cubes_lambda,new_files)
        
        dist = np.abs(lambda_ap - PSFs[0][4])
        base_pixel_scale = PSFs[0][6] #the first cube pixel scale / shortest lambda        
        base_l = PSFs[0][4] #the shortest lambda/ starting lambda
        for i in range(1,len(PSFs)):
            
            if dist > np.abs(lambda_ap - PSFs[i][4]):
                base_pixel_scale = PSFs[i][6] #the first cube pixel scale / shortest lambda        
                base_l = PSFs[i][4] #the shortest lambda/ starting lambda
                dist =  np.abs(lambda_ap - PSFs[i][4])
                
        base_r = user_r_arcsec #user defined radius in arcsec
        pixel_scale = []
        base_l_list = []
        # base_r_list = []
        #add first cubes pixel scale, lambda, radius
        pixel_scale.append(base_pixel_scale) 
        base_l_list.append(base_l)
        # base_r_list.append(base_r)
        
        for i in range(len(PSFs)):
            # if PSFs[i][4] != base_l:
                base_l_list.append(PSFs[i][4]) # add i-th cubes's starting lambda
                pixel_scale.append(PSFs[i][6]) # add i-th cube's pixel scale
                # base_r_list.append( (np.array(base_r_list[i-1])*pixel_scale[i-1]) / pixel_scale[i])
         
        #create the list with ubCubes elements    
        res = []
        for i in range(len(PSFs)):

            res.append( SubCube(path,files[i],base_r,base_l,PSFs[i][6],base_pixel_scale,point_source,  PSF, centering, background,r_in,width,aperture_type, convolve))
    
        return [res, pixel_scale, base_l_list,[],cube[11], files]
    #%%    
    def sortCubesByLambda(self,cubes,lambdas, files):
         lambdas_cp = lambdas.copy()
         cubes_cp = cubes.copy()
         files_cp = files.copy()
         res = []
         res_files = []
         for i in range(len(lambdas)):
             minLambdaIndex = lambdas_cp.index(min(lambdas_cp ))
             res.append(cubes_cp [minLambdaIndex])
             res_files.append(files_cp [minLambdaIndex])
             lambdas_cp.remove(min(lambdas_cp ))
             # print(minLambdaIndex, cubes[minLambdaIndex])
             del cubes_cp[minLambdaIndex]
             # print(lambdas,cubes)
             del files_cp[minLambdaIndex]
             
         return [res, res_files]    
                 
      #%%
    def loadUserParams(self, filename):
        
        f = open(filename, "r")
        res = []
        for x in f:
            x = x.replace("\n"," ")
            [key,value] = x.split('=')
            res.append(value)

        return res
    
    #%%Write centroids to file
    def writeCubeCentroids(self,cube,j):
            print('Writting Centroids to file')
        
        
        # for j in range(len(cubes)):
            f = open("Centroids//xys_"+cube.name_band+".csv", "w")
            for i in range(len(cube.ls)):
                line = str(cube.ls[i]) +"," +str(cube.xys[i][0])+","+str(cube.xys[i][1])+"\n"
                f.write(line)
            f.close()
            
    # Read Centroids from file        
    def readCubeCentroids(self,file):
             res = []       
             f = open(file, "r")
             for line in f:
                 # print(line)
                 [l,x,y] = line.split(",")
                 res.append([float(x),float(y)])

             f.close()    
             return res
         #%%PSF INF FLUX
    def writePSFInfFlux(self,PSFs):
            print('Writting Centroids to file')
        
        
            for j in range(len(PSFs)):
                f = open("PSF_INF//inf_"+PSFs[j].name_band+".csv", "w")
                inf_flux = preprocess.PSFInfFlux(PSFs[j].image_before, PSFs[j].CDELT1_pix*PSFs[j].CDELT2_pix)
                # print('INF FLUX IS '+str(inf_flux))
                PSFs[j].PSF_inf_flux = inf_flux
                line = str(inf_flux) +'\n'
                # print(line)
                f.write(line)
            f.close()
            
    # Read Centroids from file        
    def readPSFInfFlux(self,file):
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
         
            
         #%%PSF INF FLUX
    def writeCentroidSky(self,PSFs):
            print('Writting Centroids to file')
        
        
            for j in range(len(PSFs)):
                f = open("Centroid_Sky//sky_"+PSFs[j].name_band+".csv", "w")
                res = []
                for i in range(len(PSFs[j].ls)):
                    [jj,kk] = PSFs[j].xys[i]
                    sky = PSFs[j].wcs.pixel_to_world(jj,kk,PSFs[j].ls[i])
                    res.append(sky)
                    line = sky[0].to_string() +'\n'
                    f.write(line)
            f.close()
            
    # Read Centroids from file        
    def readCentroidSky(self,file):
             # print(file)
             res = []       
             f = open(file, "r")
             for line in f:

                 lines = line.split(' ')
                 ra = float(lines[0])
                 dec = float(lines[1])
                 c = SkyCoord(ra,dec, unit="deg") 
                 res.append(c)
             f.close()    
             return res                
#%%
    def writeResultsFile(self, filename,user_params, df,final_ratio,output_path,new_ra,new_dec, ap_l_iss,grid_extraction, grid_NX, grid_NY, grid_distance , PSF_path, Data_path):
        from astropy import units as u
        df.to_csv(output_path+filename,index=False)

        print("Writting output file...")
        
        if ap_l_iss != -1:
            warning_message = "######################################## WARNING/ERRORS \n \
            r_ap > annulus r_in from wavelength: "+str(ap_l_iss)
        else:
            warning_message=""
            
        if grid_extraction == 1:
            grid_txt = "######################################## Grid  EXtraction: {NX: "+str(grid_NX)+"(steps), NY: "+str(grid_NY) + "(steps), distance: "+str(grid_distance )+ "(pix)}\n"
        else:
             grid_txt = ""

        
        line = '########################################'\
        +'# Output file of spectrum extraction \n' + '# PSF flies path: '+PSF_path+'\n# Daata files path: '+Data_path \
             +'\n######################################### User Paramaeters\n# r_ap: '+user_params[0]+'(arcsec)\n# Input [RA,dec]: ['+user_params[1]+','+user_params[2]+'](degrees)\n'  \
             +'# Point source: '+user_params[3]+'\n# Lambda aperture: '+user_params[4]+'\n# Aperture correction: '+  user_params[5]+'\n# Centering:' +user_params[6]\
             +'\n# Centering lambda: '+user_params[7]+' (microns)'\
            +'\n# New [RA,dec] = ['+new_ra.to_string(unit=u.hour, sep=('h', 'm', 's'))+','+str(new_dec)+']'\
            +'\n# Background Subtraction: '+user_params[8]+'\n# Background Inner Radious, Annulus Width: '\
             +user_params[9]+','+user_params[10]+'(arcsec,arcsec) \n ########################################'\
             +'# Output File description \n# COLUMN_NAME : DESCRIPTION (UNIT) \n# Wave: wavelength \n# Band_name (ch_CHANNEL_BAND) :  channel and band source of information'\
             +'\n# Flux_ap: Aperture flux density (Jy)\n# Flux_err_ap: Aperture flux density error (Jy)'\
             +'\n# R_ap: Aperture radius (arcsec)\n# Background: Background flux surface brightness (MJ/sr) \n# Flux_ap_PSC: Flux density after point source correction (Jy)'\
             +'\n# Flux_Err_ap_PSC: Flux density error after point source correction (Jy) \n# PSC: Point-source aperture correction factor'\
             +'\n# Flux_ap_scl: Flux density after band scaling (Jy)'\
             +'\n# Flux_err_ap_scl: Flux density error after band scaling (Jy)'\
             +'\n# DQ: Data Quality Flag. 0: OK\n'\
             + warning_message+grid_txt\
             +'######################################## Results'
                 
        with open(output_path+filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\n') + '\n' )
            f.write('Stitching Ratio: '+str(final_ratio).rstrip('\n') + '\n' + content) 
            f.close()
        # print(str(final_ratio))
        

             