# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 11:09:13 2022

@author: roub
"""


 

def write_grid_multipart(filename, aperture_correction,output_name): 
       import matplotlib.pyplot as plt
       import numpy as np
       from astropy import units as u
       from specutils import Spectrum1D
       import astropy
       from astropy.nddata import StdDevUncertainty
       from specutils import SpectrumList
       import pandas as pd  
       from astropy.io import fits  
       hdu_list = fits.open(filename)
       res_spec1d = []
       for i in range(len(hdu_list[1].data)):
           table = hdu_list[1].data
           # print(table)
           wave = table["Wave"] * u.um 
           Flux_ap = table["Flux_ap"] * u.Jy
           Flux_err_ap = table["Flux_err_ap"] * u.Jy
           Flux_ap_st = table["Flux_ap_st"] * u.Jy
           Flux_err_ap_st = table["Flux_err_ap_st"] * u.Jy
           if aperture_correction:
               Flux_ap_PSC = table['Flux_ap_PSC'] * u.Jy
               Flux_Err_ap_PCS = table['Flux_Err_ap_PCS'] * u.Jy
           DQ = table["DQ"]
           Band_Names = table["Band_Name"]
           metad =  hdu_list[1].header[str(i)]
           dict_list = metad.split(",")
           dct={}
           for j in range( len(dict_list)):
                    line = dict_list[j]
                    key = line.split(":")[0]
                    value = line.split(":")[1]
                    dct[key] = value
                    # print(line, "   ")
           print('')        
           dct['band_name'] = Band_Names

           fluxes = [Flux_ap[i], Flux_ap_st[i], DQ[i]]
           
           errors = [Flux_err_ap[i],Flux_err_ap_st[i]]
           if aperture_correction:
               fluxes.append(Flux_ap_PSC[i])
               errors.append(Flux_Err_ap_PCS[i])
           errors.append(len(DQ[i]) * [0])
           q = astropy.units.Quantity(np.array(fluxes), unit=u.Jy) 
           unc = StdDevUncertainty(np.array(errors))
           pec1d = Spectrum1D(spectral_axis=wave[i].T, flux=q ,uncertainty = unc, meta = dct) 
           res_spec1d.append(pec1d)
       res = SpectrumList(res_spec1d)    
       dct = {}   
       for i in range(len(res)):
            w = list(res[i].meta.values())
            xx = int(w[12])
            yy = int(w[13])
            dct[xx] = {}
           
       import matplotlib.pyplot as plt  
       plt.ion() 
            
       for i in range(len(res)):
            w = list(res[i].meta.values())
            xx = int(w[12])
            yy = int(w[13])
            dct[xx][yy] = res[i]    
            
       fluxes = np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)])  
       fluxes_stitched = np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)]) 
       ERR =np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)])
       ERR_stitched  = np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)]) 
       DQ =   np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)])
        
       if aperture_correction:
            fluxes_PSC = np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)]) 
            ERR_PSC= np.empty([ len(res[0].flux[0]) ,    len(dct[0]) ,len(dct)]) 
        
       for i in range(len(dct)):
            for j in range(len(dct[i])):
                fluxes [:,j,i] = dct[i][j].flux[0,:]
                fluxes_stitched[:,j,i] = dct[i][j].flux[1,:]
                DQ[:,j,i] =   dct[i][j].flux[2,:]  
                ERR[:,j,i]  = dct[i][j].uncertainty.array[0,:]
                ERR_stitched[:,j,i] = dct[i][j].uncertainty.array[1,:]
                if aperture_correction:
                    fluxes_PSC[:,j,i] = dct[i][j].flux[2,:]
                    ERR_PSC= dct[i][j].uncertainty.array[2,:]
                    DQ[:,j,i] =   dct[i][j].flux[2,:] 
                    
                
         
        
        #%%write FITS multicard
       keys = res[0].meta.keys()
       values = list(res[0].meta.values())
       NAXIS1, NAXIS2, NAXIS3 =  fluxes.shape
        
       import astropy.io.fits    as fits
       hdu = fits.PrimaryHDU()
        
       h_flux = fits.ImageHDU(fluxes, name='FLUX')
       header = h_flux.header  
       dictionary =res[i].meta
       header['PCOUNT'] = 0
       header['GCOUNT'] = 1
       header['EXTNAME'] = 'FLUX'
       header['SRCTYPE'] = 'EXTENDED'
       header['BUNIT'] = 'MJy/sr'
       header['WCSAXES'] = 3 
       header['CRPIX1'] = (int(dictionary[" 'grid_x'"]) + 1) #CRPIX1 starts from 1 
       header['CRPIX2'] = (int(dictionary[" 'grid_y'"]) + 1) #CRPIX2 starts from 1 
       header['CRPIX3'] = (1) #CRPIX2 starts from 1 
       header['CRVAL1'] = float( dictionary["'extraction_RA'"].split(" ")[2]) #Extraction RA
       header['CRVAL2'] = float( dictionary[" 'extraction_DEC'"].split(" ")[2])  #Extraction DEC
       header['CRVAL3'] = float( dictionary[" 'first_wave'"])
       header['CTYPE1'] = 'RA---TAN'
       header['CTYPE2'] = 'DEC---TAN'
       header['CTYPE3'] = 'WAVE'
       header['CUNIT1'] = 'deg'
       header['CUNIT2'] = 'deg'
       header['CUNIT3'] = 'um '
       header['CDELT1'] = float(dictionary[" 'CDELT1'"] )/ 3600 #in degrees
       header['CDELT2'] = float(dictionary[" 'CDELT1'"] )/ 3600 #in degrees
       header['CDELT3'] = float(0)
       header['PC1_1']  = -1
       header['PC1_2']  = 0.0
       header['PC1_3']  =  0
       header['PC2_1']  = 0.0
       header['PC2_2']  = 1.0
       header['PC2_3']  =  0
       header['PC3_1']  = 0
       header['PC3_2']  = 0
       header['PC3_3']  =  1
       values[1] = values[1].replace("'", "")
       header['EXTRTYPE'] =dictionary[" 'exrtaction_type'"]
       header['REXTR'] = float(dictionary[" 'r'"])
       #add GRCNTRA , dec
       from astropy.table import Table
       df_names = pd.DataFrame(res[0].meta['band_name'][0])
       df_names.columns = ['Band_name']
       t_names = Table.from_pandas(df_names)
       
       h_err= fits.ImageHDU(ERR, name='ERR')
       h_flux_stitched = fits.ImageHDU(fluxes_stitched, name='FLUX_STCHT')
       h_wave = fits.ImageHDU(res[0].spectral_axis, name='Wave')
       h_err_stitched = fits.ImageHDU(ERR_stitched, name='ERR_STCHT')
       h_dq = fits.ImageHDU(DQ, name='DQ')
        
       if aperture_correction:
            h_flux_PSC= fits.ImageHDU(fluxes_PSC, name='FLUX_PSC')
            h_ERR_PSC = fits.ImageHDU(ERR_PSC, name='ERR_PSC')
        
        
       names_array = np.array(list(df_names['Band_name']))
       col1 = fits.Column(name='BandName', format='20A', array=names_array)
       coldefs = fits.ColDefs([col1])
       h_names = fits.BinTableHDU.from_columns(coldefs, name = "BANDNAMES")
        
       if aperture_correction:
            hdulist = fits.HDUList([hdu,  h_flux,h_err, h_flux_PSC,h_ERR_PSC,  h_flux_stitched, h_err_stitched, h_dq,h_wave, h_names])
       else:
                    hdulist = fits.HDUList([hdu,  h_flux,h_err, h_flux_stitched, h_err_stitched, h_dq,h_wave, h_names]) 
            
       hdulist.writeto(output_name, overwrite=True)
        
      
       
  
