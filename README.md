# CRETA Spectrum Extraction Tool
---
# Introduction #
CRETA is a spectra extraction tool for astronomical images. It has two different modes: the single point and the grid extraction. CRETA also provides a rich set of preprocessing options and extractions can highly parameterized. 

 ![picture alt](https://github.com/roumpakis/CRETA/blob/main/Images/22.png?raw=true "CRETA")

# Files & Folders Description #
| Filename    |Description |
|--------------|:-----:|
| cube_cp.py   |  Contains main functionality of CRETA, the only objects that user will interact with|   
| SubCube.py       |  Represents the information associated with each sub-band |
| UserAPI.py   | Class for file handling |   
| MIRIPreproc.py       |Class that contains functionality in order to create a modular system  |
| convolve_miricubes.py      | aaa |
| mylmfit2dfun.py   | aaa |   
| params.txt     |Text file that contains parameters for single point extraction |
---

| Filename    |Description |
|--------------|:-----:|
| Data   |  Contains data .fits files  |   
| PSF    | Contains PSF fits files|
| PSF_INF   | Infinite Flux per PSF sub band  |   
| Centroid_Sky  |Centroids of PSF sub bands in sky coordinates 
| Centroids |  Centroids of PSF sub bands in pixel coordinates  |
| Sigma_Eff |Sigma Effective per PSF sub band |
| Results |  Folder that CRETA saves output files |




#### How to run it
##### 1. Create a ```cube_cp``` object that gives access to both extraction options
```python 
c = cube_cp()
```
##### 2. Extraction with default options. Run a single point extraction with default parameters set with ```singlePointExtraction``` method and a grid  extraction with ```gridExtraction ``` method.
```python 
c.singlePointExtraction(parameters_file = False)
c.gridExtraction()
```
```gridExtraction ``` with default parameters  defines a grid that covers the whole FoV of the sub-band with at  the larger larger wavelength. Then, CRETA performs the extraction based on that grid. 

![image](https://user-images.githubusercontent.com/60132957/177275175-442c3ca5-2a65-48ab-a4ca-38b121aceff2.png)


##### 3. Parameterization of both spectrum extraction methods can be applied.

Single Extraction parameters <br/>
```aperture_type:```  Aperture type: 0 for Circular, 1 for Rectangular. (int)<br/>
```convolve:``` Fix resolution option. (Boolean)<br/>
```parameters_file:``` Use the parameters file or the command execution option. (boolean)<br/>
```user_ra:``` Center RA in degrees. (float)<br/>
```user_dec:``` Center Dec in degrees. (float)<br/>
```user_r_ap:``` user defined radius in arcsec. (float)<br/>
```point_source:``` Point or extended source extraction option. (boolean)<br/>
```lambda_ap:``` Wavelength that aperture is defined, only for point source (float)<br/>
```apperture_correction:``` Apperture correction option (boolean)<br/>
```centering:``` Center user input with a 11x11 box (boolean)<br/>
```lambda_cent:``` Wavelength of centering (float)<br/>
```background:``` Background subtraction option (boolean)<br/>
```r_ann_in:``` Inner annulus radius (float)<br/>
```width:``` Width of annulus (float)  <br/>

---
gridExtraction( first_subband = 'G140H', last_subband = 'ch_4_LONG', x_steps = -1, y_steps = -1, r = -1, distance = -1,  user_centroid=False,

Grid Extraction parameters <br/>
```convolve:``` Fix resolution option. (Boolean)<br/>
```parameters_file:``` Use the parameters file or the command execution option. (boolean)<br/>
```user_ra:``` Center RA in degrees. (float)<br/>
```user_dec:``` Center Dec in degrees. (float)<br/>
```user_r_ap:``` user defined radius in arcsec. (float)<br/>
```point_source:``` Point or extended source extraction option. (boolean)<br/>
```lambda_ap:``` Wavelength that aperture is defined, only for point source (float)<br/>
```apperture_correction:``` Apperture correction option (boolean)<br/>
```centering:``` Center user input with a 11x11 box (boolean)<br/>
```lambda_cent:``` Wavelength of centering (float)<br/>
```first_subband:``` Point or extended source extraction option. (boolean)<br/>
```last_subband:``` Wavelength that aperture is defined, only for point source (float)<br/>
```x_steps:``` Apperture correction option (boolean)<br/>
```y_steps:``` Center user input with a 11x11 box (boolean)<br/>
```r:``` Wavelength of centering (float)<br/>
```distance:``` Center user input with a 11x11 box (boolean)<br/>
```user_centroid:``` Wavelength of centering (float)<br/>

###### 3.1 Single Point Extraction parameters file
For single point extraction parametrization can be applied by changing the values of the parameters that ```params.txt``` file contains. 
```python 
c.singlePointExtraction(parameters_file = True)
```
---
---
#### Single Point Extraction Parameters
| Parameter    | Default Value | Data Type |
|--------------|:-----:|-----------:|
| aperture_type   |  1 |        int |
| convolve        |  False |    bool|
| parameters_file |  False |    bool|
| user_ra         |  0 |        float|
| user_dec        |  0 |        float|
| user_r_ap      | [0.25]|       list|
| point_source   | True |       bool|
| parameters_file |  False |    bool|
| lambda_ap       |  5 |        float|
| apperture_correction        |  False |       bool|
| centering |  False |    bool|
| lambda_cent       |  5 |        float|
| background        |  False |       bool|
| r_ann_in         |  1.23 |        float|
| ann_width        |  0.2 |        float|


#### Grid Extraction Parameters
| Parameter    | Default Value | Data Type |
|--------------|:-----:|-----------:|
| convolve        |  False |    bool|
| parameters_file |  False |    bool|
| user_ra         |  0 |        float|
| user_dec        |  0 |        float|
| r     | -1|       float|
| point_source   | False |       bool|
| parameters_file |  False |    bool|
| lambda_ap       |  0 |        float|
| apperture_correction        |  False |       bool|
| centering |  False |    bool|
| lambda_cent       |  5 |        float|
| first_subband | 'G140H' | string |
| last_subband  | 'ch_4_LONG' | string |
| x_steps | -1 | int |
|y_steps|-1|int|
|distance|-1|float|
|user_centroid|False|bool|
---
