# CRETA Spectrum Extraction Tool
---
# Introduction #
CRETA is a spectra extraction tool for astronomical images. It has two different modes: the single point and the grid extraction. CRETA also provides a rich set of preprocessing options and extractions can highly parameterized. 

 ![picture alt](https://github.com/roumpakis/CRETA/blob/main/Images/22.png?raw=true "CRETA")
---
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
#### How to run it
##### 1. Create a ```cube_cp``` object that gives access to both extraction options
```python 
c = cube_cp()
```
##### 2. Extraction with default options. Run a single point extraction with default parameters set with ```singlePointExtraction``` method and a grid point extraction with ```gridExtraction ``` method.
```python 
c.singlePointExtraction(parameters_file = False)
c.gridExtraction()
```

##### 3. Parameterization of both spectrum extraction methods can be applied.
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
