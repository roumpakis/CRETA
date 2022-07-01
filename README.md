# CRETA Spectrum Extraction Tool
---
# Introduction #
CRETA is a spectra extraction tool for astronomical images. It has two different modes: the single point and the grid extraction. CRETA also provides a rich set of preprocessing options and extractions can highly parameterized. 

 ![picture alt](https://github.com/roumpakis/CRETA/blob/main/Images/22.png?raw=true "CRETA")
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
