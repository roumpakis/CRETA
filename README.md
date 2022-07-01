# CRETA Spectrum Extraction Tool
---
# Introduction #


Markup : ![picture alt](http://via.placeholder.com/200x150 "CRETA")
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

##### 3. Parameterization
 Parameterization of both spectrum extraction methods can be applied.
###### 3.1 Single Point Extraction parameters file
For single point extraction parametrization can be applied by changing the values of the parameters that ```params.txt``` file contains. 
```python 
c.singlePointExtraction(parameters_file = True)
```
---
