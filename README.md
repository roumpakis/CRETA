# CRETA

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
---
