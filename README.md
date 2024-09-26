# Signed-Distance-Field
Fortran + CUDA based software to compute the signed-distance-field for OBJ file format geometry over cartesian grid. The code assumes regular grid along the streamwise and spanwise directions, while a non-uniform grid can be used in the vertical ($z$) direction.

## How to compile the code

### For small problems - serial CPU version of the code can be used
```
cd src
make ISCUDA="False"
```
### For large problems - single GPU version of the code can be used
```
cd src
make ISCUDA="True"
```

## How to change parameters
Follow the instructions in `parameters.in` to change the input parameters as per requirements

## How to run the code
Edit `parameters.in`
```
# Run the compiled program 
./gensdf
```
Example output

```
*** File: 'data/small_armadillo.obj' sucessfully read in  
   0.1900181770324707      seconds ***
 Geometry is bounded by (minimum)   0.4328440100000000      
   2.0000000000000000E-002   3.8979180000000002E-002
 Geometry is bounded by (maximum)   0.5671559900000001      
   0.1800000000000000        0.1610208200000000     
 Geometry has         99976 number of faces..
 Geometry has         49990 number of vertices..
 *** Min-Max Index-Value pair ***
 Min-Max x:          216   0.4218750000000000      |          295 
   0.5761718750000000     
 Min-Max y:            8   1.1718750000000000E-002 |          120 
   0.1867187500000000     
 Min-Max z:           20   3.0468750000000003E-002 |          108 
   0.1679687500000000     
 - - - - - - - - Finished Pre-Processing in             0 seconds...
 *** 1st - Pass: Distance calculation ***
||||||||||||||||||||||||||||||||||||||||||||||||||||100.00%    49990/   49990
 *** 2nd - Pass: Narrow band tagging - Ld <   3.1250000000000002E-003   ***
 Narrow band point:        164305 of       778624 total points.
||||||||||||||||||||||||||||||||||||||||||||||||||||100.00%       88/      88
 *** 3rd - Pass: Computing the sign for the distance | Computed on the GPU ***
 Device Name: NVIDIA GeForce RTX 3060
 Max Threads Per Block:          1024
   X-dimension:    2147483647
   Y-dimension:         65535
   Z-dimension:         65535
 GPU - Threads:            64            1            1
 GPU - Blocks:          2568            1            1
 *** Signed-Distance-Field (SDF) computed in            14 seconds ***
 *** Writing SDF to file ***
 *** Finished in            14 seconds ***
```
<center><img src="armadillo.png" height=200></center>

<center> 
Figure: Comparison of the calculated SDF around the Stanford Armadillo geometry scaled to a smaller size
</center>

## How to visualise the results
1. Use `python/numpy2vtk.py` to convert the `mask.bin` binary file to VTK file to read the data in Paraview/VisIT or other post-processing software
2. Use `python/readmask.py` to directly read the `mask.bin` binary file and plot slices

