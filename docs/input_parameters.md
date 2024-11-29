## Input Parameters

The input parameters file `parameters.in` looks as shown below.

```
! Name and location of the inputfile
'data/armadillo_withnormals.obj'
! Scalar value (real), sdfresolution (int), progressbarsize (int)  
100.0 10 20
! nx, ny, nz (Computational Grid)
512 128 128
! r0 (Origin of the computational grid)
0.0 0.0 0.0
! non_uniform_grid
0
```

All lines preceeding with a `!` are comments and not used in the code
1. `'data/armadillo_withnormals.obj'` - Geometry file used to compute the SDF
2. `100.0 10 20` 
   `100.0` - Large positive value to be set for points outside the bounding box where the SDF is not calculated
   `10` - stencil width around the bounding box (number of grid points)
   `20` - Based on the size of your terminal you can adjust how wide the progress bar should be, usually a value of 20 to 30 is best.
3. `512 128 128` - Number of grid points in x, y, and z directions respectively
4. `0.0 0.0 0.0` - Origin of the computational domain as defined in CaNS.
5. `0` - If the grid is non-uniform in the vertical direction then the value should be `1`, for example, when using grid stretching.

When using a non-uniform-grid it is important to keep in mind that the grid will be read from the `data/grid.bin` file, so it is pertinent to replace the file to avoid conflict in the SDF.