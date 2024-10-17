## Portions of this code taken from CaNS-world
# Credits - P. Costa
import numpy as np
import os

# USER INPUT PARAMETERS
iseek = 0                                               # Number of bytes to skip relative to the origin of the binary file (0 for Fortran binary)
iprecision = 8                                          # Precision of real-valued data (8 bytes for double precision)
r0 = np.array([0.0, 0.0, 0.0])                          # Domain origin
non_uniform_grid = False                                # Boolean flag for non_uniform_grid
geofile = "data/geometry.out"                           # Name of the geometry file
grid_bin_file = "data/grid.bin"                         # Grid binary data

# Preliminary checks (Change location of file write if needed below)
xgridfile = "data/xgrid.bin"
ygridfile = "data/ygrid.bin"
zgridfile = "data/zgrid.bin"
my_dtype = 'float64' if iprecision == 8 else 'float32'

# Query for geometry file
if not os.path.exists(geofile):
    raise FileNotFoundError(f"Geometry file '{geofile}' not found!")
if(non_uniform_grid):
    if not os.path.exists(grid_bin_file):
        raise FileNotFoundError(f"Grid binary file '{grid_bin_file}' not found!")

# Read the grid data and dimensions
geo_data = np.loadtxt(geofile, comments="!", max_rows=2)
ng = geo_data[0].astype(int)        
l = geo_data[1]               
dl = l / ng   

# Define grid points
x = np.arange(r0[0] + dl[0] / 2.0, r0[0] + l[0], dl[0])
y = np.arange(r0[1] + dl[1] / 2.0, r0[1] + l[1], dl[1])
z = np.arange(r0[2] + dl[2] / 2.0, r0[2] + l[2], dl[2])

# Create grid files for x, y, z if non-uniform grid, update z based on the actual grid
if non_uniform_grid:
    with open(grid_bin_file, 'rb') as f:
        grid_z = np.fromfile(f, dtype=my_dtype).reshape((ng[2], 4), order='F')
        z = r0[2] + grid_z[:, 2]

# Write x, y, z grid files
x.astype(my_dtype).tofile(xgridfile)
y.astype(my_dtype).tofile(ygridfile)
z.astype(my_dtype).tofile(zgridfile)
