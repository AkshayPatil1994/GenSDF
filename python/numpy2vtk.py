
import numpy as np
import os 
import sys

def read_grid(loc='data/',iprecision=8,ng=[10,10,10],r0=[0.,0.,0.],non_uniform_grid = False):
    '''
        This function reads the grid information generated within CaNS
    INPUT
        loc:                [string] Location where the grid information is saved
        iprecision          [integer] Precision used for the arrays
        ng:                 [3 x 1 -- list] Size of the grid in x, y, and z
        r0:                 [3 x 1 -- list] Location of the origin
        non_uniform_grid:   [Boolean] Is the grid non uniform
    OUTPUT
        xp, yp, zp:         [Numpy arrays] Cell-Center grid
        xu, yv, zw:         [Numpy arrays] Cell-Face grid
    '''
    #
    # Check if the file directory exists
    #
    filestat = os.path.exists(loc)
    if(filestat == False):
        sys.exit("The input file at %s does not exist!"%(loc))
    #
    # setting up some parameters
    #
    r0 = np.array(r0) # domain origin
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read geometry file
    #
    geofile  = loc+"geometry.out"
    geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = geo[0,:].astype('int')
    l  = geo[1,:]
    dl = l/(1.*ng)
    #
    # read and generate grid
    #
    xp = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0]) # centered  x grid
    yp = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1]) # centered  y grid
    zp = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2]) # centered  z grid
    xu = xp + dl[0]/2.                              # staggered x grid
    yv = yp + dl[1]/2.                              # staggered y grid
    zw = zp + dl[2]/2.                              # staggered z grid
    if(non_uniform_grid):
        grdfile = loc+'grid.bin'                    # Specify the location of the grid file
        f   = open(grdfile,'rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + np.transpose(grid_z[:,2]) # centered  z grid
        zw = r0[2] + np.transpose(grid_z[:,3]) # staggered z grid

    return xp, yp, zp, xu, yv, zw  

def write_vtk(array, filename, x, y, z):
    """
    Write a 3D NumPy array to a VTK file in structured grid format with arbitrary x, y, z coordinates.

    Parameters:
    array (numpy.ndarray): The 3D array of scalar values (nx x ny x nz).
    filename (str): The name of the output VTK file.
    x (numpy.ndarray): The 1D array of x-coordinates (size nx).
    y (numpy.ndarray): The 1D array of y-coordinates (size ny).
    z (numpy.ndarray): The 1D array of z-coordinates (size nz).
    """
    nx, ny, nz = array.shape
    assert len(x) == nx, "Length of x must match the first dimension of array"
    assert len(y) == ny, "Length of y must match the second dimension of array"
    assert len(z) == nz, "Length of z must match the third dimension of array"

    with open(filename, 'w') as f:
        # Write VTK header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK file containing 3D scalar data\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")

        # Write the coordinates of the grid points
        f.write(f"POINTS {nx * ny * nz} float\n")
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f"{x[i]} {y[j]} {z[k]}\n")

        # Write the scalar field data
        f.write(f"POINT_DATA {nx * ny * nz}\n")
        f.write("SCALARS scalar_data float 1\n")
        f.write("LOOKUP_TABLE default\n")
        
        # Flatten the array and write the scalar values
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f"{array[i, j, k]}\n")

# Example usage
nx, ny, nz = 512, 128, 128
filename1 = '../src/data/mask.bin'
# Load the grid for accurate location vtk write
[xp,yp,zp,xf,yf,zf] = read_grid(loc='data/',iprecision=8,ng=[10,10,10],r0=[0.,0.,0.],non_uniform_grid = False)
# Load the data from binary file 
with open(filename1, 'rb') as f:
    f.read(4)  # Skip the 4-byte Fortran marker
    sdf = np.fromfile(f, count=nx * ny * nz)
sdf = np.reshape(sdf, (nx, ny, nz), order='F')
# Write data to file
write_vtk(sdf, "data/sdf.vtk", xf, yp, zp)
