import numpy as np
import matplotlib.pyplot as plt
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

# Example usage
nx, ny, nz = 512, 128, 128
filename1 = 'data/sdfu.bin'
# Load the grid for accurate location vtk write
[xp,yp,zp,xf,yf,zf] = read_grid(loc='data/',iprecision=8,ng=[nx,ny,nz],r0=[0.,0.,0.],non_uniform_grid = True)
# Load the data from binary file
with open(filename1, 'rb') as f:
    #f.read(4)  # Skip the 4-byte Fortran marker
    sdf = np.fromfile(f, count=nx * ny * nz)
sdf = np.reshape(sdf, (nx, ny, nz), order='F')

plt.pcolor(xf,yp,sdf[:,:,64].T,cmap='turbo')
plt.show()
