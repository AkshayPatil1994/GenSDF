
import numpy as np

def write_vtk(array, filename, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    """
    Write a 3D NumPy array to a VTK file in structured points format.

    Parameters:
    array (numpy.ndarray): The 3D array of scalar values.
    filename (str): The name of the output VTK file.
    spacing (tuple): The spacing between points in the grid.
    origin (tuple): The origin of the grid.

    """
    nx, ny, nz = array.shape
    with open(filename, 'w') as f:
        # Write VTK header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("VTK file containing 3D scalar data\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        f.write(f"ORIGIN {origin[0]} {origin[1]} {origin[2]}\n")
        f.write(f"SPACING {spacing[0]} {spacing[1]} {spacing[2]}\n")
        f.write(f"POINT_DATA {nx*ny*nz}\n")
        f.write("SCALARS scalar_data float 1\n")
        f.write("LOOKUP_TABLE default\n")

        # Write the array data
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    f.write(f"{array[x, y, z]}\n")

# Example usage
nx, ny, nz = 512, 128, 128
filename1 = '../src/data/mask.bin'
filename2 = '../src/data/narrowbandpoints.bin'

x = np.linspace(32/nx,32,nx)
y = np.linspace(8/nx,8,ny)
z = np.linspace(8/nx,8,nz)

# Load the data from binary file 
with open(filename1, 'rb') as f:
    f.read(4)  # Skip the 4-byte Fortran marker
    sdf = np.fromfile(f, count=nx * ny * nz)
sdf = np.reshape(sdf, (nx, ny, nz), order='F')

# Load the data from binary file 
with open(filename2, 'rb') as f:
    f.read(4)  # Skip the 4-byte Fortran marker
    # First read the number of points
    numpoints = np.fromfile(f,count=1,dtype=np.int32)
    print("Total number of points:",numpoints)
    # Now read the points
    points = np.zeros((3,numpoints[0]),dtype=np.int32)
    for i in range(0,numpoints[0]):
        dummy = np.fromfile(f, count=3,dtype=np.int32)
        points[0,i] = dummy[0]
        points[1,i] = dummy[1]
        points[2,i] = dummy[2]        
points[1,0] = points[1,1]
# Write the array to VTK file
point_indices = np.zeros((nx,ny,nz),dtype=np.int32)
for i in range(0,numpoints[0]):
    point_indices[points[2,i],points[1,i],points[0,i]] = 1

# Write data to file
write_vtk(sdf, "data/sdf.vtk", spacing=(32/nx, 8/ny, 8/nz))
write_vtk(point_indices, "data/points.vtk", spacing=(32/nx, 8/ny, 8/nz))
