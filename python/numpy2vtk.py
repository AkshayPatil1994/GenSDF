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
nx, ny, nz = 512, 256, 128
filename = 'mask.bin'

x = np.linspace(0,30,nx)
y = np.linspace(0,10,ny)
z = np.linspace(0,8,nz)

# Load the data from binary file (assuming float32 precision)
with open(filename, 'rb') as f:
    f.read(4)  # Skip the 4-byte Fortran marker
    sdf = np.fromfile(f, count=nx * ny * nz)
sdf = np.reshape(sdf, (nx, ny, nz), order='F')  # Reshape with Fortran ordering
# inside = np.zeros((nx,ny,nz),dtype=np.float64)
# inside[sdf%2>0] = 1.0
# Write the array to VTK file
write_vtk(sdf, "output.vtk", spacing=(x[-1]/nx, y[-1]/ny, z[-1]/nz))
