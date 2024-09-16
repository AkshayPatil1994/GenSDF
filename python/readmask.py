import numpy as np
import matplotlib.pyplot as plt

nx, ny, nz = 512, 256, 128
filename = 'mask.bin'
x = np.linspace(0,30,nx)
y = np.linspace(0,10,ny)
z = np.linspace(0,8,nz)

with open(filename, 'rb') as f:
    f.read(4) # skips 4 byte marker
    sdf = np.fromfile(f, count=nx*ny*nz)
sdf = np.reshape(sdf,(nx,ny,nz),order='F')

inside = np.zeros((nx,ny,nz),dtype=bool)

inside[sdf%2>0] = True
plt.pcolor(x,y,inside[:,:,64].T,cmap='turbo'); plt.colorbar()
plt.show()
