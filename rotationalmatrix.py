# rotational matrix
from numpy import *
z = array([0, 0, 1])
nz = array([1, 1, 1])

#axis = cross (z, nz)
print(nz)
imax=0

for i in range(1,len(nz)): 
    if abs(nz[i])>abs(nz[imax]) :
        imax=i

nx=nz.copy()        
nx[imax]=0
ny = cross(nz,nx)
nx = cross(ny,nz)
### normalize
nz=nz/linalg.norm(nz)
nx=nx/linalg.norm(nx)
ny=ny/linalg.norm(ny)
# create rotational matrix from n
R = transpose(array([nx, ny, nz]))
print(R)
print(R[-1,-1])
#R = np.array(())
