# scipy

from scipy import *
import scipy.special as special
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.mplot3d.axes3d as axes3d
#import scipy as sp
#import numpy as np
n=3
m=-2
theta = linspace(0, pi, 50)
phi = linspace(0, 2*pi, 100)
Theta, Phi = meshgrid(theta, phi)

Ylm = special.sph_harm(m,n, Phi, Theta)
if m>0 :
    Ylm = real(Ylm)
else :
    if m < 0 :
        Ylm = imag(Ylm)

    
        

R = (abs(Ylm))
Z = R*cos(Theta)
X = R*sin(Theta)*cos(Phi)
Y = R*sin(Theta)*sin(Phi)
C = (sign(Ylm))
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
norm = cm.colors.Normalize()
my_col = cm.jet(norm(Ylm))
#my_col = cm.jet((Ylm-amin(Ylm))/(amax(Ylm)-amin(Ylm)))
#my_col = cm.winter((C-amin(C))/(amax(C)-amin(C)))

plot = ax.plot_surface(
    X, Y, Z, facecolors=my_col, rstride=1, cstride=1,
    linewidth=0, antialiased=True, alpha=0.8, shade=False)

#, cmap=plt.get_cmap('jet')
#ax.set_aspect(aspect='equal')
#ax.pbaspect = [1.0, 1, 1]
#fig.colorbar(plot,shrink=0.5, aspect=5)
ax.auto_scale_xyz([-1, 1], [-1, 1], [-1, 1])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
