# rotational matrix
from numpy import *
eps=1e-9

def genR2():
    R2 = zeros((5,5))
            
    R2[0,0] = (3*R1(0,0)**2-1)/2
    R2[0,1] = sqrt(3)*R1(0,1)*R1(0,0)
    R2[0,-1] = sqrt(3)*R1(0,-1)*R1(0,0)
    R2[0,2] = sqrt(3)*(R1(0,1)**2-R1(0,-1)**2)/2
    R2[0,-2] = sqrt(3)*R1(0,1)*R1(0,-1)

    R2[1,0] = sqrt(3)*R1(1,0)*R1(0,0)
    R2[1,1] = R1(1,1)*R1(0,0)+R1(1,0)*R1(0,1)
    R2[1,-1] = R1(1,-1)*R1(0,0)+R1(1,0)*R1(0,-1)
    R2[1,2] = R1(1,1)*R1(0,1)-R1(1,-1)*R1(0,-1)
    R2[1,-2] = R1(1,1)*R1(0,-1)+R1(1,-1)*R1(0,1)

    R2[-1,0] = sqrt(3)*R1(-1,0)*R1(0,0)
    R2[-1,1] = R1(-1,1)*R1(0,0)+R1(-1,0)*R1(0,1)
    R2[-1,-1] = R1(-1,-1)*R1(0,0)+R1(-1,0)*R1(0,-1)
    R2[-1,2] = R1(-1,1)*R1(0,1)-R1(-1,-1)*R1(0,-1)
    R2[-1,-2] = R1(-1,1)*R1(0,-1)+R1(-1,-1)*R1(0,1)

    R2[2,0] = sqrt(3)*(R1(1,0)**2-R1(-1,0)**2)/2
    R2[2,1] = R1(1,1)*R1(1,0)-R1(-1,1)*R1(-1,0)
    R2[2,-1] = R1(1,-1)*R1(1,0)-R1(-1,-1)*R1(-1,0)
    R2[2,2] = (R1(1,1)**2-R1(1,-1)**2-R1(-1,1)**2+R1(-1,-1)**2)/2
    R2[2,-2] = R1(1,1)*R1(1,-1)-R1(-1,1)*R1(-1,-1)

    R2[-2,0] = sqrt(3)*R1(1,0)*R1(-1,0)
    R2[-2,1] = R1(1,1)*R1(-1,0)+R1(1,0)*R1(-1,1)
    R2[-2,-1] = R1(1,-1)*R1(-1,0)+R1(1,0)*R1(-1,-1)
    R2[-2,2] = R1(1,1)*R1(-1,1)-R1(1,-1)*R1(-1,-1)
    R2[-2,-2] = R1(1,1)*R1(-1,-1)+R1(1,-1)*R1(-1,1)

    return R2

def genR(nz):
#    imax=0
#    for i in range(1,len(nz)): 
#        if abs(nz[i])>abs(nz[imax]) :
#            imax=i
#    nx=nz.copy()        
#    nx[imax]=0

    imin=0
    for i in range(1,len(nz)): 
        if abs(nz[i])<abs(nz[imin]) :
            imin=i
    nx = zeros(3)
    nx[imin]=1
    ny = cross(nz,nx)
    nx = cross(ny,nz)
### normalize
    nz=nz/linalg.norm(nz)
    nx=nx/linalg.norm(nx)
    ny=ny/linalg.norm(ny)
# create rotational matrix from n
    R = transpose(array([nx, ny, nz]))
    return R

def delta(a,b):
    if a==b:
        return 1
    else :
        return 0
    
def R1(m, mp):
    i=0
    ip=0
    if m==0:
        i=2
    elif m==1:
        i=0
    elif m==-1:
        i=1
    else :
        print('Wrong m value: ', m)
    if mp==0:
        ip=2
    elif mp==1:
        ip=0
    elif mp==-1:
        ip=1
    else :
        print('Wrong mp value: ', mp)
    return R[i, ip]    

def Rl(l, m, mp):
    if l==1:
        return R1(m, mp)
    ul = u(l, m, mp)
    vl = v(l, m, mp)
    wl = w(l, m, mp)
    Ul = 0
    Vl = 0
    Wl = 0
    if ul!=0 : Ul = U(l, m, mp)
    if vl!=0 : Vl = V(l, m, mp)
    if wl!=0 : Wl = W(l, m, mp)
#    if abs(ul)>eps : Ul = U(l, m, mp)
#    if abs(vl)>eps : Vl = V(l, m, mp)
#    if abs(wl)>eps : Wl = W(l, m, mp)
#    print(ul, vl, wl)
#    print(ul==0, vl==0, wl==0)
#    print(Ul, Vl, Wl)
    return ul*Ul+vl*Vl+wl*Wl

def u(l, m, mp):
    if abs(mp)==l :
        return sqrt((l+m)*(l-m)/(2*l)/(2*l-1))
    else :
        return sqrt((l+m)*(l-m)/(l+mp)/(l-mp))

def v(l, m, mp):
    if abs(mp)==l :
        return sqrt((1+delta(m,0))*(l+abs(m)-1)*(l+abs(m))/(2*l)/(2*l-1))*(1-2*delta(m,0))/2
    else :
        return sqrt((1+delta(m,0))*(l+abs(m)-1)*(l+abs(m))/(l+mp)/(l-mp))*(1-2*delta(m,0))/2

def w(l, m, mp):
    if abs(mp)==l :
        return sqrt((l-abs(m)-1)*(l-abs(m))/(2*l)/(2*l-1))*(1-delta(m,0))/-2
    else :
        return sqrt((l-abs(m)-1)*(l-abs(m))/(l+mp)/(l-mp))*(1-delta(m,0))/-2

def P(l, i, mu, mp):
    if mp==l:
        return R1(i,1)*Rl(l-1, mu, l-1)-R1(i,-1)*Rl(l-1,mu,-l+1)
    elif mp==-l:
        return R1(i,1)*Rl(l-1, mu, -l+1)+R1(i,-1)*Rl(l-1,mu,l-1)
    else :
        return R1(i,0)*Rl(l-1, mu, mp)

def U(l, m, mp):
    return P(l, 0, m, mp)
#    if m==0:
#        return P(l,0,0,mp)
#    elif m>0:
#        return P(l,0,m,mp)
#    else :
#        return P(l,0,m,mp)

def V(l, m, mp):
    if m==0:
        return P(l,1,1,mp) + P(l,-1,-1,mp)
    elif m>0:
        return P(l,1,m-1,mp)*sqrt(1+delta(m,1))-P(l,-1,-m+1,mp)*(1-delta(m,1))
    else:
        return P(l,1,m+1,mp)*(1-delta(m,-1))+P(l,-1,-m-1,mp)*sqrt(1+delta(m,-1))

def W(l, m, mp):
    if m==0:
        return 0
    elif m>0:
        return P(l,1, m+1, mp) + P(l, -1, -m-1, mp)
    else:
        return P(l,1, m-1, mp) - P(l, -1, -m+1, mp)
    
z = array([0, 0, 1])
nz = array([1, 1, 1])
#nz = random.random(3)-0.5
#nz=array([ 0.02940479,  0.11130159, -0.47031108])
#axis = cross (z, nz)
print(nz)

R = genR(nz)
print(R)
nzInv = matmul(transpose(R), z)
print(nzInv)
#print (Rl(2,0,0))
#R = np.array(())

Rlmat = zeros((5,5))
for m in range(-2,3):
    for mp in range(-2,3):
        Rlmat[m,mp] = Rl(2,m,mp)
#        print (m,mp,Rl(2,m,mp))

R2 = genR2()
set_printoptions(suppress=True)
print(Rlmat)
print (' ')
print(R2)

#print(Rl(2,-2,-2))

print (' ')
print(abs(Rlmat-R2)<eps)

#print(Rl(2,1,2), R2[1,2])

## plot
from scipy import *
import scipy.special as special
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.mplot3d.axes3d as axes3d

n=2
mp=0
theta = linspace(0, pi, 150)
phi = linspace(0, 2*pi, 300)
Theta, Phi = meshgrid(theta, phi)

Yloc = zeros(shape(Theta)) # local rotated
for m in range (-n, n+1):
    Ylm = special.sph_harm(m,n, Phi, Theta)
    print(m, Rl(n,m,mp))
    if m==0 :
        Ylm = real(Ylm)
    elif m>0:
        Ylm = real(Ylm)*sqrt(2)*(-1)**m
    elif m<0:
        Ylm = imag(Ylm)*sqrt(2)*(-1)**m
    else :
        continue
        #Ylm = Ylm
    Yloc = Yloc + Rl(n,m,mp)*Ylm
        
#Yloc = imag(special.sph_harm(0,2, Phi, Theta))
RY = (abs(Yloc))
Z = RY*cos(Theta)
X = RY*sin(Theta)*cos(Phi)
Y = RY*sin(Theta)*sin(Phi)
C = (sign(Yloc))
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
norm = cm.colors.Normalize()
my_col = cm.jet(sign(Yloc))
#my_col = cm.jet((Ylm-amin(Ylm))/(amax(Ylm)-amin(Ylm)))
#my_col = cm.winter((C-amin(C))/(amax(C)-amin(C)))

plot = ax.plot_surface(
    X, Y, Z, facecolors=my_col, rstride=1, cstride=1,
    linewidth=0, antialiased=True, alpha=0.8, shade=False)
#, cmap=plt.get_cmap('jet')
#ax.set_aspect(aspect='equal')
#ax.pbaspect = [1.0, 1, 1]
#fig.colorbar(plot,shrink=0.5, aspect=5)

plotline = ax.plot([0,nz[0]], [0,nz[1]], [0,nz[2]]) 
#plotline = ax.plot([0,nzInv[0]], [0,nzInv[1]], [0,nzInv[2]]) 
#plotline = ax.plot([0,nz[0]], [0,nz[1]], [0,nz[2]]) 
ax.auto_scale_xyz([-1, 1], [-1, 1], [-1, 1])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

