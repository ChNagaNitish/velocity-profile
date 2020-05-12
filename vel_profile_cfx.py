import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

y_len = 0.025695
z_len = 0.035
delta = 0.0015  # Boundary layer
grid = 0.0001   # Spacing between points

#   Inlet conditions 
u_max = 1570.19
gamma = 1.32
p_tot = 1.3e6
t_tot = 2000
m_max = 2.5
u_sonic = u_max/m_max

Ny = round(y_len/grid)-1
Nz = round(z_len/grid)-1
N = (Ny+1)*(Nz+1)

x = np.zeros((N))
y = np.zeros((Ny+1))
z = np.zeros((Nz+1))

for i in range(0,Ny):
    y[i] = i*grid
y[Ny] = y_len

for i in range(0,Nz):
    z[i] = i*grid
z[Nz] = z_len

vel = np.zeros((N,3))
vel[:,0] = u_max
co_ord = np.zeros((N,3))

for i in range(0,Nz+1):
    if z[i]<=(z_len-delta):
        for j in range(0,Ny+1):
            if y[j]>=(y_len-delta):
                vel[i*(Ny+1)+j][0] = u_max*((y_len-y[j])/delta)**(1/7)
            co_ord[i*(Ny+1)+j][1] = y[j]
            co_ord[i*(Ny+1)+j][2] = z[i]
    else:
        for j in range(0,Ny+1):
            vel[i*(Ny+1)+j][0] = u_max*((z_len-z[i])/delta)**(1/7)
            if (y_len-y[j])<(z_len-z[i]):
                vel[i*(Ny+1)+j][0] = u_max*((y_len-y[j])/delta)**(1/7)
            co_ord[i*(Ny+1)+j][1] = y[j]
            co_ord[i*(Ny+1)+j][2] = z[i]


result = np.zeros((N,8))
result[:,0] = co_ord[:,0]
result[:,1] = co_ord[:,1]
result[:,2] = co_ord[:,2]
result[:,3] = vel[:,0]
result[:,4] = vel[:,1]
result[:,5] = vel[:,2]

for i in range(0,N):
    m = vel[i][0]/u_sonic
    denom = 1+(gamma-1)*0.5*m**2
    result[i][6] = p_tot/denom**(gamma/(gamma-1))
    result[i][7] = t_tot/denom

np.savetxt("vel_profile.dat",result, delimiter=",")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(result[:,1],result[:,2],result[:,3],c=result[:,3])
ax.set_xlabel('y-direction')
ax.set_ylabel('z-direction')
ax.set_zlabel('Velocity - x component')
plt.savefig('vel_profile.png')