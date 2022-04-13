import matplotlib.pyplot as plt
import AstroTools as AT
from scipy.integrate import odeint
import numpy as np
import PIL
#from mpl_toolkits.mplot3d import Axes3D

mu = 398600 # Earth gravitational parameter [km^3/s^2]
Re = 6371   # Earth radius [km]
# Initial conditions (x[km], y[km], z[km], vx[km/s], vy[km/s], vz[km/s])
x0 = [8000, 0, 0, 0, 7.120, 3.022]

# Defining the EOM of the unperturbed 2BP
mu = 398600
def tbp(x, t):
    r = np.linalg.norm(x[:3])
    drdt = x[3:]
    dvdt = -mu/np.power(r, 3)*x[:3]
    dxdt = np.concatenate((drdt,dvdt))
    return dxdt

# Specific energ
x0 = [8000, 0, 0, 0, 7.120, 3.022]
E0=np.power(np.linalg.norm(x0[3:]),2)/2-mu/np.linalg.norm(x0[:3]);
# Semi-major axis
a0=-mu/2/E0;
# Orbital period
T = 2*np.pi*np.sqrt(np.power(a0,3)/mu)
# Timespan vector
t = np.linspace(0, T, 10000)
sol = odeint(tbp, x0, t)

# xtving the ODE
xt = odeint(tbp, x0, t)

# load Earth texture with PIL
# load bluemarble with PIL
bm = PIL.Image.open('map.jpg')
# Convert the image to array, and divide by 256 to get RGB values that matplotlib accept
# it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept
bm = np.array(bm)/256.

# coordinates of the image
# coordinates of the image - don't know if this is entirely accurate, but probably close
lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180
lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180

# Plot
plt.figure()
# Grid and black background
plt.rcParams['axes.grid'] = True
plt.style.use('dark_background')
ax = plt.axes(projection='3d') # 3D plot
# NON ho capito
ax = plt.axes(projection='3d')
x = np.outer(np.cos(lons), np.cos(lats)).T
y = np.outer(np.sin(lons), np.cos(lats)).T
z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
# Plot Earth
ax.plot_surface(x*6371, y*6371, z*6371, rstride=4, cstride=4, facecolors = bm)
# Plot orbit
ax.plot(xt[:,0], xt[:,1], xt[:,2],'r')
# Rescale the plot
ax.auto_scale_xyz([-9000, 9000], [-9000, 9000], [-9000, 9000])

lat, lon = AT.lat_long(xt[:,:3], t)
plt.figure()
plt.imshow(bm, extent=[-180,180,-90,90])
plt.plot(np.rad2deg(lon),np.rad2deg(lat),'r')
ax.plot(sol[:,0], sol[:,1], sol[:,2],'r')
plt.show()
