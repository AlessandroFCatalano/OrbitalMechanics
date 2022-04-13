import numpy as np
import numpy.ma as ma

om_E = np.deg2rad(15.04/3600)
theta_G0 = 0

def wrapToPi(angle):
    """
    Wrap the angle between -pi and pi.
    Args:
        angle (float): angle to wrap.
    Returns:
         The wrapped angle.
    """
    angle=np.mod(angle, 2*np.pi)
    for i in range(len(angle)):
        if angle[i] > np.pi:
            angle[i] -= 2*np.pi
        if angle[i] < -np.pi:
            angle[i] += 2*np.pi
    return angle

def lat_long(x,t):
    delta = [None]*len(x[:,0])
    alpha = [None]*len(x[:,0])
    for i in range(len(x[:,0])):
        delta[i]=np.arcsin(x[i,2]/np.linalg.norm(x[i,:])) # Declination
    lat = delta;
    for i in range(len(x[:,0])):
        if x[i,1] > 0:
            alpha[i]=np.arccos(x[i,0]/np.linalg.norm(x[i,:])/np.cos(delta[i]))
        else:
            alpha[i]=2*np.pi-np.arccos(x[i,0]/np.linalg.norm(x[i,:])/np.cos(delta[i]))

    theta_G = theta_G0+om_E*t
    lon2pi = alpha-theta_G                  # longhitude expressed in [0; 2pi]
    lon = wrapToPi(lon2pi)            # longhitude expressed in [-pi; pi]
    # Add NaN to avoid plot discontinuities at the boundaries:
    for i in range(len(lon)):
        if lon[i-1]-lon[i]>np.pi:
            lon[i]=ma.masked
    return lat, lon
