# SMCEF_P2
import numpy as np
import scipy
import matplotlib.pyplot as plt

####Constants####

m_shuttle = 12000
a_shuttle = 4 * scipy.pi * m_shuttle**2
cd_shuttle = 1.2
cl_shuttle = 1

a_parachute = 301
cd_parachute = 1

g = 10
earth_radius = 6371000

heights = [-1000, 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
             15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000]

densities = [1.347, 1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900,
             0.5258, 0.4671, 0.4135, 0.1948, 0.08891, 0.04008, 0.01841, 0.003996,
             0.001027, 0.0003097, 0.00008283, 0.00001846]

def interpolation():
    xs = np.arange(-1000, 140000)
    cs = scipy.interpolate.CubicSpline(heights, densities)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(heights, densities, 'o', label='data')
    ax.plot(xs, cs(xs))




    ax.legend(loc='upper right', ncol=1)
    plt.show()
    
    return cs


def air_resistance(density, a, cd, velocity):
    
    resistance = - 0.5 * density * a * cd * velocity**2
    
    return resistance

def proj_distance(angle):
    
    dist = earth_radius * angle
    
    return dist

def eulerFalling(v, dh, m, resistance):
    
    

    return g*dh - resistance/m * dh + v

def reentry(v0, dh, htotal):
      
    cs = interpolation()

    size = int(htotal / dh + 1)
    #Number of steps to take plus one to have the initial value
    
    out = np.zeros((size, 2), dtype=float) #Output numpy array
    out[0] = [130000, v0]

    v = v0
    for i in range(size):
        h = size - i
        density = cs(h)
        if density < 0:
            density = 0
        
        r = air_resistance(density, a_shuttle, cd_shuttle, v) + air_resistance(density, a_shuttle, cl_shuttle, v)
        
        v = eulerFalling(v, dh, m_shuttle, r)
        
        out[i] = [h, v]
      
    return out


result = reentry(10000, 1, 130000)
fig, ax = plt.subplots()
ax.plot(result[:, 0], result[:, 1])
ax.set_xlabel('x/m')
ax.set_ylabel('v/$ms^{-1}$')
pass
