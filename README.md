
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

####Constants####

m_shuttle = 12000
a_shuttle = 4 * np.pi
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
    cs = CubicSpline(heights, densities)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(heights, densities, 'o', label='data')
    ax.plot(xs, cs(xs))




    ax.legend(loc='upper right', ncol=1)
    plt.show()
    
    return cs


def air_resistance(density, a, cd, velocity):
    
    resistance =  0.5 * density * a * cd * velocity ** 2
    
    return resistance

def proj_distance(angle):
    
    dist = earth_radius * angle
    
    return dist


def forces(v, angle, density):
    
    vertical_forces = m_shuttle * g - air_resistance(density, a_shuttle, cd_shuttle, v)* np.cos(angle) - air_resistance(density, a_shuttle, cl_shuttle, v)
    horizontal_forces = - air_resistance(density, a_shuttle, cd_shuttle, v) * np.sin(angle)
    
    return vertical_forces, horizontal_forces
    

def forces_chute(v, angle, density):
    
    forces = m_shuttle * g - air_resistance(density, a_shuttle, cd_shuttle, v)* np.cos(angle) - air_resistance(density, a_parachute, cd_parachute, v)
    
    return forces



def eulerFalling(v, dh, m, resistance):
    
    

    return (g*dh) - ((resistance/m) *  dh) + v










def eulerFalling_newton(v0, dh, m, density):
    """
    Método de Euler Regressivo utilizando o Método de Newton-Raphson.
    
    Args:
        v0: Velocidade inicial.
        dh: Passo de tempo.
        m: Massa do objeto.
        density: Densidade do ar.
    
    Returns:
        v: Nova velocidade após o passo de tempo dh.
    """
    
    def f(v):
        resistance = air_resistance(density, a_shuttle, cd_shuttle, v) + air_resistance(density, a_shuttle, cl_shuttle, v)
        return v - v0 - dh * (g - resistance / m)
    
    def df_dv(v):
        d_resistance_dv = (density * a_shuttle * cd_shuttle * v + density * a_shuttle * cl_shuttle * v)
        return 1 - dh * (-d_resistance_dv / m)
    
    v_curr = v0
    tol = dh*0.01
    max_iter = 3000
    for i in range(max_iter):
        f_val = f(v_curr)
        df_val = df_dv(v_curr)
        v_next = v_curr - f_val / df_val
        if abs(v_next - v_curr) < tol:
            break
        v_curr = v_next
    
    return v_curr

def reentry(v0, dh, htotal, method):
    cs = interpolation()

    size = int(htotal / dh + 1)
    out = np.zeros((size, 2), dtype=float) # Output numpy array
    out[0] = [htotal, v0]

    v = v0
    for i in range(1, size):
        h = htotal - i * dh
        density = cs(h)
        if density < 0:
            density = 0
        
        if method == 'euler':
            r = air_resistance(density, a_shuttle, cd_shuttle, v) + air_resistance(density, a_shuttle, cl_shuttle, v)
            print(r)
            v = eulerFalling(v, dh, m_shuttle, r)
        elif method == 'newton':
            v = eulerFalling_newton(v, dh, m_shuttle, density)
        
        out[i] = [h, v]
      
    return out

# Comparação dos dois métodos
result_euler = reentry(1000, 0.05, 1300, method='euler')
result_newton = reentry(1000, 0.05, 1300, method='newton')

fig, ax = plt.subplots()
ax.plot(result_euler[:, 0], result_euler[:, 1], linewidth= 3, label='Euler Progressivo')
ax.plot(result_newton[:, 0], result_newton[:, 1], label='Euler Regressivo')
ax.set_xlabel('x/m')
ax.set_ylabel('v/$ms^{-1}$')
ax.legend()
plt.title('Comparação de Métodos de Reentrada')
plt.show()
