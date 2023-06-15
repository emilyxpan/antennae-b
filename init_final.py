import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 1
M = 20/3 # total mass of one galaxy
e = 0.5
Rp = 6/2
Ra = Rp * (1+e) / (1-e)
a = Rp/(1-e) # semi major long axis
Tfactor = 18
t_start = 1.0 * Tfactor
t_end = 0.0 * Tfactor 
t = t_start
dt = -1e-3
nsteps = int((t_end-t_start)/dt)

x1 = np.array([0.0, -Rp, 0.0])
x2 = np.array([0.0, Rp, 0.0])

v1 = np.array([np.sqrt((G*M/(4*a))*((1+e)/(1-e))), 0.0, 0.0])
v2 = np.array([-np.sqrt((G*M/(4*a))*((1+e)/(1-e))), 0.0, 0.0])

def accel(x1, x2):
    r = np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
    core_vector = np.array([x1[0]-x2[0], x1[1]-x2[1], x1[2]-x2[2]])
    acceleration = - G * M * core_vector/r**3
    
    return acceleration

def printstate(x1core, x2core, v1core, v2core, t):
    # Print header
    print("2 " + str(t))
    
    # Print cores
    print('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1, x1core[0], x1core[1], x1core[2], v1core[0], v1core[1], v1core[2]))
    print('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1, x2core[0], x2core[1], x2core[2], v2core[0], v2core[1], v2core[2]))

for i in range(nsteps):
    # Next timestep
    t += dt
    
    # Velocity half time step
    core_accel = accel(x1, x2)
    v1 += core_accel * dt / 2
    v2 -= core_accel * dt / 2
    
    # Position full time step
    x1 += v1 * dt
    x2 += v2 * dt
    
    # Velocity half time step
    core_accel = accel(x1, x2)
    v1 += core_accel * dt / 2
    v2 -= core_accel * dt / 2
        

print("galaxy 1 position: " + str(x1)) 
print("galaxy 1 velocity: " + str(v1))
print("galaxy 2 position: " + str(x2))
print("galaxy 2 velocity: " + str(v2))

# Results: 
'''
galaxy 1 position: [-4.6591016   5.65650751  0.        ]
galaxy 1 velocity: [-0.16545881 -0.38691876  0.        ]
galaxy 2 position: [ 4.6591016  -5.65650751  0.        ]
galaxy 2 velocity: [0.16545881 0.38691876 0.        ]
'''