import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from numpy import cross, eye, dot
from scipy.linalg import expm, norm
import os

# Constants
G = 4522.022248
M = 1.0 # [1e11 solar mass]
e = 0.5
Rmin = 25.0 # min distance between galaxtic cores [kpc]
Rmax = Rmin / (1 - e) * (1 + e) # max distance between galaxtic cores [kpc]
Ra = Rmin / 2 / (1 - e) # semi long axis [kpc]
epsilon = 0.2*Rmin
nstep = 2668
dt = 1e-2
t = -11.68
nparticles = 342

# Rotation function
def RotMatrix(axis, theta):
    return expm(cross(eye(3), axis/norm(axis)*theta))

# Galaxy 1 Position Disk Setup
r1 = []
x1 = np.zeros((3, nparticles))
x1core = np.array([0.0, -Rmax / 2, 0.0])
a = 0.2
for i in range(12):
    r1 = np.append(r1, a*Rmin)
    a += 0.05

n_particles = 12
k = 0
for i in range(12):
    for j in range(n_particles):
        theta = j*2*np.pi/n_particles
        x1[0][k] = r1[i]*np.cos(theta)
        x1[1][k] = r1[i]*np.sin(theta)
        k += 1
    n_particles += 3

# Galaxy 1 Disk Particles Velocity Setup
v1 = np.zeros((3, nparticles))
v1_r = np.sqrt((G*M*r1)/(r1**2+epsilon**2))

n_particles = 12
k = 0
for i in range(12):
    for j in range(n_particles):
        theta = np.arctan2(x1[1][k], x1[0][k])
        v1[0][k] = v1_r[i] * np.cos(theta+np.pi/2)
        v1[1][k] = v1_r[i] * np.sin(theta+np.pi/2)
        k += 1
    n_particles += 3

# Galaxy 1 Peri Angle (w) Rotation
n1, axis1, w1 = [0, 1, 0], [0, 0, -1], -np.pi/6
M1 = RotMatrix(axis1, w1)
n1 = dot(M1, n1)

# Galaxy 1 Inclination Angle (i) Rotation
i1 = np.pi/3
rot_matrix1 = RotMatrix(-n1, i1)
x1 = dot(rot_matrix1, x1)
v1 = dot(rot_matrix1, v1)

# Galaxy 1 Core Velocity
v1core = np.array([np.sqrt((G*M/(4*Ra))*((1-e)/(1+e))), 0.0, 0.0])
v1[0] += v1core[0]

# Galaxy 2 Position Disk Setup
r2 = []
x2 = np.zeros((3, nparticles))
x2core = np.array([0.0, Rmax / 2, 0.0])
a = 0.2
for i in range(12):
    r2 = np.append(r2, a*Rmin)
    a += 0.05

n_particles = 12
k = 0
for i in range(12):
    for j in range(n_particles):
        theta = j*2*np.pi/n_particles
        x2[0][k] = r2[i]*np.cos(theta)
        x2[1][k] = r2[i]*np.sin(theta)
        k += 1
    n_particles += 3

# Galaxy 2 Disk Particles Velocity Setup
v2 = np.zeros((3, nparticles))
v2_r = np.sqrt((G*M*r2)/(r2**2+epsilon**2))

n_particles = 12
k = 0
for i in range(12):
    for j in range(n_particles):
        theta = np.arctan2(x2[1][k], x2[0][k])
        v2[0][k] = v2_r[i] * np.cos(theta+np.pi/2)
        v2[1][k] = v2_r[i] * np.sin(theta+np.pi/2)
        k += 1
    n_particles += 3
    
# Galaxy 2 Peri Angle (w) Rotation
n2, axis2, w2 = [0, -1, 0], [0, 0, -1], -np.pi/6
M2 = RotMatrix(axis2, w2)
n2 = dot(M2, n2)

# Galaxy 2 Inclination Angle (i) Rotation
i2 = np.pi/3
rot_matrix2 = RotMatrix(-n2, i2)
x2 = dot(rot_matrix2, x2)
v2 = dot(rot_matrix2, v2)

# Galaxy 2 Core Velocity
v2core = np.array([-np.sqrt((G*M/(4*Ra))*((1-e)/(1+e))), 0.0, 0.0])
v2[0] += v2core[0]

# Positioning the galaxies (starting from the apocenter)
# Galaxy 1 and Galaxy 2 are Rmax apart in the y directions
for i in range(nparticles):
    x1[:,i] += x1core
    x2[:,i] += x2core

# Leapfrog Integration
def accel_core(x1, x2):
    r = np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
    core_vector = np.array([x1[0]-x2[0], x1[1]-x2[1], x1[2]-x2[2]])
    a = - G * M * core_vector/r**3
    
    return a
    
def accel_particle(x1, x2, x1core, x2core):
    # Particles in Galaxy 1
    r1_1 = np.sqrt((x1[0]-x1core[0])**2 + (x1[1]-x1core[1])**2 + (x1[2]-x1core[2])**2 + epsilon**2)
    r1_1_vector = np.array([x1[0]-x1core[0], x1[1]-x1core[1], x1[2]-x1core[2]])
    a1_1 = - G * M * r1_1_vector/r1_1**3
    r1_2 = np.sqrt((x1[0]-x2core[0])**2 + (x1[1]-x2core[1])**2 + (x1[2]-x2core[2])**2 + epsilon**2)
    r1_2_vector = np.array([x1[0]-x2core[0], x1[1]-x2core[1], x1[2]-x2core[2]])
    a1_2 = - G * M * r1_2_vector/r1_2**3
    a1_final = a1_1 + a1_2
    
    # Particles in Galaxy 2
    r2_1 = np.sqrt((x2[0]-x1core[0])**2 + (x2[1]-x1core[1])**2 + (x2[2]-x1core[2])**2 + epsilon**2)
    r2_1_vector = np.array([x2[0]-x1core[0], x2[1]-x1core[1], x2[2]-x1core[2]])
    a2_1 = - G * M * r2_1_vector/r2_1**3
    r2_2 = np.sqrt((x2[0]-x2core[0])**2 + (x2[1]-x2core[1])**2 + (x2[2]-x2core[2])**2 + epsilon**2)
    r2_2_vector = np.array([x2[0]-x2core[0], x2[1]-x2core[1], x2[2]-x2core[2]])
    a2_2 = - G * M * r2_2_vector/r2_2**3
    a2_final = a2_1 + a2_2
        
    return a1_final, a2_final

def printstate(x1, x2, x1core, x2core, v1, v2, v1core, v2core, t):
    # Print header
    print("686 " + str(t))
    
    # Print cores
    print('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1, x1core[0], x1core[1], x1core[2], v1core[0], v1core[1], v1core[2]))
    print('{:d}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1, x2core[0], x2core[1], x2core[2], v2core[0], v2core[1], v2core[2]))
    
    # Print particles
    for i in range(nparticles):
        print('{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1e-20, x1[0][i], x1[1][i], x1[2][i], v1[0][i], v1[1][i], v1[2][i]))
    for i in range(nparticles):
        print('{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}'.format(1e-20, x2[0][i], x2[1][i], x2[2][i], v2[0][i], v2[1][i], v2[2][i]))
    
# Print initial state
printstate(x1, x2, x1core, x2core, v1, v2, v1core, v2core, t)

for i in range(nstep):
    # Next timestep
    t += dt
    
    # Velocity half time step
    core_accel = accel_core(x1core, x2core)
    v1core += core_accel * dt / 2
    v2core += -core_accel * dt / 2
    
    a1_particle_accel, a2_particle_accel = accel_particle(x1, x2, x1core, x2core)
    v1 += a1_particle_accel * dt / 2
    v2 += a2_particle_accel * dt / 2
    
    # Position full time step
    x1core += v1core * dt
    x2core += v2core * dt
    x1 += v1 * dt
    x2 += v2 * dt
    
    # Velocity half time step
    core_accel = accel_core(x1core, x2core)
    v1core += core_accel * dt / 2
    v2core += -core_accel * dt / 2
    
    a1_particle_accel, a2_particle_accel = accel_particle(x1, x2, x1core, x2core)
    v1 += a1_particle_accel * dt / 2
    v2 += a2_particle_accel * dt / 2
    
    # Print state
    if i % 3 == 0:
        printstate(x1, x2, x1core, x2core, v1, v2, v1core, v2core, t)
    
# Print final state
printstate(x1, x2, x1core, x2core, v1, v2, v1core, v2core, t)