import math
from vpython import *

n = int(input("How many objects? "))

G = 6.674e-11  # Gravitational constant

m, d, radi = {}, {}, {}  # mass, density, radius
x, y, z = {}, {}, {}  # positions
rx, ry, rz = {}, {}, {}  # displacement components
r = {}  # displacement magnitude
ox, oy, oz = {}, {}, {}  # unit vectors
F, Fx, Fy, Fz = {}, {}, {}, {}  # forces
vx, vy, vz = {}, {}, {}  # velocities
ax, ay, az = {}, {}, {}  # accelerations
Fnetx, Fnety, Fnetz = {}, {}, {}  # net forces
t = 0  # time
dt = 0.1  # time step
objects = {}
running = True  # Simulation state

# Input
for i in range(n):
    m[i] = float(input(f"Enter the mass of object: "))
    d[i] = float(input(f"Enter the density of object: "))
    radi[i] = ((3 * m[i]) / (4 * math.pi * d[i])) ** (1 / 3)
    x[i] = float(input(f"Enter its x coordinate: "))
    y[i] = float(input(f"Enter its y coordinate: "))
    z[i] = float(input(f"Enter its z coordinate: "))
    vx[i] = float(input(f"Enter the x component of its velocity: "))
    vy[i] = float(input(f"Enter the y component of its velocity: "))
    vz[i] = float(input(f"Enter the z component of its velocity: "))

# Spheres
for i in range(n):
    objects[i] = sphere(pos=vector(x[i], y[i], z[i]), radius=radi[i], color=color.red, make_trail=True)

# Pause function
def Pause():
    global running  
    running = not running  
    if running == True:
        pausebut.text = "Pause"
    elif running == False:
        pausebut.text = "Play"

# Pause button
pausebut = button(bind=Pause, text="Pause")

# Time loop
while True:
    rate(100)  # FPS control
    if running:
        for i in m:
            for j in m:
                if i != j:
                    # Displacement vector
                    rx[(i, j)] = x[j] - x[i]
                    ry[(i, j)] = y[j] - y[i]
                    rz[(i, j)] = z[j] - z[i]
                    r[(i, j)] = math.sqrt(rx[(i, j)] ** 2 + ry[(i, j)] ** 2 + rz[(i, j)] ** 2)

                    # Unit vector
                    ox[(i, j)] = rx[(i, j)] / r[(i, j)]
                    oy[(i, j)] = ry[(i, j)] / r[(i, j)]
                    oz[(i, j)] = rz[(i, j)] / r[(i, j)]

                    # Force vector
                    F[(i, j)] = (G * m[i] * m[j]) / (r[(i, j)] ** 2)

                    Fx[(i, j)] = F[(i, j)] * ox[(i, j)]
                    Fy[(i, j)] = F[(i, j)] * oy[(i, j)]
                    Fz[(i, j)] = F[(i, j)] * oz[(i, j)]

                    Fx[(j, i)] = -Fx[(i, j)]
                    Fy[(j, i)] = -Fy[(i, j)]
                    Fz[(j, i)] = -Fz[(i, j)]

        # Net force
        for i in m:
            Fnetx[i] = sum(Fx[(i, j)] for j in m if i != j)
            Fnety[i] = sum(Fy[(i, j)] for j in m if i != j)
            Fnetz[i] = sum(Fz[(i, j)] for j in m if i != j)

            # Acceleration
            ax[i] = Fnetx[i] / m[i]
            ay[i] = Fnety[i] / m[i]
            az[i] = Fnetz[i] / m[i]

            # Velocity
            vx[i] += ax[i] * dt
            vy[i] += ay[i] * dt
            vz[i] += az[i] * dt

            # Position
            x[i] += vx[i] * dt
            y[i] += vy[i] * dt
            z[i] += vz[i] * dt

            objects[i].pos = vector(x[i], y[i], z[i])

        t += dt

        print(f"Time: {t:.1f}s")
        for i in m:
            print(f"Object {i}: x = {x[i]:.3f}, y = {y[i]:.3f}, z = {z[i]:.3f}")
        print()
