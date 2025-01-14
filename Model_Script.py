## In this section I am importing all the libraries I will need
import numpy as np
import matplotlib.pyplot as pl

## In this section I am setting the domain of solution and the discretised grid

# Dimensions of room in metres
L = 5
H = 5
Z = 5

#Discretised grid intervals
dx = 0.2
dy = 0.2
dz = 0.2
dt = 60  # Time in seconds
tend = 60*60*1  # Number of hours (has been set as 1 now) (multiply by 24 for number of days)

# Defining the location of the diffuser from the origin in metres
xlocation = 2
ylocation = 1.5
zlocation = 1

mass = 75 # Mass of diffusing liquid in grams
volume = 4*4*4 # Dimensions of liquid in diffuser (4cm sided cube)

# Grid (NOTE: All x, y and z axis have been defined separately despite being identical such that this code can also be used to model differently sized rooms)
x = np.arange(0, L+dx, dx) # Array along x direction
y = np.arange(0, H+dy, dy) # Array along y
zr = np.arange(0, Z+dz, dz)  # Array along z
t = np.arange(0, tend+dt, dt)

## In this section I am setting the initial values
# Creating matrix containing initial solutions (0s) for the Gauss Seidel iterations
phi = np.zeros((len(t), len(x), len(y), len(zr)))

# Initial Condition: Setting the density value at the location of the diffuser
phi[0, int(xlocation//dx)+1, int(ylocation//dy)+1, int(zlocation//dz)+1] = mass/volume

# Function of Gauss Seidel for this model
def GaussSeidel(phi, k): 
    tol = 0.001 # Tolerence of the iterations convergence
    D = 10**(-6) # Diffusivity
    r = dt/(dx**2)
    
    ## In this section I am setting the boundary conditions
    phi[k+1, :, :, 0] = phi[k+1, :, :, 1] # Derivativex = 0 at x = 0  
    phi[k+1,:, :, -1] = phi[k+1, :, :, -2] # Derivativex = 0 at x = L
    
    phi[k+1, :, 0, :] = phi[k+1, :, 1, :] # Derivativey = 0 at y = 0
    phi[k+1, :, -1, :] = phi[k+1,:, -2, :] # Derivativey = 0 at y = H

    phi[k+1, 0, : , :] = phi[k+1, 1, :, :] # Derivativez = 0 at z = 0
    phi[k+1, -1, : , :] = phi[k+1, -2, :, :] # Derivativez = 0 at z = Z
    
    A = np.copy(phi[k+1]) #Making a copy of the matrix, to check for tolerences 
    
    # Carrying out Gauss Seidel Iterations along each z, y and x directions
    for n in range(1, len(x)-1):
        for j in range(1, len(y)-1):
            for i in range(1, len(zr)-1): # First estimating values in the xy-planes (z) and then travelling along y and then along x
                # Crank-Nicholson Numerical Expression
                phi[k+1, n, j, i] = ((2/D - 6*r)*phi[k, n, j, i] + r*(phi[k+1, n, j, i-1] + phi[k, n, j, i+1] + phi[k+1, n+1, j, i] + phi[k+1, n-1, j, i] + phi[k+1, n, j-1, i] + phi[k+1, n, j+1, i] + phi[k, n, j, i-1] + phi[k, n, j, i+1] + phi[k, n, j-1, i] + phi[k, n, j+1, i] + phi[k, n+1, j, i] + phi[k, n-1, j, i]))/(2/D + 6*r)

    if (A-phi[k+1]).any() > tol : # Checking for tolerences: if tolerence not reached, Gauss Seidel iteration is repeated
        return GaussSeidel(phi, k)

    else: # if tolerences are satisfied, the array with the solutions is returned
        return phi   

## In this section I am implementing the numerical method
# Carrying Gauss Seidel to solve for the density throughout the room starting from time = 0 up until end time tend
for k in range(0, len(t)-1): 
    phi = GaussSeidel(phi, k)


## In this section I am showing the results
# 3D Plot in space x,y,z at the final time tend
Xg, Yg, Zg = np.meshgrid(x, y, zr)
fig = pl.figure()
ax = pl.axes(projection="3d")
img = ax.scatter3D(Xg, Yg, Zg, c=phi[-1], alpha=0.05, marker='.')
fig.colorbar(img)
pl.show()

# 2D plot in the z-plane where the diffuser is located
x2d, y2d = np.meshgrid(x, y)
pl.contourf(y2d, x2d, phi[-1, :,  :,int(zlocation//dz)+1])
pl.colorbar()
pl.show()
pl.contour(y2d, x2d, phi[-1, :, :,int(zlocation//dz)+1]) #2d contour without colourbar
pl.xlabel('X-axis') 
pl.ylabel('Y-axis')
pl.show()

# Graph of mass of liquid remaining at the diffuser container over time
# Calculation of mass at the location of diffuser
m = phi[:, int(xlocation//dx)+1, int(ylocation//dy)+1, int(zlocation//dz)+1]*volume
pl.plot(t/3600, m) #Time array is converted from seconds to hours
pl.ylabel('Mass of Fluid in Diffuser (g)')
pl.xlabel('Time (hours)')