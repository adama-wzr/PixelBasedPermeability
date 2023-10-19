import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# simulation properties

grid_size = 128*2

# Read image
imgName = "00010.jpg"
img = np.uint8(mpimg.imread(imgName))

# df1 = pd.read_csv("ExpUV_mod.csv")
# df1 = pd.read_csv("UV.csv")
df1 = pd.read_csv("Circles00010_sim.csv")
# get raw data from csv files

raw_data = df1[["P", "U", "V", "x", "y"]].to_numpy()

P = raw_data[:,0]
P = np.reshape(P, [grid_size,grid_size])

U = raw_data[:,1]
U = np.reshape(U, [grid_size,grid_size])

V = raw_data[:,2]
V = np.reshape(V, [grid_size,grid_size])

# Post-Process P, add a mask where P = 0 (solid)

(m,n) = P.shape

mask = np.zeros_like(P, dtype=bool)

for i in range(m):
	for j in range(n):
		if P[i][j] < 1:
			mask[i][j] = True		



P = np.ma.array(P, mask=mask)
U = np.ma.array(U, mask=mask)
V = np.ma.array(V, mask=mask)

# Create the mesh grid

Xp, Yp = np.meshgrid(np.linspace(0, 0.1, grid_size), np.linspace(0.1, 0, grid_size))

# plotting

fig1, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, constrained_layout=True)

fig1.set_dpi(100)
fig1.set_size_inches(15, 7)

# First axis is just the image

ax1.imshow(img)
ax1.set_title(imgName)

# Second axis is U-velocity contour

CS2 = ax2.contourf(Xp, Yp, U, 40, cmap=plt.cm.viridis)
cbar2 = fig1.colorbar(CS2, ax=ax2)
ax2.set_title("U-velocity Contour")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

# Third axis is the V-velocity contour

CS3 = ax3.contourf(Xp, Yp, V, 40, cmap=plt.cm.viridis)
cbar3 = fig1.colorbar(CS3, ax=ax3)
ax3.set_title("V-velocity Contour")
ax3.set_xlabel("x")
ax3.set_ylabel("y")

# Fourth axis is Pressure contour

CS4 = ax4.contourf(Xp, Yp, P, 40, cmap=plt.cm.viridis)
cbar = fig1.colorbar(CS4, ax=ax4)
ax4.set_title("Pressure Contour")
ax4.set_xlabel("x")
ax4.set_ylabel("y")

# Fifth axis is a vector field

velocityMag = np.sqrt(U**2 + V**2)

step = 12

CS5 = ax5.pcolor(Xp, Yp, velocityMag, cmap='rainbow')
ax5.quiver(Xp[::step, ::step], Yp[::step, ::step], U[::step, ::step], V[::step, ::step])
cbar = fig1.colorbar(CS5, ax=ax5)
ax5.set_title("Velocity Vector Field")

# Sixth Axis is a contour of the velocity magnitudes

CP6 = ax6.contour(Xp, Yp, velocityMag, 10)
ax6.clabel(CP6, inline=True, fontsize=12)
CS6 = ax6.pcolor(Xp, Yp, velocityMag, cmap='rainbow')
cbar = fig1.colorbar(CS6, ax=ax6)
ax6.set_title("Velocity Magnitude Contour")

plt.show()
