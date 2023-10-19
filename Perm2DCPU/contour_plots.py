import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# simulation properties

grid_size = 128*2

# df1 = pd.read_csv("ExpUV_mod.csv")
# df1 = pd.read_csv("UV.csv")
df1 = pd.read_csv("test.csv")
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

fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, constrained_layout=True)

fig1.set_dpi(100)
fig1.set_size_inches(15, 4)

CS1 = ax1.contourf(Xp, Yp, P, 20, cmap=plt.cm.viridis)
cbar = fig1.colorbar(CS1, ax=ax1)
ax1.set_title("Pressure Contour")
ax1.set_xlabel("x")
ax1.set_ylabel("y")


CS2 = ax2.contourf(Xp, Yp, U, 20, cmap=plt.cm.viridis)
cbar2 = fig1.colorbar(CS2, ax=ax2)
ax2.set_title("U-velocity Contour")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

CS3 = ax3.contourf(Xp, Yp, V, 20, cmap=plt.cm.viridis)
cbar3 = fig1.colorbar(CS3, ax=ax3)
ax3.set_title("V-velocity Contour")
ax3.set_xlabel("x")
ax3.set_ylabel("y")


plt.show()
