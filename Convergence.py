import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams["font.family"] = "Times New Roman"

df1 = pd.read_csv("ConvergenceData.csv")
# get raw data from csv files

raw_data = df1[["iter", "K", "R", "alpha", "mesh"]].to_numpy()

# Get final permeability

perm = f"Perm={raw_data[-1,1]:09f}"

fig1, (ax1, ax2) = plt.subplots(1, 2, tight_layout=True)

fig1.set_dpi(100)
fig1.set_size_inches(6, 3)

ax1.plot(raw_data[:, 0], raw_data[:,1], 'b*')
ax1.set_title('File 00001.jpg', fontsize=16)
ax1.set_xlabel('Iterations', fontsize=16)
ax1.set_ylabel(r'k$^{*}[/]$', fontsize=16)

ax2.plot(raw_data[:, 0], raw_data[:,2], 'k--')
ax2.set_xlabel('Iterations', fontsize=16)
ax2.set_ylabel('Residual', fontsize=16)
ax2.set_title(perm, fontsize=16)

plt.show()