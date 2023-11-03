import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams["font.family"] = "Times New Roman"

df1 = pd.read_csv("ConvergenceData.csv")
# get raw data from csv files

raw_data = df1[["iter", "K", "R", "alpha", "mesh"]].to_numpy()

fig1, (ax1, ax2) = plt.subplots(1, 2, tight_layout=True)

fig1.set_dpi(100)
fig1.set_size_inches(6, 3)

ax1.plot(raw_data[:, 0], raw_data[:,1], 'b*')
ax1.set_xlabel('Iterations', fontsize=14)
ax1.set_ylabel(r'k^{\*}[/]', fontsize=14)

ax2.plot(raw_data[:, 0], raw_data[:,2], 'k--')
ax2.set_xlabel('Iterations', fontsize=14)
ax2.set_ylabel('Residual', fontsize=14)

plt.show()