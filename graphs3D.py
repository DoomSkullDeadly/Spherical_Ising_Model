import matplotlib.pyplot as plt
import numpy as np
import os

data = np.loadtxt("cmake-build-debug/output/MvTvB.txt", skiprows=1).T

fig = plt.figure(figsize=plt.figaspect(0.5))

ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax1.scatter(data[0], data[1], data[2])
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('B Field (T)')
ax1.set_zlabel('Average Magnetisation per Spin')

ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.scatter(data[0], data[1], data[3])
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('B Field (T)')
ax2.set_zlabel('Average Energy per Spin (J)')

plt.show()
