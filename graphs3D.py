import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("cmake-build-debug/output/MvTvB.txt", skiprows=1).T

fig = plt.figure(figsize=plt.figaspect(0.5))

ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.plot_trisurf(data[0], data[1], data[2], cmap='coolwarm')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('B Field (T)')
ax1.set_zlabel('Average Magnetisation per Spin')

ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.plot_trisurf(data[0], data[1], data[3], cmap='coolwarm')
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('B Field (T)')
ax2.set_zlabel('Average Energy per Spin (J)')

ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax3.plot_trisurf(data[0], data[1], data[4], cmap='coolwarm')
ax3.set_xlabel('Temperature (K)')
ax3.set_ylabel('B Field (T)')
ax3.set_zlabel('Magnetic Susceptibility per Spin')

ax4 = fig.add_subplot(2, 2, 4, projection='3d')
ax4.plot_trisurf(data[0], data[1], data[5], cmap='coolwarm')
ax4.set_xlabel('Temperature (K)')
ax4.set_ylabel('B Field (T)')
ax4.set_zlabel('Heat Capacity per Spin')

plt.show()
