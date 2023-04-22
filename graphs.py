import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("results.txt").T

fig, ax = plt.subplots(1, 2, figsize=(12, 4))

ax[0].scatter(data[0], data[1], marker='.')
ax[0].set_title("Average Magnetisation vs Temperature")
ax[0].set_xlabel("Temperature (K)")
ax[0].set_ylabel("Average Magnetisation per Spin")

ax[1].scatter(data[0], data[2], marker='.')
ax[1].set_title("Average Energy vs Temperature")
ax[1].set_xlabel("Temperature (K)")
ax[1].set_ylabel("Average Energy per Spin (J)")

plt.show()
