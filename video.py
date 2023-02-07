import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

data_dir = "cmake-build-debug/output/"

files = os.listdir(data_dir)

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')


def animate(i):
    print(i)  # because it's slow af and i wanna know where it is
    ax.clear()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    data = np.loadtxt(data_dir+str(i)+'.txt')
    for point in data:
        ax.scatter(point[0], point[1], point[2], color='red' if point[3] == 1 else 'blue', s=50)
    return ax,


ani = animation.FuncAnimation(fig, animate, frames=len(files), interval=33.333333, blit=False)
writervideo = animation.FFMpegWriter(fps=30)
ani.save("lattice evolution.mp4", writer=writervideo)
plt.close()
