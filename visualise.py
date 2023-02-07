import matplotlib.pyplot as plt
import numpy as np
import os


output_dir = "cmake-build-debug/output/"

files = os.listdir(output_dir)


data = np.loadtxt(output_dir+f'{len(files)-1}.txt')
print(data)


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

for point in data:
    ax.scatter(point[0], point[1], point[2], color='red' if point[3] == 1 else 'blue', s=50)

plt.show()
