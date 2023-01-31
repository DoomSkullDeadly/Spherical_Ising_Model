import matplotlib.pyplot as plt
import numpy as np


output_dir = "cmake-build-debug/output/"

data = np.loadtxt(output_dir+'50.txt')
print(data)


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

for point in data:
    ax.scatter(point[0], point[1], point[2], color='red' if point[3] == 1 else 'blue', s=50)

plt.show()
