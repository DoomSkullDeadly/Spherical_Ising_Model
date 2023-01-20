import matplotlib.pyplot as plt
import numpy as np


output_dir = "cmake-build-debug/output/"

data = np.loadtxt(output_dir+'0.txt')
print(data)


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

