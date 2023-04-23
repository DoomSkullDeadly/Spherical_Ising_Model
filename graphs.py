import matplotlib.pyplot as plt
import numpy as np
import os

dir_path = os.getcwd()
files = os.listdir(dir_path+"//cmake-build-debug//output//")

good_files = [f for f in files if 'MvB' in f or 'MvT' in f]
if len(good_files) == 0:
    print("No files!")
    exit()

fig, ax = plt.subplots(1, 2, figsize=(12, 4))

for f in good_files:
    data = np.loadtxt(dir_path+"//cmake-build-debug//output//"+f, skiprows=1).T
    ax[0].scatter(data[0], data[1], marker='.', s=10, label=f'{f[f.find("-")+1:f.find(".txt")]}')
    ax[1].scatter(data[0], data[2], marker='.', s=10, label=f'{f[f.find("-")+1:f.find(".txt")]}')


variable = 'Temperature' if good_files[0][2] == 'T' else 'B Field'

ax[0].set_title(f"Average Magnetisation vs {variable}")
ax[0].set_xlabel(f"{variable} ({'K' if variable == 'Temperature' else 'T'})")
ax[0].set_ylabel("Average Magnetisation per Spin")
ax[0].legend()

ax[1].set_title(f"Average Energy vs {variable}")
ax[1].set_xlabel(f"{variable} ({'K' if variable == 'Temperature' else 'T'})")
ax[1].set_ylabel("Average Energy per Spin (J)")
ax[1].legend()

plt.savefig("graphs.png")
plt.show()
