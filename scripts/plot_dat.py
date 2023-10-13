import numpy as np
import matplotlib.pylab as plt


data = np.loadtxt("000001_TM.dat")

fig, ax = plt.subplots(3, 1, figsize=(10, 5))
ax[0].plot(data[:, 0], data[:, 1])
ax[1].plot(data[:, 0], data[:, 2])
ax[2].plot(data[:, 0], data[:, 3])
plt.show()
