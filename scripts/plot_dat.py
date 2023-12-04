import numpy as np
import matplotlib.pylab as plt


N = 81
for i in range(N):
    data = np.loadtxt("../%06d_TM.dat"%i)

    fig, ax = plt.subplots(3, 1, figsize=(10, 9))
    plt.subplots_adjust(hspace=0.0)

    ax[0].plot(data[:, 0], data[:, 1])
    ax[0].set(ylabel=r"$\rho$")
    ax[0].set(xscale="log")
    ax[0].set(yscale="log")

    ax[1].plot(data[:, 0], data[:, 2] * 1.e-6)
    ax[1].set(ylabel=r"$U/c [10^6]$")
    ax[1].set(xscale="log")

    ax[2].plot(data[:, 0], data[:, 3])
    ax[2].set(ylabel=r"$P$")
    ax[2].set(yscale="log")
    ax[2].set(xscale="log")

    plt.suptitle("ID = %06d"%i)
    plt.savefig("%06d.png"%i)
    plt.close()
    #plt.show()
