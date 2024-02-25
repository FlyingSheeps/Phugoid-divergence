import numpy as np
import matplotlib.pyplot as plt

diag1 = np.load("eigen-divergence.npy")
diag2 = np.load("eigen-phugpid-divergence-0.13.npy")
diag3 = np.load("eigen-phugpid-divergence-0.2.npy")

U = np.linspace(0,20,100)
plt.figure(1)
plt.plot(U,diag1*100/6,'o')
plt.plot(U,diag2,'o')
plt.plot(U,diag3,'o')
plt.legend(["divergence","phugoid-divergence1","phugoid-divergence2"])
plt.ylim([-200,200])
plt.hlines(0,0,20)
plt.grid()
plt.savefig("eigenplot.png")
plt.show()


