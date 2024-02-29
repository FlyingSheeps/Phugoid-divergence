import numpy as np
import matplotlib.pyplot as plt

diag0 = np.load("eigen-phugpid-divergence0.npy")
diag1 = np.load("eigen-phugpid-divergence1.npy")
diag2 = np.load("eigen-phugpid-divergence2.npy")
diag3 = np.load("eigen-phugpid-divergence3.npy")

U = np.load("U-phugpid-divergence.npy")
plt.figure(1)
plt.plot(U,diag0,'o')
plt.plot(U,diag1,'o')
plt.plot(U,diag2,'o')
plt.plot(U,diag3,'o')
plt.legend(["$C_m=0.0$","$C_m=-0.1$","$C_m=-0.2$","$C_m=-0.3$"])
plt.hlines(0,0,20,linestyles='--')
plt.grid()
plt.savefig("eigenplot.pdf")
plt.show()


