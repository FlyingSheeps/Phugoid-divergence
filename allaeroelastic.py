import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.linalg import eig
import sympy as sp

L = 15
N = 101  # 奇数であること
x_list = np.linspace(0, L, N)
dx = x_list[1] - x_list[0]

sigma = 12 / L
c = 1
U = 8.0
rho = 1.2
hr = 0.25
h = 0.4
rx = c * (hr - h)
b = c / 2
a = (hr - 0.5) * 2
J = 1 / 3 * sigma * (c / 2) ** 2 + rx ** 2 * sigma
EI = 85000
GIp = 16000

Mzz = np.zeros((N, N))
Mzt = np.zeros((N, N))
Mtz = np.zeros((N, N))
Mtt = np.zeros((N, N))

Czz = np.zeros((N, N))
Czt = np.zeros((N, N))
Ctz = np.zeros((N, N))
Ctt = np.zeros((N, N))

Kzz = np.zeros((N, N))
Kzt = np.zeros((N, N))
Ktz = np.zeros((N, N))
Ktt = np.zeros((N, N))

# 形状関数は二次バーなので全部同じの3*3
x = sp.Symbol('x')
Ni = (1 - x / dx) * (1 - x / (2 * dx))
Nj = 2 * x / dx * (1 - x / (2 * dx))
Nk = -x / (2 * dx) * (1 - x / dx)

NN = np.array([[quad(sp.lambdify(x, Ni * Ni), 0, 2 * dx)[0], quad(sp.lambdify(x, Ni * Nj), 0, 2 * dx)[0], quad(sp.lambdify(x, Ni * Nk), 0, 2 * dx)[0]],
               [quad(sp.lambdify(x, Nj * Ni), 0, 2 * dx)[0], quad(sp.lambdify(x, Nj * Nj), 0, 2 * dx)[0], quad(sp.lambdify(x, Nj * Nk), 0, 2 * dx)[0]],
               [quad(sp.lambdify(x, Nk * Ni), 0, 2 * dx)[0], quad(sp.lambdify(x, Nk * Nj), 0, 2 * dx)[0], quad(sp.lambdify(x, Nk * Nk), 0, 2 * dx)[0]]])
NN = np.array(NN, dtype=float)
print(NN)

NxNx = np.array([[quad(sp.lambdify(x, sp.diff(Ni) * sp.diff(Ni)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Ni) * sp.diff(Nj)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Ni) * sp.diff(Nk)), 0, 2 * dx)[0]],
                 [quad(sp.lambdify(x, sp.diff(Nj) * sp.diff(Ni)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nj) * sp.diff(Nj)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nj) * sp.diff(Nk)), 0, 2 * dx)[0]],
                 [quad(sp.lambdify(x, sp.diff(Nk) * sp.diff(Ni)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nk) * sp.diff(Nj)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nk) * sp.diff(Nk)), 0, 2 * dx)[0]]])
NxNx = np.array(NxNx, dtype=float)
print(NxNx)

NxxNxx = np.array([[quad(sp.lambdify(x, sp.diff(Ni, x, 2) * sp.diff(Ni, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Ni, x, 2) * sp.diff(Nj, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Ni, x, 2) * sp.diff(Nk, x, 2)), 0, 2 * dx)[0]],
                   [quad(sp.lambdify(x, sp.diff(Nj, x, 2) * sp.diff(Ni, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nj, x, 2) * sp.diff(Nj, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nj, x, 2) * sp.diff(Nk, x, 2)), 0, 2 * dx)[0]],
                   [quad(sp.lambdify(x, sp.diff(Nk, x, 2) * sp.diff(Ni, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nk, x, 2) * sp.diff(Nj, x, 2)), 0, 2 * dx)[0], quad(sp.lambdify(x, sp.diff(Nk, x, 2) * sp.diff(Nk, x, 2)), 0, 2 * dx)[0]]])
NxxNxx = np.array(NxxNxx, dtype=float)
print(NxxNxx)

print((sigma + np.pi * rho * b ** 2) * NN)
print(Mzz[0:3, 0:3])

for i in range(0, N - 2, 2):

    Mzz[i:i+3, i:i+3] += (sigma + np.pi * rho * b ** 2) * NN
    Mzt[i:i+3, i:i+3] += (-sigma * rx + np.pi * rho * b ** 3 * a) * NN
    Mtz[i:i+3, i:i+3] += (-sigma * rx + np.pi * rho * b ** 3 * a) * NN
    Mtt[i:i+3, i:i+3] += (J + np.pi * rho * b ** 4 * (1 / 8 + a ** 2)) * NN

    Czz[i:i+3, i:i+3] += 2 * np.pi * rho * U * b * NN
    Czt[i:i+3, i:i+3] += -(np.pi * rho * U * b ** 2 + 2 * np.pi * rho * U * b ** 2 * (0.5 - a)) * NN
    Ctz[i:i+3, i:i+3] += 2 * np.pi * rho * U * b ** 2 * (a + 0.5) * NN
    Ctt[i:i+3, i:i+3] += -2 * np.pi * rho * U * b ** 3 * (a + 0.5) * (0.5 - a) * NN

    Kzz[i:i+3, i:i+3] += EI * NxxNxx
    Kzt[i:i+3, i:i+3] += -2 * np.pi * rho * U ** 2 * b * NN
    Ktz[i:i+3, i:i+3] += 0 * NN
    Ktt[i:i+3, i:i+3] += (np.pi * rho * b ** 3 * U * (0.5 - a) - 2 * np.pi * rho * U ** 2 * b ** 2 * (a + 0.5)) * NN + GIp * NxNx

Mzz = Mzz[2:,2:]
Mzt = Mzt[2:,1:]
Mtz = Mtz[1:,2:]
Mtt = Mtt[1:,1:]

Kzz = Kzz[2:,2:]
Kzt = Kzt[2:,1:]
Ktz = Ktz[1:,2:]
Ktt = Ktt[1:,1:]

Czz = Czz[2:,2:]
Czt = Czt[2:,1:]
Ctz = Ctz[1:,2:]
Ctt = Ctt[1:,1:]

M = np.block([[Mzz, Mzt], [Mtz, Mtt]])
C = np.block([[Czz, Czt], [Ctz, Ctt]])
K = np.block([[Kzz, Kzt], [Ktz, Ktt]])

leftMatrix = np.block([[np.eye(2*N-3),np.zeros((2*N-3,2*N-3))],[np.zeros((2*N-3,2*N-3)),M]])
rightMatrix = np.block([[np.zeros((2*N-3,2*N-3)),np.eye(2*N-3)],[-K,-C]])


D,V = eig(a=rightMatrix,b=leftMatrix)


plt.figure(1)
plt.plot(np.real(D[0:5]),'o')
plt.grid()
plt.show()