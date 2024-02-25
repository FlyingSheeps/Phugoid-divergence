import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy import interpolate
import sympy as sp
import pandas as pd

#import mass and elasitcity
df = pd.read_csv("wing.csv")
print(df)

# 値の取得
xold = df['span'].values/1000
GIp = df['torsion'].values
mass = df['mass'].values
f_GIp = interpolate.interp1d(xold, GIp, kind='linear')
f_mass = interpolate.interp1d(xold, mass, kind='linear')


#subdivision
N = 150
L = max(xold)
x_list = np.linspace(0,L,N)
dx = x_list[1]-x_list[0]

#subdivisionに従って値を補間
GIp = f_GIp(x_list)
mass = f_mass(x_list)

# 形状関数は1次バー
x = sp.Symbol('x')
Ni = (1 - x / dx)
Nj = x / dx

NN = np.array([[quad(sp.lambdify(x, Ni * Ni), 0, dx)[0], quad(sp.lambdify(x, Ni * Nj), 0, dx)[0]],
               [quad(sp.lambdify(x, Nj * Ni), 0, dx)[0], quad(sp.lambdify(x, Nj * Nj), 0, dx)[0]]])
NN = np.array(NN, dtype=float)
print(NN)

NxNx = np.array([[quad(sp.lambdify(x, sp.diff(Ni) * sp.diff(Ni)), 0, dx)[0], quad(sp.lambdify(x, sp.diff(Ni) * sp.diff(Nj)), 0, dx)[0]],
                 [quad(sp.lambdify(x, sp.diff(Nj) * sp.diff(Ni)), 0, dx)[0], quad(sp.lambdify(x, sp.diff(Nj) * sp.diff(Nj)), 0, dx)[0]]])
NxNx = np.array(NxNx, dtype=float)
print(NxNx)

K = np.zeros((N, N))
M = np.zeros((N, N))

#行列の組み立て
for i in range(N-1):
    K[i:i+2, i:i+2] += GIp[i] * NxNx #弾性効果
    M[i:i+2, i:i+2] += mass[i] * NN

K = K[1:,1:] #固定端条件
M = M[1:,1:]

D,V = eigh(a=-K,b=M)
D = D[::-1]
V = V[:,::-1]
theta_mode = np.vstack((np.zeros((1,N-1)),V))


np.save("theta_mode.npy",theta_mode[:,0:5])
np.save("theta_x_list.npy",x_list)

print(np.sqrt(abs(D[0:5])))
plt.figure(1)
plt.plot(x_list,theta_mode[:,0:5])
plt.grid()
plt.savefig("theta_mode.png")
plt.show()
