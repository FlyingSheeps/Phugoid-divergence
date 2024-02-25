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
he = df['T.C.'].values
c = df['chord'].values/1000
f_GIp = interpolate.interp1d(xold, GIp, kind='linear')
f_he = interpolate.interp1d(xold, he, kind='linear')
f_c = interpolate.interp1d(xold, c, kind='linear')

#subdivision
N = 150
L = max(xold)
x = np.linspace(0,L,N)
dx = x[1]-x[0]

#subdivisionに従って値を補間
GIp = f_GIp(x)
he = f_he(x)
c = f_c(x)

# 形状関数は二次バーなので全部同じの3*3
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

#ダイバージェンス速度の探索範囲
Utry = 20
Udiv = 100
U = np.linspace(0,Utry,Udiv)
rho = 1.2
diag_list = np.zeros(Udiv)

for n in range(Udiv):
    print(U[n])
    Kelastic = np.zeros((N, N))
    Kaero = np.zeros((N, N))
    #行列の組み立て
    for i in range(N-1):
        Kelastic[i:i+2, i:i+2] += -GIp[i] * NxNx #弾性効果
        Kaero[i:i+2, i:i+2] += 0.5*rho*U[n]**2*c[i]**2*(he[i]-0.25)*2*np.pi * NN #空力効果

    K = Kelastic + Kaero
    K = K[1:,1:] #固定端条件

    D,V = eigh(K)
    diag_list[n] = max(D)

np.save("eigen-divergence.npy",diag_list)
np.save("U-divergence.npy",U)

plt.figure(1)
plt.plot(U,diag_list,'o')
plt.title('Eigenvalues of divergence')
plt.xlabel('Air speed (m/s)')
plt.ylabel('Maximum real part of eigenvalues')
plt.grid()
plt.savefig('Eigenvalues_of_divergence.png')
plt.show()