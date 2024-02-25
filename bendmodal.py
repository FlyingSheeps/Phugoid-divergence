import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate



#import mass and elasitcity
df = pd.read_csv("spar.csv")
print(df)

# ヘッダーから列の名前を取得
columns = df.columns

# x座標とy座標のデータを取得
xold = df[columns[0]].values/1000
mass = df[columns[1]].values
EI = df[columns[2]].values

#subdivision
N = 150
L = max(xold)
x = np.linspace(0,L,N)
dx = x[1]-x[0]

# 線形補間関数を作成
f_mass = interpolate.interp1d(xold, mass, kind='linear')
f_EI = interpolate.interp1d(xold, EI, kind='linear')

mass = f_mass(x)
mass = mass*13/sum(mass)
EI = f_EI(x)

M = np.zeros((2*N,2*N))
K = np.zeros((2*N,2*N))

for i in range(N-1):
    n = 2*i

    M[n:n+4,n:n+4] += np.array([[156,    22*dx,    54,    -13*dx    ],
                                [22*dx,  4*dx**2,  13*dx,  -3**dx**2],
                                [54,     13*dx,    156,    -22*dx   ],
                                [-13*dx, -3*dx**2, -22*dx, 4*dx**2  ]])*mass[i]/420
    K[n:n+4,n:n+4] += np.array([[12,   6*dx,    -12,   6*dx   ],
                                [6*dx, 4*dx**2, -6*dx, 2*dx**2],
                                [-12,  -6*dx,   12,    -6*dx  ],
                                [6*dx, 2*dx**2, -6*dx, 4*dx**2]])*EI[i]/dx**3

M = M[2:,2:]
K = K[2:,2:]

D,V = eigh(a=K, b=M)
print(V.shape)
V = np.vstack((np.zeros((2,2*N-2)),V))

#print(M)
#print(K)
omega = abs(np.sqrt(D))
print(omega[0:5])
print(3.516/L**2*np.sqrt(np.mean(EI)/np.mean(mass)*dx))

count = 5

u_mode = V[0::2,0:count]
ux_mode = V[1::2,0:count]
uxx_mode = np.gradient(ux_mode,dx,axis=0)
np.save("u_mode.npy",u_mode)
np.save("ux_mode.npy",ux_mode)
np.save("uxx_mode.npy",uxx_mode)

plt.figure(1)
plt.plot(omega[:count])
plt.hlines(3.516/L**2*np.sqrt(np.mean(EI)/np.mean(mass)*dx), 0, count)
plt.grid()

plt.figure(2)
plt.plot(x,u_mode)
plt.grid()

plt.figure(3)
plt.plot(x,ux_mode)
plt.grid()

plt.figure(4)
plt.plot(x,uxx_mode)
plt.grid()


plt.show()


