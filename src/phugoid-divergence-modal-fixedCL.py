import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig, det
from scipy import interpolate
import sympy as sp
import pandas as pd

#import mass and elasitcity
df = pd.read_csv("wing.csv")
print(df)

# 値の取得
xold = df['span'].values/1000
GIp = df['GIp'].values
he = df['T.C.'].values
c = df['c'].values/1000
CL = df['CL'].values
Cm = df['Cm'].values
U0 = df['U0'].values[0]
f_GIp = interpolate.interp1d(xold, GIp, kind='linear')
f_he = interpolate.interp1d(xold, he, kind='linear')
f_c = interpolate.interp1d(xold, c, kind='linear')
f_CL = interpolate.interp1d(xold, CL, kind='linear')
f_Cm = interpolate.interp1d(xold, Cm, kind='linear')

#モードの読み込み
temp_mode = np.load("theta_mode.npy")
xold_mode = np.load("theta_x_list.npy")
print(xold_mode)
theta_mode1 = temp_mode[:,0]
f_theta_mode1 = interpolate.interp1d(xold_mode, theta_mode1, kind='linear')
theta_mode2 = temp_mode[:,1]
f_theta_mode2 = interpolate.interp1d(xold_mode, theta_mode2, kind='linear')
theta_mode3 = temp_mode[:,2]
f_theta_mode3 = interpolate.interp1d(xold_mode, theta_mode3, kind='linear')


#subdivision
N = 150
L = max(xold)
x = np.linspace(0,L,N)
dx = x[1]-x[0]

#subdivisionに従って値を補間
GIp = f_GIp(x)
he = f_he(x)
c = f_c(x)
Cm = f_Cm(x)
CL0 = f_CL(x)
theta_mode1 = abs(f_theta_mode1(x))
thetax_mode1 =  np.gradient(theta_mode1,dx,axis=0)
theta_mode2 = abs(f_theta_mode2(x))
thetax_mode2 =  np.gradient(theta_mode2,dx,axis=0)
theta_mode3 = abs(f_theta_mode3(x))
thetax_mode3 =  np.gradient(theta_mode3,dx,axis=0)


#ダイバージェンス速度の探索範囲
Utry = 20
Udiv = 100
U = np.linspace(1,Utry,Udiv)
rho = 1.2
diag_list = np.zeros(Udiv)

for n in range(Udiv):
    print(U[n])
    Ttheta = np.zeros((3,3))
    Ztheta = np.zeros(3)
    Tu = np.zeros(3)
    #行列の組み立て
    Ttheta[0,0] = sum(0.5*rho*U[n]**2*c**2*(he-0.25)*2*np.pi * theta_mode1**2 - GIp * thetax_mode1**2) * dx 
    Ttheta[1,1] = sum(0.5*rho*U[n]**2*c**2*(he-0.25)*2*np.pi * theta_mode2**2 - GIp * thetax_mode2**2) * dx 
    Ttheta[2,2] = sum(0.5*rho*U[n]**2*c**2*(he-0.25)*2*np.pi * theta_mode3**2 - GIp * thetax_mode3**2) * dx
    Ztheta[0] = sum(-0.5*rho*U[n]**2*c*2*np.pi * theta_mode1) * dx
    Ztheta[1] = sum(-0.5*rho*U[n]**2*c*2*np.pi * theta_mode2) * dx
    Ztheta[2] = sum(-0.5*rho*U[n]**2*c*2*np.pi * theta_mode3) * dx
    Tu[0] = sum(rho*U[n]*c**2*(Cm+(he-0.25)*CL0)*theta_mode1) * dx
    Tu[1] = sum(rho*U[n]*c**2*(Cm+(he-0.25)*CL0)*theta_mode2) * dx
    Tu[2] = sum(rho*U[n]*c**2*(Cm+(he-0.25)*CL0)*theta_mode3) * dx
    Zu = sum(-rho*U[n]*c*CL0) * dx
    
    K = np.zeros((4,4))
    K[0,0] = 9.8/U[n]*Zu
    K[0,1:] = 9.8/U[n]*Ztheta
    K[1:,0] = Tu
    K[1:,1:] = Ttheta
    print(K)
    D,V = eig(K)
    diag_list[n] = max(D)

np.save("eigen-phugpid-divergence3.npy",diag_list)
np.save("U-phugpid-divergence.npy",U)

plt.figure(1)
plt.plot(U,diag_list,'o')
plt.title('Eigenvalues of phugoid-divergence')
plt.xlabel('Air speed (m/s)')
plt.ylabel('Maximum real part of eigenvalues')
plt.grid()
plt.savefig('Eigenvalues_of_phugoid-divergence-modal'+str(Cm[0])+'.pdf')
plt.show()