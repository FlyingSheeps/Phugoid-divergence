import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#import mass and elasitcity
df = pd.read_csv("wing.csv")
print(df)

# 値の取得とミラーの作成
xold = df['span'].values/1000
GIp = df['GIp'].values
he = df['T.C.'].values
c = df['c'].values/1000
CL = df['CL'].values
Cm = df['Cm'].values
U0 = df['U0'].values[0]
#原点が被らないようにミラー処理する
if xold[0] == 0:
    xold = np.concatenate((-xold[::-1], xold[1:]))
    GIp = np.concatenate((GIp[::-1], GIp[1:]))
    c = np.concatenate((c[::-1], c[1:]))
    CL = np.concatenate((CL[::-1], CL[1:]))
    Cm = np.concatenate((Cm[::-1], Cm[1:]))
    f_GIp = interpolate.interp1d(xold, GIp, kind='linear')
    f_he = interpolate.interp1d(xold, he, kind='linear')
    f_c = interpolate.interp1d(xold, c, kind='linear')
    f_CL = interpolate.interp1d(xold, CL, kind='linear')
    f_Cm = interpolate.interp1d(xold, Cm, kind='linear')
else:
    xold = np.concatenate((-xold[::-1], xold[:]))
    GIp = np.concatenate((GIp[::-1], GIp[:]))
    c = np.concatenate((c[::-1], c[:]))
    CL = np.concatenate((CL[::-1], CL[:]))
    Cm = np.concatenate((Cm[::-1], Cm[:]))
    f_GIp = interpolate.interp1d(xold, GIp, kind='linear')
    f_he = interpolate.interp1d(xold, he, kind='linear')
    f_c = interpolate.interp1d(xold, c, kind='linear')
    f_CL = interpolate.interp1d(xold, CL, kind='linear')
    f_Cm = interpolate.interp1d(xold, Cm, kind='linear')


#subdivision
N = 150
L = max(xold)
x = np.linspace(-L,L,N)
dx = x[1]-x[0]

#subdivisionに従って値を補間
GIp = f_GIp(x)
he = f_he(x)
c = f_c(x)
Cm = f_Cm(x)
CL = f_CL(x)
CD = 0.05

aw = 2.0*np.pi
at = 2.0*np.pi

#全体のパラメータ
mass = 100
grav = 9.8
rho = 1.2
hnw = 0.25

S = sum(c)*dx
b = 2*L
cbar = S/b
St = 0.25
Vh = 0.47

h = 0.42

#トリム探索の分割数
trimDivision = 100
#書き出したいデータ列or変数列or固定する値
u_list = np.linspace(1,15,trimDivision) #速度を振る
gamma_list = np.ones(trimDivision)*0 #水平飛行
theta_list = np.zeros(trimDivision)
alpha_list = np.zeros(trimDivision)
balance_list = np.zeros(trimDivision)
thrust_list = np.zeros(trimDivision)

def TR797(x):
    u = x[0]
    alpha = x[1]
    theta = x[2]
    z = x[3:3+N]
    vartheta = x[3+N:3+2*N]



    return Lift, Drag, Moment

for i in range(trimDivision):

    # Define the objective function
    def objective_function(x):
        u = u_list[i]
        gamma = x[1]
        theta = x[2]
        alpha = theta - gamma
        trust = x[4]
        z = x[4:4+N]
        vartheta = x[4+N:4+2*N]




        X = -Drag -mass*grav*np.sin(gamma) + thrust
        Z = Lift - mass*grav*np.cos(gamma)

        return X**2 + Z**2 + Moment**2

    # Initial guess
    initial_guess = np.array([5, 0, 0])

    # Bounds for variables (optional)
    bounds = ((0, 20), (-np.pi/2, np.pi/2), (-np.pi/2, np.pi/2))  # bounds for each variable (x, y)

    # Minimize the objective function
    result = minimize(objective_function, initial_guess, bounds=bounds)

    # Print results
    lt = 0.27 + (1-h)*cbar
    Vh = St*lt/S/cbar
    
    u_list[i] = result.x[0]
    gamma_list[i] = result.x[1]*180/np.pi
    theta_list[i] = result.x[2]*180/np.pi
    alpha_list[i] = (result.x[2]-result.x[1])*180/np.pi
    LP_list[i] = -1/(result.x[0]*np.tan(result.x[1]))
    margin_list[i] = (h-hnw)-Vh
    balance_list[i] = result.fun

plt.figure(0)
plt.plot(h_list,alpha_list)
plt.plot(h_list,theta_list)
plt.plot(h_list,gamma_list)
plt.plot(h_list,u_list)
plt.legend(["angle of attack","pitch angle","path angle","airspeed"])
plt.xlabel("C.G. position (MAC)")
plt.ylabel("angle (deg) / airspeed (m/s)")
plt.grid()

plt.figure(1)
plt.semilogy(h_list,balance_list)
plt.grid()

plt.show()

