import matplotlib.pyplot as plt
from math import *

def read_res(filename):
    pres = []
    vel1 = []
    vel2 = []
    with open(filename,"r") as fres:
        for line in fres:
            vals = line.strip().split()
            pres.append(float(vals[0]))
            vel1.append(float(vals[1]))
            vel2.append(float(vals[2]))
    return pres, vel1, vel2

D = 0.1
L = 20.0
r = 0.5
kB = 1.38e-16
AUM = 1.6e-24
sigma = 2.0e-15
sigma_ph = 9e-18
I1 = 5e-2 / 4.6 / 1.6e-19
I2 = 5e-1 / 4.6 / 1.6e-19
molmass = 32
V_T = sqrt(8*kB*300/(pi*AUM*molmass))
Diff0 = 200.0   # cm2/s at 1 Torr

pres, vel1, vel2 = read_res("res_2D.dat")

def diff(pres):
    diff_g = Diff0/pres
    diff_kn = D*V_T/3
    diff_tot = diff_g * diff_kn / (diff_kn + diff_g)
    return diff_tot

def tau_diff(pres):
    diff_coef = diff(pres)
    tau_d = (L/2)**2/(2*diff_coef)
    return tau_d


tau_g = [L/u for u in vel1]
tau_d = [tau_diff(p) for p in pres]
tau_p1 = [1/(sigma_ph * I1) for p in pres]
tau_p2 = [1/(sigma_ph * I2) for p in pres]

plt.plot(pres, tau_g, label=r"$\tau_g$")
plt.plot(pres, tau_d, label=r"$\tau_d$")
plt.plot(pres, tau_p1, label=r"$\tau_p$, I=50 mW/cm2")
plt.plot(pres, tau_p2, label=r"$\tau_p$, I=500 mW/cm2")
plt.xscale('log')
plt.xlim((1e-1, 1e1))
plt.xlabel("Pressure, Torr")
plt.ylabel("Characteristic time, s")
plt.legend()
plt.show()