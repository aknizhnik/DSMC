import matplotlib.pyplot as plt
from math import *

def read_distr(filename):
    vx_vals = []
    z_vals = []
    with open(filename, "r") as fdistr:
        for line in fdistr:
            vals = line.strip().split()
            x, y, z = [int(vals[0]), int(vals[1]), int(vals[2])]
            z_vals.append(z)
            vx_vals.append(abs(float(vals[4])))
    return vx_vals


D = 0.9
L = 1.0
r = 0.5
kB = 1.38e-16
AUM = 1.6e-24
sigma = 1.0e-15
molmass = 32
Pgas = 1e6 / 760
grad_P = 0.01*Pgas/L
Ngas = Pgas / (kB*300)
rho = Ngas * molmass * AUM
lambda0 = 1.0/(Ngas * sigma)
V_T = sqrt(8*kB*300/(pi*AUM*molmass))
nu = lambda0 * V_T / 3
eta = nu * rho
ksi = (2.0 - r)/r * lambda0

def u_anal(y):
    u = 1/(8*eta) * (D**2 - 4*y**2 + 4 * D * ksi) * grad_P
    return u

vx_vals = read_distr("distr_30.dat")
dy = D/len(vx_vals)
y_vals = [(i+0.5)*dy for i in range(len(vx_vals))]
vx_avg = sum(vx_vals)/len(vx_vals)

vx_anal_vals = [u_anal(y - D/2) for y in y_vals]

U_anal_avg = D**2/(12*eta) * (1 + 6 * ksi/D) * grad_P

print("Ngas =", Ngas)
print("Lambda =", lambda0)
print("V_T =", V_T)
print("nu =", nu)
print("eta =", eta)
print("U anal avg =",U_anal_avg)
print("U avg =",vx_avg)


fig, ax = plt.subplots()
ax.plot(vx_vals, y_vals, linewidth=3, label="DSMC")
ax.plot(vx_anal_vals, y_vals, linewidth=3, label="Analytic")
ax.set_xlabel("Velocity, cm/s")
ax.set_ylabel("Y coordinate, cm")
ax.set_xlim([0.0, 3000])
plt.legend()

plt.show()