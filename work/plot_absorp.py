import matplotlib.pyplot as plt
from math import *

D = 0.1
L = 20.0
r = 0.5
Temp = 300
kB = 1.38e-16
AUM = 1.6e-24
sigma = 2.0e-15
sigma_ph = 9e-18
frac_O3 = 0.06
I1 = 5e-2 / 4.6 / 1.6e-19
I2 = 5e-1 / 4.6 / 1.6e-19
molmass = 32
V_T = sqrt(8*kB*300/(pi*AUM*molmass))

def get_den(pres):
    den = pres*1e6/760/(kB*Temp)
    return den

def get_inten(pres, x):
    den = get_den(pres)
    l_optic = sigma_ph * den * frac_O3 * x
    inten = exp(-l_optic)
    return inten

Np = 50
xvals = [L/(Np-1)*i for i in range(Np)]
pvals = [0.5, 1.0, 2.0, 4.0]

for p in pvals:
    ints = [get_inten(p,x) for x in xvals]
    label = F"P = {p} Torr"
    plt.plot(xvals, ints, label=label, linewidth=3)

plt.xlabel("Distance along channel, cm")
plt.ylabel("Intensity, a.u.")
plt.legend()
plt.show()
