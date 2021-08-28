import matplotlib.pyplot as plt

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

pres, vel1, vel2 = read_res("res_2D.dat")

plt.plot(pres, vel1, label="r=1")
plt.plot(pres, vel2, label="r=0.5")
plt.xscale('log')
plt.xlim((1e-1, 1e1))
plt.xlabel("Pressure, Torr")
plt.ylabel("Mean velocity, cm/s")
plt.legend()
plt.show()
