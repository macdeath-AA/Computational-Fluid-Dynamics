import numpy as np
import matplotlib.pyplot as plt

x, y, data = np.loadtxt('T_xy_numer_060_020_0001.dat', unpack=True)
x1, y1, data1 = np.loadtxt('T_xy_exact_060_020_0001.dat', unpack=True)

x= np.reshape(x,(60,20))
y = np.reshape(y, (60, 20))
data = np.reshape(data, (60, 20))

x1 = np.reshape(x1, (60, 20))
y1 = np.reshape(y1, (60, 20))
data1 = np.reshape(data1, (60, 20))

fig, ax = plt.subplots(2,1)
fig.suptitle('Numerical vs exact solution for 60x20')

ax[0].contourf(x, y, data)
ax[1].contourf(x1, y1, data1)
ax[0].set_xlabel('Numerical')
ax[1].set_xlabel('Exact')

fig.tight_layout()
fig.savefig('60x20.png', dpi=300)
