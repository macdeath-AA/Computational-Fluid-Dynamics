import numpy as np
import matplotlib.pyplot as plt

nx, ny = 20,20
x, y, data = np.loadtxt('T_xy_numer_020_020_0000.dat', unpack=True)
# x1, y1, data1 = np.loadtxt('T_xy_exact_020_020_0000.dat', unpack=True)

x = np.reshape(x, (nx, ny))
y = np.reshape(y, (nx, ny))
data = np.reshape(data, (nx, ny))

# x1 = np.reshape(x1, (nx, ny))
# y1 = np.reshape(y1, (nx, ny))
# data1 = np.reshape(data1, (nx, ny))

fig, ax = plt.subplots(1, 1)
fig.suptitle(f'Central Difference Solution PE = 4')

ax.contourf(x, y, data)
# ax[1].contourf(x1, y1, data1)
ax.set_xlabel('x')
ax.set_ylabel('y')
# ax[1].set_xlabel('Exact')

fig.tight_layout()
fig.savefig('b_cd.png',dpi=300)
