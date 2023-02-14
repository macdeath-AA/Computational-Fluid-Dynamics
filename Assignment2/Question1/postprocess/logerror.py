import numpy as np
import matplotlib.pyplot as plt


x, y, err = np.loadtxt('l2norm.dat', unpack=True)
plt.figure()
plt.plot(np.log(x), np.log(err), '-.', label='L2 Norm')

plt.axline((4, -9), slope=-2, c='orange', label='Slope = -2')
plt.xlabel(' Log x')
plt.ylabel('Log L2 Norm')
plt.legend()
# plt.ylim(-20,50)
plt.title('L2 norm order of accuracy')
plt.savefig('logerror_1e.png')
