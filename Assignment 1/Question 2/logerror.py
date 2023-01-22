import numpy as np
import matplotlib.pyplot as plt


x, y = np.loadtxt('norm_2b.txt', unpack=True)
plt.figure()
plt.plot(np.log(x),np.log(y), '-.', label = 'L2 Norm')
# plt.xscale("log")
# plt.yscale("log")

plt.axline((6.2,-9.25), slope=-2,c='orange', label='Slope = -2')
plt.xlabel(' Log nx')
plt.ylabel('Log L2 Norm')
plt.legend()
# plt.ylim(-20,50)
plt.title('Log L2 error norm at final time')
plt.savefig('logerror_2b.png')
