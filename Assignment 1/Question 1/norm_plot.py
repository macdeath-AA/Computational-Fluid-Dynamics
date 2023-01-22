import numpy as np
import matplotlib.pyplot as plt

x1, y1 = np.loadtxt('von1b.txt', unpack=True)
# x2, y2 = np.loadtxt('norm1c.txt', unpack=True)
plt.figure()
plt.plot(x1, y1, '-.', label = 'b')
# plt.plot(x2, y2, '-.', label='a')
plt.xlabel('nx')
plt.ylabel('L2 Norm')
plt.ylim(-0.25,1.25)
plt.legend()
plt.title('L2 error norm at final time')
plt.savefig('von1b.png')
