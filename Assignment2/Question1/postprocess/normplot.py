import numpy as np
import matplotlib.pyplot as plt

x,y,err = np.loadtxt('l2norm.dat', unpack=True)

plt.figure()
plt.title('L2 norm vs grid points')

plt.plot(x,err*10e4, '-o')

plt.xlabel('nx')
plt.ylabel('ny')
plt.savefig('error.png')
