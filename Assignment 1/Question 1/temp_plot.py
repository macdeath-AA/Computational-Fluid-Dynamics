import numpy as np
import matplotlib.pyplot as plt

plt.figure()
for i in range(0, 11, 3):
    x, y = np.loadtxt(f'T_x_{i*512}.txt', unpack=True)

    plt.plot(x, y, marker='.', label=f'time step={512*i}')
    plt.xlabel('x')
    plt.ylabel('T')
    plt.title('Temp profiles for NX = 1000')
    plt.legend()
    plt.savefig('nx1000.png')
