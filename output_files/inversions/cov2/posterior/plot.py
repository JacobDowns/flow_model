import numpy as np
import matplotlib.pyplot as plt

Pxx = np.loadtxt('Pxx_new.txt')

for i in range(len(Pxx)):
    if i % 3 == 0:
        plt.plot(Pxx[i], marker = 'o')
    
plt.show()
