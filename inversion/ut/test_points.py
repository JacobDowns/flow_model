import numpy as np
from sigma_points import *
import matplotlib.pyplot as plt

N = 10
x = np.zeros(N)
Pxx = np.eye(N)

points = SigmaPoints(x, Pxx)
X, w, w = points.__get_mysovskikh_set__(1.)
print(w.sum())

quit()
print(X)


def m(indexes, powers):
    m = np.ones(X.shape[0])
    k = 0
    for i in indexes:
        m *= X[:, i]**(powers[k])
        k += 1
        
    return m@w
    #print(m)

for i in range(len(X)):
    plt.plot(X[i])
    print(np.linalg.norm(X[i]))
    print(np.sqrt(N + 2.))
plt.show()

#print(m([1], [2]))
quit()

def partition(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))

    return answer

print(partition(5))
