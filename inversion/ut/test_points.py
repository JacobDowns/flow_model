import numpy as np
from sigma_points import *

N = 5
x = np.zeros(N)
Pxx = np.eye(N)

points = SigmaPoints(x, Pxx)
X, w, w = points.__get_mysovskikh_set__(1.)
print(w.sum())


def m(indexes, powers):
    m = np.ones(X.shape[0])
    k = 0
    for i in indexes:
        m *= X[:, i]**(powers[k])
        k += 1
        
    return m@w
    #print(m)

print(4.*m([0,4], [2,2]))
quit()

def partition(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))

    return answer

print(partition(5))
