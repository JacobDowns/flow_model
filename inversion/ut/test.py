import numpy as np
from sigma_points import *

N = 5
x = np.zeros(N)
Pxx = np.eye(N)
points = SigmaPoints(x, Pxx)

X, wm, wc = points.get_set('fifth_order', 1.)

print(X)
#quit()
def m(cols, pows):
    print(np.prod((X[:, cols]@wm)**np.array(pows), axis = 1).sum()) 

for i in range(N):
    for j in range(N):
        for k in range(N):
            m([i,j,k], [2,1,1])  
