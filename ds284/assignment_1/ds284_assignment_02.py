from numpy.linalg import norm
from numpy.random import randn
import numpy
x = []
Ax_norms = [0,0,0,0,0,0]
A_norms = []
A = randn(100, 2)
for i in range(1000):
    temp = randn(2, 1)
    temp_normed = temp/norm(temp)
    x.append(temp_normed)
    # calculate matrix norms
    for p in range(6):
        norm_of_Ax = 0
        if p == 0:
            norm_of_Ax = norm(A.dot(x[i]), numpy.inf)
        else:
            norm_of_Ax = norm(A.dot(x[i]), p)
        if(Ax_norms[p] < norm_of_Ax):
            Ax_norms[p] = norm_of_Ax
    # calculate vector norms

for p in range(3):
    if p == 0:
        norm_of_A = norm(A, numpy.inf)
    else:
        norm_of_A = norm(A, p)
    A_norms.append(norm_of_A)

print("NORM \t\tA_norm Ax_norm")
for i in range(len(Ax_norms)):
    if i == 0:
        which_norm = "inf-norm"
    else:
        which_norm = "  " + str(i) + "-norm"
    print(which_norm + "\t" + str(round(A_norms[i],3)) + "\t" + str(round(Ax_norms[i],3)))
