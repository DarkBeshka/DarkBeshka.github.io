#Task: По начальным данным(матрицам Р и Q, вектору собственных чисел mu) построить программное управление - выяснить, является ли система полностью управляемой,
#найти матрицу А и вектор nu


import numpy as np
from numpy.linalg import solve
import sympy as sym
from sympy.abc import k

# Ищем матрицу А====================================
P=sym.matrices.Matrix([[0,1,0,0], [0,0,1,0],[0,0,0,1],[1,0,0,0]])
Q = sym.matrices.Matrix([[1,0],[1,1],[1,0],[0,1]])
f = sym.matrices.Matrix([0, 0,0,0])
x0 = sym.matrices.Matrix([0,1,1,0])
x1 = sym.matrices.Matrix([1, 0, 0, 1])
m=2

temp=1
for j in range(m-1,0,-1):
    temp=temp*P.subs(k, j)
A=temp*Q.subs(k,0)
temp=1
for i in range(1,m-1):
    print('i=',i)
    for j in range (m-1,i,-1):
        temp=temp*P.subs(k, j)
        print(j)
    temp=temp*Q.subs(k, i)
    A=A.row_join(temp)
    temp=1
A=A.row_join(Q.subs(k, m-1))
print("A=", A)
if A.rank()==sym.shape(P)[1]:
    print("Система полностью управляема")
else:
    print("Система НЕ полностью управляема")

# Ищем ню==================================
nu=x1
temp=1
for j in range(m-1,-1,-1):
        temp=temp*P.subs(k, j)
        #print(j)

nu=nu-temp*x0
for i in range(1,m):
    #print("i=",i)
    for j in range(m-1,i-1,-1):
        temp=temp*P.subs(k, j)
        #print(j)
    temp=temp*f.subs(k, i-1)
    nu=nu-temp

nu=nu-f.subs(k,m-1)
print("nu=",nu)



#U = np.linalg.lstsq(A, eta)

#x = np.dot(P, x0) + np.dot(Q, U[0:2, :])
#x = np.dot(P, x) + np.dot(Q, U[2:4, :])
#print('U', U)
#print(x)
