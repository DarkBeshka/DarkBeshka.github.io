"""
По заданным начальным данным(матрицам Р, Q) построить стабилизирующее управление - проверить систему на полную управляемость, найти вектор С
"""

import numpy as np
import sympy as sym
from sympy.abc import z

print('n')
n = int(input())
print('r')
r = int(input())

P = np.zeros((n, n))
Q = np.zeros((n, r))

print('P')
for i in range(n):
    s = input()
    s = s.split()
    for j in range(n):
        P[i][j] = float(s[j])

print('Q')
for i in range(n):
    s = input()
    s = s.split()
    for j in range(r):
        Q[i][j] = float(s[j])


#Проверяем на полную управляемость
S = Q
for i in range(1, n):
    S = np.concatenate((S, np.dot(np.linalg.matrix_power(P, i), Q)), axis=1)
print('S=', S)

rS = np.linalg.matrix_rank(S)  #ранг матрицы S
if rS == n:
    print('Система полностью управляемая')
    T = S
else:
    print('Система не полностью управляемая, rankS =', rS)
    T = S[:, :rS]   #выбираем первые rS столбцов матрицы Т


#Дополняем до полного базиса матрицу Т
k = np.linalg.matrix_rank(T)
#print('k=', k)
for i in range(n):
    x = np.zeros((n, 1))
    x[i][0] = 1
    if np.linalg.matrix_rank(np.concatenate((T, x), axis = 1)) > np.linalg.matrix_rank(T): #склейка по горизонтали
        T = np.concatenate((T, x), axis = 1)


print('\nmu=')
s = input()
mu = np.zeros((k, 1))
s = s.split()
for i in range(k):
    mu[i][0] = float(s[i])

print('\nT=', T)

#Строим матрицы Р, Q c волнами
Pt = np.dot(np.linalg.matrix_power(T, -1), P)
Pt = np.dot(Pt, T)
Qt = np.dot(np.linalg.matrix_power(T, -1), Q)

print('\nPt=', Pt)
print('\nQt=', Qt)

P11 = Pt[:k, :k]  #Выбираем матрицу Р11 с волной размерности kxk
print('\nP11=', P11)
P22=Pt[k:, k:]
print('\nP22=', P22)
harpoly=np.poly(P22)
check=True
for i in range(len(harpoly)):
    if harpoly[i]<0:
        check=False
if check==False:
    print("Cистема не стабилизируема")
else:
    p = np.zeros(k + 1)
    phi1 = np.poly(P11)[1:]
    phi1 = np.transpose(phi1)
    print('\nphi1=', phi1)
    K = np.eye(n)   #Матрица К2 - единичная соответствующей размерности, так что считать заново не надо
    for i in range(k):
        for j in range(i + 1, k):
            K[i][j] = phi1[j - i - 1]

    print('\nK=', K)

    #Вычисляем коэффициенты пси
    y = z - int(mu[0])
    for i in range(1, k):
        y *= z - int(mu[i])
        y = sym.expand(y)

    psi = np.zeros(k)
    for i in range(k):
        psi[i] = float(sym.core.Expr.coeff(y, z, k - 1 - i))
    print('\npsi=', psi)

    y1 = phi1 - psi #вектор гамма1

    print('\ny1=', y1)

    G = np.zeros(n)
    for i in range(k):
        G[i] = y1[i]
    print('\nG=', G)  #вектор гамма
    C = np.dot(G, np.linalg.matrix_power((np.dot(T, K)), -1))
    print('\nC=', C)
