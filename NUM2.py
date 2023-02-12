"""
Для заданной функции y = sqrt(x) + cos(x) выполнить следующее:
1. Выбрать из области её определения интервал непрерывности.
2. Назначить некоторое число равноотстоящих узлов и построить интерполяционный полином Лагранжа. Проделать то же самое, взяв
в качестве узлов то же количество узлов, но определяемых по формуле Чебышевского. Сравнить полученные результаты(вывести методическую погрешность).
3. Повторить п.2, увеличивая число узлов.
4. Повторить пункты 2 и 3, строя полином Ньютона.
5. Построить графики
"""


import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *


N = 20
a = 0.2
b = 5
plt.rcParams['figure.figsize'] = [10, 10]

#ищем функцию Лагранжа
"""
         Получить лагранжеву интерполяцию
         : param x: список значений x
         : param y: список значений y
         : param a: число для интерполяции(в формуле - это вместо х - мы строим по точкам, а не символьно)
         : return: вернуть результат интерполяции
"""
def Lagrange(x, y, a):
    res = 0
    for i in range(len(y)):  #sum i=0...n
        t1 = 1
        t2 = 1
        for j in range(len(x)):
            if i != j:
                t1 *= (a - x[j]) #ищем числитель(наш узловой многочлен)
                t2 *= (x[i] - x[j]) #знаменатель
        res += y[i] * t1 / t2  #не забываем домножить на f(xi)
    #print(simplify(t_t))
    return res


def f(x):
    return cos(x)+sqrt(x)

def h(x):
    return abs(x)*(cos(x)+sqrt(x))


def diff_f(X, N):
    x = Symbol('x')
    y_d = cos(x)+sqrt(x)
    for i in range(N+1):
        y_d = y_d.diff(x)
    y_d = lambdify(x,y_d)
    return(y_d(X))


#для выражения погрешности Лагранжа
def R(x_, x,N):
    w = 1
    for i in range(N):
        w *= (x_ - x[i])
    #print(diff_f(ksi, N))
    r = 1/math.factorial(N+1)*w
    return(r)

def w(t, x):
    W = 1
    for i in range(len(x)):
        W *= (t-x[i])
    return(W)

#разделенные разности
def div_dif(x):
    ans = 0
    for i in range(len(x)):
        temp = 1
        for j in range(len(x)):
            if j == i:
                continue
            temp *= (x[i] - x[j])
        ans += (f(x[i]) / temp)
    return ans

def Newton(t, x):
    ans = f(x[0])
    temp = 1
    for i in range(len(x)):
        temp *= (t - x[i])  #w
        ans += div_dif(x[: i + 2]) * temp    #+2, так как начинаем с первого слагаемого, а цикл с 0-го
    return ans

def err(x_,y_):    #обычное условие метод.погрешности - модель разницы рез-та и функции
    if (abs(y_- f(x_)) == 0):
        return 0
    else:
        return(math.log10(abs(y_- f(x_))))  #будем выводить десятичный логарифм от погрешности(для наглядности)


#Лагранж==================================================================
x_points = np.linspace(a, b, num = N)
y_points = [f(i) for i in x_points]

xnew = np.linspace(np.min(x_points), np.max(x_points), 50)
ynew = [Lagrange(x_points, y_points, i) for i in xnew]

x = np.linspace(a, b, num = 50)
y = [f(i) for i in x]

plt.subplot(2,2,1)
plt.grid()
plt.title('Лагранж')
#plt.plot(x_points,y_points,'o', label = 'Равномерное - узлы')
plt.plot(x, y, 'r', label = 'Функция')
plt.plot(xnew, ynew, 'b', label = 'Равномерное')
#==========================================================================


plt.subplot(2,2,3)
plt.grid()
#plt.ylim([-1,1])
plt.title('Узловой многочлен')
plt.plot(x, [w(i, x_points) for i in x], 'b', label = "Равномерное")

plt.subplot(2,2,4)
plt.grid()
plt.title('Многочлен Ньютона')
plt.plot(x, [f(i) for i in x], 'r', label = "Функция")
plt.plot(x, [Newton(i, x_points) for i in x], 'b', label = "Равномерное")

#погрешность=============================================================
ksi = abs(b-a)/2
diff = diff_f(ksi,N)
y_test=[f(i) for i in range(0,50)]
for i in range(50):
    y_test[i] = abs(diff*R(x[i],x_points,N))
    if y_test[i]!=0:
        y_test[i] = math.log10(y_test[i])
    else:
        y_test[i] = 0
    #y[i] = err(xnew[i], ynew[i])


plt.subplot(2,2,2)
#plt.ylim([-0.01,0.02])
#plt.xlim([0, math.pi])
plt.grid()
plt.title('Методическая погрешность, ln10')
#plt.plot(x, y_test, 'b', label = "Равномерное, формула")
plt.plot(x, [err(xnew[i], ynew[i]) for i in range(50)], 'b', label = "Равномерное") 
#plt.plot(x, [err(xnew[i], [Newton(i, x_points) for i in x][i]) for i in range(50)], 'b', label = 'Равномерное')
#======================================================================

#Чебышевское, Лагранж==================================================
for i in range(N):
    x_points[i] = 1/2*((b-a) * math.cos((2*i+1)*math.pi/(2*(N+1)))+(b+a))
    y_points[i] = f(x_points[i])

xnew = np.linspace(np.min(x_points), np.max(x_points), 50)
ynew = [Lagrange(x_points, y_points, i) for i in xnew]

plt.subplot(2,2,1)
plt.plot(xnew, ynew, 'orange', label = 'Чебышевское')
#plt.plot(x_points, y_points, '*', label = 'Чебышевские узлы')
plt.legend()
#======================================================================

#Чебышевское, узловой многочлен========================================
plt.subplot(2,2,3)
plt.plot(x, [w(i, x_points) for i in x], 'orange', label = 'Чебышевское')
plt.legend()
#======================================================================
#Ньютон, чебышевское===================================================
plt.subplot(2,2,4)
plt.plot(x, [Newton(i, x_points) for i in x], 'orange', label = 'Чебышевское')
plt.legend()
#=======================================================================

X = Symbol('X')
res = 0
for i in range(N):
    t1 = 1
    t2 = 1
    for j in range(N):
        if i != j:
            t1 *= (X - x_points[j])
            t2 *= (x_points[i] - x_points[j])
    res += y_points[i] * t1 / t2
print("L=", simplify(res))


ksi = abs(b-a)/2
diff = diff_f(ksi,N)
y_test=[f(i) for i in range(0,50)]
for i in range(0, 50):
    y_test[i] = abs(diff*R(x[i],x_points,N))
    if y_test[i]!=0:
        y_test[i] = math.log10(y_test[i])
    else:
        y_test[i] = 0

plt.subplot(2,2,2)
plt.plot(x, [err(xnew[i], ynew[i]) for i in range(50)], 'orange', label = 'Чебышевское')
#plt.plot(x, y_test, 'orange', label = 'Чебышевское')

#plt.plot(x, [err(xnew[i], [Newton(i, x_points) for i in x][i]) for i in range(50)], 'orange', label = 'Чебышевское')
plt.legend()
plt.show()


