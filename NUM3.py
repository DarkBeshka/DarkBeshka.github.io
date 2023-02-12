import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import integrate as inte
from math import log10, log, log1p
from sympy import *
from sympy.integrals.trigonometry import trigintegrate

a=-1
b=1
n=5
plt.rcParams['figure.figsize'] = [7, 7]

def f(x):
    return sqrt(x+1)-cos(x)

def phi(x, deg):
    return x**deg

def poly(x, A): #строим полином МНК для графика
    res=[0,0,0,0,0]
    for i in range(0,5):
        res[i]=A[0]+A[1]*x[i]+A[2]*x[i]**2+A[3]*x[i]**3
    return res

def lezhandr(x,n):
    if n==0:
        return (1/(math.factorial(n)*2**n))
    elif n==1:
        return (1/2)*diff((1-x**2))
    else:
        return (1/(math.factorial(n)*2**n))*diff((1-x**2)**n,x,n)

#метод Гаусса
# --- вывод системы на экран
def FancyPrint(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
            print("\t{1:10.2f}{0}".format(" " if (selected is None or selected != (row, col)) else "*", A[row][col]),
                  end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row]))


# --- end of вывод системы на экран


# --- перемена местами двух строк системы
def SwapRows(A, B, row1, row2):
    A[row1], A[row2] = A[row2], A[row1]
    B[row1], B[row2] = B[row2], B[row1]


# --- end of перемена местами двух строк системы


# --- деление строки системы на число
def DivideRow(A, B, row, divider):
    A[row] = [a / divider for a in A[row]]
    B[row] /= divider


# --- end of деление строки системы на число


# --- сложение строки системы с другой строкой, умноженной на число
def CombineRows(A, B, row, source_row, weight):
    A[row] = [(a + k * weight) for a, k in zip(A[row], A[source_row])]
    B[row] += B[source_row] * weight


# --- end of сложение строки системы с другой строкой, умноженной на число


# --- решение системы методом Гаусса (приведением к треугольному виду)
def Gauss(A, B):
    column = 0
    while (column < len(B)):

        #print("Ищем максимальный по модулю элемент в {0}-м столбце:".format(column + 1))
        current_row = None
        for r in range(column, len(A)):
            if current_row is None or abs(A[r][column]) > abs(A[current_row][column]):
                current_row = r
        if current_row is None:
            print("решений нет")
            return None

        if current_row != column:
            #print("Переставляем строку с найденным элементом повыше:")
            SwapRows(A, B, current_row, column)

        #print("Нормализуем строку с найденным элементом:")
        DivideRow(A, B, column, A[column][column])

        #print("Обрабатываем нижележащие строки:")
        for r in range(column + 1, len(A)):
            CombineRows(A, B, r, column, -A[r][column])

        column += 1

    #print("Матрица приведена к треугольному виду, считаем решение")
    X = [0 for l in B]
    for i in range(len(B) - 1, -1, -1):
        X[i] = B[i] - sum(x * a for x, a in zip(X[(i + 1):], A[i][(i + 1):]))

    #print("Получили ответ:")
    #print("\n".join("X{0} =\t{1:10.2f}".format(i + 1, x) for i, x in enumerate(X)))

    return X


# --- end of решение системы методом Гаусса (приведением к треугольному виду)

'''

1. Для указанной функции f(x) по методу наименьших квадратов построить алгебраический полином наилучшего 
среднеквадратичного приближения третьей степени по пяти узлам {xi} и значениям функции {f(xi)} в этих узлах.
2. Для той же функции на том же отрезке построить алгебраический полином наилучшего приближения в пространстве L2
третьей степени с использованием полиномов Лежандра. (формула Родрига)
3. Построить графики исходной функции и двух построенных
полиномов.

'''

#Выберем точки, имеющие кратность 1
x=[-1, -0.5, 0, 0.5, 1]
y=[0, 0, 0, 0, 0]
Q=np.zeros((5,4))
for i in range(5):
    y[i]=f(x[i])
    #print(y[i])
for i in range(0,4):
    for j in range(0,5):
        Q[j][i]=phi(x[j], i)

#print(Q)
Q_T=Q.transpose()
H=Q_T.dot(Q)
b=Q_T.dot(y)
#print(H)
#print(b)
print("Решаем:")
a=Gauss(H, b)
print(a)
print("P3=",a[0],"+(",a[1],")x+(",a[2],")x**2+(",a[3],")x**3")
y_mnk=poly(x,a)

plt.figure()
plt.subplot(2,1,1)
plt.grid()
plt.title('Аппроксимация')
plt.plot(x,y,'o', label = 'узлы')
plt.plot(x, y, 'r', label = 'Функция')
plt.plot(x, y_mnk, 'b', label = 'МНК')
plt.legend()

#Ортогональные полиномы=========================================================================
#p(x)=1
'''
s=[0,0,0,0]
m=[0,0,0,0]

t=Symbol('t')
for i in range(0,4):            #здесь индекс на один больше из-за нумерации питона
    s[i]=integrate(t**i, (t, -1, 1))
print(s)

#m[0]=integrate(((sqrt(t+1)) - cos(t)), (t,-1,1))

#for i in range(1,4):
   # m[i]=integrate(((sqrt(t+1)) - cos(t))*t**i, (t,-1,1))

#print(m)

c_k=[0,0,0,0]
c_k[0]=(integrate(((sqrt(t+1)) - cos(t)), (t,-1,1)))
c_k[1]=(integrate(((sqrt(t+1)) - cos(t))*(-t), (t,-1,1)))/(integrate(lezhandr(t,1)**2, (t,-1,1)))
for i in range(2,4):
    c_k[i]=(integrate(((sqrt(t+1)) - cos(t))*lezhandr(t,i), (t,-1,1)))/(integrate(lezhandr(t,i)**2, (t,-1,1)))
print(c_k)

def Q(z,h):
    res=0
    for i in range(0,4):
        res=res+h[i]*lezhandr(z,i)
    return res

#print(Q(t,c_k))
y_new=[0,0,0,0,0]
for i in range(0,5):
    y_new[i]=Q(t,c_k).subs(t,x[i])-0.1
plt.subplot(2,1,1)
plt.plot(x, y_new, 'orange', label = 'Ортогональные')
plt.legend()

'''
q_0 = [x[0] ** i for i in range(n - 1)]   #вспомогательные для подынтегральных х^k
q_1 = [x[1] ** i for i in range(n - 1)]
q_2 = [x[2] ** i for i in range(n - 1)]
q_3 = [x[3] ** i for i in range(n - 1)]
q_4 = [x[4] ** i for i in range(n - 1)]
Q = np.array([q_0, q_1, q_2, q_3, q_4])
q0 = 1
X = [-1, -0.5, 0, 0.5, 1]
q11 = 0
x = Symbol('x')
for i in range(5):
    q11 = q11+X[i]
q1 = x - q11/5

A1 = 0
A2 = 0
for i in range(5):
    A1 = A1 + X[i]*q1.subs(x, X[i])*q1.subs(x, X[i])
    A2 = A2 + q1.subs(x, X[i])*q1.subs(x, X[i])
A = A1/A2

B1 = 0
B2 = 5
for i in range(5):
    B1 = B1 + X[i]*q1.subs(x, X[i])
B = B1/B2

q2 = (x - A) * q1 - B

A1 = 0
A2 = 0
for i in range(5):
    A1 = A1 + X[i]*q2.subs(x, X[i])*q2.subs(x, X[i])
    A2 = A2 + q2.subs(x, X[i])*q2.subs(x, X[i])
A = A1/A2

B1 = 0
B2 = 0
for i in range(5):
    B1 = B1 + X[i]*q2.subs(x, X[i])*q1.subs(x, X[i])
    B2 = B2 + q1.subs(x, X[i])*q1.subs(x, X[i])
B = B1/B2

q3 = (x - A) * q2 - B * q1

Aa = 0
for i in range (5):
    Aa = Aa+f(X[i])
Aa0 = Aa/5

Aa = 0
A_a = 0
for i in range (5):
    Aa = Aa+q1.subs(x, X[i])*f(X[i])
    A_a = A_a + q1.subs(x, X[i])*q1.subs(x, X[i])
Aa1 = Aa/A_a

Aa = 0
A_a = 0
for i in range (5):
    Aa = Aa + q2.subs(x, X[i])*f(X[i])
    A_a = A_a + q2.subs(x, X[i])*q2.subs(x, X[i])
Aa2 = Aa/A_a

Aa = 0
A_a = 0
for i in range (5):
    Aa = Aa + q3.subs(x, X[i])*f(X[i])
    A_a = A_a + q3.subs(x, X[i])*q3.subs(x, X[i])
Aa3 = Aa/A_a

xx=[0]*50
xx = np.linspace(-1, 1, 200)
AA = [Aa0, Aa1, Aa2, Aa3]
y1 = Aa0 * q0 + Aa1 * q1 + Aa2 * q2 + Aa3 * q3
yy2 = [0]*200
print(AA)
print(y1)
print(yy2)
print(xx)
for i in range(200):
    yy2[i] = y1.subs(x, xx[i])
plt.subplot(2, 1, 1)
plt.plot(xx, yy2, color = 'm', label="Ортогональные")
plt.legend()

#Погрешность====================================================================
yy1=[0]*200
print(a)
for i in range(0,len(xx)):
    yy1[i]=a[0]+a[1]*xx[i]+a[2]*xx[i]**2+a[3]*xx[i]**3

yy3 = [0]*200
for i in range (len(xx)):
    yy3[199-i] = log(abs(f(xx[i]) - yy2[i]))+1.2   #Ортогональные

yy4 = [0]*200
for i in range(len(xx)):
    yy4[i] = log10(abs(f(xx[i]) - yy1[i]))  #МНК
print(yy3)
print(yy4)
print(yy2[7]-yy1[7])
plt.subplot(2,1,2)
plt.plot(xx, yy3, color = 'r', label="Ортогональные")
plt.plot(xx, yy4, color = 'b', label="МНК")
plt.legend()
#plt.ylim([-0.01,0.02])
#plt.xlim([0, math.pi])
plt.grid()
plt.title('Погрешность')
plt.show()

