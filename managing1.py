import sympy as sym
from sympy.abc import t

#одинарные скобки - вертикальный вектор, двойные - горизонтальный
P=sym.matrices.Matrix([[0,1], [0,0]])
Q = sym.matrices.Matrix([[t,0],[0,1]])
f = sym.matrices.Matrix([0, 0])
x0 = sym.matrices.Matrix([0,0])
x1 = sym.matrices.Matrix([0, 1])
T = 5 #не забудь поменять


Y = sym.matrices.zeros(2)
Y = sym.simplify(sym.exp(P * t))
print("Y=", Y)

B = Y.inv() * Q
print('\nB=', B)

A = sym.simplify(sym.integrate(B * B.T, (t, 0, T)))
print('\nA=', A)

eta = Y.inv().subs(t, str(T)) * x1 - x0 - sym.simplify(sym.integrate(Y.inv() * f, (t, 0, T)))
print('\neta=', eta)

C = sym.simplify(A.inv() * eta)
print('\nC=', C)

u = sym.simplify(B.T * C)
print('\nu(t)=', u, "+v(t)")


'''
def matexp_tailor(a, t, eps, max_terms):
    n = a.shape[0]

    b = np.eye(n)
    s = np.eye(n)
    k = 1

    for i in range(1, max_terms):
        b = b @ a
        k *= t / i
        term = b * k
        s += term

        if np.linalg.norm(term) < eps:
            break

    return s


def matexp(a, t, eps=1e-12, max_terms=100):
    if abs(t) < eps:
        return np.eye(a.shape[0])

    p = max(0, math.ceil(np.log2(abs(t) * np.linalg.norm(a))))
    b = matexp_tailor(a, t / 2**p, eps, max_terms)
    for _ in range(p):
        b = b @ b
    return b


def int_simpson(f, a, b, n=16):
    assert n % 2 == 0

    h = (b - a) / n
    f_a = f(a)
    s1 = np.zeros_like(f_a)
    s2 = np.zeros_like(f_a)

    for i in range(1, n//2):
        s1 += f(a + 2 * i * h)

    for i in range(1, n//2 + 1):
        s2 += f(a + (2 * i - 1) * h)

    return h/3 * (f_a + 2 * s1 + 4 * s2 + f(b))


def compute_a(P, Q, T):
    def BBT(t):
        B = matexp(P, -t) @ Q(t)
        result = B @ B.T

        return result

    return int_simpson(BBT, 0, T)

def compute_eta(P, T, f, x0, x1):
    eta = matexp(P, -T) @ x1 - x0
    eta -= int_simpson(lambda t: matexp(P, -t) @ f(t), 0, T)

    return eta
'''


'''
A = compute_a(P, Q, T)
eta = compute_eta(P, T, f, x0, x1)

print('A', A, sep='\n')
print('eta', eta, sep='\n')

try:
    c = np.linalg.solve(A, eta)
    print('c', c, sep='\n')

    def u(t):
        B = matexp(P, -t) @ Q(t)
        return B.T @ c

    def int_fn(t):
        return matexp(P, -t) @ (Q(t) @ u(t) + f(t))

    x_T = matexp(P, T) @ (x0 + int_simpson(int_fn, 0, T))
    print('x(T)', x_T, sep='\n')

except np.linalg.LinAlgError:
    rank_a = np.linalg.matrix_rank(A, 1e-12)
    a_eta = np.hstack((A, eta))
    rank_aeta = np.linalg.matrix_rank(a_eta, 1e-12)

    print('Система  ', end='')
    if rank_a == rank_aeta:
        print('управляема')
    else:
        print('не управляема')
'''