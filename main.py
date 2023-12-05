from functions import f, fd, g1, g2, g3, g4, g5
from rootClosedMethods import bisection, regulaFalsi
from rootOpenMethods import fixedPoint, newraph, secant, muller, muller1
from plot import roots


if __name__ == '__main__':
    a = 3.141592654/2.0
    b = a*2.0
    n = 400
    epsilon = 1.0e-7

    method = input("Select method: Bisection(b), Regula Falsi(r), Fixed Point(f), Newton-Raphson(n), Secant(s), Muller(m), Muller1(m1): ")

    if (method == 'b'):
        x = bisection(f, a, b, n, epsilon)
        print(x)
        roots(x)
    elif (method == 'r'):
        x = regulaFalsi(f, a, b, n, epsilon)
        print(x)
        roots(x)
    elif (method == 'f'):
        print("Choose 1 initial x:")
        # x_init = [2, 2.5, 3.1]
        x_init = float(input("Initial x: "))
        #change second parameter: g1,g2,g3,g4
        x = fixedPoint(f, g3, x_init, n, epsilon)
        print(x)
        roots(x)
    elif (method == 'n'):
        print("Choose 1 initial x.")
        x_init = float(input("Initial x: "))
        x = newraph(f, fd, x_init, n, epsilon)
        print(x)
        roots(x)
    elif (method == 's'):
        print("Choose 2 initial x.")
        x1 = float(input("Initial x1: "))
        x2 = float(input("Initial x2: "))
        x = secant(f, x1, x2, n, epsilon)
        print(x)  
        roots(x)
    elif (method == 'm'):
        print("Choose 3 initial x.")
        x1 = float(input("Initial x1: "))
        x2 = float(input("Initial x2: "))
        x3 = float(input("Initial x3: "))
        x = muller(f, x1, x2, x3, n, epsilon)
        print(x)  
        roots(x)
    elif (method == 'm1'):
        print("Choose 3 initial x.")
        x1 = float(input("Initial x1: "))
        x2 = float(input("Initial x2: "))
        x3 = float(input("Initial x3: "))
        x = muller1(f, x1, x2, x3, n, epsilon)
        print(x)  
        roots(x)
