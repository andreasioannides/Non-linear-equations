"""
    Documentation:
    The "rootClosedMethods" module contains two methods:
    a)Bisection
    b)Regula Falsi
    to calculate the root of the given function in a 
    specified closed interval.

    Inputs:
    f        : function
    (a,b)    : interval which is searched to calculate the root
    n        : number of iterations 
    epsilon  : a threshold for accepted deviation between x_new 
               and x_old which defines if the method converges
"""


from math import ceil, log, tan, exp
import numpy as np


def validInterval(fa, fb):
    if (fa*fb > 0.0e0):
        raise Exception("Invalid interval.")


def result(method, converged, epsilon, x, fx, iterations):
    res = {  
            "Method" : method,
            "Converged" : converged,
            "Epsilon" : epsilon,
            "x" : x,
            "f(x)" : fx,
            "Iterations" : iterations
           } 
    return res


def bisectionFunc(a, b, fa=None, fb=None):
    '''Divide the interval according to the Bisection method.
       For the rest code is referred as "div" function
       (division criterion)'''

    return (a + b)/2


def regulaFalsiFunc(a, b, fa, fb):
    '''Divide the interval according to the Bisection method.
       For the rest code is referred as "div" function
       (division criterion)'''
    
    return a - ((fa*(b-a))/(fb-fa))


def search(f, a, b, fa, fb, n, div, epsilon):
    conv = False   #converged
    x_old = a
    roots = np.array([[x_old, f(x_old)]], dtype=np.float64)

    for i in range(n):
        x = div(a, b, fa, fb)
        fx = f(x)
        roots = np.append(roots, [[x, fx]], axis=0)

        print(f"Iter={i+1}, x={x}, f(x)={fx}")

        # if (abs((x - x_old)/x) <= epsilon):
        if (abs((x - x_old)/x) <= epsilon):
            conv = True
            break
        
        x_old = x

        if (fx*fa > 0):
            a = x
            fa = fx
        elif (fx*fb > 0):
            b = x
            fb = fx

    return conv, x, fx, i+1, roots


def bisection(f, a, b, n, epsilon):  
    # iters = ceil(log(abs(b - a)/epsilon)/log(2.0)) 
    print("Bisection Method:\n")     
    
    fa = f(a)
    fb = f(b)

    validInterval(fa, fb)

    div = bisectionFunc
    conv, x, fx, iters, roots = search(f, a, b, fa, fb, n, div, epsilon)

    print(result('Bisection', conv, epsilon, x, fx, iters))

    return roots


def regulaFalsi(f, a, b, n, epsilon):  
    # iters = ceil(log(abs(b - a)/epsilon)/log(2.0))
    print("Regula Falsi Method:\n")      

    fa = f(a)
    fb = f(b)

    validInterval(fa, fb)

    div = regulaFalsiFunc
    conv, x, fx, iters, roots = search(f, a, b, fa, fb, n, div, epsilon)

    print(result('Regula Falsi', conv, epsilon, x, fx, iters))

    return roots