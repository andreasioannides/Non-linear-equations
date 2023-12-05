"""
    Documentation:
    The "rootOpenMethods" module contains four methods:
    a)Fixed Point 
    b)Newton&Raphson
    c)Secant
    d)Muller
    to calculate the root of the given function. Only for Muller 
    method two different approaches are proposed, muller and muller1.

    Inputs:
    f          : function
    x1, x2, x3 : initial points 
    n          : number of iterations 
    epsilon    : a threshold for accepted deviation between x_new 
                 and x_old which defines if the method converges
"""


from rootClosedMethods import result
from math import sqrt
import numpy as np


def fixedPoint(f, g, x, n, epsilon):
    print("Fixed Point Method:\n")     
    
    fx = f(x)
    roots = np.array([[x, fx]], dtype=np.float64)

    conv = False
    x_old = x

    for i in range(n):
        x_new = g(x_old)

        roots = np.append(roots, [[x_new, f(x_new)]], axis=0)
        print(f"Iter={i+1}, x={x_new}")

        if (abs((x_new - x_old)/x_new) <= epsilon):
            conv = True
            x = x_new
            fx = f(x_new)
            break

        x_old = x_new
            
    iters = i + 1

    print(result('Fixed Point', conv, epsilon, x, fx, iters))

    return roots

def newraph(f, fd, x, n, epsilon):  # fd:derivative of f 
    print("Newton-Raphson:\n") 

    fx = f(x)
    roots = np.array([[x, fx]], dtype=np.float64)

    conv = False
    x_old = x

    for i in range(n):
        fdx = fd(x_old)

        if (fdx == 0):
            iters = i+1
            x = x_old
            fx = f(x)
            break
        
        x_new = x_old - (fx/fdx)

        roots = np.append(roots, [[x_new, f(x_new)]], axis=0)
        print(f"Iter={i+1}, x={x_new}")

        if (abs((x_new - x_old)/x_new) <= epsilon):
            conv = True
            x = x_new
            fx = f(x_new)
            break

        x_old = x_new
        fx = f(x_new)

    iters = i + 1

    print(result('Newton-Raphson', conv, epsilon, x, fx, iters))

    return roots


def secant(f, x1, x2, n, epsilon):
    print("Secant:\n") 

    fx1 = f(x1)
    fx2 = f(x2)

    roots = np.array([[x1, fx1]], dtype=np.float64)

    conv = False
    x = None
    fx = None
    
    for i in range(n):        
        x_new = x1 - fx1*((x1 - x2)/(fx1 - fx2))
        roots = np.append(roots, [[x_new, f(x_new)]], axis=0)

        print(f"Iter={i+1}, x={x_new}")

        if (abs((x_new - x1)/x_new) <= epsilon):
            conv = True
            fx = f(x_new)
            break
        
        x2 = x1
        x1 = x_new
        fx2 = fx1
        fx1 = f(x_new)
  
    iters = i + 1

    print(result('Secant', conv, epsilon, x, fx, iters))

    return roots


def sort(x1, x2, x3):
    nums = [x1, x2, x3]
    nums.sort()

    return nums


def muller(f, x0, x1, x2, n, epsilon):
    print("Muller:\n") 

    nums = sort(x0, x1, x2)  #x0=nums[1], x1=nums[2], x2=nums[0] 

    roots = np.array([[nums[2], f(nums[1])]], dtype=np.float64)

    conv = False
    x = None
    fx = None

    x_old = nums[1]

    for i in range(n):
        fx0 = f(nums[1])
        fx1 = f(nums[2])
        fx2 = f(nums[0])

        h1 = abs(nums[2] - nums[1])
        h2 = abs(nums[1] - nums[0])

        a0 = fx0
        a2 = (h2*fx1 - (h1+h2)*fx0 + h1*fx2)/(h1*h2*(h1+h2))
        a1 = (fx1 - fx0 - a2*h1**2)/h1

        den1 = a1 + sqrt(a1**2 - 4*a2*a0)  #denominator
        den2 = a1 - sqrt(a1**2 - 4*a2*a0)

        if (abs(den1) < abs(den2)):   #larger den -> x3 closer to x0
            x3 = nums[1] + (-2*a0)/den2
        else:
            x3 = nums[1] + (-2*a0)/den1

        fx3 = f(x3)
        roots = np.append(roots, [[x3, fx3]], axis=0)

        print(f"Iter={i+1}, x={x3}")

        if (abs((x3 - x_old)/x3) <= epsilon):
            conv = True
            fx = fx3
            break

        if (fx3*f(nums[2]) < 0):
            nums[0] = nums[1]
            nums[1] = x3
        else:
            nums[2] = nums[1]
            nums[1] = x3
        
        x_old = x3

    iters = i + 1

    print(result('Muller', conv, epsilon, x, fx, iters))

    return roots


def muller1(f, x0, x1, x2, n, epsilon):
    print("Muller:\n") 

    nums = [x0, x1, x2]
    y = [f(nums[0]), f(nums[1]), f(nums[2])]

    roots = np.array([[nums[2], f(nums[1])]], dtype=np.float64)

    conv = False
    x = None
    fx = None

    x_old = nums[2]

    for i in range(n):
        coeffs = np.array([[nums[0]**2, nums[0], 1], 
                           [nums[1]**2, nums[1], 1], 
                           [nums[2]**2, nums[2], 1]])

        b = np.array([y[0], y[1], y[2]])

        a2, a1, a0 = np.linalg.solve(coeffs, b)

        dis = a1**2 - 4*a2*a0

        if (dis < 0):
            print("Discriminant negative.")
            x = None
            fx = None
            break

        x31 = (-a1 + sqrt(dis))/(2*a2)
        x32 = (-a1 - sqrt(dis))/(2*a2)

        d1 = abs(nums[0] - x31) + abs(nums[2] - x31)
        d2 = abs(nums[0] - x32) + abs(nums[2] - x32)

        c = 0
        for j in y:
            if (j < 0):
                c+=1

        if (c == 0 or c == 1):
            idx = y.index(max(y))
        elif (c == 2 or c == 3):
            idx = y.index(min(y))

        if (d1 < d2):
            nums[idx] = x31
            y[idx] = f(x31)
        else:
            nums[idx] = x32
            y[idx] = f(x32)
        
        x_new = nums[idx]
                 
        roots = np.append(roots, [[x_new, y[idx]]], axis=0)

        print(f"Iter={i+1}, x={x_new}")

        if (abs((x_new - x_old)/x_new) <= epsilon):
            conv = True
            x = x_new
            fx = f(x)
            break

        x_old = x_new        

    iters = i + 1

    print(result('Muller', conv, epsilon, x, fx, iters))

    return roots  