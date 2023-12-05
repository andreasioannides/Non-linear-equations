from math import tan, cos, atan, exp, log, e, tanh, sqrt


def sec(x):
    return 1/cos(x)

def sech(x):
    return 2.0/(exp(x) + exp(-x))

def f(x):
    return tan(x)*tanh(x) + 1.0

def fd(x):  # derivative of f 
    return (sec(x)**2)*tanh(x) + (sech(x)**2)*tan(x)

'''The rest g(x) functions are used in Fixed Point 
   Method.'''

def g1(x):
    return x - tan(x)*tanh(x) - 1.0

def g2(x):
    return x + tan(x)*tanh(x) + 1.0

def g3(x):
    # The logarithmic base by default is 'e'
    return 0.5*log(1 - ((e**(2*x) + 1.0)/tan(x)))  

def g4(x):
    return sqrt(x**2 - tan(x)*tanh(x) - 1)

def g5(x):  
    return atan(-1.0/tanh(x)) 