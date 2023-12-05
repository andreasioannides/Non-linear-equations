import matplotlib.pyplot as plt
import numpy as np


'''Plot iterations-x Diagram''' 

def roots(roots):  
    iters = roots.shape[0]
    n = np.arange(0, iters)

    plt.plot(n, roots[:,0], color='blue', linewidth=1, label='Iteration roots')
    plt.scatter(n, roots[:,0], color='blue', s=4)
    plt.scatter(iters-1, 2.34704557e+00, color='red', s=20, label='root')
    
    plt.title('Convergence Diagram')
    plt.xlabel('n')
    plt.ylabel('x')
    plt.xlim(0, iters)
    plt.legend()
    plt.show()

