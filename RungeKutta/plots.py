import numpy as np
import matplotlib.pyplot as plt
    
t = np.loadtxt('timeStamps.txt')
numeric = np.loadtxt('heatSolution.txt', delimiter=',')

for step in range(200, 1200, 200):
    t_fixed = t[step]
    x = np.arange(0.1, 1.0, 0.1)
    n = np.arange(1, 10000).reshape(-1, 1)
    exact = 2/np.pi*np.sum((-1)**(n+1)/n*np.exp(-0.01*np.pi**2*n**2*t_fixed)*np.sin(np.pi*n*x), axis=0)
    #plt.plot(x, exact, label='аналитическое', marker='.')
    #plt.plot(x, numeric[step,:], label='численное', marker='.')
    plt.plot(x, np.abs(exact - numeric[step,:]), label=f't = {t_fixed}', marker='.')
plt.legend(fontsize='small')
plt.title(f'погрешности при разных t')
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.show()