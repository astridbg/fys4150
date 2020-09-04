import numpy as np
import matplotlib.pyplot as plt

def f(x):
    f = 100*np.exp(-10*x)
    return f

def u(x):
    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return u

n = int(input("Number of grid points: n = "))
h = 1./(n+1)
x_start = 0
x_end = 1
x = np.linspace(x_start,x_end,n+2)
v_start = 0
v_end = 0
v = np.zeros(n+2)

a = np.ones(n)*(-1.)
b = np.ones(n)*2.
c = np.ones(n)*(-1.)
d = h**2*f(x[1:n+1])
c_new = np.zeros(n)
d_new = np.zeros(n)

c_new[0] = c[0]/b[0]
d_new[0] = d[0]/b[0]

for i in range(1,n-1):
    c_new[i] = c[i]/(b[i] - a[i-1]*c_new[i-1])
    d_new[i] = (d[i] - a[i-1]*d_new[i-1])/(b[i] - a[i-1]*c_new[i-1])

d_new[n-1] = (d[n-1] - a[n-2]*d_new[n-2])/(b[n-1] - a[n-2]*c_new[n-2])
v[n-1] = d_new[n-1]

for i in range(n-2,0,-1):
    v[i] = d_new[i] - c_new[i]*v[i+1]

plt.figure()
plt.plot(x,u(x),label="closed-form solution")
plt.plot(x,v,label="results")
plt.legend()
plt.show()
