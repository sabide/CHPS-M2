import matplotlib
import numpy as np
import matplotlib.pyplot as plt

data=np.fromfile('out.dat',sep=" ")
print data
x=data[0]
nelems=len(data)
data=data.reshape(nelems/3,3)
nelems=nelems/3
n=200
x=data[:,0].reshape(n,n)
y=data[:,1].reshape(n,n)
u=data[:,2].reshape(n,n)

# trcer des isovaleurs
fig=plt.figure()
ax1 = fig.add_subplot(111)
lev = np.linspace(np.min(u),np.max(u),num=10)
CF = ax1.contourf(x,y,u,levels=lev)
CS = ax1.contour(x,y,u,levels=lev, colors = 'k',)

cbar = plt.colorbar(CF, format='%.4f')
plt.show()
