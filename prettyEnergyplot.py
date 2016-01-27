import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pyl
plt.rc('font', family='serif',size = 24)
#read data
coeff="coeff.txt"
bas="basis.txt"
coefficents=np.loadtxt(coeff,dtype=np.float_)
basis=np.loadtxt(bas,dtype=np.float_)
numpts =1000
xminimum= 0
xmaximum = 5
yminimum= -1
ymaximum = 1
xvec = np.linspace(xminimum,xmaximum,numpts)
def f(r,a,c):
	temp = 0.0
	for i in range(4):
		temp+= np.abs(c[i]) * np.exp(-a[i]*(r**2))
	return temp

fig = plt.figure(facecolor='white',figsize=(8,8))
ax = fig.add_subplot(111)
plt.xlim(xminimum, xmaximum)
plt.ylim(yminimum, ymaximum)
yvec = np.copy(xvec)
#cm = pyl.get_cmap('jet')
cm = ['red','orange','green','blue']
for i in range(1):
	for j in range(numpts):
		yvec[j] = f(xvec[j],basis,coefficents[i])#*f(xvec[j],basis,coefficents[i])*xvec[j]
    xvec2 = xvec + ((xvec[1] +  xvec[0])/2.0)
    temp = 0.0
    for i in range(len(yvec)-1):
       temp += xvec2[i]*((yvec[2]+yvec[1])/2.0) 
	plt.plot(xvec,yvec/temp,color='blue',linewidth=2.1)
	plt.plot(xvec,np.exp(-xvec),color='red',linewidth=2.1)
	plt.plot([xminimum, xmaximum], [0, 0], color='black', linestyle='--', linewidth=2)
for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(3)
#plt.xticks([0,np.pi,2*np.pi],['0',r'$\pi$',r'$2\pi$'])
ax.xaxis.set_tick_params(width=2.5,length=7)
ax.yaxis.set_tick_params(width=2.5,length=7)
plt.xlabel(r'Separation($\AA$)')
plt.ylabel(r'Energy($\frac{\mathrm{kcal}}{\mathrm{mol}}$)')
plt.title(r'Dimer Energy',y=1.08)
fig.tight_layout()
fig.savefig('prettyplot.png')
#plt.show()
