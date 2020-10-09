import numpy as np
import matplotlib.pylab as plt


data = np.loadtxt('FMS_IC.out')
p = data[:,1] / 1e5
T = data[:,2]
prc = data[:,3] /1e5

data = np.loadtxt('FMS_IC_2.out')
p2 = data[:,1] / 1e5
T2 = data[:,2]
prc2 = data[:,3] /1e5


fig = plt.figure()

plt.plot(T,p)
plt.plot(T2,p2)

plt.hlines(prc[0],100,4000,ls='solid')
plt.hlines(prc2[0],100,4000,ls='dashed')


print(prc[0])


plt.ylim(1e-4,1e4)
plt.xlim(100,4000)
plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()


plt.show()
