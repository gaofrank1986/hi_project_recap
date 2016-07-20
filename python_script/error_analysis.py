import numpy as np
a = np.loadtxt("./fort.6000")
b=np.arange(17)
b=b*0.05
import matplotlib.pyplot as pl
c = np.loadtxt("./fort.6001")
d=np.arange(172)
d=d*0.005
pl.plot(d,c,'r',label=r'$\delta t=0.05$')
pl.plot(b,a,'g',label=r'$\delta t=0.5$')
pl.xlabel('time(s)')
pl.ylabel('error')
pl.title('HOBEM w/ hyper singular')
pl.legend(loc='upper left')
pl.show()
