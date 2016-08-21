from numpy import *
cycle=3
tper=1.2
total_time=tper*cycle
x = linspace(0,total_time,10000)
# print(f(8.))
y = zeros_like(x)
for i in range(len(x)):
    if x[i]<=2*tper:
        y[i]= 0.5*(1-cos(pi*x[i]/2/tper))
    else:
        y[i]=1.


from matplotlib import pyplot



pyplot.plot(x,y)
pyplot.axis('equal')
# pyplot.axis([0, 9, -0.5, 2])
pyplot.grid()
pyplot.show()

