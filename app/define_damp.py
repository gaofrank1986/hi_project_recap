from dutwav.mesh import Mesh

# file1="./INPUT/DATBDMS.txt"
file2="./INPUT/outer_sf.txt"
surface_type=1
maxn=999
prefix="fort.7"

m1 = Mesh()

m1.read_mesh(file2,surface_type)
from dutwav.util import def_damp_func
r = 6.
L = 1.9
alpha = 1.1 
f = def_damp_func(r,alpha,L)

from numpy import *
x = linspace(0,8,1000)
print(f(8.))
y = zeros_like(x)
for i in range(len(x)):
    if x[i]<=8:
        y[i]=f(x[i])

from matplotlib import pyplot

pyplot.plot(x,y)
pyplot.axis('equal')
pyplot.axis([0, 9, -0.5, 2])
pyplot.grid()
pyplot.show()
        
m1._generate_damp_info(f)
assert(len(m1.nodes)==len(m1.damp_info))
m1.tecplt_surface("./damp_info.dat",[m1.damp_info],1)
m1.output_mesh("./damped_surface.txt",3)

# from dutwav.util import create_animation
# create_animation("./",m1,maxn,prefix)
