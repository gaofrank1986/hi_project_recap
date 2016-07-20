from dutwav.mesh import Mesh

file1="./INPUT/DATBDMS.txt"
file2="./INPUT/DATWPMS.txt"
surface_type=1
maxn=999
prefix="fort.7"

m1 = Mesh()
m2 = Mesh()

m1.read_mesh(file2,surface_type)
# m1.read_mesh(file1,0)
# m1.draw_model()

m2.read_mesh(file2,surface_type)
m2.read_mesh(file1,0)
m2.tecplt_nrml("./model_nrml.dat")

from dutwav.util import create_animation
create_animation("./",m1,maxn,prefix)
