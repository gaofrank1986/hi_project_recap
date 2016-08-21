from dutwav.mesh import Mesh

file1="./INPUT/DATBDMS.txt"
file2="./INPUT/DATWPMS.txt"
surface_type=1
maxn=999
prefix="fort.9"

m1 = Mesh()
m2 = Mesh()

m1.read_mesh(file2,surface_type)
# m1.read_mesh(file1,0)
# m1.draw_model()

m2.read_mesh(file2,surface_type)
m2.read_mesh(file1,0)
# m2.tecplt_nrml("./model_nrml.dat")

# from dutwav.util import create_animation
# create_animation("./",m1,maxn,prefix)
from dutwav.util import display_value
display_value("./fort.9999","./solid_angle1.dat",m2)
# display_value("./test_solid_angle.txt","./solid_angle.dat",m2)
# display_value("./out_sla.txt","./solid_angle.dat",m2)
