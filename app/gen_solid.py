from dutwav.mesh import Mesh

file1="./INPUT/DATBDMS.txt"
file2="./INPUT/DATWPMS.txt"
path="./out_sla.txt"
surface_type=1


m1 = Mesh()
m2 = Mesh()


m2.read_mesh(file2,surface_type)
nnf=len(m2.nodes)
m2.read_mesh(file1,0)
m2._mark_elems_at_z(-1,'bottom')
# m2._count_elem()
ws = m2._get_waterline_node()
bs = m2._get_btm_node()
sol={}

for i in range(len(m2.elems)):
    nl=m2.elems[i+1][2]
    for j in range(len(nl)):
        if(m2.elems[i+1][0]=='free surface'):
            sol[nl[j]]=0.5
        if(m2.elems[i+1][0]=='bottom'):
            sol[nl[j]]=0.5
        if(m2.elems[i+1][0]=='body'):
            sol[nl[j]]=0.391



for i in range(len(m2.nodes)):
    if(i in ws):
        sol[i] = 0.391/2
    if(i in bs):
        sol[i] = 0.391/2

with open(path,"wb") as f:
    for i in m2.nodes:
        f.write('{:6d}    '.format(i))
        f.write('{0:<7.4f}\n'.format(float(sol[i])))
        # f.write('{:6d} \n'.format(i)
        # f.write('{:6d}    {0:<7.5f}\n'.format(i,sol[i]))

