import sys
from itertools import product
nv =  int(sys.argv[1])
f = open('omega_list_path','w')
for v1,v2 in product(range(-nv,nv),range(-nv,nv)) :
    f.write('{} {} \n'.format(2*v1+1,2*v2+1))
f.close()
