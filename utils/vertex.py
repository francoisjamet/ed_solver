import h5py
import numpy as np
import sys

nomb = int(sys.argv[1])
beta = float(sys.argv[2])
chiloc = np.loadtxt('chiloc_spin_1_1111')
nomv = int(np.sqrt(chiloc.shape[0]))



g = np.loadtxt('g1.inp').T
g = g[0::2] +1j * g[1::2]

g = g[:,:int(nomv/2)]

norb = int(g.shape[0]**0.5)

gg = np.zeros((norb**2, nomv),dtype = complex)
gg[:,int(nomv/2):]  = g[:]
gg[:,:int(nomv/2)] = np.conjugate(g[:,::-1])
gg = gg.reshape(norb,norb,nomv)
g2 = np.zeros((nomv,nomv),dtype = complex)

chiloc_spin = np.zeros((nomv, nomv, norb,norb,norb,norb, nomb),dtype = complex)
chiloc_charge = np.zeros((nomv, nomv, norb,norb,norb,norb ,nomb),dtype = complex)


from itertools import product

for iom,i1,i2,i3,i4 in product(range(nomb),range(norb),range(norb),range(norb),range(norb)) :
    chi = np.loadtxt('chiloc_spin_{}_{}{}{}{}'.format(iom+1,i1+1,i2+1,i3+1,i4+1)).T
    chi = chi[0] +1j* chi[1]
    k = 0

    for i in range(nomv) :
        for j in range(nomv) :
            chiloc_spin[i,j,i1,i2,i3,i4,iom] = chi[k]/beta
            k+=1
    chi = np.loadtxt('chiloc_charge_{}_{}{}{}{}'.format(iom+1,i1+1,i2+1,i3+1,i4+1)).T
    chi = chi[0] +1j * chi[1]
    k = 0
    for i in range(nomv) :
        for j in range(nomv) :
            chiloc_charge[i,j,i1,i2,i3,i4,iom] = chi[k]/beta
            k+=1
    # if iom == 0 :
    #     print('read green_output_matsu and add for chi_charge')
    #     for v1, v2 in product(range(nomv),range(nomv)) :
    #         g2[v1,v2] = -gg[i1,i4,v1]*gg[i2,i3,v2]*2
    #     chiloc_charge[:,:,i1,i2,i3,i4,iom] += g2

f = h5py.File('chiS.h5', 'w')
f.create_dataset('chis',data=chiloc_spin)
f.create_dataset('chic',data=chiloc_charge)
f.create_dataset('beta',data=[beta])
f.close()
