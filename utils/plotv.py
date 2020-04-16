import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('chiS.h5','r')
c = f['chis']
print(c.shape)
plt.pcolormesh((c[:,:,0,0,0,0,0].real))
plt.show()
