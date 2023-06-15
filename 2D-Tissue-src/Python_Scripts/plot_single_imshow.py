import numpy as np
import numpy as np
import matplotlib.pyplot as plt







file = 'Bin_kmf.0.25.ISO.1.0'
freq = np.loadtxt(file+'.freq.csv',delimiter=',')

freq = np.reshape(freq, (120,125)) # C-like index ordering


plt.figure()
plt.imshow(freq, aspect='equal',origin='lower', interpolation='bicubic',vmin=0, vmax=10,cmap='jet') #cmap=plt.get_cmap('jet'),

plt.show()