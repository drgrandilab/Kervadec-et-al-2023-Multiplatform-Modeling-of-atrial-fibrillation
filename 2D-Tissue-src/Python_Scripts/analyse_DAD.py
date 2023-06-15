import numpy as np
import numpy as np
import matplotlib.pyplot as plt



#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy as np
import matplotlib.pyplot as plt

	
import scipy
	
from scipy import signal

def process(time, newAP,cellid=0, plot_fig=False):

	min_x = min(newAP[:,cellid])
	

	min_x = min(newAP[2000:,cellid]);
	# print(min_x)
	
	vec = newAP[2000:,cellid];
	t = np.array(time[2000:])  # make sure this excludes the first paced beat.
	threshold = min(vec)+4
	index,pks=signal.find_peaks(vec, height=threshold,distance=100)

	Amp = max(vec) - min(vec);
	# print(index)
	
	# print(t[index])

	if plot_fig:
		plt.plot(time, newAP[:,cellid])
		plt.show()
		plt.plot(time[:], newAP[:,cellid])
		plt.plot(t[index],pks['peak_heights'],'x')
	
		plt.show()
	
		print(len(pks['peak_heights']))
	return len(pks['peak_heights']),Amp




folder=['Bin_kmf.1.0.ISO.0.0','Bin_kmf.1.0.ISO.1.0','Bin_kmf.0.25.ISO.0.0','Bin_kmf.0.25.ISO.1.0']


for file in folder:
	time = range(10000,20000,1)
	AP = []
	newAP=[]
	for i in time:
		data = np.fromfile(file+'/'+'v_%d.bin'%i, dtype=np.float32);
		# print(data[0])
		AP.append(data)
	
	
	# In[3]:
	
	
	# len(AP)
	
	
	# In[7]:
	
	
	newAP= np.array(AP)
	
	
	# In[23]:
	
	
	# len(newAP[:,0])
	# cellid=5000
	
	freq = []
	amplitude = []
	for i in range(15000):
		d,amp = process(time, newAP,cellid=i, plot_fig=False);
		# d=i
		freq.append(d)
		amplitude.append(amp)
	
	# print(time[2000])
	# freq = newAP[2000,:];
	print (max(freq),max(amplitude))

	freq = np.array(freq);
	freq = np.reshape(freq, (120,125)) # C-like index ordering
	amplitude = np.array(amplitude);
	amplitude = np.reshape(amplitude, (120,125)) # C-like index ordering


	plt.figure()
	
	# plt.pcolor(freq)
	plt.imshow(freq, aspect='equal',origin='lower', interpolation='bicubic',vmin=0, vmax=10,cmap='jet') #cmap=plt.get_cmap('jet'),
	
	# plt.axes().set_aspect('equal')
	plt.colorbar()
	plt.xticks([])  # Disable xticks.
	plt.yticks([])  # Disable xticks.
	plt.savefig(file+'.freq'+'.pdf',dpi=400)

	plt.figure()
	
	# plt.pcolor(freq)
	plt.imshow(amplitude, aspect='equal',origin='lower', interpolation='bicubic',vmin=0, vmax=120,cmap='jet') #cmap=plt.get_cmap('jet'),
	plt.xticks([])  # Disable xticks.
	plt.yticks([])  # Disable xticks.
	# plt.axes().set_aspect('equal')
	plt.colorbar()
	plt.savefig(file+'.amp'+'.pdf',dpi=400)


	np.savetxt(file+'.freq.csv', freq, delimiter=",")
	

	# plt.show()
# In[26]:





# In[27]:


# print(min_x)


# In[57]:

