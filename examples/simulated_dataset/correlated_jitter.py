import numpy as np
import matplotlib.pyplot as plt
x= np.arange(0.,200.,1.)

y = -25.*np.sin(2*np.pi*x/4.5678)
z = 5.*np.sin(2*np.pi*(x-2.5245)/14.5678)

j = 1. + (z + 5)

y_obs = y*0
for ii in range(0, np.size(y_obs)):
    y_obs[ii] = np.random.normal(y[ii], j[ii])
z_obs = np.random.normal(z, 1.)


plt.scatter(z,j)

plt.show()




plt.scatter(x,y_obs, c='C0', alpha=0.5)
plt.scatter(x,z_obs, c='C1', alpha=0.5)
plt.show()

fileout_Y = open('Y_jitter_dataset.dat','w')
fileout_Z = open('Z_jitter_dataset.dat','w')

for xi, yo, zo in zip(x, y_obs, z_obs):
    fileout_Y.write('{0:f} {1:f} 1.0 0 0 -1 \n'.format(xi,yo))
    fileout_Z.write('{0:f} {1:f} 1.0 0 0 -1 \n'.format(xi,zo))
fileout_Y.close()
fileout_Z.close()

p = np.polyfit(z_obs+200, y_obs, 1)
print(p)
