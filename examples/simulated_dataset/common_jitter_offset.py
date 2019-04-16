import numpy as np
import matplotlib.pyplot as plt
x1 = np.arange(0.,100.,1.)
x2 = np.arange(50.25,150.25,1.)
y = 1000 + -25.*np.sin(2*np.pi*x1/4.5678)
z = 1000 + -25.*np.sin(2*np.pi*x2/4.5678)


z_obs = np.random.normal(z, 5.)
y_obs = np.random.normal(y, 5.)

fileout_Y = open('Y_common_jitter_offset.dat','w')
fileout_Z = open('Z_common_jitter_offset.dat','w')

for xi, yo in zip(x1, y_obs):
    fileout_Y.write('{0:f} {1:f} 1.0 0 0 -1 \n'.format(xi,yo))

for xi, zo in zip(x2, z_obs):
    fileout_Z.write('{0:f} {1:f} 1.0 0 0 -1 \n'.format(xi,zo))

fileout_Y.close()
fileout_Z.close()

p = np.polyfit(z_obs+200, y_obs, 1)
print(p)
