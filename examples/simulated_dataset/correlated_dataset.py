import numpy as np
import matplotlib.pyplot as plt
x= np.arange(0.,200.,1.)

y = 50.*np.sin(2*np.pi*(x-2.5245)/14.5678)
z = 3*y - 200.

y_obs = np.random.normal(y, 5)
z_obs = np.random.normal(z, 5)

print(np.average(z_obs))

plt.scatter(z_obs,y_obs

 )
p = np.polyval([0.33291117, 0.000], z_obs+200)
plt.scatter(z_obs, p)
plt.show()


y_obs -= 25.*np.sin(2*np.pi*x/4.5678)


plt.scatter(x,y_obs, c='C0', alpha=0.5)
plt.scatter(x,z_obs, c='C1', alpha=0.5)
plt.show()


fileout_YZ = open('YZ_dataset.dat','w')

for xi, yo, zo in zip(x, y_obs, z_obs):
    fileout_YZ.write('{0:f} {1:f} 5.0 {2:f} 5.00 0 -1 \n'.format(xi,yo,zo))
fileout_YZ.close()


#fileout_Y = open('Y_dataset.dat','w')
#fileout_Z = open('Z_dataset.dat','w')
#
#for xi, yo, zo in zip(x, y_obs, z_obs)
#    fileout_Y.write('{0:f} {1:f} 5.0 0 0 -1 \n'.format(xi,yo))
#    fileout_Z.write('{0:f} {1:f} 5.0 0 0 -1 \n'.format(xi,zo))
#fileout_Y.close()
#fileout_Z.close()

#p = np.polyfit(z_obs+200, y_obs, 1)
#print(p)
