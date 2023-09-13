from george import kernels, GP
from george.metrics import Metric
import numpy as np

x = np.random.uniform(-1, 1., 20)
y = np.random.uniform(-1, 1., 20)
z = np.random.uniform(-1, 1., 20)

lx = 1.0
ly = 1.0
lz = 1.0



r = np.asarray([x, y, z]).T
print(r.shape)

env = np.ones_like(x)*0.001
res = np.ones_like(x)*0.001 

#k = 2.0 * kernels.Matern32Kernel(metric=[1./np.exp(lx**2), 1./np.exp(ly**2), 1./np.exp(lz**2)], ndim=3)
k = 2.0 * kernels.Matern32Kernel(metric=[lx**2, ly**2, lz**2], ndim=3)
print(k.get_parameter_vector())



lx = 3.0
ly = 4.0
lz = 2.0

l2_lx = np.log(lx**2)
l2_ly = np.log(ly**2)
l2_lz = np.log(lz**2)

l_amp = np.log(2.00)
k.set_parameter_vector([l_amp, l2_lx, l2_ly, l2_lz])


gp = GP(k)
gp.compute(r, env)
mu = gp.sample_conditional(res, r)
mat = gp.get_matrix(r)
print(mat.shape)

import matplotlib.pyplot as plt
plt.matshow(mat)
plt.show()

x1, x2 = np.meshgrid(x,x)
y1, y2 = np.meshgrid(y,y)
z1, z2 = np.meshgrid(z,z)

d_mat = np.sqrt(((x1-x2)/lx)**2 + ((y1-y2)/ly)**2 + ((z1-z2)/lz)**2)
k_mat = 2.0 * (1 + np.sqrt(3)*d_mat)*np.exp(-np.sqrt(3)*d_mat)
plt.matshow(k_mat)
plt.show()

print(mat)

diff = mat-k_mat
print(diff)
plt.matshow(diff)
plt.show()