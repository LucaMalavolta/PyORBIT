from multiprocessing import Pool
import numpy as np 
import time


time_0 = time.time()
def cube(x):
    return np.arange(0, 25, 0.1)**x



pool = Pool(processes=16)
results = pool.map(cube, range(1,20))

pool.close()
print(np.shape(results))


time_1 = time.time()
print("Time taken: ", time_1 - time_0)


# output:
# [1, 8, 27, 64, 125, 216]


x= np.ones([10,3]) * 5.5
y = np.ones([10,3]) * 1.1
k = np.arange(0, 10)
print(x)

for i in zip(x):
    print(i)



def multiple(theta):
    return np.sum(theta)

pool = Pool(processes=16)
results = pool.map(multiple, zip(x))
print(results)