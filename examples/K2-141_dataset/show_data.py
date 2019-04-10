import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt('K2-141_detrended_bonly.csv', delimiter=',', skip_header=1, dtype=np.double)

planetb_period = 0.2803234429
sel = (data[:,0] < 0.00000)
TIME = data[:,0]+planetb_period*sel
PHA = TIME/planetb_period
VAL = data[:,2]

sel = (PHA> 0.2) & (PHA<0.8)
K2_error_estimate = np.std(VAL[sel])  #

K2_error_estimate = 0.00002 # Halved and rounded in order to allow for variation of the jitter parameter

#create file with transit of only planet b
TIME = data[:,1] + 4833.00
fileout = open('K2-141_detrended_bonly.dat','w')
for t,f in zip(TIME,VAL):
    fileout.write('{0:f} {1:f} {2:f} 0 -1 -1 \n'.format(t,f,K2_error_estimate))
fileout.close()




data = np.genfromtxt('K2-141_raw_noplanets.dat')
TIME = data[:,0]
RAW = data[:,1]

size_full = np.size(data[:,0])
size_tobereshaped = size_full - size_full % 5

BJD = data[:size_tobereshaped,0]
VAL = data[:size_tobereshaped,1]

BJD_avg = np.mean(BJD.reshape(-1, 5), axis=1)
VAL_avg = np.mean(VAL.reshape(-1, 5), axis=1)
ERR_avg = np.std(VAL.reshape(-1, 5), axis=1)

fileout = open('K2-141_raw_noplanets_5sampled.dat','w')
for t,f,e in zip(BJD_avg,VAL_avg,ERR_avg):
    fileout.write('{0:f} {1:f} {2:f} 0 0 -1 \n'.format(t,f,e))
fileout.close()

size_full = np.size(data[:,0])
size_tobereshaped = size_full - size_full % 10

BJD = data[:size_tobereshaped,0]
VAL = data[:size_tobereshaped,1]

BJD_avg = np.mean(BJD.reshape(-1, 10), axis=1)
VAL_avg = np.mean(VAL.reshape(-1, 10), axis=1)
ERR_avg = np.std(VAL.reshape(-1, 10), axis=1)

fileout = open('K2-141_raw_noplanets_10sampled.dat','w')
for t,f,e in zip(BJD_avg,VAL_avg,ERR_avg):
    fileout.write('{0:f} {1:f} {2:f} 0 0 -1 \n'.format(t,f,e))
fileout.close()

#create rebinned file

data = np.genfromtxt('K2-141_raw_detrended.csv', delimiter=',', skip_header=1, dtype=np.double)
TIME = data[:,0] + 4833.00000
RAW = data[:,1]
DET = data[:,2]
#create file with raw and normalized curve
fileout = open('K2-141_raw.dat','w')
for t,f in zip(TIME,RAW):
    fileout.write('{0:f} {1:f} {2:f} 0 -1 -1 \n'.format(t,f,K2_error_estimate))
fileout.close()

fileout = open('K2-141_detrended.dat','w')
for t,f in zip(TIME,DET):
    fileout.write('{0:f} {1:f} {2:f} 0 -1 -1 \n'.format(t,f,K2_error_estimate))
fileout.close()
