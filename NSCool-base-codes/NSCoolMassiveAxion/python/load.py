import time, sys, os, h5py
import os.path
from os import path

import numpy as np
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy import interpolate


cid = int(sys.argv[1])


arr = np.logspace(-11,-7,32)
res = []

res.append(0)
for i in range(len(arr)):
    res.append( str(arr[i]))
    res.append( str(-arr[i]))
res = np.array(res)
 
s = res
CSnum = np.array(Cs,dtype='float')
argsort = np.argsort(CSnum)
Cs = Cs[argsort]

EOSs = np.array(['0','1','2','3','4'])
Pairings = np.array(['0','1','2','3','4','5'])
Masses = np.array(['1.0','1.2','1.4','1.6','1.8','2.0'])
logBini = np.array(['11','12','13','14','15','16'])
Names = np.array(['J1856','J1308','TopHatAgeJ0720','J1605'])

LLall = np.zeros((len(Names),len(Cs),len(EOSs),len(Pairings),len(Masses),len(logBini))) + 1e99
Found = np.zeros(len(Cs))

for c in [cid]:
    for e in range(len(EOSs)):
        for p in range(len(Pairings)):
            for m in range(len(Masses)):
                for b in range(len(logBini)):
                    LLsingle = 0
                    for n in range(len(Names)):
                        if n==0:
                            pid = '1048577'
                        if n==1 or n==2:
                            pid = '524289'
                        if n==3:
                            pid = '2097153'
                        file = 'Runs/'+pid+'/'+EOSs[e]+'_'+Pairings[p]+'_'+Masses[m]+'_'+logBini[b]+'/0_0_0_'+Cs[c]+'/'\
                               + 'Results'+Names[n]+'_Lum.npz'
                             # +'Results30PctSystematicErrors'+Names[n]+'.npz'
                        if path.exists(file):
                            data = np.load(file)
                            LLall[n,c,e,p,m,b] = data['MinLikelihood']
                            Found[c] += 1./len(Names)
                        #else:
                        #    print(c,e,p,m,b,file)
    print(Cs[c],Found[c])

np.save('data/cid_'+str(cid)+'.npy',LLall)
