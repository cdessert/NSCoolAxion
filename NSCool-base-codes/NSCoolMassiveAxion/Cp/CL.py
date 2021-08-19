import time, sys, os, h5py
import os.path
from os import path
import numpy as np
from scipy import stats
from scipy import interpolate
from multiprocessing import Pool

C_p = 1
C_n = 0

m_p = 0.938272088 # GeV
m_n = 0.939565413

#f_a = np.logspace(7,12,64)
f_a = np.logspace(5.8,10.8,64)

fs = []
for i in range(len(f_a)):
    fs.append(f_a[i])
#for i in range(len(f_a)):
#    fs.append(-f_a[len(f_a)-i-1])
fs = np.array(fs) 
ma = 5.691e9 / fs * 1e-3 # eV

gapp = C_p * m_p / f_a
gann = C_n * m_n / f_a

EOSs = np.array(['0','1','2','3','4'])
Pairings = np.array(['0'])
Masses = np.array(['1.0','1.2','1.4','1.6','1.8','2.0'])
logBini = np.array(['11','12','13','14','15','16'])
#Names = np.array(['J1856','J1308','TopHatAgeJ0720','J1605'])#,'J0720'])
Names = np.array(['J1856','J1308','J0720','J1605'])
logDeltaM = np.array(['-20','-19','-18','-17','-16','-15','-14','-13','-12','-11','-10','-9','-8','-7','-6'])
pid = np.array(['65504','4259808'])


Coup = sys.argv[1]
ID = sys.argv[2]
low = int(sys.argv[3])
high = int(sys.argv[4])

NewAges,NewLumi = np.load('SampData/'+Coup+'.npy')

NewAges = NewAges[:,low:high]
NewLumi = NewLumi[:,low:high]

LumiCurves_All = np.zeros(( len(gann),2,len(pid),len(EOSs),len(Pairings),len(Masses),len(logBini),len(logDeltaM),1000 ))

def load(c):
#for c in range(len(gann)):
    fname = '../Runs/Cp/run_'+str(gann[c])+'_'+str(gapp[c])+'_0_0_0.h5'
    print('loading curves: ',c,fname)
    sys.stdout.flush()
    f = h5py.File(fname, 'r')
    
    for s in range(len(pid)):
        for e in range(len(EOSs)):
            for p in range(len(Pairings)):
                for m in range(len(Masses)):
                    for b in range(len(logBini)):
                        for d in range(len(logDeltaM)):
                            key = pid[s]+'_'+EOSs[e]+'_'+Pairings[p]+'_'\
                                +Masses[m]+'_'+logBini[b]+'_'+logDeltaM[d]+'_Teff'
                            if key in list(f.keys()):
                                rawdata = f[key][:]
                                length = np.shape(rawdata)[0]
                                LumiCurves_All[c,0,s,e,p,m,b,d,:length] = f[key][:,1]
                                LumiCurves_All[c,1,s,e,p,m,b,d,:length] = f[key][:,3]
    f.close()
    return LumiCurves_All[c]

print('done loading!\n')
sys.stdout.flush()

start  = time.time()
with Pool(28) as p:
    LumiCurves_All = p.map(load, range(len(gann)))
end = time.time()
print('loading time: ',end-start)



print('average: ',np.mean(LumiCurves_All))

cgs2Lsol = 2.599e-34 # 1 erg/s = 2.599 10^-34 solar luminosities
K2keV = 8.617e-8 # 1 Kelvin = 8.617 10^-8 keV
obj_names_short = ['J1856','J1308','J0720','J1605']
Bfields = ["2.9d13","6.8d13","6.8d13","2.0d13"]
Bfield_dict = {o:b for o,b in zip(obj_names_short,Bfields)}
Tsurfs = np.array([50,70,92,78])*1e-3
Terrs = np.array([14,20,10,42])*1e-3
Lums = np.array([0.065,0.33,0.19,0.2535])*1e33 # erg/s
Lerrs = np.array([0.015,0.06,0.095,0.2465])*1e33 # erg/s
Ages = np.array([0.42,0.55,0.69,0.44])*1e6
Agerrs = np.array([0.08,0.25,0.21,0.07])*1e6

# Computes the likelihood as a function of NS cooling time
def likelihood_along_evolution(data,data_err,model_evolution):
    likelihood_arr = -1./2.*(data-model_evolution)**2/data_err**2
    LL2 = -2*np.sum(likelihood_arr,axis=1)
    return LL2

# Same but for 0720 which has a top hat prior on age
def likelihood_along_evolution_for_0720(data,data_err,model_evolution):
    #not_possible_age_idxs = np.where((model_evolution[:,0] > age0720_hi)|(model_evolution[:,0] < age0720_lo))
    likelihood_arr = -1./2.*(data[1]-model_evolution[:,1])**2/data_err[1]**2
    #likelihood_arr[not_possible_age_idxs] = -1e9
    LL2 = -2*likelihood_arr
    return LL2

age0720_lo = 0.20*1e6
age0720_hi = 1.18*1e6

#age_arr = np.logspace(1,8,7001)

SampBFc = np.zeros((len(NewAges[0]),3)) - 1

def fit(i):
    LL_s_all = np.zeros((len(gann),len(pid),len(EOSs),len(Pairings),len(logBini),len(Masses),len(logDeltaM),len(Names)))\
         + 1e99

    N = 512
    age_arr0= np.logspace(np.log10(np.max([1,NewAges[0][i] - 4*Agerrs[0]])),
                          np.log10(NewAges[0][i] + 4*Agerrs[0])
                          ,512)
     
    age_arr1= np.logspace(np.log10(np.max([1,NewAges[1][i] - 4*Agerrs[1]])),
                          np.log10(NewAges[1][i] + 4*Agerrs[1])
                          ,512)
     
    age_arr2= np.logspace(np.log10(np.max([1,NewAges[2][i] - 4*Agerrs[2]])),
                          np.log10(NewAges[2][i] + 4*Agerrs[2])
                          ,512)

    age_arr3= np.logspace(np.log10(np.max([1,NewAges[3][i] - 4*Agerrs[3]])),
                          np.log10(NewAges[3][i] + 4*Agerrs[3])
                          ,512)
    
    age_arr = np.array([age_arr0,age_arr1,age_arr2,age_arr3]) 
   
    for c in range(len(gann)):
        print(i,c)
        sys.stdout.flush()
        for s in range(len(pid)):
            for e in range(len(EOSs)):
                for p in range(len(Pairings)):
                    for m in range(len(Masses)):
                        for b in range(len(logBini)):
                            for d in range(len(logDeltaM)):
                                samp_age_NSCool = LumiCurves_All[c][0,s,e,p,m,b,d,:]
                                samp_lum_NSCool = LumiCurves_All[c][1,s,e,p,m,b,d,:]
                                locs = np.where(samp_age_NSCool > 1e4)
                                   
                                if len(locs[0])==0:
                                    continue
 
                                for iNS,NSname in enumerate(Names):
                                    # Restrict to the relevant data
                                    NS_data_L = np.array([NewAges[iNS][i],NewLumi[iNS][i]])
                                    NS_L_err = np.array([Agerrs[iNS],Lerrs[iNS]])

                                    LumCoolingCurve = np.interp(age_arr[iNS],samp_age_NSCool[locs],samp_lum_NSCool[locs])

                                    # Compute the likelihood
                                    model_evolution_L = np.vstack((age_arr[iNS], LumCoolingCurve)).T
                                    LL2_evolution_L = likelihood_along_evolution(NS_data_L,NS_L_err,
                                                                                 model_evolution_L)
                                    # Store the minimum along the array
                                    LL2_min_L = np.min(LL2_evolution_L)
                                    LL_s_all[c,s,e,p,b,m,d,iNS] = LL2_min_L
                  
    LLjoint = np.zeros((len(gann),len(pid),len(EOSs),len(Pairings)))
    for c in range(len(gann)):
        for s in range(len(pid)):
            for e in range(len(EOSs)):
                for p in range(len(Pairings)):
                    for n in range(len(Names)):
                        dat = LL_s_all[c,s,e,p,:,:,:,n]
                        LLjoint[c,s,e,p] += np.min(dat)

    LL = np.zeros((len(gann),len(pid)))
    for c in range(len(gann)):
        for s in range(len(pid)):
            dat = LLjoint[c,s,:,:]
            LL[c,s] = np.min(dat)

    loc = np.where(LL == np.min(LL))
    BFc = loc[0][0]
    BFs = loc[1][0]
    SampBFc[i,0] = BFc
    SampBFc[i,1] = BFs
    SampBFc[i,2] = np.min(LL)

    return SampBFc[i]
    
with Pool(28) as p:
    SampBFc = p.map(fit, range(len(NewAges[0])))

np.save('CLs/BFs_'+Coup+'_'+ID+'.npy',SampBFc)


