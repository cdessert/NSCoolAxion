import numpy as np
from glob import glob
import os
import sys
#sys.path.append('/clusterfs/heptheory/dessert/axion-bfield-production/python/')
import NSCool_process as NSC


# Take in these inputs
# gamm10=sys.argv[1] # Axion coupling
# NSmass=sys.argv[2] # This one specifies the mass and actually also the EOS
# Pairing=sys.argv[3] # Superfluidity model
# Target=sys.argv[5] # Which NS was being modeled (matters because they have different Bfields)
NSCool_res = sys.argv[1]#'/Path/To/NSCool/Sim/Model_1/'

# Temp,Age,Bfield data
# Note 1308 and 0720 have the same Bfield so they have the same cooling at fixed nuisance parameters
K2keV = 8.617e-8
obj_names_short = ['J1856','J1308','J0720','J1605']
Bfields = ["2.9d13","6.8d13","6.8d13","2.0d13"]
Bfield_dict = {o:b for o,b in zip(obj_names_short,Bfields)}
Tsurfs = np.array([50,70,92,78])*1e-3
Terrs = np.array([14,20,10,42])*1e-3
Ages = np.array([0.42,0.55,0.69,0.44])*1e6
Agerrs = np.array([0.08,0.25,0.21,0.07])*1e6

# Computes the likelihood as a function of NS cooling time
def likelihood_along_evolution(data,data_err,model_evolution):
    likelihood_arr = -1./2.*(data-model_evolution)**2/data_err**2
    LL2 = -2*np.sum(likelihood_arr,axis=1)
    return LL2

# Same but for 0720 which has a top hat prior on age
def likelihood_along_evolution_for_0720(data,data_err,model_evolution):
    not_possible_age_idxs = np.where((model_evolution[:,0] > age0720_hi)|(model_evolution[:,0] < age0720_lo))
    likelihood_arr = -1./2.*(data[1]-model_evolution[:,1])**2/data_err[1]**2
    likelihood_arr[not_possible_age_idxs] = -1e9
    LL2 = -2*likelihood_arr
    return LL2

age0720_lo = 0.20*1e6
age0720_hi = 1.18*1e6

# Computes these likelihoods for all NSs
# I do this for 0720 just in case I want to compare later
for iNS,NSname in enumerate(obj_names_short):
    print(NSname)
    # Restrict to the relevant data
    NS_data = np.array([Ages[iNS],Tsurfs[iNS]])
    NS_err = np.array([Agerrs[iNS],np.sqrt(Terrs[iNS]**2+(0.3*Tsurfs[iNS])**2)])
    Bfield_NS = Bfield_dict[NSname]
    # Retrieve NSCool data, so this goes to whatever path you put it in
    data = NSC.NScool(NSCool_res)
    age_NSCool = data.Teff_times
    temp_NSCool = data.Teff_Teffs*K2keV
    # Compute the likelihood
    model_evolution = np.vstack((age_NSCool, temp_NSCool)).T
    LL2_evolution = likelihood_along_evolution(NS_data,NS_err,model_evolution)
    # Store the minimum along the array
    LL2_min = np.min(LL2_evolution)
    # And also where that minimum is located
    LL2_argmin = np.argmin(LL2_evolution)
    
    np.savez(NSCool_res + 'Results30PctSystematicErrors'+NSname+'.npz',AgeYr=age_NSCool,TempkeV=temp_NSCool,CoolingLikelihood=LL2_evolution,MinLikelihood=LL2_min,ArgMinLikelihood=LL2_argmin)
    
# Then just 0720, do the same for top hat prior
print('Tophat')
iNS = 2
NSname = obj_names_short[iNS]
NS_data = np.array([Ages[iNS],Tsurfs[iNS]])
NS_err = np.array([Agerrs[iNS],np.sqrt(Terrs[iNS]**2+(0.3*Tsurfs[iNS])**2)])
Bfield_NS = Bfield_dict[NSname]
data = NSC.NScool(NSCool_res)
age_NSCool = data.Teff_times
temp_NSCool = data.Teff_Teffs*K2keV
model_evolution = np.vstack((age_NSCool, temp_NSCool)).T
LL2_evolution = likelihood_along_evolution_for_0720(NS_data,NS_err,model_evolution)
LL2_min = np.min(LL2_evolution)
LL2_argmin = np.argmin(LL2_evolution)
    
np.savez(NSCool_res + 'Results30PctSystematicErrorsTopHatAge'+NSname+'.npz',AgeYr=age_NSCool,TempkeV=temp_NSCool,CoolingLikelihood=LL2_evolution,MinLikelihood=LL2_min,ArgMinLikelihood=LL2_argmin)

