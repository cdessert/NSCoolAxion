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
    not_possible_age_idxs = np.where((model_evolution[:,0] > age0720_hi)|(model_evolution[:,0] < age0720_lo))
    likelihood_arr = -1./2.*(data[1]-model_evolution[:,1])**2/data_err[1]**2
    likelihood_arr[not_possible_age_idxs] = -1e9
    LL2 = -2*likelihood_arr
    return LL2

age0720_lo = 0.20*1e6
age0720_hi = 1.18*1e6

age_arr = np.logspace(1,8,7001)

# Computes these likelihoods for all NSs
# I do this for 0720 just in case I want to compare later
for iNS,NSname in enumerate(obj_names_short):
    print(NSname)
    # Restrict to the relevant data
    NS_data_T = np.array([Ages[iNS],Tsurfs[iNS]])
    NS_T_err = np.array([Agerrs[iNS],np.sqrt(Terrs[iNS]**2+(0.3*Tsurfs[iNS])**2)])
    NS_data_L = np.array([Ages[iNS],Lums[iNS]])
    NS_L_err = np.array([Agerrs[iNS],Lerrs[iNS]])
    Bfield_NS = Bfield_dict[NSname]
    # Retrieve NSCool data, so this goes to whatever path you put it in
    data = NSC.NScool(NSCool_res)
    age_NSCool = data.Teff_times
    temp_NSCool = data.Teff_Teffs*K2keV
    lum_NSCool = data.Teff_Lphots
    B_NSCool = data.Teff_Bfields
    Tcentral_NSCool = data.Teff_Tcentrals
    last_age = age_NSCool[-1]
    if last_age < 1.0e6:
        sys.exit('NSCool stopped before 10^6 years. Last age: {:.2e} years.'.format(last_age))
    TempCoolingCurve = np.interp(age_arr,age_NSCool,temp_NSCool)
    LumCoolingCurve = np.interp(age_arr,age_NSCool,lum_NSCool)
    # Compute the likelihood
    model_evolution_T = np.vstack((age_arr, TempCoolingCurve)).T
    LL2_evolution_T = likelihood_along_evolution(NS_data_T,NS_T_err,model_evolution_T)
    model_evolution_L = np.vstack((age_arr, LumCoolingCurve)).T
    LL2_evolution_L = likelihood_along_evolution(NS_data_L,NS_L_err,model_evolution_L)
    # Store the minimum along the array
    LL2_min_T = np.min(LL2_evolution_T)
    LL2_min_L = np.min(LL2_evolution_L)
    # And also where that minimum is located
    LL2_argmin_T = np.argmin(LL2_evolution_T)
    LL2_argmin_L = np.argmin(LL2_evolution_L)
    
    np.savez(NSCool_res + 'Results30PctSystematicErrors'+NSname+'_Temp.npz',AgeYr=age_arr,TempkeV=TempCoolingCurve,CoolingLikelihood=LL2_evolution_T,\
             MinLikelihood=LL2_min_T,ArgMinLikelihood=LL2_argmin_T,BGauss=B_NSCool,TcentralkeV=Tcentral_NSCool)
    np.savez(NSCool_res + 'Results'+NSname+'_Lum.npz',AgeYr=age_arr,LumErgs=LumCoolingCurve,CoolingLikelihood=LL2_evolution_L,\
             MinLikelihood=LL2_min_L,ArgMinLikelihood=LL2_argmin_L,BGauss=B_NSCool,TcentralkeV=Tcentral_NSCool,MinLikelihoodT=LL2_min_T)
   

