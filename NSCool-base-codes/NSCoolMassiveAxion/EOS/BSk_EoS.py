import numpy as np
import pandas as pd
import scipy
import shutil
from units import *

u = my_units()

eos_names = ['22', '24', '25', '26']
row_names = ['rho','press','nbar','Ye','Ymu','Yn','Yp','Yla','Ysm','Ys0','Ysp','mstp','mstn','mstla','mstsm','msts0','mstsp']

e_eq_data =  pd.read_csv('./data_fits/e_eq.dat', delimiter=r'\s+', names = eos_names, engine='python') 
P_data =     pd.read_csv('./data_fits/P.dat', delimiter=r'\s+', names = eos_names, engine='python')
Ye_data =    pd.read_csv('./data_fits/Ye.dat', delimiter=r'\s+', names = eos_names, engine='python')
mu_n_data =  pd.read_csv('./data_fits/mu_n.dat', delimiter=r'\s+', names = eos_names, engine='python')
m_eff_data = pd.read_csv('./data_fits/eff_mass_calc.dat', delimiter=r'\s+', names = ['19','20','21','22','23','24','25','26'], engine='python')

# table 14 in Pearson et al 
core_crust_boundary = {}
for eos in eos_names :
    core_crust_boundary[eos] = {}

# neutron density fm^-3
core_crust_boundary['22']['n'] = 0.0716068
core_crust_boundary['24']['n'] = 0.0807555
core_crust_boundary['25']['n'] = 0.0855534
core_crust_boundary['26']['n'] = 0.0849477

# electron abundance
core_crust_boundary['22']['Ye'] = 0.028294
core_crust_boundary['24']['Ye'] = 0.033671
core_crust_boundary['25']['Ye'] = 0.035829
core_crust_boundary['26']['Ye'] = 0.035721

# pressure in MeV fm^-3
core_crust_boundary['22']['P'] = 0.290934
core_crust_boundary['24']['P'] = 0.267902
core_crust_boundary['25']['P'] = 0.210878
core_crust_boundary['26']['P'] = 0.363049



def e_eq(n, eos) :
    ''' Eq. C1 in Pearson et al
    input: 
        - n    average neutron density fm**-3
        - eos  switches between different NS EoS

    output: 
        - e_eq equilibrium energy per baryon
    '''

    if 3. < n < 1E-9 :
        raise ValueError('n-bar outside range of validity')
    if eos not in eos_names : 
        raise ValueError('Unknown EoS')

    p = np.concatenate(([np.nan],np.array(e_eq_data[eos])))

    def w1(n) :
        return 1./(1.+p[9]*n)
    def w2(n) :
        return 1./(1. + (p[13]*n)**p[14])

    res = -9.1536 + (p[1]*n)**(7./6.) /(1. + np.sqrt(p[2]*n) ) * (1. + np.sqrt(p[4]*n)) / ( (1. + np.sqrt(p[3]*n) ) * (1. + np.sqrt(p[5]*n))) * w1(n)
    res += p[6] * n**p[7] * (1. + p[8]*n) * (1-w1(n))*w2(n)
    res += (p[10] * n)**p[11] / (1. + p[12]*n) * (1. - w2(n))

    return res

e_eq = np.vectorize(e_eq)


def rho(n, eos) :
    ''' 
    Based on Eq. 4 of Pearson et al
    '''

    n_cm3 = n / (u.fm**3.) * u.cm**3.
    e_eq_grams = e_eq(n, eos) * u.MeV/u.grams 

    return n_cm3 * ( e_eq_grams + u.m_n/u.grams )

def press(n, eos) :
    ''' 
    Based on Eq. C4 of Pearson et al 
    '''
    xi = np.log10( rho(n,eos) )
    K = 0.         # use this for output in dyn cm^-2
    # K = -33.2047 # use this for output in MeV fm^-3

    p = np.concatenate(([np.nan],np.array(P_data[eos])))

    res = K + (p[1] + p[2] * xi + p[3] * xi**3.)/(1. + p[4]*xi) * 1. / ( np.exp(p[5]*(xi-p[6])) + 1.)
    res += (p[7]  + p[8]*xi)   * 1. / ( np.exp(p[9] *(p[6] -xi)) + 1.)
    res += (p[10] + p[11]*xi)  * 1. / ( np.exp(p[12]*(p[13]-xi)) + 1.)
    res += (p[14]  + p[15]*xi) * 1. / ( np.exp(p[16]*(p[17]-xi)) + 1.)
    res += p[18]/ (1. + (p[20] * (xi - p[19]))**2.)
    res += p[21]/ (1. + (p[23] * (xi - p[22]))**2.)

    return 10**res

def Ye(n, eos) :
    ''' 
    Based on Eq. C17 of Pearson et al 
    '''

    p = np.concatenate(([np.nan],np.array(Ye_data[eos])))

    res = (p[1] + p[2]*n + p[6]*n**(3./2.) + p[3]*n**p[7])/( 1. + p[4]*n**(3./2.) + p[5]*n**p[7] )

    return res

def Ymu(n, eos) :
    ''' 
    Based on Eqs. C14 - 16 of Pearson et al 
    '''

    n_e = Ye(n, eos) * n                                            # in units of fm^-3
    x_e = (3. * np.pi**2. * n_e * u.fm**(-3.))**(1./3.) / u.m_e     # this should be dimensionless

    if 1. + x_e**2. - (u.m_mu/u.m_e)**2. < 0. :
        return 0.

    else :
        n_mu = u.m_e**3. * u.fm**3. / (3.*np.pi**2) * \
        (1. + x_e**2. - (u.m_mu/u.m_e)**2. )**(3./2.)             # in units of fm^-3

        return n_mu/n

Ymu = np.vectorize(Ymu)

def Yp(n, eos) :
    '''
    Simple relationship from beta equilibrium
    '''
    return Ye(n, eos) + Ymu(n, eos)

def Yn(n, eos) :
    '''
    Looking at the NSCool eos Yn + Yp = 1. 
    We therefore also enforce this using this calculated Yp!
    '''
    return 1. - Yp(n,eos)

def mu_n(n, eos) :
    ''' 
    Based on Eqs. C20 of Pearson et al 
    '''

    p = np.concatenate(([np.nan],np.array(mu_n_data[eos])))
    ncc = core_crust_boundary[eos]['n']

    res = p[1]*n**p[2] * (1. + (p[3]*n)**6. )**p[4] / ( (1.+(p[5]*n)**7.)**p[6] * (1. + 1.5 * (n/ncc - 1.))  )

    return res

def mst(n, eos, p_or_n) :
    '''
    Calculate effective mass using using underlying parameters for 
    the Skyrme parametrization of the interactions

    Equations based on Eq. A10 of https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
    Parameters from 
     - Table II of "Further explorations of Skyrme-Hartree-Fock-Bogoliubov mass formulas. XIII. The 2012 atomic
       mass evaluation and the symmetry coefficient"
        S. Goriely, N. Chamel, and J. M. Pearson 
        https://journals.aps.org/prc/pdf/10.1103/PhysRevC.88.024308
     - Table 1 of "Further explorations of Skyrme-Hartree-Fock-Bogoliubov mass
        formulas. XII: Stiffness and stability of neutron-star matter"
        S. Goriely, N. Chamel, and J. M. Pearson
        https://arxiv.org/pdf/1009.3840.pdf
    '''

    t1   = m_eff_data[eos][0]
    t2x2 = m_eff_data[eos][1]
    t4   = m_eff_data[eos][2]
    t5   = m_eff_data[eos][3]

    x1 = m_eff_data[eos][4]
    x4 = m_eff_data[eos][5]
    x5 = m_eff_data[eos][6]

    beta  = m_eff_data[eos][7]
    gamma = 1./12.

    Y_p = Yp(n, eos)
    Y_n = Yn(n, eos)
    # Y_n = 0.8

    Df_n =  t1*((1.+.5*x1)-(.5+x1)*Y_n)
    Df_n += t2x2*(.5+Y_n)
    Df_n += t4*((1.+.5*x4)-(.5+x4)*Y_n)*n**beta
    Df_n += t5*((1.+.5*x5)+(.5+x5)*Y_n)*n**gamma

    Df_p =  t1*((1.+.5*x1)-(.5+x1)*Y_p)
    Df_p += t2x2*(.5+Y_p)
    Df_p += t4*((1.+.5*x4)-(.5+x4)*Y_p)*n**beta
    Df_p += t5*((1.+.5*x5)+(.5+x5)*Y_p)*n**gamma

    # print(Df)

    # Df has dimensions MeV fm^5
    # n  has dimesnions fm^-3
    # print(0.5 * (Df*u.MeV*u.fm**5.) * (n/u.fm**3.) * u.m_n )
    res_n = 1. / (1. + 0.5 * (Df_n*u.MeV*u.fm**5.) * (n/(u.fm**3.)) * u.m_n)
    res_p = 1. / (1. + 0.5 * (Df_p*u.MeV*u.fm**5.) * (n/(u.fm**3.)) * u.m_p)

    # The code below was for a comparison to 
    # "Thermal properties of supernova matter: The bulk homogeneous phase"
    # Constantinou, Muccioli, Prakash and Lattimer
    # https://journals.aps.org/prc/pdf/10.1103/PhysRevC.89.065802
    
    # t1 =  570.88
    # t2 =  -67.7 

    # x1 = 0.
    # x2 = 0.


    # Y_n = 0.1

    # Df =  t1*((1.+.5*x1)-(.5+x1)*Y_n)
    # Df += t2*((1+.5+Y_n) + x2*(0.5 + Y_n))
    # # Df += t4*((1.+.5*x4)-(.5+x4)*Y_n)*n**beta
    # # Df += t5*((1.+.5*x5)+(.5+x5)*Y_n)*n**gamma

    # # Df has dimensions MeV fm^5
    # # n  has dimesnions fm^-3
    # # print(0.5 * (Df*u.MeV*u.fm**5.) * (n/u.fm**3.) * u.m_n )
    # res = 1. / (1. + 0.5 * (Df*u.MeV*u.fm**5.) * (n/(u.fm**3.)) * u.m_n)

    # return res

    if p_or_n == 'p' :
        return res_p
    if p_or_n == 'n' :
        return res_n


mst = np.vectorize(mst)
Yn = np.vectorize(Yn)
Yp = np.vectorize(Yp)
Ye = np.vectorize(Ye)
Ymu = np.vectorize(Ymu)

def output_eos(eos, nbar_max, steps) :

    core_boundary_Ye = 3.1606E-02
    nbar_cc_boundry = scipy.optimize.fsolve(
        lambda x: Yp(x, eos) - core_boundary_Ye, [.1]
        )[0]
    print('nbar at crust-core boundary =', nbar_cc_boundry)
    lg10_nbar_max        = np.log10(nbar_max)
    lg10_nbar_cc_boundry = np.log10(nbar_cc_boundry)
    nbar_vals            = np.logspace(lg10_nbar_max, lg10_nbar_cc_boundry, steps)

    rho_vals   = rho(nbar_vals, eos)
    press_vals = press(nbar_vals, eos)

    Ye_vals    =  Ye(nbar_vals, eos)
    Ymu_vals   = Ymu(nbar_vals, eos)
    Yn_vals    =  Yn(nbar_vals, eos)
    Yp_vals    =  Yp(nbar_vals, eos)

    Yla_vals   = np.zeros(steps)
    Ysm_vals   = np.zeros(steps) 
    Ys0_vals   = np.zeros(steps) 
    Ysp_vals   = np.zeros(steps) 

    mstp_vals  = mst(nbar_vals, eos, 'p')
    mstn_vals  = mst(nbar_vals, eos, 'n')
    mstla_vals = np.zeros(steps)
    mstsm_vals = np.zeros(steps)
    msts0_vals = np.zeros(steps)
    mstsp_vals = np.zeros(steps)

    full_array = np.transpose(np.array([
        rho_vals, press_vals, nbar_vals, 
        Ye_vals, Ymu_vals, Yn_vals, Yp_vals,
        Yla_vals, Ysm_vals, Ys0_vals,
        Ysp_vals, mstp_vals, mstn_vals,
        mstla_vals, mstsm_vals, msts0_vals,
        mstsp_vals
        ]))

    pd_df = pd.DataFrame(full_array,columns=row_names)

    filename = '/Users/tobyopferkuch/BSk'+eos+'_test.dat'
    og_eos = '/Users/tobyopferkuch/Dropbox/CurrentProjects/MuonicNSCooling/ns-muons/to/packages/nscool/EOS/APR_EOS_Cat_backup.dat'

    pd_df.to_csv(filename, 
        sep=' ', index=False, header=False, float_format='%8.6E')

    with open(filename, 'r+') as f:
        readcontent = f.read()  
        f.seek(0, 0)
        f.write('6    '+str(steps+1+62)+'    '+str(steps+1)+'          Itext Imax Icore \n')
        f.write('\n')
        f.write('      BSk'+str(eos)+' with crust from APR \n')
        f.write('\n')
        f.write('   Rho        Press       nbar       Ye         Ymu        Yn         Yp         Yla        Ysm        Ys0        Ysp       mstp       mstn       mstla      mstsm      msts0      mstsp \n')
        f.write('  g/cm3     dyne/cm2     #/fm3     [A_cell]    [A_ion]     [Z] \n')
        f.write('\n')
        f.write(readcontent) 

        with open(og_eos, 'r+') as orig_EOS:
            for line_no, line in enumerate(orig_EOS):
                if line_no > 185 :
                    f.write(line)


    dest_eos = '/Users/tobyopferkuch/Dropbox/CurrentProjects/MuonicNSCooling/ns-muons/to/packages/nscool/EOS'
    shutil.copy(filename, dest_eos+'/BSk'+eos+'_test.dat')

    return 'printing complete'
    



    














