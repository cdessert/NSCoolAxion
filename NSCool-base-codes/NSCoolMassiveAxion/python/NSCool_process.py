import numpy as np
from scipy.interpolate import interp1d

import copy

Ifm_to_GeV = 1./5.068 #GeV
Icm_to_GeV = 1.98e-14 #Icm to GeV
g_to_GeV = 5.62e23 #gram to GeV
rho_s = 2.8e14*g_to_GeV*Icm_to_GeV**3 #GeV^4
m_n = 0.939 #GeV
m_pi = 0.140# GeV
m_e = 0.511e-3 #GeV
m_mu = 105.7e-3 #GeV

Kelvin_to_keV = 8.617e-8

def fexp(x):
	return np.exp(x)#np.exp(max(x,-700))

def return_Rs(T,Tc):
	'''
	Returns superfluid suppression factors Rpp, Rnp
	'''
	Rpp = np.ones(len(Tc))
	Rnp = np.ones(len(Tc))
	whs = np.where(T<Tc)[0]
	if len(whs)>0:
		tau = T/Tc
		v = np.sqrt(1-tau[whs])*(1.456-0.157/np.sqrt(tau[whs])+1.764/tau[whs])
		app = 0.1747 + np.sqrt(0.8253**2+0.07933**2*v**2)
		bpp = 0.7333 + np.sqrt(0.2667**2+0.1678**2*v**2)
		anp = 0.9982 + np.sqrt(0.0018**2+0.3815**2*v**2)
		bnp = 0.3949 + np.sqrt(0.6051**2+0.2666**2*v**2)

		Rpp[whs] = 1/2.*(app**2*fexp(4.228-np.sqrt(4.228**2+4*v**2)) + bpp**7.5*fexp(7.762-np.sqrt(7.762**2+9*v**2)) )

		Rnp[whs] = 1/2.732*(anp*np.exp(1.306-np.sqrt(1.306**2+v**2)) + 1.732*bnp**7*np.exp(3.303-np.sqrt(3.303**2+4*v**2)))

	return Rpp, Rnp


def F(x):
	return 1 - 3/2.*x*np.arctan(1/x) + x**2/2./(1.+x**2)

class NScool:
	def __init__(self,mod_folder="/nfs/turbo/bsafdi/bsafdi/github/white-dwarf-axion/local/NSCool/Model_1/"):
		self._mod_folder = mod_folder
		self._Teff_file = mod_folder + "Teff_Try.dat"
		self._star_file = mod_folder + "Star_Try.dat"
		self._Temp_file = mod_folder + "Temp_Try.dat"
		self._EOS_file = mod_folder + "../EOS/APR_EOS_Cat.dat"
		
		self._read_files()
# 		self._process_star_crust()
# 		self._process_star_core()
# 		self._process_EOS()
# 		self._compute_electron_star_factor()
# 		self._compute_nucelon_star_factors()

		self._process_Temp()
		
	def _read_files(self):
		Teff_open = open(self._Teff_file, "r") 
		star_open = open(self._star_file, "r") 
		Temp_open = open(self._Temp_file, "r") 
#		EOS_open = open(self._EOS_file, "r") 

		# Teff
		it = 0
		lines_Teff = []
		for line in Teff_open:
			lines_Teff += [line]
			
		self.lines_Teff = np.array(lines_Teff[26:])
		self.Teff_times = np.zeros(len(self.lines_Teff))
		self.Teff_Teffs = np.zeros(len(self.lines_Teff))
		for i in range(len(self.lines_Teff)):
			lt = self.lines_Teff[i].split()
			self.Teff_times[i] = float(lt[1])
			self.Teff_Teffs[i] = float(lt[2])

		# star
		it = 0
		lines_star = []
		for line in star_open:
			lines_star += [line]

		self.lines_star = np.array(lines_star[6:])

		# Temp
		it = 0
		lines_Temp = []
		for line in Temp_open:
			lines_Temp += [line]

		self.lines_Temp = np.array(lines_Temp[4:])
		
#		it = 0
#		lines_EOS = []
#		for line in EOS_open:
#			lines_EOS += [line]

#		self.lines_EOS = np.array(lines_EOS[7:186])       
		

	def _process_star_crust(self):
		len_arr = np.zeros(len(self.lines_star))
		for i in range(len(self.lines_star)):
			lt = self.lines_star[i]
			len_arr[i] = len(lt.split())
		iout_s = np.where(len_arr==13)[0][0]
		iout_e = np.where(len_arr==13)[0][-1]

		self.line_inner_crust = self.lines_star[iout_s:iout_e]


		lines = self.line_inner_crust
		Ns = len(lines)
		self.star_rad = np.zeros(Ns)
		self.star_emas = np.zeros(Ns)
		self.star_rho = np.zeros(Ns)
		self.star_pres = np.zeros(Ns)
		self.star_nb = np.zeros(Ns)
		self.star_kfe = np.zeros(Ns)
		self.star_kfn = np.zeros(Ns)
		self.star_Tcn = np.zeros(Ns)
		self.star_what = np.empty( Ns,dtype="S10")
		self.star_acell = np.zeros(Ns)
		self.star_aion = np.zeros(Ns)
		self.star_zion = np.zeros(Ns)
		# i,rad(i)/100.,emas(i),rrho(i),pres(i),bar(i),kfe(i),kfn(i),tcn(i),what,a_cell(i),a_ion(i),z_ion(i)
		beg = 1
		for i in range(Ns):
			lt = lines[i]
			ar = lt.split()
			if len(ar) != 13:
				print('len(ar) != 13')
			self.star_rad[i] = float(ar[beg+0])
		#     print star_rad[i]
			self.star_emas[i] = float(ar[beg+1])
			self.star_rho[i] = float(ar[beg+2])
			self.star_pres[i] = float(ar[beg+3])
			self.star_nb[i] = float(ar[beg+4])*Ifm_to_GeV**3 #GeV^3
			self.star_kfe[i] = float(ar[beg+5])*Ifm_to_GeV
			self.star_kfn[i] = float(ar[beg+6])*Ifm_to_GeV
			self.star_Tcn[i] = float(ar[beg+7])
			self.star_what[i] = str(ar[beg+8])
			#print str(ar[beg+8])
			self.star_acell[i] = float(ar[beg+9])
			self.star_aion[i] = float(ar[beg+10])
			self.star_zion[i] = float(ar[beg+11])

	def _process_star_core(self):
		len_arr = np.zeros(len(self.lines_star))
		for i in range(len(self.lines_star)):
			lt = self.lines_star[i]
			len_arr[i] = len(lt.split())
		iout_s = np.where(len_arr==22)[0][0]
		iout_e = np.where(len_arr==22)[0][-1]

		self.line_core = self.lines_star[iout_s:iout_e]


		lines = self.line_core
		Ns = len(lines)
		self.star_rad_core = np.zeros(Ns)
		self.star_emas_core = np.zeros(Ns)
		self.star_rho_core = np.zeros(Ns)
		self.star_pres_core = np.zeros(Ns)
		self.star_nb_core = np.zeros(Ns)
		self.star_kfe_core = np.zeros(Ns)
		self.star_kfmu_core = np.zeros(Ns)
		self.star_kfp_core = np.zeros(Ns)
		self.star_kfn_core = np.zeros(Ns)
		self.star_Tcn_core = np.zeros(Ns)
		self.star_Tcp_core = np.zeros(Ns)
		self.star_nsf_core = np.empty( Ns,dtype="S10")

		# i,rad(i)/100.,emas(i),rrho(i),pres(i),bar(i),kfe(i),kfn(i),tcn(i),what,a_cell(i),a_ion(i),z_ion(i)
		beg = 1
		for i in range(Ns):
			lt = lines[i]
			ar = lt.split()
			if len(ar) != 22:
				print(len(ar) != 22)
			self.star_rad_core[i] = float(ar[beg+0])
		#     print star_rad[i]
			self.star_emas_core[i] = float(ar[beg+1])
			self.star_rho_core[i] = float(ar[beg+2])
			self.star_pres_core[i] = float(ar[beg+3])
			self.star_nb_core[i] = float(ar[beg+4])*Ifm_to_GeV**3 #GeV^3
			self.star_kfe_core[i] = float(ar[beg+5])*Ifm_to_GeV
			self.star_kfmu_core[i] = float(ar[beg+6])*Ifm_to_GeV
			self.star_kfp_core[i] = float(ar[beg+7])*Ifm_to_GeV
			self.star_kfn_core[i] = float(ar[beg+8])*Ifm_to_GeV
			self.star_Tcn_core[i] = float(ar[beg+13])
			self.star_Tcp_core[i] = float(ar[beg+14])
			self.star_nsf_core[i] = str(ar[beg+17])

	def _process_EOS(self):
		self.EOS_rho  = np.zeros(len(self.lines_EOS))        
		self.EOS_mstn = np.zeros(len(self.lines_EOS))
		self.EOS_mstp = np.zeros(len(self.lines_EOS))
		for i in range(len(self.lines_EOS)):
			lt = self.lines_EOS[i]
			self.EOS_rho[i]  = float(lt.split()[0])
			self.EOS_mstp[i] = float(lt.split()[11])
			self.EOS_mstn[i] = float(lt.split()[12])
						
			
#	def return_mstn(self,rho):            

	def _compute_electron_star_factor(self):
		self._electron_star_factor = (self.star_nb/10./rho_s/m_n)*(0.5e-3/self.star_kfe/1e-2)**2*self.star_zion**2/self.star_aion #fix Z and A here
		self._neutron_crust_star_factor = (self.star_kfn/0.337)**6

	def _compute_nucelon_star_factors(self):
		x = m_pi / (2*self.star_kfn_core)
		y = m_pi / (2*self.star_kfp_core)
		#self.Fx_test = F(x)
		#self.Fy_test = F(y)
		Fx = F(x)
		Fy = F(y)
		self.Fxyp  = F(2*x*y/(x+y))
		self.Fxym  = F(2*x*y/(y-x))
		self._eann_star_factor = (self.star_kfn_core/0.337)*Fx/0.607211 #all in GeV
		self._eapp_star_factor = (self.star_kfp_core/0.337)*Fy/0.607211

 ##
		Fxyp = self.Fxyp
		Fxym = self.Fxym #

		Gg = 1/2.*Fy+( (Fxyp + Fxym) +y/x*(Fxyp-Fxym) )+(1-y*np.arctan(1/y))
		Gh = 1/2.*Fy+1/2.*( (Fxyp + Fxym) +y/x*(Fxyp-Fxym) )+(1-y*np.arctan(1/y))

		self._eanp_star_factor_g = (self.star_kfp_core/0.337)*Gg
		self._eanp_star_factor_h = (self.star_kfp_core/0.337)*Gh

		self._PBF_s_n_star_factor = (self.star_kfn/0.337)**3 
		self._PBF_s_p_star_factor = (self.star_kfp_core/0.337)**3 
		self._PBF_p_star_factor = (self.star_kfn_core/0.337)

		gamma_e = self.star_kfe_core/m_e
		gamma_mu = self.star_kfmu_core/m_mu
		self._pe_star_factor = (self.star_kfp_core/0.337)**3*(100*m_e/self.star_kfe_core)**2*(2*np.log(2*gamma_e)+np.log(1/137./np.pi))/4.53
		self._pmu_star_factor = (self.star_kfp_core/0.337)**3*(100*m_e/np.sqrt(self.star_kfmu_core**2+m_mu**2))**2*(2*np.log(2*gamma_mu)+np.log(1/137./np.pi))/4.53

		whs = np.where(self.star_kfmu_core<=1e-10)[0]
		self._pmu_star_factor[whs] = 0.0

	def _process_Temp(self):
		len_arr = np.zeros(len(self.lines_Temp))
		Time_arr = []
		Te_inf_arr = []

		zone_arr = []
		Rad_arr = []
		Rho_arr = []
		ephi_arr = []
		dvol_arr = []
		Temp_arr = []

		zt = []
		rt = []
		rhot = []
		ephit = []
		dvolt = []
		Tt = []

		in_sesh=False
		go=False
		start=False
		for i in range(len(self.lines_Temp)):
			lt = self.lines_Temp[i].split()
			len_arr[i] = len(self.lines_Temp[i].split())
			if len(lt)>0:
				#print in_sesh,go,start
				if in_sesh and len(lt)==7:
					go = True
					start = True
				else:
					go = False
					#print "go!"

				if in_sesh and start and not go:
					#print "adding!"
					in_sesh=False
					go = False
					start = False
					zone_arr += [np.array(zt)]
					Rad_arr += [np.array(rt)]
					Rho_arr += [np.array(rhot)]
					ephi_arr += [np.array(ephit)]
					dvol_arr += [np.array(dvolt)]
					Temp_arr += [np.array(Tt)]
					
					zt = []
					rt = []
					rhot = []
					ephit = []
					dvolt = []
					Tt = []

				if lt[0] == "Time=":
					Time_arr += [float(lt[1])]
					Te_inf_arr += [float(lt[3])]
					in_sesh = True
					#start = False
					#print "Time!!!"
					#print lt[1]


				if in_sesh and go:
					#print "adding!!"
					#print lt
					zt += [float(lt[0])]
					rt += [float(lt[1])]
					rhot += [float(lt[2])]
					ephit += [float(lt[3])]
					dvolt += [float(lt[4])]
					Tt += [float(lt[5])]

		zone_arr += [np.array(zt)]
		Rad_arr += [np.array(rt)]
		Rho_arr += [np.array(rhot)]
		ephi_arr += [np.array(ephit)]
		dvol_arr += [np.array(dvolt)]
		Temp_arr += [np.array(Tt)]
					
		self.zone_arr = np.array(zone_arr)
		self.Rad_arr = np.array(Rad_arr)
		self.Rho_arr = np.array(Rho_arr)
		self.ephi_arr = np.array(ephi_arr)
		self.dvol_arr = np.array(dvol_arr)
		self.Temp_arr = np.array(Temp_arr)
					
		self.Time_arr =np.array(Time_arr)
		self.Te_inf_arr = np.array(Te_inf_arr)

	def _find_time(self,time):
		self._temp_arg = np.argmin(np.abs(self.Time_arr-time))

	def _find_Temp(self,Temp):
		'''
		Temp in keV
		'''
		Temp_K = Temp/Kelvin_to_keV
		#print Temp_K
		self._temp_arg = np.argmin(np.abs(self.Te_inf_arr-Temp_K))

				#print self._temp_arg
				#print self.Rad_arr[self._temp_arg]
				#print self.star_rad_core[1:]
				#print self.Temp_arr[self._temp_arg]
							  
		Ts_local_interp = interp1d(self.Rad_arr[self._temp_arg],self.Temp_arr[self._temp_arg],bounds_error=False,fill_value=self.Temp_arr[self._temp_arg][-1])(self.star_rad_core[1:])
		Ts_local_interp = np.insert(Ts_local_interp,0,Ts_local_interp[0])
		self.Ts_local_interp_core = Ts_local_interp

		Ts_local_interp = interp1d(self.Rad_arr[self._temp_arg],self.Temp_arr[self._temp_arg])(self.star_rad[1:])
		Ts_local_interp = np.insert(Ts_local_interp,0,Ts_local_interp[0])
		self.Ts_local_interp_crust = Ts_local_interp


		#z,valuesA = np.load('/nfs/turbo/bsafdi/buschman/NSaxion-master/python/IpnA.npy')
		#z,valuesB = np.load('/nfs/turbo/bsafdi/buschman/NSaxion-master/python/IpnB.npy')
		#self.IpnA_interp = interp1d(z,valuesA)
		#self.IpnB_interp = interp1d(z,valuesB)



	def _do_electron(self,gaee,gann,Ts_local_interp):
		Gamma = 22.73 * self.star_zion**2 / self.star_aion**(1./3.) * (1e6/Ts_local_interp) * (self.star_rho/1e6)**(1./3.)
		whr_solid = np.where( Gamma > 180. )		
		whr_liquid = np.where( Gamma <= 180. )		
	
		#print Gamma
		#print 'rho:', self.star_rho
		#print 'nb:', self.star_nb
	
		x = np.log10( self.star_rho/1e6 )  ## check log!!
		whr_xs = np.where( x <= 11.4 ) 
		whr_xl = np.where( x > 11.4 ) 

		uSL = np.empty( len(Gamma) )
		uSL[whr_solid]  = 0.488049 + 1.25585*Gamma[whr_solid]/1e3 - 0.743902*(Gamma[whr_solid]/1e3)**2
		uSL[whr_liquid] = 0.672409 + 0.182774*Gamma[whr_liquid]/1e3 + 0.144817*(Gamma[whr_liquid]/1e3)**2

		logFSL = np.empty( len(Gamma) )
		logFSL[whr_xs] = 0.21946 + 0.00287263*x[whr_xs]**2 - 0.000142016*x[whr_xs]**4 - (1.-uSL[whr_xs])
		logFSL[whr_xl] =-6.47808 + 0.068645*x[whr_xl]**2 - 0.000252677*x[whr_xl]**4 - (1.-uSL[whr_xl])

		self.electron_epsilon = 2.39e4 * self.star_zion**2 / self.star_aion * (gaee/1e-14)**2 * (Ts_local_interp/1e6)**4 * (self.star_nb/1.292e-3) * 10**( logFSL ) # check log!!
		#print 10**logFSL


		self.neutron_crust_epsilon = 1.2e15*self._neutron_crust_star_factor*(Ts_local_interp/1e8)**(11/2.)*(gann/1e-10)**2 #erg/cm^3/s
		drt = (self.Rad_arr[self._temp_arg][:-1]-self.Rad_arr[self._temp_arg][1:])*1e2 #cm
#print drt[-1]
		drt =np.append(drt,drt[-1])
		dvdr_interp = interp1d(self.Rad_arr[self._temp_arg],self.dvol_arr[self._temp_arg]/drt)(self.star_rad) #cm^2
		ephi_interp = interp1d(self.Rad_arr[self._temp_arg],self.ephi_arr[self._temp_arg])(self.star_rad)
		self.ephi_interp = ephi_interp
		self.dvdr_interp = dvdr_interp
		self.drt = drt
		dr_star =(self.star_rad[1:] -self.star_rad[:-1])*1e2 #cm
		dr_star = np.append(dr_star,dr_star[-1]) 
		self.dr_star = dr_star
		self.dv_star = dr_star*dvdr_interp #cm^3
		self.electron_flux_binned = self.electron_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		self.neutron_crust_flux_binned = self.neutron_crust_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		self.total_electron_flux = np.sum(self.electron_flux_binned)
		self.total_neutron_crust_flux = np.sum(self.neutron_crust_flux_binned)
		self.Teff_infty = np.mean(Ts_local_interp*self.ephi_interp)

	def _do_nucleon(self,gann,gapp,gaee,gamumu,Ts_local_interp,Ts_local_interp_crust,Es,superfluid=False,PBF=False):
		mstn_interp = interp1d(self.EOS_rho,self.EOS_mstn)(self.star_rho_core)
		mstp_interp = interp1d(self.EOS_rho,self.EOS_mstp)(self.star_rho_core)

		# aNN
		self.ann_epsilon = 1.8274e12*self._eann_star_factor*(Ts_local_interp/1e8)**6*(gann/1e-10)**2 * mstn_interp**2#erg/cm^3/s
		self.app_epsilon = 1.8274e12*self._eapp_star_factor*(Ts_local_interp/1e8)**6*(gapp/1e-10)**2 * mstn_interp**2#erg/cm^3/s
		g = gapp+gann
		h = gapp-gann
		self.anp_epsilon = 2.008e12*(self._eanp_star_factor_h*h**2+self._eanp_star_factor_g*g**2)/(1e-10)**2*(Ts_local_interp/1e8)**6 * mstn_interp**2
	
		#print self.ann_epsilon
		#print self.app_epsilon
		#print self.anp_epsilon
	#print (self._eanp_star_factor_h*h**2+self._eanp_star_factor_g*g**2)
		if PBF: 
			# 1s0 n
			whst_s_n = np.where(Ts_local_interp_crust<self.star_Tcn)[0]
			Delta_T_s_n = np.zeros(len(Ts_local_interp_crust))
		if len(whst_s_n)>0:
			tau = Ts_local_interp_crust[whst_s_n]/self.star_Tcn[whst_s_n]
			Delta_T_s_n[whst_s_n] = Ts_local_interp_crust[whst_s_n] * np.sqrt(1.-tau)*( 1.456 - 0.157/np.sqrt(tau) + 1.764/tau )
						#3.06*self.star_Tcn[whst_s_n]*np.sqrt(1-Ts_local_interp_crust[whst_s_n]/self.star_Tcn[whst_s_n])
		zn = Delta_T_s_n / Ts_local_interp_crust
		Ias_n = (0.158151*zn**2+0.543166*zn**4)*np.sqrt(1+np.pi*zn/4./0.543166**2)*np.exp(0.0535359-np.sqrt(4*zn**2+0.0535359**2))
		self.PBS_s_n_epsilon = 7.01e14*(gann/1e-10)**2*self._PBF_s_n_star_factor*(Ts_local_interp_crust/3e8)**5   * ( 1./m_n )**2   *  Ias_n/2.2e-2

		# 1s0 p
		whst_s_p = np.where(Ts_local_interp <self.star_Tcp_core)[0]
		Delta_T_s_p = np.zeros(len(Ts_local_interp))
		if len(whst_s_p)>0:
			tau = Ts_local_interp[whst_s_p]/self.star_Tcn_core[whst_s_p]
			Delta_T_s_p[whst_s_p] = Ts_local_interp[whst_s_p] * np.sqrt(1.-tau)*( 1.456 - 0.157/np.sqrt(tau) + 1.764/tau )	
			#Delta_T_s_p[whst_s_p] = 3.06*self.star_Tcp_core[whst_s_p]*np.sqrt(1-Ts_local_interp[whst_s_p]/self.star_Tcp_core[whst_s_p])
		zn = Delta_T_s_p / Ts_local_interp
		Ias_p = (0.158151*zn**2+0.543166*zn**4)*np.sqrt(1+np.pi*zn/4./0.543166**2)*np.exp(0.0535359-np.sqrt(4*zn**2+0.0535359**2))
		self.PBS_s_p_epsilon = 7.01e14*(gapp/1e-10)**2*self._PBF_s_p_star_factor*(Ts_local_interp/3e8)**5   * ( 1./mstp_interp )**2   *  Ias_p/2.2e02

		# 3p2 A
		whst_3p2A = np.where(Ts_local_interp < self.star_Tcn_core)[0]
		Delta_T_3p2A = np.zeros(len(Ts_local_interp))
		if len(whst_3p2A)>0:
			tau = Ts_local_interp[whst_3p2A]/self.star_Tcn_core[whst_3p2A]
			Delta_T_3p2A[whst_3p2A] = Ts_local_interp[whst_3p2A] * np.sqrt(1.-tau)*( 0.7893 + 1.118/tau ) 
				#3.06*self.star_Tcn_core[whst_3p2]*np.sqrt(1-Ts_local_interp[whst_3p2]/self.star_Tcn_core[whst_3p2])
		zn = Delta_T_3p2A / Ts_local_interp
		IanPA = self.IpnA_interp(zn)
		self.PBS_pA_epsilon = 1.67e15*(gann/1e-10)**2*self._PBF_p_star_factor*(Ts_local_interp/3e8)**5*(IanPA/5.96e-3)

		# 3p2 B
		whst_3p2B = np.where(Ts_local_interp < self.star_Tcn_core)[0]
		Delta_T_3p2B = np.zeros(len(Ts_local_interp))
		if len(whst_3p2B)>0:
			tau = Ts_local_interp[whst_3p2B]/self.star_Tcn_core[whst_3p2B]
			Delta_T_3p2B[whst_3p2B] = Ts_local_interp[whst_3p2B] * np.sqrt(1.-tau**4)/tau*( 2.030 - 0.4903*tau**4 + 0.1727*tau**8)
				#3.06*self.star_Tcn_core[whst_3p2]*np.sqrt(1-Ts_local_interp[whst_3p2]/self.star_Tcn_core[whst_3p2])
		zn = Delta_T_3p2B / Ts_local_interp
		IanPB = self.IpnB_interp(zn)
		self.PBS_pB_epsilon = 1.67e15*(gann/1e-10)**2*self._PBF_p_star_factor*(Ts_local_interp/3e8)**5*(IanPB/5.96e-3)


		self.Delta_T_s_n = Delta_T_s_n
		self.Delta_T_s_p = Delta_T_s_p
		self.Delta_T_3p2A = Delta_T_3p2A
		self.Delta_T_3p2B = Delta_T_3p2B
		self.Ts_local_interp_crust = Ts_local_interp_crust
		self.Ts_local_interp = Ts_local_interp

		self.ann_epsilon_old = copy.deepcopy(self.ann_epsilon)
		if superfluid:
			Rnn,Rnp_n = return_Rs(Ts_local_interp,self.star_Tcn_core)
			Rpp,Rnp_p = return_Rs(Ts_local_interp,self.star_Tcp_core)
			Rnp = np.minimum(Rnp_n,Rnp_p)
			#print Rnn
			self.ann_epsilon*=Rnn
			self.app_epsilon*=Rpp
			self.anp_epsilon*=Rnp


		drt = (self.Rad_arr[self._temp_arg][:-1]-self.Rad_arr[self._temp_arg][1:])*1e2 #cm
#print drt[-1]
		drt =np.append(drt,drt[-1])
		# I'm here
		dvdr_interp = interp1d(self.Rad_arr[self._temp_arg],self.dvol_arr[self._temp_arg]/drt,bounds_error=False,fill_value=0)(self.star_rad_core[1:]) #cm^2
		ephi_interp = interp1d(self.Rad_arr[self._temp_arg],self.ephi_arr[self._temp_arg],bounds_error=False,fill_value=0)(self.star_rad_core[1:])
		dvdr_interp = np.insert(dvdr_interp,0,dvdr_interp[0])
		ephi_interp = np.insert(ephi_interp,0,ephi_interp[0])
		self.ephi_interp_core = ephi_interp
		self.dvdr_interp_core = dvdr_interp
		self.drt_core = drt
		dr_star =(self.star_rad_core[1:] -self.star_rad_core[:-1])*1e2 #cm
		dr_star = np.append(dr_star,dr_star[-1]) 
		self.dr_star_core = dr_star
		self.dv_star_core = dr_star*dvdr_interp #cm^3

		dvdr_interp_crust = interp1d(self.Rad_arr[self._temp_arg],self.dvol_arr[self._temp_arg]/drt)(self.star_rad[1:]) #cm^2
		ephi_interp_crust = interp1d(self.Rad_arr[self._temp_arg],self.ephi_arr[self._temp_arg])(self.star_rad[1:])
		dvdr_interp_crust = np.insert(dvdr_interp_crust,0,dvdr_interp_crust[0])
		ephi_interp_crust = np.insert(ephi_interp_crust,0,ephi_interp_crust[0])
		self.ephi_interp_crust = ephi_interp_crust
		self.dvdr_interp_crust = dvdr_interp_crust
		self.drt_crust = drt
		dr_star_crust =(self.star_rad[1:] -self.star_rad[:-1])*1e2 #cm
		dr_star_crust = np.append(dr_star_crust,dr_star_crust[-1]) 
		self.dr_star_crust = dr_star_crust
		self.dv_star_crust = dr_star_crust*dvdr_interp_crust #cm^3

		self.ann_flux_binned = self.ann_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		self.ann_flux_binned_old = self.ann_epsilon_old*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		self.app_flux_binned = self.app_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		self.anp_flux_binned = self.anp_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
		if PBF:
			self.PBS_s_n_flux_binned = self.PBS_s_n_epsilon*dr_star_crust*dvdr_interp_crust*ephi_interp_crust #erg/s #*volume_interp 
			self.PBS_s_p_flux_binned = self.PBS_s_p_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
			self.PBS_pA_flux_binned = self.PBS_pA_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 
			self.PBS_pB_flux_binned = self.PBS_pB_epsilon*dr_star*dvdr_interp*ephi_interp #erg/s #*volume_interp 

		self.ephi_interp = ephi_interp
		self.ephi_interp_crust = ephi_interp_crust

		if PBF:
			self.PBS_s_n_spectrum = np.zeros((len(Es),len(Ts_local_interp_crust)))
			self.PBS_s_n_total = np.zeros((len(Es),len(Ts_local_interp_crust)))
			for i in range(len(Es)):
				omega = Es[i]/8.61833e-8 / ephi_interp_crust[whst_s_n] # K
				argwt = omega / (2.*Delta_T_s_n[whst_s_n])
				self.PBS_s_n_spectrum[i,whst_s_n] = 1./(2.*Delta_T_s_n[whst_s_n] * 8.61833e-8) * argwt**3 / np.sqrt(argwt**2-1) / (np.exp(omega/(2.*Ts_local_interp_crust[whst_s_n]))+1.)**2
				self.PBS_s_n_spectrum[i,np.where( np.isnan(self.PBS_s_n_spectrum[i] ))] = 0
				zn = Delta_T_s_n / Ts_local_interp_crust
				self.N_PBS_s_n = Ias_n / zn**5
				self.PBS_s_n_total[i] = self.PBS_s_n_flux_binned * self.PBS_s_n_spectrum[i] / self.N_PBS_s_n
				self.PBS_s_n_total[i,np.where( zn == 0 )] = 0

			self.PBS_s_p_spectrum = np.zeros((len(Es),len(Ts_local_interp)))
			self.PBS_s_p_total = np.zeros((len(Es),len(Ts_local_interp)))
			for i in range(len(Es)):
				omega = Es[i]/8.61833e-8 / ephi_interp_crust[whst_s_p] # K
				argwt = omega / (2.*Delta_T_s_p[whst_s_p])
				self.PBS_s_p_spectrum[i,whst_s_p] = 1./(2.*Delta_T_s_p[whst_s_p] * 8.61833e-8) * argwt**3 / np.sqrt(argwt**2-1) / (np.exp(omega/(2.*Ts_local_interp[whst_s_p]))+1.)**2
				self.PBS_s_p_spectrum[i,np.where( np.isnan(self.PBS_s_p_spectrum[i] ))] = 0
				zn = Delta_T_s_p / Ts_local_interp
				self.N_PBS_s_p = Ias_p / zn**5
				self.PBS_s_p_total[i] = self.PBS_s_p_flux_binned * self.PBS_s_p_spectrum[i] / self.N_PBS_s_p
				self.PBS_s_p_total[i,np.where( zn == 0 )] = 0


			self.PBS_pA_spectrum = np.zeros((len(Es),len(Ts_local_interp)))
			self.PBS_pA_total = np.zeros((len(Es),len(Ts_local_interp)))
			#print Delta_T_3p2A,Delta_T_3p2B
			for i in range(len(Es)):
				omega = Es[i]/8.61833e-8 / ephi_interp[whst_s_p] # K
				argwt = omega / (2.*Delta_T_s_p[whst_s_p])
				exis = np.where( argwt > 1. )
			#print np.max(Es[-1]/8.61833e-8 / ephi_interp),np.min(Es[0]/8.61833e-8 / ephi_interp)
			#print Delta_T_s_p	
			#print Ts_local_interp


		self.total_ann_flux = np.sum(self.ann_flux_binned)
		self.total_ann_flux_old = np.sum(self.ann_flux_binned_old)
		self.total_app_flux = np.sum(self.app_flux_binned)
		self.total_anp_flux = np.sum(self.anp_flux_binned)
		if PBF:
			self.total_PBS_s_n_flux = np.sum(self.PBS_s_n_total,axis=1)
			self.total_PBS_s_p_flux = np.sum(self.PBS_s_p_total,axis=1)
			self.total_PBS_pA_flux = np.sum(self.PBS_pA_total,axis=1)
			#self.total_PBS_pA_flux = np.sum(self.PBS_pA_total,axis=1)
		else:
			self.total_PBS_s_n_flux = 0
			self.total_PBS_s_p_flux = 0 
			self.total_PBS_pA_flux = 0 
			#self.total_PBS_pA_flux = np.sum(self.PBS_pA_total,axis=1)
	
		#print self.total_ann_flux	
		#print self.total_app_flux	
		#print self.total_anp_flux	
		self.total_aNN_flux = self.total_ann_flux + self.total_app_flux+self.total_anp_flux#+self.total_PBS_s_p_flux+self.total_PBS_s_n_flux+self.total_PBS_p_flux
		self.Teff_infty_core = np.mean(Ts_local_interp*self.ephi_interp_core) ####

	def do_nucleon_Temp(self,alpha,Es,gann=1e-10,gapp=1e-10,gaee=5e-14,gamumu=1e-11,superfluid=False,PBF=False):
		self._do_nucleon(gann,gapp,gaee,gamumu,self.Ts_local_interp_core*alpha,self.Ts_local_interp_crust*alpha,Es,superfluid=superfluid,PBF=PBF)

	def do_electron_Temp(self,alpha,gaee=5e-14,gann=1e-10):
		self._do_electron(gaee,gann,self.Ts_local_interp_crust*alpha)

	def do_electron_time(self,time,gaee=5e-14):
		self._find_time(time)
		self._do_electron(gaee)







