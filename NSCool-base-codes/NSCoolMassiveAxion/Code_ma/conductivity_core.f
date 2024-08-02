c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c ******************** CORE CONDUCTIVITY ******************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine con_core(icon_core,debug,
     1                    Temp,kf_e,kf_mu,
     2                    kf_p ,mst_p ,Tc_p ,        ! proton
     3                    kf_n ,mst_n ,Tc_n ,isfn,   ! neutron
     4                    kf_la,mst_la,Tc_la,        ! lambda
     5                    kf_sm,mst_sm,Tc_sm,        ! sigma-
     6                    kf_s0,mst_s0,Tc_s0,        ! sigma0
     7                    kf_sp,mst_sp,Tc_sp,        ! sigma+
     8                    f_had,
     9                    sigma,lambda,
     x                    nu_e_s,nu_e_l)
c *********************************************************************
c ***** Checked on Sept. 7, 2009                                      *
c *********************************************************************
       implicit real*8(a-h,k-z)
c Use simple Flowers & Itoh (1981) formula: ***************************
        if (icon_core.eq.1) then
         lambda=1.d23*(kf_n**3/1.6d0) / (Temp/1.d8)
         sigma=0.d0
c Use Yakovlev et al. calculations: ***********************************
        else if (icon_core.eq.2) then
         call con_core_lep(Temp,kf_e,kf_mu,
     1                     kf_p ,mst_p ,Tc_p ,
     4                     kf_sm,mst_sm,Tc_sm,
     6                     kf_sp,mst_sp,Tc_sp,
     8                     sigma_lep,lambda_e,lambda_mu,debug,
     9                     nu_e_s1,nu_e_l1)
         icontrol=1         ! =1 uses the full thing
         call con_core_bar(Temp,kf_e,kf_mu,
     1                     kf_p ,mst_p ,Tc_p ,
     2                     kf_n ,mst_n ,Tc_n ,isfn,
     3                     kf_la,mst_la,Tc_la,
     4                     kf_sm,mst_sm,Tc_sm,
     5                     kf_s0,mst_s0,Tc_s0,
     6                     kf_sp,mst_sp,Tc_sp,
     8                     sigma_bar,lambda_bar,debug,
     9                     nu_e_s2,nu_e_l2,icontrol)
c        Quark conductivity NOT defined !
         lambda_qrk=0.d0
         sigma_qrk =0.d0
c ***
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         lambda=lambda_e+lambda_mu+
c     1          lambda_bar*f_had+
c     2          lambda_qrk*(1.d0-f_had)
cccccccccccccc
         lambda=lambda_e+lambda_mu
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         sigma =0.d0
c Use simple 1/T formula: *********************************************
        else if (icon_core.ge.20) then
         lambda=float(icon_core) / (Temp/1.d8)
         sigma=0.d0
c Use simple T-independent formula: ***********************************
        else if (icon_core.le.-20) then
         lambda=abs(float(icon_core))
c That's All ! ********************************************************
        end if
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine con_core_lep(Temp,kf_e,kf_m,
     1                    kf_p ,mst_p0 ,Tc_p ,        ! proton
     2                    kf_sm,mst_sm0,Tc_sm,        ! sigma-
     3                    kf_sp,mst_sp0,Tc_sp,        ! sigma+
     4                    sigma_lep,lambda_e,lambda_m,debug,
     5                    nu_e_s,nu_e_l)
c *********************************************************************
c Checked on Sept 4-7, 2009
c *********************************************************************
c Calculates the thermal conductivity in the core from
c Shternin & Yakovlev, Phys. Rev. D75, 103004, 2007
c for electrons and muons contributions.
c No electrical conductivity is calculated !
c Sigma- and Sigma+ scattering is NOT implemented !
c *********************************************************************
       implicit real*8(a-h,k-z)
       parameter(pi=3.1415926535d0)
       parameter(hbar=1.0546d-27,c=2.99792d10,kb=1.3806d-16)
       parameter(mp=1.6726d-24,me=9.1095d-28,e=4.803d-10)
c *********************************************************************
c Pairing correction factors: [t=T/Tc y=u(t)]
c 1S0:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
c *********************************************************************
       if (debug.eq.1.2d0) print *,'Entering con_core_lep:',
     1                   ' T, kfeo=',T,kf_e
c *** In case there are no leptons:
       if (kf_e.eq.0.d0) then
        lambda_e  =0.d0
        lambda_m  =0.d0
        lambda_lep=0.d0
        sigma_lep =0.d0
        if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1           ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
        return
       end if
c *** In case there are no muons:
      if (kf_m.eq.0.d0) then
       muons=0.d0
      else
       muons=1.d0
      end if
c ***
       T8=Temp/1.d8
c THERMAL CONDUCTIVITY: ***********************************************
c Trick to automatically eliminate absent baryons:
c That's because they show up only through the phase space integrals
c and this comes to zero if mst=0 !
       mst_p=mst_p0
       if (kf_p .eq.0.d0) mst_p =0.d0
       mst_sm=mst_sm0
       if (kf_sm.eq.0.d0) mst_sm=0.d0
       mst_sp=mst_sp0
       if (kf_sp.eq.0.d0) mst_sp=0.d0
c **** Define Fermi momenta ratios (massively used):
       rkf_0 =1.68d0/kf_e
       rkf_m =kf_m  /kf_e
       rkf_p =kf_p  /kf_e
c **** Screening momenta ratios:
c      (kf_e/q_l)**3:
       rkf_e_ql3=1./(0.00929d0*
     x          (1.d0+rkf_m+2.83d0*mst_p*rkf_0*rkf_p))**1.5
c      (kf_e/q_t)**2:
       rkf_e_qt2=1./(0.00929d0*
     x          (1.d0+rkf_m**2+rkf_p**2))
c **** Longitudinal collisional frequencies:
       nu_ee_par=1.43d11 *          rkf_0   *rkf_e_ql3 * T8**2
       nu_em_par=nu_ee_par *muons
       nu_ep_par=1.15d12 * mst_p**2*rkf_0**2*rkf_e_ql3 * T8**2
       if (muons.eq.1.d0) then
        nu_mm_par=nu_ee_par /rkf_m
        nu_me_par=nu_mm_par
        nu_mp_par=nu_ep_par /rkf_m
       else
        nu_mm_par=0.d0
        nu_me_par=0.d0
        nu_mp_par=0.d0
       end if
c **** Transverse ("perpendicular") collisional frequencies:
       nu_ee_per=6.49d14 * rkf_e_qt2 * T8
       nu_em_per=nu_ee_per * rkf_m**2
       nu_ep_per=nu_ee_per * rkf_p**2
       nu_mm_per=nu_ee_per * rkf_m**3
       nu_me_per=nu_ee_per * rkf_m
       nu_mp_per=nu_ee_per * rkf_m*rkf_p**2
c **** Cross ("prime") collisional frequencies:
       nu_ee_pri=4.38d12 * rkf_0**(2./3.)*rkf_e_qt2**(1./3.)*
     x                     rkf_e_ql3**(2./3.) * T8**(5./3.)
       if (muons.ne.0.d0) then
        nu_em_pri=nu_ee_pri *rkf_m**2
        nu_mm_pri=nu_ee_pri *rkf_m**3
        nu_me_pri=nu_ee_pri /rkf_m
       else
        nu_em_pri=0.d0
        nu_mm_pri=0.d0
        nu_me_pri=0.d0
       end if
c **** Effect of pairing:
c Proton 1S0 pairing:
       if ((Temp.le.Tc_p).and.(kf_p.gt.0.d0)) then
        y=u_1s0(Temp/Tc_p)
        r=(kf_e**2+kf_m**2)/kf_p**2
        R_l_pri  =(r+1.d0)**(1./3.) /
     x            ((r+1.d0)**2-0.757d0*y+(0.50651d0*y)**2)**(1./6.)
         p1=0.48d0-0.17d0*r
         p3=((1.d0-p1)*54.d0/4.d0/pi**2/r)**2
        R_tot_per=p1*dexp(-0.14d0*y**2)+(1.d0-p1)/dsqrt(1.d0+p3*y**2)
c       This is R_p_par fron the paper:
c
c        R_p_par  =(1.d0
c     x   + (26.33d0*y**2+0.376d0*y**4)*dexp(-dsqrt((3.675d0)**2+y**2))
c     x   + 0.742d0*(dexp((1.673d0)**2-dsqrt((1.673d0)**2+y**2))-1.d0)
c     y    ) *dexp((1.361d0)**2-dsqrt((1.361d0)**2+y**2))
c
c            which is obviously wrong !
c       Here is my guess for the correct value:
c
c        R_p_par  =(1.d0
c     x   + (26.33d0*y**2+0.376d0*y**4)*dexp(-dsqrt((3.675d0)**2+y**2))
c     x   + 0.742d0*(dexp((1.673d0)-dsqrt((1.673d0)**2+y**2))-1.d0)
c     y    ) *dexp((1.361d0)-dsqrt((1.361d0)**2+y**2))
c
c       Here is Peter Shternin's new version of the fit 
c       (it gives essentially the same result as my guess):
c
        R_p_par=(0.998d0 + 
     x           (2.04d0 + 0.68d0*dsqrt(y) + 5.7d0*y**2 + 1.71d0*y**4 ) 
     y           * dexp(-1.04d0*y) ) * dexp(-dsqrt(1.23d0+y**2))
       else
        R_tot_per=1.d0
        R_p_par  =1.d0
        R_l_pri  =1.d0
       end if
c **** Adjust collisional frequencies for pairing and add them
       nu_e_par=nu_ee_par+nu_em_par+nu_ep_par*R_p_par
       nu_e_per=(nu_ee_per+nu_em_per+nu_ep_per)*R_tot_per
       nu_e_pri=nu_ee_pri*R_l_pri
       nu_m_par=nu_mm_par+nu_me_par+nu_mp_par*R_p_par
       nu_m_per=(nu_mm_per+nu_me_per+nu_mp_per)*R_tot_per
       nu_m_pri=nu_mm_pri*R_l_pri
       nu_em_pri=nu_em_pri*R_l_pri
       nu_me_pri=nu_me_pri*R_l_pri
c **** Total collisional frecuencies:
       nu_e=nu_e_par+nu_e_per+nu_e_pri
       nu_m=nu_m_par+nu_m_per+nu_m_pri
c **** Relaxation times:
       if (muons.ne.0.d0) then
        tau_e = (nu_m-nu_em_pri) / (nu_e*nu_m-nu_em_pri*nu_me_pri)
        tau_m = (nu_e-nu_me_pri) / (nu_e*nu_m-nu_em_pri*nu_me_pri)
       else
        tau_e=1.d0/nu_e
        tau_m=0.d0
       end if
c **** thermal conductivies:
       lambda_e = 1.70d24 * T8 *(1.d15*tau_e) * (kf_e/1.68d0)**2
       lambda_m = 1.70d24 * T8 *(1.d15*tau_m) * (kf_m/1.68d0)**2
       lambda_lep=lambda_e+lambda_m

c ELECTRICAL CONDUCTIVITY: ********************************************
        sigma_lep=0.d0
c *********************************************************************
        nu_e_l=nu_e
        nu_e_s=0.d0
c ****
       if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1          ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
c ****
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c     THIS SUBROUTINE SHOULD NOT BE USED:
c     wrong screening for transverse plasmons.
c     Use subroutine con_core_lep instead
      subroutine con_core_lep_OLD(Temp,kf_e,kf_mu,
     1                    kf_p ,mst_p ,Tc_p ,        ! proton
     4                    kf_sm,mst_sm,Tc_sm,        ! sigma-
     6                    kf_sp,mst_sp,Tc_sp,        ! sigma+
     8                    sigma_lep,lambda_lep,debug,
     9                    nu_e_s,nu_e_l)
c *********************************************************************
c Checked on December 2, 2000 and March 16, 2001 and Sept 2, 2009
c *********************************************************************
c This corrects a error in the paper and uses all nu' = 0
c It does NOT correct the error that zl_1S0=1 (but that's easy to do !)
c *********************************************************************
c Calculates the thermal conductivity in the core from
c Gnedin & Yakovlev, Nucl. Phys. A582 (1995), p. 697
c for electrons and muons contributions.
c Calculates the electrical conductivity
c from Baym, Pethick & Pine, Nature 224 (1969), p. 674
c CHECK IF NOT BETTER TO USE WIEDEMANN-FRANZ LAW !
c *********************************************************************
       implicit real*8(a-h,k-z)
       parameter(pi=3.1415926535d0)
       parameter(hp=1.0546d-27,c=2.99792d10,kb=1.3806d-16)
       parameter(mp=1.6726d-24,me=9.1095d-28,e=4.803d-10)
c *********************************************************************
c Pairing correction factors: [t=T/Tc y=u(t)]
c 1S0:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
      rl_1S0(y)=(0.7694d0+dsqrt((0.2306d0)**2+(0.07207d0*y)**2)+
     1    (27.d0*y**2+0.1476d0*y**4)*dexp(-dsqrt((4.273d0)**2+y**2))+
     2     0.5051d0*(dexp(4.273d0-dsqrt((4.273d0)**2+y**2))-1.d0))*    
     3     dexp(1.187d0-dsqrt((1.187d0)**2+y**2))
      zl_1S0(y)=dsqrt(0.9443d0+dsqrt((0.0557d0)**2+(0.1886d0*y)**2))*
     1         dexp(1.753d0-dsqrt((1.753d0)**2+y**2))
c *********************************************************************
       if (debug.eq.1.2d0) print *,'Entering con_core_lep:',
     1                   ' T, kfe=',T,kf_e
c *** In case there are no leptons:
       if (kf_e.eq.0.) then
        lambda_lep=0.d0
        sigma_lep =0.d0
        if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1           ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
        return
       end if
c THERMAL CONDUCTIVITY:
c Trick to automatically eliminate absent baryons:
c That's because they show up only through the phase space integrals
c and this comes to zero if mst=0 !
       if (kf_p .eq.0.) mst_p =0.
       if (kf_sm.eq.0.) mst_sm=0.
       if (kf_sp.eq.0.) mst_sp=0.
c **** Effect of pairing:
c Proton 1S0 pairing:
       if (Temp.le.Tc_p) then
        y=u_1s0(Temp/Tc_p)
        z_p=zl_1S0(y)
        r_p=rl_1S0(y)
       else
        z_p=1.d0
        r_p=1.d0
       end if
c Sigma- 1S0 pairing:
       if (Temp.le.Tc_sm) then
        y=u_1s0(Temp/Tc_sm)
        z_sm=zl_1S0(y)
        r_sm=rl_1S0(y)
       else
        z_sm=1.d0
        r_sm=1.d0
       end if
c Sigma+ 1S0 pairing:
       if (Temp.le.Tc_sp) then
        y=u_1s0(Temp/Tc_sp)
        z_sp=zl_1S0(y)
        r_sp=rl_1S0(y)
       else
        z_sp=1.d0
        r_sp=1.d0
       end if
c **** Screening momentum:
       q0_over_pfe=0.0964d0*sqrt(
     1       1.d0+abs(kf_mu/kf_e)+
     2       2.83d0*mst_p *1.68*kf_p /kf_e**2*z_p +
     3       2.83d0*mst_sm*1.68*kf_sm/kf_e**2*z_sm+
     4       2.83d0*mst_sp*1.68*kf_sp/kf_e**2*z_sp)
       pfe_over_q0_cube=1.d0/q0_over_pfe**3
c ****
       xxx=1.15d12*pfe_over_q0_cube*abs(1.68/kf_e)**3
       xnu_ep =xxx * mst_p **2*r_p
       xnu_esm=xxx * mst_sm**2*r_sm
       xnu_esp=xxx * mst_sp**2*r_sp
       xnu_ee=3.58d11*pfe_over_q0_cube*abs(1.68/kf_e)
       if (kf_mu.ne.0.d0) then
        xnu_mup =xnu_ep *abs(kf_e/kf_mu)
        xnu_musm=xnu_esm*abs(kf_e/kf_mu)
        xnu_musp=xnu_esp*abs(kf_e/kf_mu)
        xnu_emu =1.43d11*pfe_over_q0_cube*abs(1.68/kf_e)*
     1           (1.d0+0.5d0*(kf_mu/kf_e)**2)
        xnu_mue =xnu_emu*abs(kf_e/kf_mu)
        xnu_mumu=xnu_ee *abs(kf_mu/kf_e)**3*
     1      (1.d0+6./5.*(0.536d0/kf_mu)**2+2./5.*(0.536d0/kf_mu)**4)
       else
        xnu_mup =0.d0
        xnu_musm=0.d0
        xnu_musp=0.d0
        xnu_emu =0.d0
        xnu_mue =0.d0
        xnu_mumu=0.d0
       end if
       xnu_e =(xnu_ep +xnu_esm +xnu_esp +xnu_ee  +xnu_emu)
       lambda_e=1.7d24*1.2*(1.d15/xnu_e)*(kf_e/1.68)**2*(1.d8/Temp)
       if (kf_mu.ne.0.d0) then
        xnu_mu=(xnu_mup+xnu_musm+xnu_musp+xnu_mumu+xnu_mue)
        lambda_mu=lambda_e*(kf_mu/kf_e)**3*(xnu_e/xnu_mu)
       else
        lambda_mu=0.d0
       end if
       lambda_lep=lambda_e+lambda_mu
c ELECTRICAL CONDUCTIVITY:
c Core conductivity: from Baym, pethick & Pine, Nature 224 (1969), p. 674
c       nb=(kf_n**3+kf_p**3)/3./pi**2
c        sigma_lep=4.2e26*(nb/0.16)**3/(Temp/1.e9)**2
        sigma_lep=4.2d26/(Temp/1.d9)**2
c ****
       if (debug.eq.1.2d0) print *,'Exiting con_core_lep:',
     1          ' sigma_lep, lambda_lep=',sigma_lep,lambda_lep
c ****
        nu_e_l=xnu_e
        nu_e_s=nu_e_l
c ****
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine con_core_bar(Temp,kf_e,kf_mu,
     1                    kf_p ,mst_p ,Tc_p ,        ! proton
     2                    kf_n ,mst_n ,Tc_n ,isfn,   ! neutron
     3                    kf_la,mst_la,Tc_la,        ! lambda
     4                    kf_sm,mst_sm,Tc_sm,        ! sigma-
     5                    kf_s0,mst_s0,Tc_s0,        ! sigma0
     6                    kf_sp,mst_sp,Tc_sp,        ! sigma+
     7                    sigma_bar,lambda_bar,debug,
     8                    nu_e_s,nu_e_l,icontrol)
c *********************************************************************
c Checked on 2 Septembre 2009 against Fig 1, 3, & 4                   *
c *********************************************************************
c Calculates the neutron thermal conductivity in the core from        *
c Baiko, Haensel & Yakovlev, 2001, A&A 374, p 151                     *
c                                                                     *
c Chemical potentials mu's in MeV                                     *
c Fermi momenta kf's in fm^-1                                         *
c *********************************************************************
       implicit real*8(a-h,k-z)
       parameter(pi=3.1415926535d0)
       parameter(hbar=1.0546d-27,c=2.99792d10,kb=1.3806d-16)
       parameter(mu=1.6726d-24,me=9.1095d-28,e=4.803d-10)
c *********************************************************************
c Pairing correction factors: [t=T/Tc y=u(t)]
c 1S0:
       u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
       u_3p2(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c *********************************************************************
        if (debug.eq.1.2d0) print *,'Entering con_core_bar:',
     1                   ' T, kfe=',T,kf_e
c *********************************************************************
c      Neutrons contribution:
c ***
        Sn1=14.57d0/kf_n**1.5 *
     1      (1.d0-0.0788d0*kf_n+0.0883d0*kf_n**2) /
     2      (1.d0-0.1114d0*kf_n)
        Sn2=7.880d0/kf_n**2 *
     1      (1.d0-0.2241d0*kf_n+0.2006d0*kf_n**2) /
     2      (1.d0-0.1742d0*kf_n)
        Sp1=0.8007d0*kf_p/kf_n**2 *
     1      (1.d0+31.28d0*kf_p-0.0004285d0*kf_p**2+
     1       26.85d0*kf_n+0.08012d0*kf_n**2) /      
     2      (1.d0-0.5898d0*kf_n+0.2368d0*kf_n**2+
     2       0.5838d0*kf_p**2+0.884d0*kf_n*kf_p)      
        Sp2=0.3830d0*kf_p**4/kf_n**5.5 *
     1      (1.d0+102.d0*kf_p+53.91d0*kf_n) /
     2      (1.d0-0.7087d0*kf_n+0.2537d0*kf_n**2+
     2       9.404*kf_p**2-1.589d0*kf_n*kf_p)
c ***
        u=kf_n-1.665d0
        Kn1=(0.4583d0+0.892d0*u**2-0.5497d0*u**3-0.06205d0*kf_p
     1      +0.04022d0*kf_p**2+0.2122d0*u*kf_p)
     2      / mst_n**2
        u=kf_n-1.556d0
        Kn2=(0.4891d0+1.111d0*u**2-0.2283d0*u**3+0.01589d0*kf_p
     1      -0.02099*kf_p**2+0.2773*u*kf_p)
     2      / mst_n**2
        u=kf_n-2.126d0
        Kp1=(0.04377d0+1.100d0*u**2+0.1180d0*u**3+0.1626d0*kf_p 
     1       +0.3871d0*u*kf_p-0.2990d0*u**4)
     2      / mst_p**2
        u=kf_n-2.116d0
        Kp2=(0.0001313d0+1.248d0*u**2+0.2403d0*u**3+0.3257d0*kf_p
     1       +0.5536d0*u*kf_p-0.3237d0*u**4+0.09786d0*u**2*kf_p)
     2      / mst_p**2
        if (icontrol.ge.2) then
         Kn1=1.d0
         Kn2=1.d0
         Kp1=1.d0
         Kp2=1.d0
        end if
c *** Pairing effects:
        if (Temp.le.Tc_p) then
         tau=Temp/Tc_p
         yp=u_1s0(tau)
        else
         yp=0.d0
        end if
        if (Temp.le.Tc_n) then
         tau=Temp/Tc_n
         if (isfn.eq.1) then
          yn=u_1s0(tau)
         else if (isfn.eq.3) then
          yn=u_3p2(tau)
         else
          pause 'Subroutine con_core_bar: isfn badly defined !'
         end if
        else
         yn=0.d0
        end if
        call con_core_bar_pairing_supr(yn,yp,Rn1,Rn2,Rp1,Rp2,RC)
c ***
        Snn=Sn2*Kn2*Rn2 + 3.0d0*Sn1*Kn1*(1.d0*Rn1-Rn2)
        Snp=Sp2*Kp2*Rp2 + 0.5d0*Sp1*Kp1*(3.d0*Rp1-Rp2)
c       When T<<Tc_p and/or Tc_n Snn and Snp may vanish, so better use:
        Snn = max(Snn,1.d-200)
        Snp = max(Snp,1.d-200)
        if (icontrol.eq.3) Snp=0.d0
        nu_nn=3.48d15*  mst_n**3    *(Temp/1.d8)**2 *Snn
        nu_np=3.48d15*mst_n*mst_p**2*(Temp/1.d8)**2 *Snp
        tau_n=RC/(nu_nn+nu_np)
c ***
        lambda_n=7.2d23* (Temp/1.d8) *
     x           RC**2/mst_n* 1.d15/(nu_nn+nu_np) * (kf_n/1.68)**3
c *********************************************************************
c      Lambda contribution:
       lambda_la=0.d0
c *********************************************************************
c      Sigma0 contribution:
       lambda_s0=0.d0
c *********************************************************************
       lambda_bar=lambda_n+lambda_la+lambda_s0
c *********************************************************************
c ELECTRICAL CONDUCTIVITY:
        sigma_bar=0.d0
c *********************************************************************
       if (debug.eq.1.2d0) print *,'Exiting con_core_bar:',
     1            ' sigma_bar, lambda_bar=',sigma_bar,lambda_bar
c ****
        nu_e_l=0.d0
        nu_e_s=0.d0
c ****
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine con_core_bar_pairing_supr(yn,yp,Rn1,Rn2,Rp1,Rp2,RC)
       implicit real*8(a-h,k-z)
        if (yn.eq.0.d0) then
         Rn1=1.d0
         Rn2=1.d0
         RC =1.d0
        else
         Rn1=
     1       (2.d0/3.d0) * 
     1       (0.9468d0+dsqrt((0.0532d0)**2+0.5346d0*yn**2))**3 *
     1       dexp(0.377d0-dsqrt((0.377d0)**2+4.d0*yn**2)) +
     2       (1.d0/3.d0) *
     2       (1.d0+1.351d0*yn**2)**2 *
     2       dexp(0.169d0-dsqrt((0.169d0)**2+9.d0*yn**2))
         Rn2=
     1       0.5d0 *
     1       (0.6242d0+dsqrt((0.3758d0)**2+0.07198d0*yn**2))**3 *
     1       dexp(3.6724d0-dsqrt((3.6724d0)**2+4.d0*yn**2)) +
     2       0.5d0 *
     2       (1.d0+0.01211d0*yn**2)**9 *
     2       dexp(7.5351d0-dsqrt((7.5351d0)**2+9.d0*yn**2))
         RC=(0.647d0+dsqrt((0.353d0)**2+0.109d0*yn**2))**1.5 *
     x      dexp(1.39d0-dsqrt((1.39d0)**2+yn**2))
        end if
        if ((yn.eq.0.d0).and.(yp.eq.0.d0)) then
         Rp1=1.d0
         Rp2=1.d0
        else if ((yn.gt.0.d0).and.(yp.eq.0.d0)) then
         Rp1=(0.4459d0+dsqrt((0.5541d0)**2+0.03016d0*yn**2))**2 *
     1       dexp(2.1178d0-dsqrt((2.1178d0)**2+yn**2))
         Rp2=(0.801d0+dsqrt((0.199d0)**2+0.04645d0*yn**2))**2 *
     1       dexp(2.3569d0-dsqrt((2.3569d0)**2+yn**2))
        else if ((yn.eq.0.d0).and.(yp.gt.0.d0)) then
         Rp1=
     1       0.5d0 *
     1       (0.3695d0+dsqrt((0.6305)**2+0.01064d0*yp**2)) *
     1       dexp(2.4451d0-dsqrt((2.4451d0)**2+yp**2)) +
     2       0.5d0 *       
     2       (1.d0+0.1917*yp**2)**1.4 *
     2       dexp(4.6627d0-dsqrt((4.6627d0)**2+4.d0*yp**2))
         Rp2=
     1       0.0436d0 *       
     1       (dsqrt((3.345d0)**2+19.55*yp**2)-3.345d0) *
     1       dexp(2.0247d0-dsqrt((2.0247d0)**2+yp**2))
     2       + 0.0654d0 * dexp(8.992d0-dsqrt((8.992d0)**2+1.5d0*yp**2))
     3       + 0.8910d0 * dexp(9.627d0-dsqrt((9.627d0)**2+9.0d0*yp**2))
        else
         y_p=max(yn,yp)
         y_m=min(yn,yp)
         u_p=dsqrt(y_p**2+(1.485d0)**2)-1.485d0
         u_m=dsqrt(y_m**2+(1.485d0)**2)-1.485d0
         up =dsqrt(yp**2 +(1.485d0)**2)-1.485d0
         un =dsqrt(yn**2 +(1.485d0)**2)-1.485d0
         Rp1=
     1       dexp(-u_p-u_m) * (0.7751d0+0.4823d0*un+0.1124d0*up+
     1       0.04991d0*un**2+0.08513d0*un*up+0.01284d0*un**2*up)       
     2       + dexp(-2.d0*u_p) * (0.2249d0+0.3539d0*u_p-0.2189d0*u_m-
     2       0.6069d0*un*u_m+0.7362d0*up*u_p)       
         u_p=dsqrt(y_p**2+(1.761d0)**2)-1.761d0
         u_m=dsqrt(y_m**2+(1.761d0)**2)-1.761d0
         up =dsqrt(yp**2 +(1.761d0)**2)-1.761d0
         un =dsqrt(yn**2 +(1.761d0)**2)-1.761d0
         Rp2=
     1       dexp(-u_p-u_m)*(1.1032d0+0.8645d0*un+0.2042d0*up+
     1       0.07937d0*un**2+0.1451d0*un*up+0.01333d0*un**2*up)       
     2       +dexp(-2.d0*u_p)*(-0.1032d0-0.2340d0*u_p+0.06152d0*un*u_p+
     2       0.7533d0*un*u_m-1.007d0*up*u_p)       
        end if
       return
      end
c *********************************************************************
c *********************************************************************
