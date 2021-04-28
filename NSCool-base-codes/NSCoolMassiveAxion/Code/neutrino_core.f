c**********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************

      subroutine numurca_nucl(i,t,qmurca_nucl)

c **** checked partially on March 30, 1993

c**********************************************************************
c Rates and suppression factors from:                                 *
c Yakovlev & Levenfish, A&A 297 (1995): p. 717.                       *
c**********************************************************************
      implicit real*8 (a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      parameter (pi = 3.14159265d0)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
c ***************************
      fexp(x)=dexp(max(x,-7.d2))
c ***************************
c SINGLET PAIRING:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
c Murca_n:  n+n -> n+p+e+nu
      rmurca_n_p1s0(u)=
     1    fexp(3.4370d0-dsqrt((3.4370d0)**2+(1.d0*u)**2))*0.5d0*
     2  ( (0.1477d0+dsqrt((0.8523d0)**2+(0.1175d0*u)**2))**7.5d0 +
     3    (0.1477d0+dsqrt((0.8523d0)**2+(0.1297d0*u)**2))**5.5d0 )
c      rmurca_n_n1s0(u)=rmurca_p_p1s0(u)
      rmurca_n_n1s0(u)=
     1    fexp(5.3390d0-dsqrt((5.3390d0)**2+(2.d0*u)**2))*
     2    (0.2414d0+dsqrt((0.7586d0)**2+(0.1318d0*u)**2))**7.0d0
c Murca_p:  n+p -> p+p+e+nu
      rmurca_p_p1s0(u)=
     1     fexp(5.3390d0-dsqrt((5.3390d0)**2+(2.d0*u)**2))*
     2     (0.2414d0+dsqrt((0.7586d0)**2+(0.1318d0*u)**2))**7.0d0
c      rmurca_p_n1s0(u)=rmurca_n_p1s0(u)
      rmurca_p_n1s0(u)=
     1     fexp(3.4370d0-dsqrt((3.4370d0)**2+(1.d0*u)**2))*0.5d0*
     2   ( (0.1477d0+dsqrt((0.8523d0)**2+(0.1175d0*u)**2))**7.5d0 +
     3     (0.1477d0+dsqrt((0.8523d0)**2+(0.1297d0*u)**2))**5.5d0 )
c ***************************
c TRIPLET PAIRING:
      u_3p2B(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c Murca_p:  n+p -> p+p+e+nu
      rmurca_p_n3p2B(u)=
     1     fexp(2.3980d0-dsqrt((2.3980d0)**2+(1.d0*u)**2))*0.5d0*
     2   ( (0.1612d0+dsqrt((0.8388d0)**2+(0.1117d0*u)**2))**7 +
     3     (0.1612d0+dsqrt((0.8388d0)**2+(0.1274d0*u)**2))**5 )
c Murca_n:  n+n -> n+p+e+nu
      rmurca_n_n3p2B(u,t)=39.1d0*t*fexp(-1.188d0/t)*rmurca_p_n3p2B(u)
c This `rmurca_p_n3p2B' is exact in the limit t<<1 and is only approximate
c when t~<1 since Yakovlev and Levenfish only calculated the t<<1 limit.

c      u_3p2C(t)=dsqrt(1.d0-t**4)/t*
c     1          (5.875d0-1.419d0*t**4+0.500d0*t**8)
c      r_3p2C(u)=1.d0/0.d0
c ******************************
c ******************************
      if (kfn(i)*kfp(i).eq.0.) then
       qmurca_nucl=0.d0
       return
      end if
c ******************************

      rmn=mstn(i)
      rmp=mstp(i)

c Murca_n:  n+n -> n+p+e+nu:
      alpha_n =1.76d0-0.63d0*(1.68d0/kfn(i))**2
      beta_n  =0.68d0
      qmurca_n=8.55d21 *rmn**3*rmp   * 
     1         (kfe(i)/1.68d0 + 
     1          kfm(i)/1.68d0) *
     2         alpha_n * beta_n * (t/1.d9)**8 
c Murca_p:  n+p -> p+p+e+nu:
      alpha_p =alpha_n
      beta_p  =beta_n
c      qmurca_p=8.55d21 *rmn   *rmp**3* 
c     1         (kfe(i)/1.68d0 * (1.d0-kfe(i)/4.d0/kfp(i))+
c     1          kfm(i)/1.68d0 * (1.d0-kfm(i)/4.d0/kfp(i)) ) *
c     2         alpha_p * beta_p * (t/1.d9)**8
      qmurca_p=8.55d21 *rmn   *rmp**3* 
     1         (kfe(i)/1.68d0) * 
     1         (kfe(i)+3.d0*kfp(i)-kfn(i))**2/(8.d0*kfe(i)*kfp(i)) *    ! From Oleg
     2         alpha_p * beta_p * (t/1.d9)**8

c **** effect of superfluidity :

      if(t .lt. tcn(i))then
       if (i.ge.isf) then
        tt=t/tcn(i)
        u=u_1s0(tt)
        rmurca_n_n=rmurca_n_n1s0(u)
        rmurca_p_n=rmurca_p_n1s0(u)
       else
        tt=t/tcn(i)
        u=u_3p2B(tt)
        rmurca_n_n=rmurca_n_n3p2B(u,tt)
        rmurca_p_n=rmurca_p_n3p2B(u)
       end if
      else
       rmurca_n_n=1.0d0
       rmurca_p_n=1.0d0
      end if

      if(t .lt. tcp(i))then
        tt=t/tcp(i)
        u=u_1s0(tt)
        rmurca_n_p=rmurca_n_p1s0(u)
        rmurca_p_p=rmurca_p_p1s0(u)
      else
       rmurca_n_p=1.0d0
       rmurca_p_p=1.0d0
      end if

      rmurca_n=min(rmurca_n_p,rmurca_n_n)
      rmurca_p=min(rmurca_p_p,rmurca_p_n)

      qmurca_n=rmurca_n*qmurca_n
      qmurca_p=rmurca_p*qmurca_p
      qmurca_nucl=qmurca_n+qmurca_p

      return

      end


c**********************************************************************
c**********************************************************************

      subroutine nubrem_crust_nn(i,t,vion,qbrem_nn,qasync,ProcessID)

c **** checked partially on March 30, 1993

c**********************************************************************
c Includes the Levenfish & Yakovlev suppression factors for DURCA,
c and only Boltzmann factors for MURCA.
c**********************************************************************

      implicit real*8 (a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'pid.inc.f'
      INCLUDE 'rho_limits.inc.f'
      parameter (pi = 3.14159265d0)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
      integer :: ProcessID
      real*8 Ias_n

c ***************************
      fexp(x)=dexp(max(x,-7.d2))
c ***************************
c SINGLET PAIRING:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
c Brem_nn:  n+n -> n+n+2nu
      rbrem_nn_n1s0(u)=
     1     (0.1747d0+dsqrt((0.8253d0)**2+(.07933d0*u)**2))**2*
     2     fexp(4.228d0-dsqrt((4.228d0)**2+(4.d0*u)**2))/2.0d0 +
     3     (0.7333d0+dsqrt((0.2667d0)**2+(0.1678d0*u)**2))**7.5d0*
     4     fexp(7.762d0-dsqrt((7.762d0)**2+(9.d0*u)**2))/2.0d0
c ***************************
c TRIPLET PAIRING:
      u_3p2B(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c Brem_nn:  n+n -> n+n+2nu
      rbrem_nn_n3p2B(u)=rbrem_nn_n1s0(u)
c **********************************************
      n_nu=3.d0    ! Number of neutrino famillies !
c Brem_nn:  n+n -> n+n+2nu
      alpha_nn =0.59d0
      beta_nn  =0.56d0
      qbrem_nn=n_nu * 7.4d19 *       mstn(i)**4     * 
     1         (kfn(i)/1.68d0) * alpha_nn * beta_nn * (t/1.d9)**8
c **** effect of superfluidity :
      if(t .lt. tcn(i))then
       if (i.ge.isf) then
        tt=t/tcn(i)
        u=u_1s0(tt)
        rbrem_nn=rbrem_nn_n1s0(u)
       else
        tt=t/tcn(i)
        u=u_3p2B(tt)
        rbrem_nn=rbrem_nn_n3p2B(u)
       end if
      else
       rbrem_nn=1.0d0
      end if
c **** Reduction from ion's volume:
      
      qbrem_nn=rbrem_nn*qbrem_nn * (1.d0-vion)

 



c *************************** _compute_nucleon_star_factors
c all in GeV
      fmG=1.973d-1
      m_pi = 0.140d0
      m_n = 0.939

      star_kfn_crust = kfn(i) * fmG

      PBF_s_n_star_factor = (star_kfn_crust/0.337d0)**3d0
      PBF_p_star_factor = (star_kfn_crust/0.337d0)

      m_x = m_pi / (2.d0*star_kfn_crust) 
      Fx    = 1.d0 - 3.d0/2.d0 * m_x * atan(1.d0/m_x) +
     &        m_x**2.d0 / 2.d0 / (1.d0+x**2)
      eann_star_factor = (star_kfn_crust/0.337d0) * Fx/0.607211d0
 
c *************************** PBF

c 1s0 n
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_s_n = t * sqrt( 1.d0 - tau ) * ( 1.456d0 
     &               - 0.157d0/sqrt(tau) + 1.764/tau ) 
       zn = Delta_T_s_n / t
       Ias_n = (0.158151d0*zn**2d0+0.543166d0*zn**4d0)
     &          *sqrt(1d0+pi*zn/4.d0/0.543166d0**2d0)
     &          *exp(0.0535359d0-sqrt(4d0*zn**2d0+0.0535359d0**2d0))
       PBF_s_n_epsilon = 7.01d14 * (gann/1d-10)**2d0
     &                   * PBF_s_n_star_factor * (t/3d8)**5d0
     &                   * (1d0/m_n)**2d0 * (Ias_n/2.2d-2)
      else
       PBF_s_n_epsilon = 0d0
      endif
      
      if( gann.lt.0 ) then
       PBF_s_n_epsilon = -PBF_s_n_epsilon
      endif

c 3p2 A
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_3p2A = t * sqrt(1d0-tau)*(0.7893d0 + 1.118d0/tau)
       zn = Delta_T_3p2A / t
       IanPA = IpnA_interp(zn)
       PBF_pA_epsilon = 1.67d15 * (gann/1d-10)**2d0 * PBF_p_star_factor
     &                * (t/3d8)**5d0 * (IanPA/5.96d-3)
      else
       PBF_pA_epsilon = 0d0
      endif

      if( gann.lt.0 ) then
       PBF_pA_epsilon = -PBF_pA_epsilon
      endif

c 3p2 B
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_3p2B = t * sqrt(1d0-tau**4d0)/tau
     &              *( 2.03d0 - 0.4903d0*tau**4d0 + 0.1727d0*tau**8d0 )
       zn = Delta_T_3p2B / t
       IanPB = IpnB_interp(zn)
       PBF_pB_epsilon = 1.67d15 * (gann/1d-10)**2d0 * PBF_p_star_factor
     &                * (t/3d8)**5d0 * (IanPB/5.96d-3)
      else
       PBF_pB_epsilon = 0d0
      endif

      if( gann.lt.0 ) then
       PBF_pB_epsilon = -PBF_pB_epsilon
      endif

c *************************** _do_nucelon
c in erg/cm^3/s
      qabrem_nn = 1.827e12 * eann_star_factor * (t/1.d8)**6d0
     &            * (gann/1d-10)**2d0 * mstn(i)**2d0
      qabrem_nn_super = qabrem_nn * rbrem_nn

      if( gann.lt.0 ) then
       qabrem_nn = -qabrem_nn
       qabrem_nn_super = -qabrem_nn_super
      endif

c ***************************

      qasync = 0d0

      if (IAND(pid_nn_inner_crust,ProcessID).gt.0) then
       qasync = qasync + qabrem_nn
      endif       
      if (IAND(pid_nn_inner_crust_super,ProcessID).gt.0) then
       qasync = qasync + qabrem_nn_super
      endif      
      if (IAND(pid_PBF_s_n_inner_crust,ProcessID).gt.0) then
       qasync = qasync + PBF_s_n_epsilon
      endif
      if (IAND(pid_PBF_pA_inner_crust,ProcessID).gt.0) then
       qasync = qasync + PBF_pA_epsilon
      endif
      if (IAND(pid_PBF_pB_inner_crust,ProcessID).gt.0) then
       qasync = qasync + PBF_pB_epsilon
      endif
 

c **** Reduction from ion's volume:
      qasync = qasync * (1.d0-vion)

      return

      end


c**********************************************************************
c**********************************************************************

      subroutine nubrem_nucl(i,t,time,qbrem_nucl,qasync,ProcessID)

c **** checked partially on March 30, 1993

c**********************************************************************
c Includes the Levenfish & Yakovlev suppression factors for DURCA,
c and only Boltzmann factors for MURCA.
c**********************************************************************

      implicit real*8 (a-h,k-z)

      INCLUDE 'size.inc.f'
      INCLUDE 'pid.inc.f'
      INCLUDE 'rho_limits.inc.f'
      parameter (pi = 3.14159265d0)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
      integer :: ProcessID
      real*8 Ias_p,Ias_n,IanPA,IanPB
      
c      common/h_tables/h_bfield,h_temp,h_pfermi,h_emissivity

c ***************************
      fexp(x)=dexp(max(x,-7.d2))
c ***************************
c SINGLET PAIRING:
      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
c Brem_nn:  n+n -> n+n+2nu
      rbrem_nn_p1s0(u)=1.0d0    ! Not affected by proton pairing !
      rbrem_nn_n1s0(u)=
     1     (0.1747d0+dsqrt((0.8253d0)**2+(.07933d0*u)**2))**2*
     2     fexp(4.228d0-dsqrt((4.228d0)**2+(4.d0*u)**2))/2.0d0 +
     3     (0.7333d0+dsqrt((0.2667d0)**2+(0.1678d0*u)**2))**7.5d0*
     4     fexp(7.762d0-dsqrt((7.762d0)**2+(9.d0*u)**2))/2.0d0
c Brem_np:  n+p -> n+p+2nu
      rbrem_np_p1s0(u)=
     1     (0.9982d0+dsqrt((0.0018d0)**2+(0.3815d0*u)**2))**1*
     2     fexp(1.306d0-dsqrt((1.306d0)**2+(1.d0*u)**2))/2.732d0 +
     3     (0.3949d0+dsqrt((0.6051d0)**2+(0.2666d0*u)**2))**7*
     4     fexp(3.303d0-dsqrt((3.303d0)**2+(4.d0*u)**2))/1.577d0
      rbrem_np_n1s0(u)=rbrem_np_p1s0(u)
c Brem_pp:  p+p -> p+p+2nu
      rbrem_pp_p1s0(u)=rbrem_nn_n1s0(u)
      rbrem_pp_n1s0(u)=1.0d0   ! Not affected by neutron pairing !
c ***************************
c TRIPLET PAIRING:
      u_3p2B(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c Brem_nn:  n+n -> n+n+2nu
      rbrem_nn_n3p2B(u)=rbrem_nn_n1s0(u)
c Brem_np:  n+p -> n+p+2nu
      rbrem_np_n3p2B(u)=rbrem_np_n1s0(u)
c Brem_pp:  p+p -> p+p+2nu
      rbrem_pp_n3p2B(u)=1.0  ! Not affected by neutron pairing !
c ***************************
c ****** murca : **********************************************
      n_nu=3.d0    ! Number of neutrino famillies !
c Brem_nn:  n+n -> n+n+2nu
      alpha_nn =0.59d0
      beta_nn  =0.56d0
      qbrem_nn=n_nu * 7.4d19 *       mstn(i)**4     * 
     1         (kfn(i)/1.68d0) * alpha_nn * beta_nn * (t/1.d9)**8
c Brem_np:  n+p -> n+p+2nu
      alpha_np =1.06d0
      beta_np  =0.66d0
      qbrem_np=n_nu * 1.5d20 * mstn(i)**2*mstp(i)**2 * 
     1         (kfp(i)/1.68d0) * alpha_np * beta_np * (t/1.d9)**8
c Brem_pp:  p+p -> p+p+2nu
      alpha_pp=0.11d0
      beta_pp=0.7d0
      qbrem_pp=n_nu * 7.4d19 *       mstp(i)**4      * 
     1         (kfp(i)/1.68d0) * alpha_pp * beta_pp * (t/1.d9)**8
   
c *************************** superfluidity
      if(t .lt. tcn(i))then
       if (i.ge.isf) then
        tt=t/tcn(i)
        u=u_1s0(tt)
        rbrem_nn_n=rbrem_nn_n1s0(u)
        rbrem_np_n=rbrem_np_n1s0(u)
        rbrem_pp_n=rbrem_pp_n1s0(u)
       else
        tt=t/tcn(i)
        u=u_3p2B(tt)
        rbrem_nn_n=rbrem_nn_n3p2B(u)
        rbrem_np_n=rbrem_np_n3p2B(u)
        rbrem_pp_n=rbrem_pp_n3p2B(u)
       end if
      else
       rbrem_nn_n=1.0d0
       rbrem_np_n=1.0d0
       rbrem_pp_n=1.0d0
      end if

      if(t .lt. tcp(i))then
        tt=t/tcp(i)
        u=u_1s0(tt)
        rbrem_nn_p=rbrem_nn_p1s0(u)
        rbrem_np_p=rbrem_np_p1s0(u)
        rbrem_pp_p=rbrem_pp_p1s0(u)
      else
       rbrem_nn_p=1.0d0
       rbrem_np_p=1.0d0
       rbrem_pp_p=1.0d0
      end if

      rbrem_nn=min(rbrem_nn_p,rbrem_nn_n)
      rbrem_np=min(rbrem_np_p,rbrem_np_n)
      rbrem_pp=min(rbrem_pp_p,rbrem_pp_n)
 

c **** effect of superfluidity :

      qbrem_nn=rbrem_nn*qbrem_nn
      qbrem_np=rbrem_np*qbrem_np
      qbrem_pp=rbrem_pp*qbrem_pp

      qbrem_nucl=qbrem_nn+qbrem_np+qbrem_pp








c *************************** 

 
c      efac=303d-3 
c      K2GeV=8.617d-14
c      gaee=1d-10
c      GeVtoem=3.155d62
c      GtoG2=1.95d-20
c      Bfie=4.41d13*GtoG2
c      invfm2GeV=1.973d-1
c      Bfield_G = 2.0d13
c      gamm = 0d-9
     

c *************************** _compute_nucleon_star_factors
c all in GeV
      g_c = gapp + gann
      h_c = gapp - gann

      fmG=1.973d-1
      m_pi = 0.140d0
      m_n = 0.939d0
      K2keV=8.617d-8
      star_kfn_core = kfn(i) * fmG
      star_kfp_core = kfp(i) * fmG

      m_x = m_pi / (2.d0*star_kfn_core) 
      m_y = m_pi / (2.d0*star_kfp_core) 
      xyp = 2.d0 * m_x * m_y /(m_x+m_y)
      xym = 2.d0 * m_x * m_y /(m_x-m_y)

      Fx    = 1.d0 - 3.d0/2.d0 * m_x * atan(1.d0/m_x) +
     &        m_x**2.d0 / 2.d0 / (1.d0+x**2)
      Fy    = 1.d0 - 3.d0/2.d0 * m_y * atan(1.d0/m_y) +
     &        m_y**2.d0 / 2.d0 / (1.d0+y**2)
      Fxyp  = 1.d0 - 3.d0/2.d0 * xyp * atan(1.d0/xyp) +
     &        m_y**2.d0 / 2.d0 / (1.d0+xyp**2)
      Fxym  = 1.d0 - 3.d0/2.d0 * xym * atan(1.d0/xym) +
     &        m_y**2.d0 / 2.d0 / (1.d0+xym**2)

      eann_star_factor = (star_kfn_core/0.337d0) * Fx/0.607211d0
      eapp_star_factor = (star_kfp_core/0.337d0) * Fy/0.607211d0
                
      gfacg = 0.5d-1*Fy+       ( (Fxyp + Fxym) +m_y/m_x*(Fxyp-Fxym) )
     &     + (1.d0-m_y*atan(1.d0/m_y))
      gfach = 0.5d-1*Fy+0.5d-1*( (Fxyp + Fxym) +m_y/m_x*(Fxyp-Fxym) )
     &     + (1.d0- m_y*atan(1.d0/m_y))

      eanp_star_factor_g = (star_kfn_core/0.337d0) * gfacg
      eanp_star_factor_h = (star_kfn_core/0.337d0) * gfach

      PBF_s_p_star_factor = (star_kfp_core/0.337d0)**3d0
      PBF_s_n_star_factor = (star_kfn_core/0.337d0)**3d0
      PBF_p_star_factor = (star_kfn_core/0.337d0)


c *************************** PBF

c 1s0 p
      if(t .lt. tcp(i))then
       tau = t/tcp(i)
       Delta_T_s_p = t * sqrt( 1.d0 - tau ) * ( 1.456d0 
     &               - 0.157d0/sqrt(tau) + 1.764/tau ) 
       zn = Delta_T_s_p / t
       Ias_p = (0.158151d0*zn**2d0+0.543166d0*zn**4d0)
     &          *sqrt(1d0+pi*zn/4.d0/0.543166d0**2d0)
     &          *exp(0.0535359d0-sqrt(4d0*zn**2d0+0.0535359d0**2d0))
       PBF_s_p_epsilon = 7.01d14 * (gapp/1d-10)**2d0
     &                   * PBF_s_p_star_factor * (t/3d8)**5d0
     &                   * (1d0/mstp(i))**2d0 * (Ias_p/2.2d-2)
      else
       PBF_s_p_epsilon = 0d0
      endif

      if( gapp.lt.0 ) then
       PBF_s_p_epsilon = -PBF_s_p_epsilon
      endif

c 1s0 n
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_s_n = t * sqrt( 1.d0 - tau ) * ( 1.456d0 
     &               - 0.157d0/sqrt(tau) + 1.764d0/tau ) 
       zn = Delta_T_s_n / t
       Ias_n = (0.158151d0*zn**2d0+0.543166d0*zn**4d0)
     &          *sqrt(1d0+pi*zn/4.d0/0.543166d0**2d0)
     &          *exp(0.0535359d0-sqrt(4d0*zn**2d0+0.0535359d0**2d0))
       PBF_s_n_epsilon = 7.01d14 * (gann/1d-10)**2d0
     &                   * PBF_s_n_star_factor * (t/3d8)**5d0
     &                   * (1d0/m_n)**2d0 * (Ias_n/2.2d-2)
      else
       PBF_s_n_epsilon = 0d0
      endif


      if( gann.lt.0 ) then
       PBF_s_n_epsilon = -PBF_s_n_epsilon
      endif

c 3p2 A
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_3p2A = t * sqrt(1d0-tau)*(0.7893d0 + 1.118d0/tau)
       zn = Delta_T_3p2A / t
       IanPA = IpnA_interp(zn)
       PBF_pA_epsilon = 1.67d15 * (gann/1d-10)**2d0 * PBF_p_star_factor
     &                * (t/3d8)**5d0 * (IanPA/5.96d-3)
      else
       PBF_pA_epsilon = 0d0
      endif

      if( gann.lt.0 ) then
       PBF_pA_epsilon = -PBF_pA_epsilon
      endif

c 3p2 B
      if(t .lt. tcn(i))then
       tau = t/tcn(i)
       Delta_T_3p2B = t * sqrt(1d0-tau**4d0)/tau
     &              *( 2.03d0 - 0.4903d0*tau**4d0 + 0.1727d0*tau**8d0 )
       zn = Delta_T_3p2B / t
       IanPB = IpnB_interp(zn)
       PBF_pB_epsilon = 1.67d15 * (gann/1d-10)**2d0 * PBF_p_star_factor
     &                * (t/3d8)**5d0 * (IanPB/5.96d-3)
      else
       PBF_pB_epsilon = 0d0
      endif

      if( gann.lt.0 ) then
       PBF_pB_epsilon = -PBF_pB_epsilon
      endif


c *************************** _do_nucelon
c in erg/cm^3/s
      qabrem_nn = 1.827e12 * eann_star_factor * (t/1.d8)**6d0
     &            * (gann/1d-10)**2d0 * mstn(i)**2d0
      qabrem_pp = 1.827e12 * eapp_star_factor * (t/1.d8)**6d0
     &            * (gapp/1d-10)**2d0 * mstp(i)**2d0
      qabrem_np = 2.008d12 * (eanp_star_factor_h * h_c**2d0 
     &            + eanp_star_factor_g * g_c**2d0 )/(1d-10)**2d0
     &            * (t/1.d8)**6d0 * mstn(i)**2d0 

      qabrem_nn_super = qabrem_nn * rbrem_nn
      qabrem_pp_super = qabrem_pp * rbrem_pp
      qabrem_np_super = qabrem_np * rbrem_np

      if( gann.lt.0 ) then
       qabrem_nn = -qabrem_nn 
       qabrem_nn_super = -qabrem_nn_super
      endif

      if( gapp.lt.0 ) then
       qabrem_pp = -qabrem_pp
       qabrem_pp_super = -qabrem_pp_super
      endif

      if( gapp.lt.0 .or. gann.lt.0) then
       qabrem_np = -qabrem_np 
       qabrem_np_super = -qabrem_np_super
      endif
c *************************** leptons

      alphaEM = 1.d0/137.d0
      alpha_e = gaee**2.d0/(4.d0*pi)
      alpha_m = gamm**2.d0/(4.d0*pi)

c GeV
      me = 5.11e-4
      mm = 0.105658d0
      mu = 0.9315d0
      rhoGeV = mu * star_kfp_core**3.d0 / 3.d0 / pi**2.d0
c              rrho(i) * 4.31013e-18
      TGeV = t * 8.61733e-14

      GammaI = 22.73*1.d6/t * ( rhoGeV/4.3103d-18 / 1.d6 )**(1.d0/3.d0)
      xState = LOG10( rhoGeV/4.3103d-18 )


      if( xState.gt.11.4d0 ) then
       a = -6.47808d0
       b = 0.068645d0
       c = -0.000252677d0
      else
       a = 0.21946d0
       b = 0.00287263d0
       c = -0.000142016d0
      endif

      if( GammaI.gt.178 ) then
       u = 0.488049d0 + 1.25585d0*GammaI/1.d3 - 
     &     0.743902d0*(GammaI/1.d3)**(2.d0)
      else
       u = 0.672409d0 + 0.182774d0*GammaI/1.d3 +
     &     0.144817d0*(GammaI/1.d3)**(2.d0)
      endif

      Fep = 10.d0**( a + b*xState**2.d0 + c*xState**4.d0 - 1.d0+u )
      Fem = Fep * mm**2.d0 / me**2.d0

      qabrem_e = 3.168d62 * pi**2.d0 / 15.d0 * alphaEM**2.d0 * alpha_e 
     &         * rhoGeV * TGeV**4.d0 / me**2.d0 / mu * Fep
      qabrem_m = 3.168d62 * pi**2.d0 / 15.d0 * alphaEM**2.d0 * alpha_m 
     &         * rhoGeV * TGeV**4.d0 / mm**2.d0 / mu * Fem

      if( gaee.lt.0 ) then
       qabrem_e = -qabrem_e 
      endif

      if( gamm.lt.0 ) then
       qabrem_m = -qabrem_m
      endif

c      write(*,*) qabrem_e, GammaI, xState, Fep, u

c *************************** synchotron
 
c      Temp_keV = t * K2keV
c      pFermi_GeV = kfm(i) * fmG
c      if ((Temp_keV>1).and.(Temp_keV<900)
c     & .and.(pFermi_GeV>1d-3).and.(pFermi_GeV<3d-1)) then
c       qasync_o = gamm*gamm*emissivity(Bfield_G, Temp_keV, pFermi_GeV)
c       write(6,*)gamm,Bfield_G,Temp_keV,pFermi_GeV,qasync
c      else
c       qasync_o=0d0
c      end if
c      if (gamm < 0) then
c       qasync_o = -qasync_o
c      end if



     

c ***************************
      qasync = 0d0
      
      if (IAND(pid_synchotron,ProcessID).gt.0) then
       qasync = qasync + qasync_0
      endif
      if (IAND(pid_nn_core,ProcessID).gt.0) then
       qasync = qasync + qabrem_nn 
      endif
      if (IAND(pid_pp_core,ProcessID).gt.0) then
       qasync = qasync + qabrem_pp
      endif
      if (IAND(pid_np_core,ProcessID).gt.0) then
       qasync = qasync + qabrem_np
      endif
      if (IAND(pid_nn_core_super,ProcessID).gt.0) then
       qasync = qasync + qabrem_nn_super
      endif
      if (IAND(pid_pp_core_super,ProcessID).gt.0) then
       qasync = qasync + qabrem_pp_super
      endif
      if (IAND(pid_np_core_super,ProcessID).gt.0) then
       qasync = qasync + qabrem_np_super
      endif
      if (IAND(pid_PBF_s_p_core,ProcessID).gt.0) then
       qasync = qasync + PBF_s_p_epsilon
      endif
      if (IAND(pid_PBF_s_n_core,ProcessID).gt.0) then
       qasync = qasync + PBF_s_n_epsilon
      endif
      if (IAND(pid_PBF_pA_core,ProcessID).gt.0) then
       qasync = qasync + PBF_pA_epsilon
      endif
      if (IAND(pid_PBF_pB_core,ProcessID).gt.0) then
       qasync = qasync + PBF_pB_epsilon
      endif
      if (IAND(pid_ep_core,ProcessID).gt.0) then
       qasync = qasync + qabrem_e
      endif
      if (IAND(pid_mp_core,ProcessID).gt.0) then
       qasync = qasync + qabrem_m
      endif
      
c      write(*,*)'coupling:',gann,gapp,gaee,gamm,qasync,qbrem_nucl

      return

      end


c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c                     H Y P E R O N S
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine numurca_hyp(i,t,qmurca_hyp)
       implicit real*8(a-h,k-z)
        qmurca_hyp=0.d0
      return
      end
c *********************************************************************
c *********************************************************************
      subroutine nubrem_hyp(i,t,qbrem_hyp)
       implicit real*8(a-h,k-z)
        qbrem_hyp=0.0d0
       return
      end
c *********************************************************************
c *********************************************************************

      subroutine nudurca_h(i,t,rho,
     1                     qdurca_np,qdurca_lap,
     2                     qdurca_smn,qdurca_smla,qdurca_sms0)

c **** checked partially on March 30, 1993

c**********************************************************************
c    From Levenfish & Yakovlev, Astron. Lett. 20 (1994), p. 43        *
c**********************************************************************

      implicit real*8 (a-h,k-z)
      parameter (pi = 3.14159265d0)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
c ***************************
      fexp(x)=dexp(max(x,-7.d2))

      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
      u_3p2B(t)=dsqrt(1.d0-t)*(0.7893d0+1.188d0/t)
c      u_3p2C(t)=dsqrt(1.d0-t**4)/t*
c     1          (5.875d0-1.419d0*t**4+0.500d0*t**8)

      r_1s0(u)=(0.2312d0+dsqrt(0.7688d0**2+(0.1438d0*u)**2))**5.5d0*
     1         fexp(3.427d0-dsqrt(3.427d0**2+u**2))
      r_3p2B(u)=(0.2546d0+dsqrt(0.7454d0**2+(0.1284d0*u)**2))**5*
     1         fexp(2.701d0-dsqrt(2.701d0**2+u**2))
      r_1s0_1s0(u1,u2)=func_r_1s0_1s0(u1,u2)
      r_1s0_3p2B(t1,t2)=func_r_1s0_3P2B(t1,t2) ! WARNING: this uses t1 & t2
                                               ! instead of u1 & u2
c      r_3p2C(u)=1.d0/0.d0

c****************************

      rmn = mstn(i)
      rmp = mstp(i)
      rmla=mstla(i)
      rmsm=mstsm(i)
      rms0=msts0(i)
      rmsp=mstsp(i)

      rate_np  =1.0000d0
      rate_lap =0.0394d0
      rate_smn =0.0125d0
      rate_smla=0.2055d0
      rate_sms0=0.6052d0

c **** n-p:
      if     (idurca_np(i).eq.1)then
       qdurca_np= rate_np * 4.24d27*rmn*rmp*(t/1.d9)**6*
     1            (abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)
      else if(idurca_np(i).eq.2)then
       qdurca_np= rate_np * 4.24d27*rmn*rmp*(t/1.d9)**6*
     1           ((abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)+
     2            (abs( ymuon(i))*bar(i)/0.16d0)**(1.d0/3.d0))
      else
       qdurca_np= 0.d0
      end if
c Pairing suppression:
      if ((t.gt.tcn(i)).and.(t.gt.tcp(i))) then
       r_np=1.d0
      else if ((t.gt.tcn(i)).and.(t.le.tcp(i))) then
       tt=t/tcp(i)
       u=u_1s0(tt)
       r_np=r_1s0(u)
      else if ((t.le.tcn(i)).and.(t.gt.tcp(i))) then
       if (i.ge.isf) then 
        tt=t/tcn(i)
        u=u_1s0(tt)
        r_np=r_1s0(u)
       else
        tt=t/tcn(i)
        u=u_3p2B(tt)
        r_np=r_3p2B(u)
       end if
      else
       tt1=t/tcp(i)
       u1=u_1s0(tt1)
       if (i.ge.isf) then 
        tt2=t/tcn(i)
        u2=u_1s0(tt2)
        r_np=r_1s0_1s0(u1,u2)
       else
        tt2=t/tcn(i)
        u2=u_3p2B(tt2)   ! Not needed for r_1s0_3p2B
        r_np=r_1s0_3p2B(tt1,tt2)
       end if
      end if
c **** la-p:
      if     (idurca_lap(i).eq.1)then
       qdurca_lap= rate_lap * 4.24d27*rmla*rmp*(t/1.d9)**6*
     1             (abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)
      else if(idurca_lap(i).eq.2)then
       qdurca_lap= rate_lap * 4.24d27*rmla*rmp*(t/1.d9)**6*
     1            ((abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)+
     2             (abs( ymuon(i))*bar(i)/0.16d0)**(1.d0/3.d0))
      else
       qdurca_lap= 0.d0
      end if
c **** sm-n:
      if     (idurca_smn(i).eq.1)then
       qdurca_smn= rate_smn * 4.24d27*rmn*rmsm*(t/1.d9)**6*
     1             (abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)
      else if(idurca_smn(i).eq.2)then
       qdurca_smn= rate_smn * 4.24d27*rmn*rmsm*(t/1.d9)**6*
     1            ((abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)+
     2             (abs( ymuon(i))*bar(i)/0.16d0)**(1.d0/3.d0))
      else
       qdurca_smn= 0.d0
      end if
c **** sm-la:
      if     (idurca_smla(i).eq.1)then
       qdurca_smla= rate_smla * 4.24d27*rmsm*rmla*(t/1.d9)**6*
     1              (abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)
      else if(idurca_smla(i).eq.2)then
       qdurca_smla= rate_smla * 4.24d27*rmsm*rmla*(t/1.d9)**6*
     1             ((abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)+
     2              (abs( ymuon(i))*bar(i)/0.16d0)**(1.d0/3.d0))
      else
       qdurca_smla= 0.d0
      end if
c **** sm-s0:
      if     (idurca_sms0(i).eq.1)then
       qdurca_sms0= rate_sms0 * 4.24d27*rmsm*rms0*(t/1.d9)**6*
     1            (abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)
      else if(idurca_sms0(i).eq.2)then
       qdurca_sms0= rate_sms0 * 4.24d27*rmsm*rms0*(t/1.d9)**6*
     1           ((abs(yelect(i))*bar(i)/0.16d0)**(1.d0/3.d0)+
     2            (abs( ymuon(i))*bar(i)/0.16d0)**(1.d0/3.d0))
      else
       qdurca_sms0= 0.d0
      end if

c **** effect of superfluidity :

c (Note : isf is the zone where neutron pairing shifts from 3P2 to 1S0)

c Neutron pairing suppression:
      if(t.lt.tcn(i)) then
       if (i.ge.isf) then 
        tt=t/tcn(i)
        u=u_1s0(tt)
        rn=r_1s0(u)
       else
        tt=t/tcn(i)
        u=u_3p2B(tt)
        rn=r_3p2B(u)
       end if
      else
       rn=1.0 d0
      end if
c Proton pairing suppression:
      if(t .lt. tcp(i))then
       tt=t/tcp(i)
       u=u_1s0(tt)
       rp=r_1s0(u)
      else
       rp=1.0d0
      end if
c Lambda pairing suppression:
      if(t .lt. tcla(i))then
       tt=t/tcla(i)
       u=u_1s0(tt)
       rla=r_1s0(u)
      else
       rla=1.0d0
      end if
c Sigma- pairing suppression:
      rsm=1.0d0
c Sigma0 pairing suppression:
      rs0=1.0d0
c Sigma+ pairing suppression:
      rsp=1.0d0
c****
c      r_np  =min(rn,rp)
      r_lap =min(rla,rp)
      r_smn =min(rsm,rn)
      r_smla=min(rsm,rla)
      r_sms0=min(rsm,rs0)

      qdurca_np  =r_np  *qdurca_np
      qdurca_lap =r_lap *qdurca_lap
      qdurca_smn =r_smn *qdurca_smn
      qdurca_smla=r_smla*qdurca_smla
      qdurca_sms0=r_sms0*qdurca_sms0

c      qdurca =qdurca_np  + qdurca_lap +
c     2        qdurca_smn + qdurca_smla+ qdurca_sms0

      return

      end

c *********************************************************************
c *********************************************************************
      function func_r_1s0_1s0(v1,v2)
c CHECKED ON MARCH 7, 2001 against Fig. 2 of L&Y Paper
       implicit real*8 (a-h,k-z)
       parameter (pi = 3.14159265d0,gamma=5040.d0/457.d0/pi**6)
        u=v1**2+v2**2
        w=v1**2-v2**2
        u1=1.8091d0+dsqrt(v1**2+2.2476d0**2)
        u2=1.8091d0+dsqrt(v2**2+2.2476d0**2)
        p=(u+12.421d0+dsqrt(w**2+16.350d0*u+45.171d0))/2.d0
        q=(u+12.421d0-dsqrt(w**2+16.350d0*u+45.171d0))/2.d0
        ps=(u+dsqrt(w**2+5524.8d0*u+6.7737d0))/2.d0
        pe=(u+0.43847d0+dsqrt(w**2+8.3680d0*u+491.32d0))/2.d0
        D=(u1*u2)**1.5/(2.d0*4.0567d0**5)*(u1**2+u2**2)*
     1    dexp(-u1-u2+8.1134d0)
        K0=dsqrt(p-q)/120.d0*(6.d0*p**2+83.d0*p*q+16.d0*q**2)-
     1     dsqrt(p)*q/8.d0*(4.d0*p+3.d0*q)*
     2     dlog((dsqrt(p)+dsqrt(p-q))/dsqrt(q))
        K1=pi**2*dsqrt(p-q)/6.d0*(p+2.d0*q)-
     1     pi**2/2.d0*q*dsqrt(p)*
     2     dlog((dsqrt(p)+dsqrt(p-q))/dsqrt(q))
        K2=7.d0*pi**4/60.d0*dsqrt(p-q)
        S=gamma*(K0+K1+0.42232d0*K2)*dsqrt(pi/2.d0)*ps**0.25*
     1    dexp(-dsqrt(pe))
        
        func_r_1s0_1s0=u/(u+0.9163d0)*S+D

       return
      end
c *********************************************************************
c *********************************************************************
c CHECKED ON MARCH 8, 2001 against Fig. 3 of L&Y Paper
      function func_r_1s0_3p2B(t1,t2)
       implicit real*8 (a-h,k-z)
c       INCLUDE 'Base_Dir.inc.f'
       dimension lgtau1(35),lgtau2(35),lgr(35,35),lgr2(35,35)
       save lgtau1,lgtau2,lgr,lgr2
c************************************************
c Notice: the table must be read backward so that lgtau are
c ordered in increasing values, for the Slpine interpolation
       if (tread.ne.123.) then
        open(unit=33,file='Code/Data_Files/SF_suppression.dat',
     x       status='old')
         read(33,*)(lgtau1(i),i=35,1,-1)
         read(33,*)(lgtau2(i),i=35,1,-1)
         do i=35,1,-1
          read(33,*)(lgr(i,j),j=35,1,-1)
         end do
        close(unit=33,status='keep')
        call spline2(lgtau1,lgtau2,lgr,35,35,lgr2)
        tread=123.
       end if
c************************************************
        lt1=log10(t1)
        lt2=log10(t2)
        call splint2(lgtau1,lgtau2,lgr,lgr2,35,35,lt1,lt2,lr)
        func_r_1s0_3p2B=10.d0**lr
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccc
        lt=sqrt(lt1**2+lt2**2)
        lt_limit=3.
        if (lt.le.lt_limit) then
         call splint2(lgtau1,lgtau2,lgr,lgr2,35,35,lt1,lt2,lr)
         func_r_1s0_3p2B=10.d0**lr
        else
         lt1=lt1/lt
         lt2=lt2/lt
         call splint2(lgtau1,lgtau2,lgr,lgr2,35,35,lt1,lt2,lr)
         func_r_1s0_3p2B=10.d0**lr * exp(-lt/lt_limit)
        end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
      end

c *********************************************************************
c *********************************************************************

      subroutine nufast(i,t,rho,qfast)

c **** checked partially on March 30, 1993

c**********************************************************************
c Includes the Levenfish & Yakovlev suppression factors for DURCA,
c**********************************************************************

      implicit real*8 (a-h,k-z)
      INCLUDE 'size.inc.f'
      INCLUDE 'rho_limits.inc.f'
      parameter (pi = 3.14159265d0)
        
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
c ***************************
      fexp(x)=dexp(max(x,-7.d2))

      u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
      r_1s0(u)=(0.2312d0+dsqrt(0.7688d0**2+(0.1438d0*u)**2))**5.5d0*
     1         fexp(3.427d0-dsqrt(3.427d0**2+u**2))

      u_3p2B(t)=dsqrt(1.-t)*(5.596d0+8.424d0/t)
      r_3p2B(u)=(0.2546d0+dsqrt(0.7454d0**2+(0.01811d0*u)**2))**5*
     1         fexp(2.701d0-dsqrt(2.701d0**2+u**2/(16.d0*pi)))

c      u_3p2C(t)=dsqrt(1.d0-t**4)/t*
c     1          (5.875d0-1.419d0*t**4+0.500d0*t**8)
c      r_3p2C(u)=d01./0.d0

c*******************************

      rmn=mstn(i)
      rmp=mstp(i)

      u=bar(i)/0.16d0
      ratio=.319d0/(abs(yelect(i))*u)**(1.d0/3.d0)
      zero=0.d0
      f=max(zero,1-ratio**2)
      f=dsqrt(f)

c **** Kaon urca: 

      if (theta_k(i).ne.0.) then
       g_a=1.0d0
       mu_el=kfe(i)*197.d0
       qkaon= 5.d0/4.d0*dsin(theta_k(i))**2 * sin(0.223d0)**2 *
     1        2.21d26 * mstn(i) * mstp(i) * (mu_el/100.d0) *
     2        (1.d0+3.d0*g_a**2) * (t/1.d9)**6
c From Thorsson et al, Phys. Rev. D52, p. 3739, 1995
      else
       qkaon=0.d0
      end if

c **** eurca :

      if(rho.ge.rhoexo)then
       qexo=cexo*(rho/2.8d14)**(2.d0/3.d0)*(t/1.d9)**pexo
c       rhoexo2=1.4d25
c       if (rho.ge.rhoexo2) then
c        qexo=101.d0*qexo
c       end if
      else 
       qexo=0.d0
      end if

c **** effect of superfluidity :

c (Note : isf is the zone where neutron pairing shifts from 3P2 to 1S0)

      if( (t.lt.tcp(i)) .and. (t.lt.tcn(i)) )then
       if (i.ge.isf) then 
        un=u_1s0(t/tcn(i))
        rn=r_1s0(un)
       else
        un=u_3p2B(t/tcn(i))
        rn=r_3p2B(un)
       end if
       up=u_1s0(t/tcp(i))
       rp=r_1s0(up)
       r=min(rn,rp)
      else if(t .lt. tcn(i)) then
       if (i.ge.isf) then 
        u=u_1s0(t/tcn(i))
        r=r_1s0(u)
       else
        u=u_3p2B(t/tcn(i))
        r=r_3p2B(u)
       end if
      else if(t .lt. tcp(i))then
       u=u_1s0(t/tcp(i))
       r=r_1s0(u)
      else
       r=1.d0
      end if

      qkaon=qkaon*r

      r_exo=r
      if ((pexosn.eq.0.).and.(pexosp.eq.0.)) r_exo=1.
      qexo =qexo *r_exo

      qfast=qkaon+qexo
      return

      end

c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine nu_1s0_pbf(T,Tc,mst,kf,q_1s0_pbf)
c     This subroutine uses only the Axial part: good for n & p !
c     From : Neutrino emission due to proton pairing in neutron stars
c            Kaminker, A. D.; Haensel, P.; Yakovlev, D. G.
c            1999A&A...345L..14K
c     Vector part is put to zero according to:
c     "Vector current conservation and neutrino emission from 
c     singlet-paired baryons in neutron stars"
c     Leinson, L.B., & Perez, A.
c     2006PhLB..638..114L
       implicit real*8 (a-h,k-z)
        if (T.le.Tc) then
         pf=kf*197.d0
         vf=pf/mst/940.d0
         a_v=0.d0
         a_a=1.60d0*vf**2*(mst**2+11.d0/42.d0)
         a=a_v+a_a
c        The vector part in "a" has been put to zero !
         tau=T/Tc
         u=dsqrt(1.d0-tau)*(1.456d0-0.157d0/dsqrt(tau)+1.764d0/tau)
         q_1s0_pbf=1.170d21*mst**2*vf*(T/1.d9)**7*3.d0*a*
     1              control_pbf_1S0(u)
        else
         q_1s0_pbf=0.0d0
        end if
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine nu_n3p2_B_pbf(T,Tc,mst,kf,q_n3p2_pbf)
c     This subroutine uses only the Axial part from
c     Neutrino emission due to Cooper pairing of nucleons in cooling neutron stars
c     Yakovlev, D. G.; Kaminker, A. D.; Levenfish, K. P.
c     1999A&A...343..650Y
c     Vector part is put to zero according to:
c     Neutrino emission due to Cooper pairing in neutron stars
c     Leinson, L.B., & Perez, A.
c     2006astro.ph..6653L
c
c     For the j=2 m=0 gap
       implicit real*8 (a-h,k-z)
       parameter (pi = 3.14159265d0,rhonucl=2.8d14)
        if (T.le.Tc) then
         pf=kf*197.d0
         vf=pf/mst/940.d0
         g_A=1.26d0
         a_v=0.d0
         a_a=2.d0*g_A**2
         a=a_v+a_a
c        The vector part in "a" has been put to zero !
         tau=T/Tc
         u=dsqrt(1.d0-tau)*(0.7893d0+1.764d0/tau)
         q_n3p2_pbf=1.170d21*mst**2*vf*(T/1.d9)**7*3.d0*a*
     1              control_pbf_3P2_B(u)
        else
         q_n3p2_pbf=0.0d0
        end if
       return
      end
c *********************************************************************
c *********************************************************************
      function control_pbf_1S0(v)
       implicit real*8 (a-h,k-z)
        x=0.602d0*v**2+0.5942d0*v**4+0.288d0*v**6
        y=dsqrt(0.5547d0+dsqrt(0.4453d0**2+0.01130*v**2))
        z=dexp(-dsqrt(4.d0*v**2+2.245d0**2)+2.245d0)
        control_pbf_1S0=x*y*z
       return
      end
c *********************************************************************
c *********************************************************************
      function control_pbf_3P2_B(v)
       implicit real*8 (a-h,k-z)
        x=(1.204d0*v**2+3.733d0*v**4+0.3191*v**6)/(1.d0+0.3511*v**2)
        y=(0.7591+dsqrt(0.2409d0**2+0.3145d0*v**2))**2
        z=dexp(-dsqrt(4.d0*v**2+0.4616d0**2)+0.4616d0)
        control_pbf_3P2_B=x*y*z
       return
      end
c *********************************************************************
c *********************************************************************
      function control_pbf_3P2_C(v)
       implicit real*8 (a-h,k-z)
        x=0.4013d0*v**2-0.043d0*v**4+0.002172d0*v**6
        y=1.d0/
     1   (1.d0-2.018d-1*v**2+2.601d-2*v**4-1.477d-3*v**6+4.34d-5*v**8)
        control_pbf_3P2_C=x*y
       return
      end
c *********************************************************************
c *********************************************************************

c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c
c                     Q U A R K   P R O C E S S E S
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************

c *********************************************************************
c *********************************************************************
      subroutine nudurca_q(i,t,rho,qdurca_q)
c *** Checked on April 28, 2004              <======  TO BE DONE !?!?!?
c**********************************************************************
c     Includes the Levenfish & Yakovlev 
c     suppression factors for nucleon DURCA
c**********************************************************************
       implicit real*8 (a-h,k-z)
       parameter (pi = 3.14159265d0)
       parameter (g_fermi=1.436d-49,theta_c=0.239,mev=1.6d-6)
       parameter (h_bar=1.054d-27,kb=1.38d-16,c_light=3.d10)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'control_nu.inc.f'
       INCLUDE 'quark.inc.f'
c ***************************
       fexp(x)=dexp(max(x,-7.d2))

       u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
       r_1s0(u)=(0.2312d0+dsqrt(0.7688d0**2+(0.1438d0*u)**2))**5.5d0*
     1          fexp(3.427d0-dsqrt(3.427d0**2+u**2))

c       u_3p2b(t)=dsqrt(1.d0-t)*(5.596d0+8.424d0/t)
c       r_3p2b(u)=(0.2546d0+dsqrt(0.7454d0**2+(0.01811d0*u)**2))**5*
c     1           fexp(2.701d0-dsqrt(2.701d0**2+u**2/(16.d0*pi)))

       r_1s0_1s0(u1,u2)=func_r_1s0_1s0(u1,u2)
c       r_1s0_3p2B(t1,t2)=func_r_1s0_3P2B(t1,t2) ! WARNING: this uses t1 & t2
c                                                ! instead of u1 & u2
c       u_3p2c(t)=dsqrt(1.d0-t**4)/t*
c     1           (5.875d0-1.419d0*t**4+0.500d0*t**8)
c       r_3p2c(u)=1.d0/0.d0
c****************************

c ******************************************
c **** u-d: for each color:
        coeff_ud= (1.d0/3.d0) * 914./315.*
     1           (g_fermi*dcos(theta_c)/(h_bar**5*c_light**3))**2*
     2           alpha_c
        if (idurca_quqd(i).eq.1)then
         qdurca_quqd=coeff_ud*1.d+39*
     1               kfqd(i)*kfqu(i)*kfe(i)*h_bar**3*(kb*t)**6
        else if(idurca_quqd(i).eq.2)then
         qdurca_quqd=coeff_ud*1.d+39*
     1               (kfqd(i)*kfqu(i)*kfe(i)*h_bar**3*(kb*t)**6
     2               +kfqd(i)*kfqu(i)*kfm(i)*h_bar**3*(kb*t)**6)
        else
         qdurca_quqd= 0.d0
        end if
c Pairing suppression: COLOR 1
        if ((t.gt.tcu1(i)).and.(t.gt.tcd1(i))) then
         r_ud1=1.d0
        else if ((t.gt.tcu1(i)).and.(t.le.tcd1(i))) then
         tt=t/tcd1(i)
         u=u_1s0(tt)
         r_ud1=r_1s0(u)
        else if ((t.le.tcu1(i)).and.(t.gt.tcd1(i))) then
         tt=t/tcu1(i)
         u=u_1s0(tt)
         r_ud1=r_1s0(u)
        else
         tt1=t/tcu1(i)
         u1=u_1s0(tt1)
         tt2=t/tcd1(i)
         u2=u_1s0(tt2)
         r_ud1=r_1s0_1s0(u1,u2)
        end if
c Pairing suppression: COLOR 2
        if ((t.gt.tcu2(i)).and.(t.gt.tcd2(i))) then
         r_ud2=1.d0
        else if ((t.gt.tcu2(i)).and.(t.le.tcd2(i))) then
         tt=t/tcd2(i)
         u=u_1s0(tt)
         r_ud2=r_1s0(u)
        else if ((t.le.tcu2(i)).and.(t.gt.tcd2(i))) then
         tt=t/tcu2(i)
         u=u_1s0(tt)
         r_ud2=r_1s0(u)
        else
         tt1=t/tcu2(i)
         u1=u_1s0(tt1)
         tt2=t/tcd2(i)
         u2=u_1s0(tt2)
         r_ud2=r_1s0_1s0(u1,u2)
        end if
c Pairing suppression: COLOR 3
        if ((t.gt.tcu3(i)).and.(t.gt.tcd3(i))) then
         r_ud3=1.d0
        else if ((t.gt.tcu3(i)).and.(t.le.tcd3(i))) then
         tt=t/tcd3(i)
         u=u_1s0(tt)
         r_ud3=r_1s0(u)
        else if ((t.le.tcu3(i)).and.(t.gt.tcd3(i))) then
         tt=t/tcu3(i)
         u=u_1s0(tt)
         r_ud3=r_1s0(u)
        else
         tt1=t/tcu3(i)
         u1=u_1s0(tt1)
         tt2=t/tcd3(i)
         u2=u_1s0(tt2)
         r_ud3=r_1s0_1s0(u1,u2)
        end if
c Putting everything together (3 colors + pairing suppression) :

         qdurca_quqd= (r_ud1+r_ud2+r_ud3) * qdurca_quqd

c ******************************************
c **** u-s: for each color
        theta_34=pi/4.     ! <----- A rough approximation !
        coeff_us=(1.d0/3.d0) * 457.*pi/840.*
     1           (g_fermi*dsin(theta_c)/(h_bar**5*c_light**3))**2*
     2           (1.d0-dcos(theta_34))
        if (idurca_quqs(i).eq.1)then
         mus=dsqrt(kfqs(i)**2*1.d+26*c_light**2*h_bar**2+
     1             (mev*strange_mass)**2)
         qdurca_quqs=coeff_us*1.e+26*
     1               mus/c_light*kfqu(i)*kfe(i)*h_bar**2*(kb*t)**6
        else if(idurca_quqs(i).eq.2)then
         mus=dsqrt(kfqs(i)**2*1.e+13*c_light**2*h_bar**2+
     1             (mev*strange_mass)**2)
         qdurca_quqs=coeff_us*1.e+26*
     1               (mus/c_light*kfqu(i)*kfe(i)*h_bar**2*(kb*t)**6
     2               +mus/c_light*kfqu(i)*kfm(i)*h_bar**2*(kb*t)**6)
        else
         qdurca_quqs= 0.d0
        end if
c Pairing suppression: COLOR 1
        if ((t.gt.tcu1(i)).and.(t.gt.tcs1(i))) then
         r_us1=1.d0
        else if ((t.gt.tcu1(i)).and.(t.le.tcs1(i))) then
         tt=t/tcs1(i)
         u=u_1s0(tt)
         r_us1=r_1s0(u)
        else if ((t.le.tcu1(i)).and.(t.gt.tcs1(i))) then
         tt=t/tcu1(i)
         u=u_1s0(tt)
         r_us1=r_1s0(u)
        else
         tt1=t/tcu1(i)
         u1=u_1s0(tt1)
         tt2=t/tcs1(i)
         u2=u_1s0(tt2)
         r_us1=r_1s0_1s0(u1,u2)
        end if
c Pairing suppression: COLOR 2
        if ((t.gt.tcu2(i)).and.(t.gt.tcs2(i))) then
         r_us2=1.d0
        else if ((t.gt.tcu2(i)).and.(t.le.tcs2(i))) then
         tt=t/tcs2(i)
         u=u_1s0(tt)
         r_us2=r_1s0(u)
        else if ((t.le.tcu2(i)).and.(t.gt.tcs2(i))) then
         tt=t/tcu2(i)
         u=u_1s0(tt)
         r_us2=r_1s0(u)
        else
         tt1=t/tcu2(i)
         u1=u_1s0(tt1)
         tt2=t/tcs2(i)
         u2=u_1s0(tt2)
         r_us2=r_1s0_1s0(u1,u2)
        end if
c Pairing suppression: COLOR 3
        if ((t.gt.tcu3(i)).and.(t.gt.tcs3(i))) then
         r_us3=1.d0
        else if ((t.gt.tcu3(i)).and.(t.le.tcs3(i))) then
         tt=t/tcs3(i)
         u=u_1s0(tt)
         r_us3=r_1s0(u)
        else if ((t.le.tcu3(i)).and.(t.gt.tcs3(i))) then
         tt=t/tcu3(i)
         u=u_1s0(tt)
         r_us3=r_1s0(u)
        else
         tt1=t/tcu3(i)
         u1=u_1s0(tt1)
         tt2=t/tcs3(i)
         u2=u_1s0(tt2)
         r_us3=r_1s0_1s0(u1,u2)
        end if
c Putting everything together (3 colors + pairing suppression) :

         qdurca_quqs= (r_us1+r_us2+r_us3) * qdurca_quqs

c All Quark Durca processes:

        qdurca_q = qdurca_quqd  + qdurca_quqs

       return
      end
c *********************************************************************
c *********************************************************************

      subroutine numurca_q(i,t,rho,qmurca_q)

c **** checked
 
c**********************************************************************
c Includes the Levenfish & Yakovlev suppression factors for MURCA,
c**********************************************************************

      implicit real*8 (a-h,k-z)
      parameter (pi = 3.14159265d0)
      parameter (g_fermi=1.436d-49,theta_c=0.239,mev=1.6d-6)
      parameter (h_bar=1.054d-27,kb=1.38d-16,c_light=3.d10)
      INCLUDE 'rho_limits.inc.f'
      INCLUDE 'size.inc.f'
      INCLUDE 'profile_star.inc.f'
      INCLUDE 'profile_comp.inc.f'
      INCLUDE 'pairing.inc.f'
      INCLUDE 'fermi.inc.f'
      INCLUDE 'control_nu.inc.f'
      INCLUDE 'quark.inc.f'

       num_coeff=1.e0     ! Remains to be calculated

       qmurca_q=num_coeff*
     1           (alpha_c*g_fermi*cos(theta_c)/
     2           h_bar**5/c_light**4)**2 *
     3           (1.e13*kfqu(i)*h_bar) * (kb*t)**8

c Up quark pairing suppression:
       if (t.lt.tcu(i)) then
        r_u=dexp(-1.76*tcu(i)/t)
       else
        r_u=1.0
       end if
c Down quark pairing suppression:
       if (t.lt.tcd(i)) then
        r_d=dexp(-1.76*tcd(i)/t)
       else
        r_d=1.0
       end if
cc Strange quark pairing suppression:
c       if (t.lt.tcs(i)) then
c        r_s=dexp(-1.76*tcs(i)/t)
c       else
c        r_s=1.0
c       end if
c****
       r_ud=r_u*r_d
       r_us=r_u*r_s

       qmurca_q = qmurca_q * r_ud

       return

      end


c *********************************************************************
c *********************************************************************
      subroutine nu_strange(i,T,qnu,debug)
c     CHECKED on March 3, 2002
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'quark.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'pairing_quark_PU2002.inc.f'
       parameter (pi=3.1415926535d0)
c ***************************
       fexp(x)=dexp(max(x,-7.d2))
       u_1s0(t)=dsqrt(1.d0-t)*(1.456d0-0.157d0/dsqrt(t)+1.764d0/t)
       r_1s0(u)=(0.2312d0+dsqrt(0.7688d0**2+(0.1438d0*u)**2))**5.5d0*
     1         fexp(3.427d0-dsqrt(3.427d0**2+u**2))
       r_1s0_1s0(u1,u2)=func_r_1s0_1s0(u1,u2)
c       u_3p2C(t)=dsqrt(1.d0-t**4)/t*
c     1           (5.875d0-1.419d0*t**4+0.500d0*t**8)
c ***************************
        if (debug.ge.2.) print *,'Entering subroutine `neutrino'' '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *,'subroutine nu_strange:'
        print *,'   no guarantee that it returns a meaningful value !'
        pause
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ***************************
c Just to be safe:
         qnu=0.d0
         qnu_one=0.d0
         qnu_12=0.d0
         qnu_3=0.d0
c ***************************
         nbar=bar(i)
         T9=T/1.d9
         alphac=alpha_c/4.d0       ! alpha_c=g^2/4pi alphac=g^2/16pi
         coeff=(1.d0+16.d0*alphac/3.d0/pi)*alphac
c Neutrino emissivity for one color d --> u+e+nubar & d+e --> u+nu :
         qnu_one=1.25d26*coeff*kfqu(i)*kfqd(i)*kfe(i)*T9**6
c Pairing suppression:
         check_normal=0.0
         check_CFL=0.0
         check_3SC=0.0
         check_2SC=0.0
         check_s  =0.0
c *******************************************************************
c NO PAIRING:
         IF (normal.eq.1.d0) THEN
          qnu=3.d0*qnu_one
          check_normal=1.0
c *******************************************************************
c CFL PHASE:
         ELSE IF (nbar.ge.nb_CFL) THEN
          qnu_qqbr=3.d19*3.d0*T9**8*(nbar/0.16d0)**(1.d0/3.d0)
          check_CFL=1.0
          If (T.ge.Tc_CFL(i)) then
           qnu=qnu_qqbr
          Else
           tt=T/Tc_CFL(i)
           u=u_1s0(tt)
           qnu_qqbr=qnu_qqbr*r_1s0(u)**2
           qnu_pbf=1.4d20*3.d0*T9**7*function_pbf_1S0(u)
     1             *(nbar/0.16d0)**(1.d0/3.d0)
           qnu=qnu_qqbr+qnu_pbf
          End if
c          qnu=1.d0                           ! No electrons present !
c *******************************************************************
c 3SC PHASE:
         ELSE IF ((nbar.ge.nb_3SC).and.
     1            (nbar.lt.nb_CFL)) THEN             
          check_3SC=1.0
          If (T.ge.Tc_3SC(i)) then
           qnu=3.d0*qnu_one
          Else
           tt=T/Tc_3SC(i)
           u=u_1s0(tt)
           if      (type_3SC.eq.0.d0) then
            r=r_1s0(u)
           else if (type_3SC.eq.1.d0) then
            r=tt**2
           else if (type_3SC.eq.2.d0) then
            r=tt
           else
            pause 'Neutrino: type_3SC wrong !'
           end if
           qnu=r*3.d0*qnu_one
          End if
c *******************************************************************
c 2SC PHASE:
         ELSE
          check_2SC=1.0
          If (T.ge.Tc_2SC(i)) then
c No quarks paired:
           qnu=3.d0*qnu_one
          Else
c u d color 1 & 2 paired:
          check_2SC=2.0
           tt=T/Tc_2SC(i)
           u=u_1s0(tt)
           if      (type_2SC.eq.0.d0) then
            r=r_1s0(u)
           else if (type_2SC.eq.1.d0) then
            r=tt**2
           else if (type_2SC.eq.2.d0) then
            r=tt
           else
            pause 'Neutrino: type_2SC wrong !'
           end if
           qnu_12=r*2.d0*qnu_one
c u d color 3 not paired
           if ((T.ge.Tc_u(i)).and.(T.ge.Tc_d(i))) then
            qnu_3=qnu_one
c u color 3 paired
           else if ((T.lt.Tc_u(i)).and.(T.ge.Tc_d(i))) then
            check_2SC=2.1
            tt=T/Tc_u(i)
            u=u_1s0(tt)
            if      (type_u.eq.0.d0) then
             r=r_1s0(u)
            else if (type_u.eq.1.d0) then
             r=tt**2
            else if (type_u.eq.2.d0) then
             r=tt
            else
             pause 'Neutrino: type_u wrong !'
            end if
            qnu_3=r*qnu_one
c d color 3 paired
           else if ((T.ge.Tc_u(i)).and.(T.lt.Tc_d(i))) then
           check_2SC=2.2
            tt=T/Tc_d(i)
            u=u_1s0(tt)
            if      (type_d.eq.0.d0) then
             r=r_1s0(u)
            else if (type_d.eq.1.d0) then
             r=tt**2
            else if (type_d.eq.2.d0) then
             r=tt
            else
             pause 'Neutrino: type_d wrong !'
            end if
            qnu_3=r*qnu_one
c u & d color 3 paired
           else
            check_2SC=2.3
            tt1=T/Tc_u(i)
            u1=u_1s0(tt1)
            tt2=T/Tc_d(i)
            u2=u_1s0(tt2)
            if      ((type_u.eq.0.d0).and.(type_d.eq.0.d0)) then
             r=r_1s0_1s0(u1,u2)
            else if ((type_u.eq.0.d0).and.(type_d.ne.0.d0)) then
             r=r_1s0(u1)
            else if ((type_u.ne.0.d0).and.(type_d.eq.0.d0)) then
             r=r_1s0(u2)
            else
             if      (type_u.eq.1.d0) then
              r1=tt1**2
             else if (type_u.eq.2.d0) then
              r1=tt1
             else
              pause 'Neutrino: type_u wrong'
             end if
             if      (type_d.eq.1.d0) then
              r2=tt2**2
             else if (type_d.eq.2.d0) then
              r2=tt2
             else
              pause 'Neutrino: type_d wrong'
             end if
             r=min(r1,r2)
            end if
            qnu_3=r*qnu_one
           end if
           qnu=qnu_12+qnu_3
          End if
         END IF
c *********************************************************************
c   Adding s-s bremstrahlung:
         qnu_ssbr=1.d20*T9**8
         if (T.le.Tc_s(i)) then
          check_s=1.0
          tt=T/Tc_s(i)
          u=u_1s0(tt)
          if (type_s.eq.0.d0) then
           rs=r_1s0_1s0(u,u)
          else if (type_s.eq.1.d0) then
           rs=tt**2
          else if (type_s.eq.2.d0) then
           rs=tt
          else
           pause 'Neutrino: type_d wrong'
          end if
         else
          rs=1.d0
         end if
         qnu_ssbr=qnu_ssbr*rs
c *********************************************************************
c   Adding e-e bremstrahlung:
c   From Kaminker & Haensel, astro-ph/9908249
c   It's wrong ! But better than nothing.
        if (nbar.lt.nb_CFL) then
         k0=1.68
         qnu_eebr=0.69d14 * (kfe(i)/k0) * (T/1.d9)**8
        else
         qnu_eebr=0.d0                          ! No electrons in CFL !
        end if
         qnu_eebr=0.d0
c *********************************************************************
        qnu=qnu+qnu_eebr+qnu_ssbr
c *********************************************************************
        if (debug.lt.-0.5d0) then
         print '(i5,1p3e12.3,3x,1p2e12.3,3x,0p5f6.2)',
     1       i,rrho(i),T,qnu-qnu_eebr-qnu_ssbr,qnu_eebr,qnu_ssbr,
     2       check_normal,check_CFL,check_3SC,check_2SC,check_s
        end if
c *********************************************************************
c OLD SUBROUTINE FROM USOV's PAPER:
cc         alpha_c=1.d0     ! It's in quark.inc.f
c         nb=nbar
c         Ye=yelect(i)
c         qnu=2.2d26*alpha_c*Ye**(1.0/3.0)*(nb/0.16)*(T/1.d9)**6
c         qnu=1.0d0*qnu
cc        Warning: alpha_c = g^2/4pi
c *********************************************************************
        if (debug.ge.2.) print *,'Exiting subroutine `neutrino'' '
       return
      end
c**********************************************************************
c**********************************************************************
      function function_pbf_1S0(v)
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415926535)
        x=0.602d0*v**2+0.5942d0*v**4+0.288d0*v**6
        y=dsqrt(0.5547d0+dsqrt(0.4453d0**2+0.01130*v**2))
        z=dexp(-dsqrt(4.d0*v**2+2.245d0**2)+2.245d0)
        function_pbf_1S0=x*y*z
       return
      end
c *********************************************************************
c *********************************************************************




