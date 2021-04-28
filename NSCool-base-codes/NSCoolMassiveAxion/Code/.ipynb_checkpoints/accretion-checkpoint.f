c-------------------------------------------------------------
c  i_acc   Accretion type
c
c   0      None
c   1      Transient FRED
c   2      Transient STEP
c-------------------------------------------------------------
      subroutine initialize_accretion_rate
c ***** checked on 
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'accretion.inc.f'
       INCLUDE 'profile_star.inc.f'
       common/gravity/gs14,compactness
c
c  FIRST: Figure out units for m_dot:
c
        if (m_dot0.gt.1.d10) then      ! m_dot in g/sec
         m_dot1=m_dot0
        else if (m_dot0.lt.1.d-6) then ! m_dot in Msun/yr
         m_dot1=m_dot0 * 2.d33 / 3.15576d7
        else                           ! m_dot in unit of M_Edd
         eta=eta_Edd                        ! Efficiency of mass-> energy
         X = X_Edd                          ! Hydrogen fraction
         opac=6.023d23*0.665d-24*(1.+X)/2.
         F_Edd=3d10*(1.d14*gs14)/opac
         L_edd=4.*3.14159*rad(imax)**2*F_Edd
         M_Edd=L_Edd/(3.d10)**2/eta
         m_dot1=m_dot0*M_Edd
        end if
c
c       INITIALIZE
c 
c No accretion: ------------------------------------------
        if (i_acc.eq.0) then
         mass_acc=0.
c Accretion in FRED transients: --------------------------------------
        else if (i_acc.eq.1) then
c        t_acc0 = beginning of accretion
c        t_acc1 = duration of each accretion cycle
c        t_acc2 = tau_acc_r = duration of exponential rise
         tau_acc_r=t_acc2
         if (m_dot0.lt.1.) then
c         In this case m_dot0 is actually DelM in Mo !
          coeff=1.d0/(exp(1.d0)-1.d0)+
     1         (1.d0-(t_acc2/t_acc1)**(alpha_acc-1.))/(alpha_acc-1.)
c          print *,'COEFF=',coeff
          m_dot0=m_dot0/t_acc2/coeff * 2.e33
         end if
         m_dot_max=m_dot0
         m_dot_ris=m_dot_max/(1.-exp(-1.))
         macc_r=t_acc2*m_dot_ris/exp(1.)
         macc_d=t_acc2*m_dot_max/(alpha_acc-1.)*
     1          (1.-(t_acc2/t_acc1)**(alpha_acc-1.))
         m_dot_av=(macc_r+macc_d)/t_acc1
c         m_dot_ini=0.0
c         m_dot_ini=0.99*m_dot_av
         m_dot_ini=m_dot_av
c         m_dot_quiet=1.e-5*m_dot_max
         m_dot_quiet=0.
         print '(7a12)',
     1      'm_dot_max','m_dot_av',
     2      'm_dot_ris','m_dot_ini',
     3      't_rec','t_rise','Del(M)'
         print '(7a12)',
     1      '[Mo/yr] ','[Mo/yr] ',
     2      '[Mo/yr] ','[Mo/yr] ',
     3      '[days]','[days]','[Mo] '
         coeff=3.15e7/2.d33
         print '(1p4e12.3,0p2f12.1,1p1e12.3)',
     1      m_dot_max*coeff,m_dot_av*coeff,
     2      m_dot_ris*coeff,m_dot_ini*coeff,
     2      t_acc1/86400.d0,t_acc2/86400.d0,m_dot_av*t_acc1/2.d33
c Accretion in STEP transients: --------------------------------------
        else if (i_acc.eq.2) then
c        t_acc0 = beginning of accretion
c        t_acc1 = duration of each accretion cycle
c        t_acc2 = duration of each accretion outburst
         m_dot_max=m_dot1
         m_dot_av =m_dot1*t_acc2/t_acc1
         m_dot_ini=m_dot_av
         m_dot_quiet=0.d0
         print '(2a12)',
     1      'm_dot_max','m_dot_av'
         print '(2a12)',
     1      '[Mo/yr] ','[Mo/yr] '
         coeff=3.15e7/2.d33
         print '(1p2e12.3)',
     1      m_dot_max*coeff,m_dot_av*coeff
c WRONG input !
        else
         pause 'Initializa_accretion_rate: i_acc badly defined'
        end if
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine accretion_rate(time,dtime,m_dot)
c ***** checked on 
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'accretion.inc.f'
c *******************************
        fr(t)=m_dot_ris*(t+t_acc2*exp(-t/t_acc2))
        fd(t)=m_dot_max*t/(1.-alpha_acc)*(t_acc2/t)**alpha_acc
c *******************************
c No accretion:
        if (i_acc.eq.0) then
         m_dot=0.d0
c Accretion in FRED transients: --------------------------------------
        else if (i_acc.eq.1) then
         if ((time-dtime).le.t_acc0) then
          m_dot=m_dot_ini
         else
          if (dtime.ge.t_acc1) then
           pause '  dtime larger than t_acc1 !'
          end if
          time0=time-dtime
          time1=time
          icycle0=int((time0-t_acc0)/t_acc1)
          t0=time0-t_acc0-float(icycle0)*t_acc1
          t1=time1-t_acc0-float(icycle0)*t_acc1
          icycle=icycle0                     ! for book-keeping only
          delt_acc=t1                        ! for book-keeping only
          if (t0.le.t_acc2) then
           if (t1.le.t_acc2) then
            macc=fr(t1)-fr(t0)
           else if (t1.le.t_acc1) then
            macc0=fr(t_acc2)-fr(t0)
            macc1=fd(t1)-fd(t_acc2)
            macc=macc0+macc1
           else
            macc0=fr(t_acc2)-fr(t0)
            macc1=fd(t_acc1)-fd(t_acc2)
            macc2=fr(t1-t_acc1)-fr(0.d0)
            macc=macc0+macc1+macc2
           end if
          else
           if (t1.le.t_acc1) then
            macc=fd(t1)-fd(t0)
           else if (t1.le.(t_acc1+t_acc2)) then
            macc0=fd(t_acc1)-fd(t0)
            macc1=fr(t1-t_acc1)-fr(0.d0)
            macc=macc0+macc1
           else
            macc0=fd(t_acc1)-fd(t0)
            macc1=fr(t_acc2)-fr(0.d0)
            macc2=fd(t1-t_acc1)-fd(t_acc2)
            macc=macc0+macc1+macc2
           end if
          end if
          m_dot=macc/dtime
         end if
c Accretion in STEP transients: --------------------------------------
        else if (i_acc.eq.2) then
         if (time+dtime.le.t_acc0) then
          m_dot=m_dot_ini
         else
c          if (dtime.ge.t_acc2) then
c           pause '  dtime larger than t_acc2 = outburst length !'
c          end if
          if (t_burst.le.t_acc2) then
           m_dot=m_dot_max
          else
           m_dot=m_dot_quiet
          end if
         end if
c Just in case:
        else
         m_dot=0.d0
        end if
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine accreted_mass(dtime,m_dot)
c ***** checked on 
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415926535)
       INCLUDE 'size.inc.f'
       INCLUDE 'accretion.inc.f'
       INCLUDE 'profile_star.inc.f'
        mass_acc=mass_acc+m_dot*dtime
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine accretion_velocity(m_dot)
c ***** checked on 
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415926535)
       INCLUDE 'size.inc.f'
       INCLUDE 'accretion.inc.f'
       INCLUDE 'profile_star.inc.f'
        do i=1,imax
         v_acc(i)= - m_dot/(4.*pi*rad(i)**2*rrho(i))
        end do
        v_acc(0)=0.
       return
      end
c *********************************************************************
c *********************************************************************


