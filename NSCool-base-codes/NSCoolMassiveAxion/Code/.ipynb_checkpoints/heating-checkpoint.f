c*********************************************************************
c*********************************************************************
c    Heating sources:
c
c    WARNING: DO NOT ASSUME ANY OF THESE SUBROUTINES ACTUALLY WORKS !
c
c    Before using any heating source, check it carefully !
c    Some may be old (or very old), and not work properly
c    in the present version of the code.
c
c    Each process is controlled by its own parameter: i_heat_*
c    If i_heat_* = 0 the process is turned off (subroutine is not called).
c
c*********************************************************************
c*********************************************************************
      subroutine heating(i,time,dtime,temp,rho,a,a1,z,tcp,
     1                   m_dot,i_eleconb,ephi,dephi,heat)
c **** Partially checked on July 25, 1996 ************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       parameter (year=3.15e7)
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'mag_field.inc.f'
       INCLUDE 'control_heat.inc.f'
       real*8 j_44,j_heat
       
c     Chris: Magnetic Field Evolution
c     Initialize parameters
      real*8 :: Bcurr,dBdt
      common/bfield_values/Bcurr,dBdt

c***** Deep crustal heating:
       if (i_heat_deep_crust.eq.1) then
        if (rrho(i).lt.rhocore) then
         call heating_deep_crust(i,m_dot,heat_deep_crust)
        else
         heat_deep_crust=0.d0
        end if
       else
        heat_deep_crust=0.d0
       end if
c***** Heat deposition:
       if (i_heat_deposit.eq.1) then
        time1=time
        time2=time+dtime
        call heating_deposit(i,time1,time2,heat_deposit)
       else
        heat_deposit=0.d0
       end if
c***** Heating by neutron -> SQM conversion
       if (i_heat_convert.eq.1) then
        call heating_convert(i,m_dot,heat_convert)
       else
        heat_convert=0.d0
       end if
c***** Crustal superfluid friccion:
       if (i_heat_vortex_creep.eq.1) then
        call heating_vortex_creep(i,time,heat_vortex_creep)
       else
        heat_vortex_creep=0.d0
       end if
c***** Joule Heating:
       if (i_heat_joule.eq.1) then
        call conduct(i,temp,rho,a,a1,z,qimp,nbfield2(i),i_eleconb,
     1               sig,lll,debug,nu_e_s,nu_e_l)
c        call elecon(i,temp,rho,a,a1,z,qimp,tcp,nbfield2(i),i_eleconb,
c     1              sig,print)
        curr=current(i)
        btheta=bf_t(i)
        call heating_joule(curr,btheta,sig,ephi,dephi,heat_joule)
       else
        heat_joule=0.d0
       end if
c***** Magnetic field decay a la GOLDREICH-REISENNEGGER
       if (i_heat_field_decay.eq.1) then
        if (rrho(i).ge.rhocore) then
         call heating_field_decay(i,bfield2(i),temp,heat_field_decay)
        else
         heat_field_decay=0.0
        end if
       else
        heat_field_decay=0.d0
       end if
c***** Chris: Magnetic Field Evolution
c      The heat gets *-1 already, so this should be positive value
c      To heat the NS. dBdt is negative.
       heat_constant_field_evol=-2.52171d-9*Bcurr*dBdt
c*****
        heat=
     1       heat_deep_crust+
     2       heat_deposit+
     3       heat_convert+
     4       heat_vortex_creep+
     5       heat_joule+
     6       heat_field_decay+
     7       heat_constant_field_evol
c*****
       return
      end
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
      subroutine heating_deep_crust(i,m_dot,heat)
       implicit real*8 (a-h,k-z)
       parameter(m_u=1.67e-24)
       INCLUDE 'size.inc.f'
       dimension q_deep_crust(0:isize)
       common/deep_crust/q_deep_crust
        heat=(m_dot/m_u) * q_deep_crust(i)
       return
      end
c*********************************************************************
c*********************************************************************
      subroutine initialize_heating_deep_crust
c
c    OLD SUBROUTINE: Probably does not work anymore !
c
c *****************************************************************
c           Haensel & Zdunik, A&A 227 (1990), p. 431
c *****************************************************************
       implicit real*8 (a-h,k-z)
       dimension rho_pycno(18),q_pycno(18)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       dimension qpycno(0:isize)
       common/pycno/qpycno
       data rho_pycno/1.494e09,1.114e10,7.848e10,
     1                2.496e11,6.110e11,9.075e11,
     2                1.131e12,1.455e12,1.766e12,
     3                2.134e12,2.634e12,3.338e12,
     4                4.379e12,5.839e12,7.041e12,
     5                8.980e12,1.127e13,1.137e13/
       data q_pycno  /  0.01  ,  0.01  ,  0.01  ,
     1                  0.01  ,  0.05  ,  0.09  ,
     2                  0.10  ,  0.47  ,  0.05  ,
     3                  0.05  ,  0.06  ,  0.07  ,
     4                  0.28  ,  0.02  ,  0.02  ,
     5                  0.03  ,  0.11  ,  0.01  /
c *****
        do i=imax,0,-2
         qpycno(i)=0
         qpycno(i-1)=0.
         do j=1,18
          if ( (rho_pycno(j).ge.rrho( i )) .and.
     1         (rho_pycno(j).lt.rrho(i-2)) ) then
           qpycno(i)=qpycno(i)+q_pycno(j)
          end if
         end do
         qpycno(i)=qpycno(i)*1.6e-6 / (dvol(i)+dvol(i+1))
        end do
       return
      end
c*********************************************************************
c*********************************************************************
      subroutine heating_deposit(i,time1,time2,heat_deposit)
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415826535)
       INCLUDE 'size.inc.f'
       INCLUDE 'control_heat.inc.f'
       dimension qdeposit(0:isize)
       common/deposit/total_heat,t_dep,del_t_dep,qdeposit,i_dep
c       common/control_heating/i_joule_heat,i_pycno_heat,
c     1                       i_deposit_heat,i_convert_heat

        t1=(time1-t_dep)/del_t_dep
        t2=(time2-t_dep)/del_t_dep
        if ((t2.le.-100.).or.(t1.ge.+100.)) then
         heat_deposit=0.d0
        else
         heat_deposit=qdeposit(i)/2./(time2-time1)*
     1                (erfcc(t1)-erfcc(t2))
        end if

       return
      end
c*********************************************************************
      function erfcc(x)
c From Numerical Recipes, `Complementary Error Fnction'
       implicit real*8 (a-h,k-z)
        z=abs(x)
        t=1.d0/(1.d0+0.5d0*z)
        erfcc=t*dexp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
     1     t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
     2     t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
        if (x.lt.0.0) erfcc=2.-erfcc
       return
      end
c*********************************************************************
c*********************************************************************
      subroutine initialize_heating_deposit
c     i_dep=zone where heat will de deposited 
c     total_heat= total amount of heat deposited
c     t_dep=time of maximum heating
c     del_t_dep=witdh of the gaussian deposition rate
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       dimension qdeposit(0:isize)
       common/deposit/total_heat,t_dep,del_t_dep,qdeposit,i_dep
        do i=0,imax
         qdeposit(i)=0
        end do
        qdeposit(i_dep)=total_heat / (dvol(i_dep)+dvol(i_dep+1))
       return
      end
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
      subroutine heating_convert(i,m_dot,heat_convert)
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415826535)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'control_heat.inc.f'
       dimension qconvert(0:isize)
       common/convert/qconvert
c       common/control_heating/i_joule_heat,i_pycno_heat,
c     1                       i_deposit_heat,i_convert_heat

        heat_convert=qconvert(i)*m_dot

       return
      end
c*********************************************************************
c*********************************************************************
c   This inject heat at surface of strange matter from baryon conversion !
      subroutine initialize_heating_convert(MeV_neutron)
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       dimension qconvert(0:isize)
       common/convert/qconvert
c HERE DANY ccccccccccccccccccccccccccccccccccccccccc
c WARNING: need i_convert_1 <= i_convert_2
        i_convert_1=icore
        i_convert_2=icore
c        i_convert_1=101
c        i_convert_2=101
c        i_convert_1=11
c        i_convert_2=51
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        print *,'Initializing deposit convert:'
        do i=0,imax
         qconvert(i)=0.d0
         if ((i.ge.i_convert_1).and.(i.le.i_convert_2)) then
          qconvert(i)=MeV_neutron * 1.6d-6 * 6.023d23
     1                / (dvol(i)+dvol(i+1))
c         ergs per cm^3 for each gram of accreted matter since m_dot is in g/sec
c          print '(a20,i5,1p1e12.3,a24)',
c     1     'i,qconvert(i) =',i,qconvert(i),'erg/cm^3/sec per gram'
c          read(5,*)
         end if
        end do
       return
      end
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
      subroutine heating_vortex_creep(i,time,heat)
c **** checked on July 25, 1996 ************
       implicit real*8 (a-h,k-z)
       real*8 j_heat
       INCLUDE 'size.inc.f'
       parameter (year=3.15e7)
       dimension mag_field(0:isize),j_heat(0:isize)
       INCLUDE 'rho_limits.inc.f'
       common/others/mag_field,j_heat
c*****
        tau=300.
c*****
        heat=1.d40*(j_heat(i)+j_heat(i+1))/(time/year+tau)**(3./2.)
       return
      end
c*********************************************************************
c*********************************************************************
      subroutine initialize_heating_vortex_creep
     1           (imax,rad,rrho,tcn,dvol,j_44)
c **** checked on July 25, 1996 ************
       implicit real*8 (a-h,k-z)
       real*8 j_44,j_heat
       INCLUDE 'size.inc.f'
       dimension rad(0:isize),rrho(0:isize),tcn(0:isize),
     1           dvol(0:isize)
       dimension mag_field(0:isize),j_heat(0:isize)
       INCLUDE 'rho_limits.inc.f'
       common/others/mag_field,j_heat

        sum=0.0
        do i=1,imax,2
         if ((rrho(i).ge.rhodrip).and.(rrho(i).lt.rhocore)) then
          j_heat(i)=tcn(i)**3*rrho(i)**(1./3.)*rad(i)**3*
     1              (rad(i+1)-rad(i))
          sum=sum+j_heat(i)
         end if
        end do
        do i=1,imax,2
         j_heat(i)=j_44*(j_heat(i)/sum)
         j_heat(i)=j_heat(i)/(dvol(i)+dvol(i+1))
        end do
cccccccccccccccccccccccccccccccc
c        sum=0.0
c        do i=1,imax,2
c         if (j_heat(i).ne.0.0) then
c          print '(1p1e12.3,0p1f10.5)',rrho(i),j_heat(i)*
c     1                                (dvol(i)+dvol(i+1))
c         end if
c         sum=sum+j_heat(i)*(dvol(i)+dvol(i+1))
c        end do
c        print '(a15,0p2f10.3)','j_44 & sum =',j_44,sum
cccccccccccccccccccccccccccccccc
       return
      end
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
c*********************************************************************
      subroutine heating_joule(current,btheta,sigma,ephi,dephi,
     1                         joule)
       implicit real*8 (a-h,k-z)
       real*8 joule
       parameter(pi=3.1415826535)
       parameter(c=2.99792e10)
       parameter (year=3.15576e7)
       INCLUDE 'control_heat.inc.f'
c       common/control_heating/i_joule_heat,i_pycno_heat,
c     1                       i_deposit_heat,i_convert_heat
       common/stuff/time,istep


        joule=(2./3.)*current**2/sigma     ! 2/3 is from spherical average
                                           ! current being at theta=90deg

        joule_gr=c/6.d0/pi/ephi*dephi*current*btheta/sigma

        joule=joule-joule_gr

       return
      end
c*********************************************************************
c*********************************************************************
c*********************************************************************
      subroutine heating_field_decay(i,field,temp,heat)
c **** checked on July 25, 1996 ************
       implicit real*8 (a-h,k-z)
       real*8 j_heat
       INCLUDE 'size.inc.f'
       parameter(pi=3.1415926535)
       dimension mag_field(0:isize),j_heat(0:isize)
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'control_nu.inc.f'
c ***************************
c      u_3p2b(t)=sqrt(1.-t)*(5.596+8.424/t)
c      u_3p2c(t)=sqrt(1.-t**4)/t*(5.875-1.419*t**4+0.500*t**8)
c      r_3p2c(u)=1./0.

c ***************************
c SINGLET PAIRING:
      u_1s0(t)=sqrt(1.-t)*(1.456-0.157/sqrt(t)+1.764/t)
      r_1s0(u)=(0.2312+sqrt(0.7688**2+(0.1438*u)**2))**5.5*
     1         exp(3.427-sqrt(3.427**2+u**2))
c Murca_n:  n+n -> n+p+e+nu
      rmurca_n_p1s0(u)=exp(3.4370-sqrt((3.4370)**2+(1.*u)**2))*
     2       0.5*( (0.1477+sqrt((0.8523)**2+(0.1175*u)**2))**7.5 +
     3             (0.1477+sqrt((0.8523)**2+(0.1297*u)**2))**5.5 )
      rmurca_n_n1s0(u)=exp(5.3390-sqrt((5.3390)**2+(2.*u)**2))*
     2             (0.2414+sqrt((0.7586)**2+(0.1318*u)**2))**7.0
c ***************************
c TRIPLET PAIRING:
      u_3p2b(t)=sqrt(1.-t)*(0.7893+1.188/t)
      r_3p2b(u)=(0.2546+sqrt(0.7454**2+(0.1284*u)**2))**5*
     1         exp(2.701-sqrt(2.701**2+u**2))
c Murca_p:  n+p -> p+p+e+nu
      rmurca_p_n3p2b(u)=exp(2.3980-sqrt((2.3980)**2+(1.*u)**2))*
     2       0.5*( (0.1612+sqrt((0.8388)**2+(0.1117*u)**2))**7 +
     3             (0.1612+sqrt((0.8388)**2+(0.1274*u)**2))**5 )
c Murca_n:  n+n -> n+p+e+nu
      rmurca_n_n3p2b(u,t)=39.1*t*exp(-1.188/t)*rmurca_p_n3p2b(u)
c This `rmurca_p_n3p2b' is exact in the limit t<<1 and is only approximate
c when t~<1 since Yakovlev and Levenfish only calculated the t<<1 limit.
c****************************
      length=1.e6      ! length scale for diffusion
c****************************


c****************************
c       m_p/tau_p:
c****************************
       movertau_p=(mstp(i)*1.672e-24)*
     1            4.7e18*(temp/1.e9)**2*(2.8e14/rrho(i))**(1./3.)
c Neutron pairing suppression:
      if(temp.lt.tcn(i)) then
       if (i.ge.isf) then 
        tt=temp/tcn(i)
        u=u_1s0(tt)
        rn=r_1s0(u)
       else
        tt=temp/tcn(i)
        u=u_3p2b(tt)
        rn=r_3p2b(u)
       end if
      else
       rn=1.0 
      end if
c Proton pairing suppression:
      if(temp .lt. tcp(i))then
       tt=temp/tcp(i)
       u=u_1s0(tt)
       rp=r_1s0(u)
      else
       rp=1.0
      end if
c ****
      rdir=min(rn,rp)
      movertau_p=movertau_p*rdir
c****************************
c         a^2:
c****************************
c FROM MODIFIED URCA PROCESSES:
      lambda=3.6e34*(yelect(i)*bar(i)/0.16)**(1./3.)*(temp/1.e9)**6
c Effect of superfluidity :
      if(temp .lt. tcn(i))then
       if (i.ge.isf) then
        tt=temp/tcn(i)
        u=u_1s0(tt)
        rmurca_n_n=rmurca_n_n1s0(u)
       else
        tt=temp/tcn(i)
        u=u_3p2b(tt)
        rmurca_n_n=rmurca_n_n3p2b(u,tt)
       end if
      else
       rmurca_n_n=1.0
      end if
      if(temp .lt. tcp(i))then
        tt=temp/tcp(i)
        u=u_1s0(tt)
        rmurca_n_p=rmurca_n_p1s0(u)
      else
       rmurca_n_p=1.0
      end if
      rmurca_n=min(rmurca_n_p,rmurca_n_n)
      lambda=lambda*rmurca_n
c FROM DIRECT URCA PROCESSES:
      if (inu_durca.eq.1) then
       if (idurca_np(i).eq.1)then
        lambda=5.1d40*(yelect(i)*bar(i)/0.16)**(1./3.)*(temp/1.e9)**4
        lambda=lambda*rdir
       end if
      end if
c****
      nc=(yelect(i)+ymuon(i))*bar(i)*1.d39
      a2=nc/lambda/movertau_p
c****************************
c      t_amb:
c****************************
      t=4.*pi*nc/field**2*movertau_p
      t_amb_s=t*(length**2)
      t_amb_i=t*(length**2+a2)
c****
      heat=field**2/4./pi/t_amb_i

       return
      end
