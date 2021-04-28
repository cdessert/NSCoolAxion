c ***********************************************************************
c ***********************************************************************
      subroutine conduct(i,T,rho,A,A1,Z,Q,magfield,
     1                   sigma,lambda,debug,
     2                   nu_e_s,nu_e_l)
c ***** CHECKED ON : ?
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'control_con.inc.f'
       INCLUDE 'profile_comp.inc.f'
c ***********************************************************************
        if (debug.ge.2.) print *,'Entering conduct:',
     1                   ' T, rho, A, A1, Z, Qimp =',T,rho,A,A1,Z,Q
c ***********************************************************************
        if  (rho.ge.rhocore) then
         if (istrange.eq.0) then
          if (i.le.isf) then
           isfn=3
          else
           isfn=1
          end if
          call con_core(icon_core,debug,
     1                  T,kfe(i),kfm(i),
     2                  kfp(i) ,mstp(i) ,Tcp(i) ,
     3                  kfn(i) ,mstn(i) ,Tcn(i) ,isfn,
     4                  kfla(i),mstla(i),Tcla(i),
     5                  kfsm(i),mstsm(i),Tcsm(i),
     6                  kfs0(i),msts0(i),Tcs0(i),
     7                  kfsp(i),mstsp(i),Tcsp(i),
     8                  fhad(i),
     9                  sigma,lambda,
     x                  nu_e_s,nu_e_l)
         else if (istrange.eq.1) then
          if (c_con_str.le.1.0) then
           call con_strange(i,T,lambda,debug)
          else
           lambda=c_con_str/(T/1.d9)**p_con_str
          end if
         else
          print *,'conduct: istrange not defined !'
          pause
         end if
        else
         call con_crust(icon_crust,debug,
     1                  T,rho,kfe(i),A,A1,Z,Q,sigma,lambda,
     2                  nu_e_s,nu_e_l)
        end if
c ***********************************************************************
        if (debug.ge.2.) print *,'Exiting conduct:',
     1                   ' sigma, lambda=',sigma,lambda
c ***********************************************************************
       return
      end
c ***********************************************************************
c ***********************************************************************

c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
c
c                QQQ    U   U   AAA   RRRR   K   K   SSSS
c               Q   Q   U   U  A   A  R   R  K  K   S
c               Q   Q   U   U  AAAAA  RRRR   KKK     SSS
c               Q   Q   U   U  A   A  R  R   K  K       S
c                QQQ Q   UUU   A   A  R   R  K   K  SSSS
c
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine con_strange(izone,T,lambda,debug)
       implicit real*8 (a-h,k-z)
       real*8 Ik
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'spec_heat.inc.f'
       INCLUDE 'pairing.inc.f'
       INCLUDE 'pairing_quark_PU2002.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'quark.inc.f'
       parameter(hbar=6.585d-22)     ! hbar in MeV sec
       parameter (pi=3.1415926535d0,z3=1.20205d0,z4=1.0823232d0)
       parameter (alpha_f=1.d0/137.d0,c=3.0d10)
       parameter (arad=7.56d-15,clight=2.99792e10)
       parameter (sigth=0.665d-24)
       dimension q(3)
       dimension nu_qq(3,3,3,3),nu_q(3,3),nu_eq(3,3),nu_qe(3,3)
       dimension cv(3,3),kf(3),red(3,3),lambda_q(3,3)
c      Indices: i=flavor j=color
c ***************************
c        fexp(x)=dexp(max(x,-5.d2))
c        u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
c        r_1s0(u)=(0.4186+dsqrt(1.007**2+(0.5010*u)**2))**2.5*
c     1            fexp(1.456-dsqrt(1.456**2+u**2))
c ***************************
        data q/0.6666666666d0,0.3333333333d0,-0.3333333333d0/
c ***************************
        if (debug.ge.2.) print *,'Entering subroutine `con_strange'' '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *,'subroutine con_strange:'
        print *,'   no guarantee that it returns a meaningful value !'
        pause
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ***************************
        nbar=bar(izone)
c INITIALIZE FERMI MOMENTA:
        kf_e=kfe(izone)
        kf(1)=kfqu(izone)
        kf(2)=kfqd(izone)
        kf(3)=kfqs(izone)
        muq=(kf(1)+kf(2)+kf(3))/3.d0 * 197.d0     ! quark    chemical potential
        mue=kf_e *  197.d0                        ! electron chemical potential
c SPECIFIC HEAT:
        cv_e=cve(izone) *T
        do j=1,3
         cv(1,j)=cvqu(izone)*T
         cv(2,j)=cvqd(izone)*T
         cv(3,j)=cvqs(izone)*T
        end do
c TEMPERATURE in MeV:
        TMeV=T/1.1604d10
c **********************************************************************
c NO PAIRING:
         if (normal.eq.1.d0) then
          cvquark=cv(1,1)+cv(1,2)+cv(1,3)+
     2            cv(2,1)+cv(2,2)+cv(2,3)+
     3            cv(3,1)+cv(3,2)+cv(3,3)
          qD=sqrt(6.d0/pi*alpha_c*muq**2)
          Ik=2.d0*z3*(TMeV/qD)**2
          nu=8.d0/pi**3 * 3.d0 * alpha_c**2 * muq**2/TMeV *Ik
          tau=hbar/nu
          lambda=1.d0/3.d0 * c**2 * cvquark * tau
c **********************************************************************
c CFL PHASE:
         else if (nbar.ge.nb_CFL) then
c          tt=T/Tc_CFL(i)
c          u=u_1s0(tt)
c          rr=exp(-2.d0*u)
c          nq=3.d0*bar(i)*1.d39 * rr
c          lambda=4.d0*arad*clight*T**3/sigth/nq
c          lambda=min(lambda,1.d30)
c          if ((izone.eq.1).or.(izone.eq.99).or.(izone.eq.199)) then
c           print '(0p1f10.5,1p2e12.3)',bar(i),T,lambda
c          end if
          lambda=1.d30*(T/1.d9)**3                    ! Photon transport 
c **********************************************************************
c 3SC PHASE:
         else if ((nbar.ge.nb_3SC).and.
     1            (nbar.lt.nb_CFL)) then  
          qD=sqrt(4.d0/pi*alpha_f*mue**2)
          Ik=2.d0*z3*(TMeV/qD)**2
          nu=12.d0/pi**3 * alpha_f**2 * mue**2/TMeV *Ik
          tau=hbar/nu
          lambda=1.d0/3.d0 * c**2 * cv_e * tau
c **********************************************************************
c 2SC PHASE:
         else
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is the correct value:
          qD=sqrt(2.d0/9.d0 * 2.d0/pi*alpha_c*muq**2)    ! Meissner mass
          Ik=6.d0*pi*z4*(TMeV/qD)**3
          nu=4.d0/pi**3 * alpha_c**2 * muq**2/TMeV *Ik
          tau=hbar/nu
          lambda=1.d0/3.d0 * c**2 * (cv(1,3)+cv(2,3)) * tau
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is from the normal phase:
c          cvquark=cv(1,1)+cv(1,2)+cv(1,3)+
c     2            cv(2,1)+cv(2,2)+cv(2,3)+
c     3            cv(3,1)+cv(3,2)+cv(3,3)
c          qD=sqrt(6.d0/pi*alpha_c*muq**2)
c          Ik=2.d0*z3*(TMeV/qD)**2
c          nu=8.d0/pi**3 * 3.d0 * alpha_c**2 * muq**2/TMeV *Ik
c          tau=hbar/nu
c          lambda=1.d0/3.d0 * c**2 * cvquark * tau
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end if
c *********************************************************************
c OLD SUBROUTINE FROM USOV's PAPER:
c        nb=bar(izone)
c        lambda_old= 6.d20/alpha_c*(nb/0.16d0)**(2.0/3.0)
c        lambda=lambda_old
c *********************************************************************
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        lambda=1.d35
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (debug.ge.2.) print *,'Exiting subroutine `con_strange'' '
       return
      end
c *********************************************************************
c *********************************************************************

