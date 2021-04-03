c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine neebrem(i,T,mu,qeebrem,qasync,ProcessID)
c *********************************************************************
c *  calculate the energy loss rate per cubic centimeter              *
c *  from the electron-electron neutrino pair bremsstrahlung.         *
c *                                                                   *  
c *  From Jaikumar et al.                                             *
c *                                                                   *
c *  Mu=mu_e in Mev    T in K                                         *
c *                                                                   *
c *********************************************************************
       implicit real*8 (a-h,k-z)
c       INCLUDE 'Base_Dir.inc.f'
       INCLUDE 'size.inc.f'
       INCLUDE 'pid.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'

       parameter (pi = 3.14159265d0)
       integer :: ProcessID
       dimension logt(56),nalpha(56),n2(56)
       save logt,nalpha,n2,i_do_it
c***************************************************
        if(i_do_it.eq.314159) goto 666
        open(unit=20,file='Code/Data_Files/Nalpha.dat',status='old')
         do i=1,56
          read(20,*)temp,nalpha(i)
          logt(i)=log10(temp)+9.d0       ! converts T9 to T
         end do
         call SPLINE(logt,nalpha,56,1.d30,1.d30,n2)
        close(unit=20)
        i_do_it=314159
 666    continue
c***************************************************
        mu10=mu/10.d0
        alpha=T/mu10                     ! converts mu to mu10
        logalpha=log10(alpha)
        call SPLINT(logt,nalpha,n2,56,logalpha,naa)

        qeebrem=2.16d14 * (T/1.d9)**7 * (mu/10.d0)**2 * naa




c *************************** e-Ion

        GammaI = 22.73 * z_ion(i)**2 * 1.d6 / t
     &       * ( rrho(i)/1.d6 / a_ion(i) )**(1.d0/3.d0)
        xState = LOG10( rrho(i) )

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
     &       0.743902d0*(GammaI/1.d3)**(2.d0)
        else
         u = 0.672409d0 + 0.182774d0*GammaI/1.d3 +
     &       0.144817d0*(GammaI/1.d3)**(2.d0)
        endif

        FeIon = 10.d0**( a + b*xState**2.d0 + c*xState**4.d0 - 1.d0+u )

        qaebremIon = 10.8d0 * rrho(i) * gaee**2.d0 / (4.d0*pi*1.d-26) 
     &           * T/1.d8 * FeIon


c      write(*,*) qaebremIon,T,FeIon,xState,GammaI

c ***************************
      qasync = 0d0
      if (IAND(pid_eI_crust,ProcessID).gt.0) then
       qasync = qasync + qaebremIon
      endif

      return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine npb_new (temp,rho,qnpb)
c *********************************************************************
c *  calculate the energy loss rate per cubic centimeter in the crust *
c *  of a neutron star from neutrino pair bremsstrahlung.             *
c *                                                                   *  
c *  From Kaminker et al, A&A 343 (1999), p. 1009, Equ. (40)          *
c *                                                                   *
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'rho_limits.inc.f'
       INCLUDE 'gamma_limits.inc.f'
        tau = log10(temp/1.d8)
        r   = log10(rho/1.d12)
        rho0= 2.8d14
        lgq=11.204+7.304*tau+0.2976*r
     1      -0.370*tau**2+0.188*tau*r-0.103*r**2+0.0547*tau**2*r
     2      -6.77*log10(1.d0+0.228*rho/rho0)
        qnpb=10.d0**lgq
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************

      subroutine npb (t,rho,a,z,qnpb)
     
c *********************************************************************
c *  calculate the energy loss rate per cubic centimeter in the crust *
c *  of a neutron star from neutrino pair bremsstrahlung.             *
c *                                                                   *  
c *                                                                   *
c *  checked on February 27, 1996 against figures of Itoh et al 1996  *
c *                                                                   *
c *********************************************************************

      implicit real*8 (a-h,k-z)
      INCLUDE 'rho_limits.inc.f'
      INCLUDE 'gamma_limits.inc.f'

      mm=(rho/1.d6*z/a)**(2.d0/3.d0)
      tf=5.930d9*((1.d0+1.018d0*mm)**.5d0-1.d0)
      gamma=2.273d5*z**2*(rho/a)**(1.d0/3.d0)/t
     
      if (rho.le.1.d4) then
       call npbpde(t,rho,a,z,qnpb)
c For rho < 10^4 only npb_pde is accurately calculated (?)
c so it is used as a default
      else if (t .gt. .35*tf) then
       call npbpde(t,rho,a,z,qnpb)
      else if (t .ge. .30d0*tf) then
       call npbpde(t,rho,a,z,qnpb1)
       call npbl  (t,rho,a,z,qnpb2)
       qnpb=(t-.25d0*tf)/(.1d0*tf)*qnpb1+(.35d0*tf-t)/(.1d0*tf)*qnpb2
      else if (gamma .lt. gammaliq) then
       call npbl  (t,rho,a,z,qnpb)
      else if (gamma .gt. gammacryst) then
       call npbc(t,rho,a,z,qnpb)
      else
       call npbl  (t,rho,a,z,qnpb1)
       call npbc  (t,rho,a,z,qnpb2)
       qnpb=(gammacryst-gamma)/(gammacryst-gammaliq)*qnpb1+
     1        (gamma-gammaliq)/(gammacryst-gammaliq)*qnpb2
      endif
c*********************************
c For extremely impure crust:
c      call npbl  (t,rho,a,z,qnpb)
c*********************************

      return
      end
     
c *********************************************************************
c *********************************************************************

      subroutine npbpde (t,rho,a,z,qnpbpde)
     
c *********************************************************************
c *  calculate the energy loss rate per cubic centimeter in the crust *
c *  of a neutron star when electrons are partially degenerate        *
c *  (T>0.3Tf) from neutrino pair bremsstrahlung.                     *
c *                                                                   *
c *  from H. Munakata, Y. Kohyama & N. Itoh,                          *
c *  Ap. J. 316 (1987): p. 708.                                       *
c *                                                                   *
c *  which is almost identical to (typing errors corrected):          *
c *                                                                   *
c *  N. Itoh, H. Hayashi, A. Nishikawa & Y. Kohyama                   *
c *  Ap. J. Suppl. XXX (1996): p. yyy                                 *
c *                                                                   *
c *  checked on February 27, 1996 against figures of Itoh et al 1996  *
c *                                                                   *
c *********************************************************************
     
      implicit real*8(a-h,k-z)
     
      n=2.
c     n = number of neutrinos, other than electron neutrino, which can be
c         considered as massless at temperature t.
     
      gamma=2.273d5*z**2*(rho/a)**(1.d0/3.d0)/t
      t8=t/1.d8
     
      eta=(rho*z/a) / (7.05d6*(t8**1.5d0)+5.12d4*(t8**3))
      f1=23.5d0+6.83d4/(t8**2)+7.81d8/(t8**5)
      f2=1.d0+1.47d0/eta+0.0329d0/(eta**2)
      f=1.d0/f1+1.26d0*(1.d0+1.d0/eta)/f2
     
      b3=7.75d5*(t8**1.5d0)+247.d0*(t8**3.85d0)
      b4=4.07d0+.0240d0*(t8**1.40d0)
      b5=4.59d-5/(t8**0.11d0)
      g1=230.d0+6.7d5/(t8**2)+7.66d9/(t8**5)
      g2=b3/rho/z*a+b4+b5*(rho*z/a)**0.656d0
      g=1.d0/(1.d0+1d-9*rho*z/a)/g1+1.d0/g2
     
      qnpbpde=.5738d0*(z**2)/a*(t8**6)*rho*
     1        (.5d0*(1.122d0+.254d0*n)*f-.5d0*(.622d0-.246d0*n)*g)

      return
     
      end

c *********************************************************************
c *********************************************************************
     
      subroutine npbl (t,rho,a,z,qnpbl)
     
c *********************************************************************
c * calculate the energy loss rate per cubic centimeter in the crust  *
c * of a neutron star in liquid or amorphous glassy state (gamma<210) *
c * from neutrino-pair bremsstrahlung.                                *
c *                                                                   *
c * from N. Ntoh & Y. Kohyama,                                        *
c * Ap. J. 275 (1983): p. 858.                                        *
c *                                                                   *
c * checked on February 27, 1996 against figures of Itoh et al 1996   *
c *********************************************************************
     
      implicit real*8(a-h,k-z)
c      INCLUDE 'files_phys.inc.f'
      INCLUDE 'files.inc.f'
      parameter (jnpbl=1050)

      dimension fu1(0:jnpbl),fu160(0:jnpbl),
     1          gu1(0:jnpbl),gu160(0:jnpbl),
     2          alpha(0:3),beta(0:3)
      save fu1,fu160,gu1,gu160,alpha,beta
      save read


      data alpha/-0.07913,0.13177,1.74940,-0.80199/
      data beta/-0.06778,0.06268,1.78740,-0.78226/
     
      if (read.eq.1.) goto 500

c ***** read the npbl2 file *******************************************

      open(unit=12,file=f_npbl,status='old')
       read(12,*) jtext,jmax
       do 120 j=1,jtext
        read(12,*)
120    continue
       do 121 j=0,jmax
        read(12,*)a1,i2,i3,fu1(j),fu160(j),gu1(j),gu160(j)
121    continue
      close(unit=12,status='keep')

      read=1.

c *********************************************************************

500   continue

      n=2.
c     n = number of neutrinos, other than electron neutrino, which can be
c         considered as massless at temperature t.
     
      gamma=2.273d5*z**2*(rho/a)**(1.d0/3.d0)/t
      jfit=int(100.d0*log10(rho/1.d4))
      if (jfit.ge.jmax) jfit=jmax

      v=0.d0
      w=0.d0
      gamthird = gamma**(1.d0/3.d0)
      do im=0,3
       v=v+alpha(im)/gamthird**im
       w=w+beta(im)/gamthird**im
      end do
     
      f=v*fu1(jfit)+(1.d0-v)*fu160(jfit)
      g=w*gu1(jfit)+(1.d0-w)*gu160(jfit)
     
      qnpbl=rho*.5738d0*(z**2)/a*(1.d-8*t)**6*
     1      (.5d0*(1.122d0+.254d0*n)*f-.5d0*(.622d0-.246d0*n)*g)
     
      return
     
      end
     
c *********************************************************************
c *********************************************************************
     
      subroutine npbc (t,rho,a,z,qnpbc)
     
c *********************************************************************
c * calculate the energy loss rate per cubic centimeter in the crust  *
c * of a neutron star in crystalline state (gamma>170) from           *
c * neutrino-pair bremsstrahlung.                                     *
c *                                                                   *
c * N. Itoh, Y. Kohyama, N. Matsumoto & M. Seki,                      *
c * Ap. J. 285 (1984): p. 304.                                        *
c * and erratum in Ap. J. 322 (1987): p. 584.                         *
c *                                                                   *
c * checked on February 27, 1996 against figures of Itoh et al 1996   *
c *********************************************************************
     
      implicit real*8(a-h,k-z)
      parameter (jnpbc=1050)

c      INCLUDE 'files_phys.inc.f'
      INCLUDE 'files.inc.f'

      dimension fu171(0:jnpbc),fu5000(0:jnpbc),f1u171(0:jnpbc),
     1          gu171(0:jnpbc),gu5000(0:jnpbc),g1u171(0:jnpbc),
     1          alpha1(0:3),beta1(0:3),alpha(0:3),beta(0:3)
      save fu171,fu5000,f1u171,gu171,gu5000,g1u171
      save alpha,beta,alpha1,beta1
      save read
     
      data alpha/0.6808,12.9514,-145.0630,289.6790/
      data beta/0.7398,10.9457,-124.3540,226.3790/
      data alpha1/.2463,-12.1573,190.2310,-552.5950/
      data beta1/0.3082,-13.2025,186.3220,-509.1340/
     
      if (read.eq.1.) goto 500

c ***** read the npbc file ********************************************
     
      open(unit=13,file=f_npbc,status='old')
       read(13,*) jtext,jmax
       do 130 j=1,jtext
        read(13,*)
130    continue
       do 131 j=0,jmax
        read(13,*)a1,i2,i3,fu171(j),fu5000(j),gu171(j),gu5000(j),
     1                      f1u171(j),g1u171(j)
131    continue
      close(unit=13,status='keep')

      read=1.

c *********************************************************************

500   continue

      n=2.
c     n= number of neutrinos, other than electron neutrinos, which can be
c        considered as massless at temperature t.
     
      gamma=2.273d5*z**2*(rho/a)**(1.d0/3.d0)/t
      jfit=int(100.d0*log10(rho/1.d4))
      if (jfit.ge.jmax) jfit=jmax

      v=0.d0
      w=0.d0
      v1=0.d0
      w1=0.d0
      gamthird = gamma**(1.d0/3.d0)
      do im=0,3
       rgamth = 1.d0/gamthird**im
       v=v+alpha(im)*rgamth
       w=w+beta(im)*rgamth
       v1=v1+alpha1(im)*rgamth
       w1=w1+beta1(im)*rgamth
      end do
     
      if (gamma.ge.5000.d0) then
       v=1.d0
       w=1.d0
      end if

      f=(1.-v)*fu171(jfit)+v*fu5000(jfit)+v1*f1u171(jfit)
      g=(1.-w)*gu171(jfit)+w*gu5000(jfit)+w1*g1u171(jfit)

      qnpbc=rho*.5738d0*(z**2)/a*(t*1d-8)**6*
     1      (.5d0*(1.122d0+.254d0*n)*f-.5d0*(.622d0-.246d0*n)*g)

      return
     
      end
     
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************
     
      subroutine npair(t,rho,a,z,qpair)
     
C     DOMAIN OF VALIDITY:
C
C     Density:      1 < rho/mu_e < 10^14 gm/cm3
C     Temperature:    10^7 K  <  T  <  10^11 K 
C
c *********************************************************************
c * Calculates the energy loss rate per cubic centimeter in the crust *
c * of aneutron star from pair neutrinos.                             *
c *                                                                   *
c * N. Itoh, T. Adachi, M. Nakagawa, Y. Kohyama & H. Munakata         *
c * Ap. J. 339 (1989): p. 354                                         *
c *                                                                   *
c * which is identical to:                                            *
c *                                                                   *
c * N. Itoh, H. Hayashi, A. Nishikawa & Y. Kohyama                    *
c * Ap. J. Suppl. XXX (1996): p. yyy                                  *
c *                                                                   *
c * checked on February 27, 1996                                      *
c * against figure 1 - 5 of Itoh et al 1996                           *
c *********************************************************************
     
      implicit real*8 (a-h,k-z)
      dimension apa(0:2),bpa_l(3),bpa_h(3)
     
      data apa  / +6.002d19 , +2.084d20 , +1.872d21 /
      data bpa_l/ +9.383d-1 , -4.141d-1 , +5.829d-2 /
      data bpa_h/ +1.2383   , -0.8141   , +0.0 /
      data cpa_l/ +5.5924  /
      data cpa_h/ +4.9924  /

      save apa,bpa,bpa_h
     
      fexp(x)=dexp(max(x,-7.d2))
c*****
      if (t.lt.1.d7) then
       qphoto=0.0d0
       return
      end if
c*****

      l=t/5.9302d9
      xi=(rho*z/a*1.d-9)**(1.d0/3.d0)/l
     
      n=2.
c     n= number of neutrinos, other than electron neutrinos, which can be
c        considered as massless at temperature t.
     
      if (t.lt.1.d10) then
       fpair=(apa(0)+apa(1)*xi+apa(2)*(xi**2))*fexp(-cpa_l*xi)/
     1       ((xi**3)+bpa_l(1)/l+bpa_l(2)/(l**2)+bpa_l(3)/(l**3))
      else
       fpair=(apa(0)+apa(1)*xi+apa(2)*(xi**2))*fexp(-cpa_h*xi)/
     1       ((xi**3)+bpa_h(1)/l+bpa_h(2)/(l**2)+bpa_h(3)/(l**3))
      end if
      g=1.d0-13.04d0*(l**2)+133.5d0*(l**4)+
     1  1534.d0*(l**6)+918.6d0*(l**8)
     
      qpa=(10.7480d0*(l**2)+0.3967d0*(l**.5d0)+1.0050d0)**(-1)*
     1    (1.+rho*z/a/(7.692d7*(l**3)+9.715d6*(l**.5d0)))**(-0.3d0)
      qpair=.5d0*(1.122d0+n*.254d0)*(1.+(.622d0-n*.246d0)/
     1      (1.122d0+n*.254d0)*qpa)*g*fexp(-2./l)*fpair
     

      return
     
      end

c *********************************************************************
c *********************************************************************
     
      subroutine nphoto(t,rho,a,z,qphoto)
     
C     DOMAIN OF VALIDITY:
C
C     Density:      1 < rho/mu_e < 10^11 gm/cm3
C     Temperature:    10^7 K  <  T  <  10^11 K 
C
c *********************************************************************
c * calculate the energy loss rate per cubic centimeter in the crust  *
c * of a neutron star from photo neutrinos.                           *
c *                                                                   *
c * N. Itoh, T. Adachi, M. Nakagawa, Y. Kohyama & H. Munakata         *
c * Ap. J. 339 (1989): p. 354                                         *
c *                                                                   *
c * which is identical to:                                            *
c *                                                                   *
c * N. Itoh, H. Hayashi, A. Nishikawa & Y. Kohyama                    *
c * Ap. J. Suppl. XXX (1996): p. yyy                                  *
c *                                                                   *
c * checked on February 27, 1996                                      *
c * against figure 1 - 5 of Itoh et al 1996                           *
c *********************************************************************
     
      implicit real*8 (a-h,l-z)
      parameter(pi=3.1415926535)
      dimension aph(0:2),bph(3),cph(0:2,0:6,3),dph(0:2,5,3)
     
      data bph/ +6.290e-3 , +7.483e-3 , +3.061e-4 /
      data (((cph(i,j,k),j=0,6),i=0,2),k=1,3) /
     1   +1.008E+11 , +0.000E+0 , +0.000E+0 , +0.000E+0 ,
     1                +0.000E+0 , +0.000E+0 , +0.000E+0 ,
     2   +8.156E+10 , +9.728E+8 , -3.806E+9 , -4.384E+9 ,
     2                -5.774E+9 , -5.249E+9 , -5.153E+9 ,
     3   +1.067E+11 , -9.782E+9 , -7.193E+9 , -6.936E+9 ,
     3                -6.893E+9 , -7.041E+9 , -7.193E+9 ,
     4   +9.889E+10 , -4.524E+8 , -6.088E+6 , +4.269E+7 ,
     4                +5.172E+7 , +4.910E+7 , +4.388E+7 ,
     5   +1.813E+11 , -7.556E+9 , -3.304E+9 , -1.031E+9 ,
     5                -1.764E+9 , -1.851E+9 , -1.928E+9 ,
     6   +9.750E+10 , +3.484E+10, +5.199E+9 , -1.695E+9 ,
     6                -2.865E+9 , -3.395E+9 , -3.418E+9 ,
     7   +9.581E+10 , +4.107E+8 , +2.305E+8 , +2.236E+8 ,
     7                +1.580E+8 , +2.165E+8 , +1.721E+8 ,
     8   +1.459E+12 , +1.314E+11, -1.169E+11, -1.765E+11, 
     8                -1.867E+11, -1.983E+11, -1.896E+11,
     9   +2.424E+11 , -3.669E+9 , -8.691E+9 , -7.967E+9 ,
     9                -7.932E+9 , -7.987E+9 , -8.333E+9 /
      data (((dph(i,j,k),j=1,5),i=0,2),k=1,3) /
     1   +0.000E+00 , +0.000E+0 , +0.000E+0 , +0.000E+0 , +0.000E+0 ,
     2   -1.879E+10 , -9.667E+9 , -5.602E+9 , -3.370E+9 , -1.825E+9 ,
     3   -2.919E+10 , -1.185E+10, -7.270E+9 , -4.222E+9 , -1.560E+9 ,
     4   -1.135E+08 , +1.256E+8 , +5.149E+7 , +3.436E+7 , +1.005E+7 ,
     5   +1.652E+09 , -3.119E+9 , -1.839E+9 , -1.458E+9 , -8.956E+8 ,
     6   -1.548E+10 , -9.338E+9 , -5.899E+9 , -3.035E+9 , -1.598E+9 ,
     7   +4.724E+08 , +2.976E+8 , +2.242E+8 , +7.937E+7 , +4.859E+7 ,
     8   -7.094E+11 , -3.697E+11, -2.189E+11, -1.273E+11, -5.705E+10,
     9   -2.254E+10 , -1.551E+10, -7.793E+9 , -4.489E+9 , -2.185E+9 /

      save aph,bph,cph,dph

      fexp(x)=dexp(max(x,-7.d2))
c*****
      if (t.lt.1.d7) then
       qphoto=0.0
       return
      end if
c*****

      l=t/5.9302e9
      xi=(rho*z/a*1.e-9)**(1./3.)/l

      if (t.lt.1.e8) then
       cphot=0.5654+dlog10(t/1.e7)
       k=1
       tau=dlog10(t/1.e7)
      else
       cphot=1.5654
       if (t.lt.1.e9) then
        k=2
        tau=dlog10(t/1.e8)
       else
        k=3
        tau=dlog10(t/1.e9)
       end if
      end if

      do i=0,2
       aph(i)=0.5*cph(i,0,k)+0.5*cph(i,6,k)*dcos(10.*pi*tau)
       do j=1,5
        aph(i)=aph(i)+cph(i,j,k)*dcos(5./3.*pi*float(j)*tau)
     1               +dph(i,j,k)*dsin(5./3.*pi*float(j)*tau)
       end do
      end do

      n=2.
c     n= number of neutrinos, other than electron neutrinos, which can be
c        considered as massless at temperature t.
     
      fphoto=(aph(0)+aph(1)*xi+aph(2)*(xi**2))/
     1      ((xi**3)+bph(1)/l+bph(2)/(l**2)+bph(3)/(l**3))
      for=1.875e8*l+1.653e8*(l**2)+8.499e8*(l**3)-1.604e8*(l**4)
      qph=0.666*(1+2.045*l)**(-2.066)/(1.+rho*z/a/for)
      qphoto=.5*(1.122+n*.254)*
     1       (1.-(.622-n*.246)/(1.122+n*.254)*qph)*
     2       (rho*z/a)*l**5*fphoto*fexp(-cphot*xi)

      return
     
      end

c *********************************************************************
c *********************************************************************
     
      subroutine nplasma(t,rho,a,z,qplasma)
     
C     DOMAIN OF VALIDITY:
C
C     Density:      1 < rho/mu_e < 10^14 gm/cm3
C     Temperature:    10^7 K  <  T  <  10^11 K 
C
c *********************************************************************
c * calculate the energy loss rate per cubic centimeter in the crust  *
c * of a neutron star from plasma neutrinos.                          *
c *                                                                   *
c * M. Haft, G. Raffelt & A. Weiss                                    *
c * Ap. J. 425 (1996): p. 222                                         *
c *                                                                   *
c * which is identical to:                                            *
c *                                                                   *
c * N. Itoh, H. Hayashi, A. Nishikawa & Y. Kohyama                    *
c * Ap. J. Suppl. XXX (1996): p. yyy                                  *
c *                                                                   *
c *  checked on February 27, 1996                                     *
c *  against figure 1 - 5 of Itoh et al 1996                          *
c *********************************************************************
     
      implicit real*8 (a-h,k-z)
     
      fexp(x)=dexp(max(x,-7.d2))

      if (z.eq.0.0d0) then
       qplasma=0.0d0
       return
      end if

      l=t/5.9302d9

      den = 1.d0+(1.019d-6*rho*z/a)**(2.d0/3.d0)
      gamma2=1.1095d11*(rho*z/a)/t**2/den**0.5d0
      gamma=dsqrt(gamma2)

      f_t = 2.4d0 + 0.6d0*gamma**0.5d0 + 
     1      0.51d0*gamma + 1.25d0*gamma**1.5d0
      f_l = (8.6d0*gamma2 + 1.35d0 * gamma**3.5d0) / 
     1      (225.d0 - 17.d0*gamma + gamma2)

      x = 1.d0/6.d0 * 
     1    (+17.5d0 + dlog10(2.d0*rho*z/a) - 3.d0*dlog10(t))
      y = 1.d0/6.d0 * 
     1    (-24.5d0 + dlog10(2.d0*rho*z/a) + 3.d0*dlog10(t))
      if ((abs(x).gt.0.7d0).or.(y.lt.0.0d0)) then
       f_xy = 1.0d0
      else
       first = 0.39d0 - 1.25d0*x - 0.35d0*dsin(4.5d0*x) 
     1         - 0.3d0*fexp(-(4.5d0*x + 0.9d0)**2)
       sec=y - 1.6d0 + 1.25d0*x
       sec=min(0.d0,sec)
       second = sec/(0.57d0 - 0.25d0*x)
     2          
       f_xy = 1.05d0 + first*fexp(-second**2)
      end if

      qplasma = 3.00d21 * l**9 * gamma**6 * fexp(-gamma) * 
     1          (f_l + f_t) * f_xy


      return
     
      end

c *********************************************************************
c *********************************************************************
     
      subroutine nplasma_old(t,rho,a,z,qplasma)
     
c ****************************************************************************
c *  calculate the energy loss rate per cubic centimeter in the crust of a   *
c *  neutron star from plasma neutrinos.                                     *
c *                                                                          *
c *  from h.munakata, y.kohyama & n.itoh, ap.j.296(1985),p.197               *
c *                                                                          *
c *  checked on oct. 15 1990                                                 *
c *                                                                          *
c ****************************************************************************
     
      implicit real*8 (a-h,k-z)
      INCLUDE 'rho_limits.inc.f'
      dimension apl(0:2),bpl(3)
      save apl,bpl,cpl
     
      data apl/2.320e-7, 8.449e-8,1.787e-8/
      data bpl/2.581e-2, 1.734e-2,6.990e-4/
      data cpl/0.56457/
     
      fexp(x)=dexp(max(x,-7.d2))

      if (z.eq.0.0) then
       qplasma=0.0
       return
      end if

      l=t/5.930e9
      xi=(rho*z/a*1.e-9)**(1./3.)/l
     
      n=2.
c     n= number of neutrinos, other than electron neutrinos, which can be
c        considered as massless at temperature t.
     
      fplasma=(apl(0)+apl(1)*xi+apl(2)*(xi**2))*fexp(-cpl*xi)/
     1      ((xi**3)+bpl(1)/l+bpl(2)/(l**2)+bpl(3)/(l**3))
      qplasma=(.872+n*.004)*(rho*z/a)**3*fplasma

      return
     
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************

      subroutine nbub (i,t,rho,a,z,qbubble)
c *********************************************************************
c *********************************************************************
c *  calculate the energy loss rate per cubic centimeter in           *
c *  bubble phase of the crust.                                       *
c *                                                                   *
c *  from L. Leinson, ApJ 415, p. 759, 1993                           *
c *                                                                   *
c *  checked on xxx xx, 1994                                          *
c *                                                                   *
c *********************************************************************
       implicit real*8(a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
       parameter (pi = 3.14159265)
       INCLUDE 'pairing.inc.f'
c ***************************
       fexp(x)=dexp(max(x,-7.d2))

       u_1s0(t)=dsqrt(1.-t)*(1.456-0.157/dsqrt(t)+1.764/t)
       r_1s0(u)=(0.2312+dsqrt(0.7688**2+(0.1438*u)**2))**5.5*
     1          fexp(3.427-dsqrt(3.427**2+u**2))

       u_3p2B(t)=dsqrt(1.-t)*(5.596+8.424/t)
       r_3p2B(u)=(0.2546+dsqrt(0.7454**2+(0.01811*u)**2))**5*
     1          fexp(2.701-dsqrt(2.701**2+u**2/(16.*pi)))
c ***************************
        rhomin=1.e14
        if ((rho.lt.rhocore).and.(rho.ge.rhomin)) then
         qbubble= 1.1e22*(t/1.e9)**6
        else
         qbubble=0.0
        end if
c **** effect of superfluidity :
c (Note : isf is the zone where neutron pairing shifts from 3P2 to 1S0)
       if(t.lt.tcn(i)) then
        if (i.ge.isf) then 
         u=u_1s0(t/tcn(i))
         r=r_1s0(u)
        else
         u=u_3p2B(t/tcn(i))
         r=r_3p2B(u)
        end if
       else
        r=1.0
       end if

c       qbubble= qbubble*r

       return
      end

c *********************************************************************
c *********************************************************************
c *********************************************************************
c *********************************************************************

      subroutine nsynch(t,bfield,kfe,qsynch)
c *********************************************************************
c *********************************************************************
c * Calculate the energy loss rate per cubic centimeter in the crust  *
c * of a neutron star from synchrotron neutrinos.                     *
c *                                                                   *
c * From Bezchastnov, Haensel, Kaminker & Yakovlev,                   *
c *      A&A 328 (1997): p. 409                                       *
c *                                                                   *
c * Checked on Nov. 6, 1997                                           *
c *                                                                   *
c**********************************************************************
      implicit real*8 (a-h,k-z)
      parameter (pi=3.1415926535)

      data a1,b1,c1/2.036e-4,7.405e-8,3.675e-4/
      data a2,b2,c2,d2,e2/3.356e-3,1.536e-5,1.436e-2,
     1                    1.024e-5,7.647e-8/
      save a1,b1,c1
      save a2,b2,c2,d2,e2


      b13=bfield/1.e13
      x=kfe*197/0.511

      tp=2.02e9*b13*x**2
      xi=tp/t
      y1=((1.+3172.*xi**(2./3.))**(2./3.)-1.)**(3./2.)
      y2=((1.+172.2*xi**(2./3.))**(2./3.)-1.)**(3./2.)
      fp=44.01*(1.+c1*y1)**2/
     1   (1.+a1*y1+b1*y1**2)**4
      fm=36.97*(1.+c2*y2+d2*y2**2+e2*y2**3)/
     1   (1.+a2*y2+b2*y2**2)**5
      s_ab=27.*xi**4/pi**2/2**9/1.037*(fp-0.175/1.675*fm)

      tb=1.34e9*b13/dsqrt(1+x**2)
      z=tb/t
      d_1=1.+0.4228*z+0.1014*z**2+0.006240*z**3
      d_2=1.+0.4535*z**(2./3.)+0.03008*z-0.05043*z**2+0.004314*z**3
      s_bc=dexp(-z/2.)*d_1/d_2

      qsynch=9.04e14*b13**2*(t/1.e9)**5*s_ab*s_bc

      WRITE(*,*) 't=',t 
      WRITE(*,*) 'bfield=',bfield 
      WRITE(*,*) 'kfe=',kfe 
      WRITE(*,*) 'qsynch=',qsynch 

      return

      end
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
c**********************************************************************
