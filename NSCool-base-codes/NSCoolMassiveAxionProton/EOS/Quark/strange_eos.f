      program strange5
      implicit real*8 (a-h,l-z)

c    
c     Alpha_c is as in Fahri & Jaffe
c
c
c      parameter(alpha_c=0.3d0,m_s=150.d0,bbag=145.d0**4,renorm=313.d0)
      parameter(epsilon=1.d-10,delta=1.d-12)
      parameter(pi=3.1415926535d0)
      parameter(MeV=197.33d0)
      parameter(qe=-1.d0,qu=2.d0/3.d0,qd=-1.d0/3.d0,qs=-1.d0/3.d0)
      parameter(ia_max=150)

      character*26 filename
      character*3 durca_ud,durca_ud1,durca_ud2
      character*3 durca_us,durca_us1,durca_us2
      dimension func(4),dfunc(4,4)
      dimension a(4,4),b(4,1)

      dimension nbar_p(0:ia_max),rho_p(0:ia_max),pres_p(0:ia_max),
     2  ne_p(0:ia_max),nu_p(0:ia_max),nd_p(0:ia_max),ns_p(0:ia_max),
     3  mue_p(0:ia_max),muu_p(0:ia_max),mud_p(0:ia_max),mus_p(0:ia_max)

c     **************************************************
      eq1(muu,mud,mus,mue)=mud-mus
      eq2(muu,mud,mus,mue)=muu+mue-mus
      eq3(muu,mud,mus,mue)=2.d0*n_u(muu,alpha_c)-n_d(mud,alpha_c)-
     1                     n_s(m_s,mus,alpha_c,renorm)-3.d0*n_el(mue)
      eq4(muu,mud,mus,mue)=n_u(muu,alpha_c)+n_d(mud,alpha_c)+
     1                     n_s(m_s,mus,alpha_c,renorm)-3.d0*nbar
      eq5(muu,mud,mus,mue)=omega_e(mue)+
     1                     omega_u(muu,alpha_c)+omega_d(mud,alpha_c)+
     2                     omega_s(m_s,mus,alpha_c,renorm)+bbag

c     ***************************************************
c     ***************************************************
      renorm=1.5862d0*MeV
c      alpha_c=0.0
c ****
c 123  write(6,*)'alpha_c, B1/4, m_s ='
c      read(5,*)alpha_c,bbag0,m_s
       alpha_c=0.900000000001d0
c 123   write(6,'(a10,0p1f3.1)')'alpha_c =',alpha_c
 123   write(6,'(a10,i4)')'alpha_c =',int(10.d0*alpha_c)
       write(6,*)'B1/4, m_s ='
       read(5,*)bbag0,m_s
       bbag=bbag0**4
c ****
c     ***************************************************
c     ***************************************************
c      renorm=1.5862d0*MeV
c      m_s=150.d0
c      alpha_c=0.3
c      bbag=(145.d0)**4
cc      bbag=64.206*MeV**3
c     ***************************************************
c     ***************************************************

c     ***************************************************
c     ***************************************************
c     ***************************************************
c     ***************************************************
c     CALCULATE NBAR FOR PRESSURE=0:

       nbar0=0.3*MeV**3
       mu0=(nbar0*pi**2/(1.-2.*alpha_c/pi))**(1./3.)
       muu=mu0
       mud=mu0
       mus=mu0
       mue=1.

 201   continue

       func(1)=eq1(muu,mud,mus,mue)
       func(2)=eq2(muu,mud,mus,mue)
 1     func(3)=eq3(muu,mud,mus,mue)
       func(4)=eq5(muu,mud,mus,mue)

c     ***************************************************
      
       dfunc(1,1)=-(eq1(     muu      ,mud,mus,mue)-
     1              eq1(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(1,2)=-(eq1(muu,     mud      ,mus,mue)-
     1              eq1(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(1,3)=-(eq1(muu,mud,     mus      ,mue)-
     1              eq1(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(1,4)=-(eq1(muu,mud,mus,     mue      )-
     1              eq1(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(2,1)=-(eq2(     muu      ,mud,mus,mue)-
     1              eq2(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(2,2)=-(eq2(muu,     mud      ,mus,mue)-
     1              eq2(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(2,3)=-(eq2(muu,mud,     mus      ,mue)-
     1              eq2(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(2,4)=-(eq2(muu,mud,mus,     mue      )-
     1              eq2(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(3,1)=-(eq3(     muu      ,mud,mus,mue)-
     1              eq3(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(3,2)=-(eq3(muu,     mud      ,mus,mue)-
     1              eq3(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(3,3)=-(eq3(muu,mud,     mus      ,mue)-
     1              eq3(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(3,4)=-(eq3(muu,mud,mus,     mue      )-
     1              eq3(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(4,1)=-(eq5(     muu      ,mud,mus,mue)-
     1              eq5(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(4,2)=-(eq5(muu,     mud      ,mus,mue)-
     1              eq5(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(4,3)=-(eq5(muu,mud,     mus      ,mue)-
     1              eq5(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(4,4)=-(eq5(muu,mud,mus,     mue      )-
     1              eq5(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

c     ***************************************************

       do i=1,4
        b(i,1)=-func(i)
        do j=1,4
         a(i,j)=dfunc(i,j)
        end do  
       end do

       call GAUSSJ(a,4,4,b,1,1)

       muu=muu+b(1,1)
       mud=mud+b(2,1)
       mus=mus+b(3,1) 
       mue=mue+b(4,1)

       bb=sqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2+b(4,1)**2)
       mm=sqrt(muu**2+mud**2+mus**2+mue**2)

c       print *,bb/mm

c      Check if solution found:
       if ((bb/mm).ge.epsilon) goto 201
c     ***************************************************
c      Solution found: get everything:
       mue=mue
       ne=n_el(mue)
       oe=omega_e(mue)

       nu=n_u(muu,alpha_c)
       ou=omega_u(muu,alpha_c)
     
       nd=n_d(mud,alpha_c)
       od=omega_d(mud,alpha_c)

       ns=n_s(m_s,mus,alpha_c,renorm)
       os=omega_s(m_s,mus,alpha_c,renorm)

       nbar=(nu+nd+ns)/3.d0

       pres=-(oe)-(ou)-(od)-(os)-bbag

       rho=(oe+mue*ne)+(ou+muu*nu)+(od+mud*nd)
     1     +(os+mus*ns)+bbag

c KEEP RESULTS & CHANGE UNITS:
       E_over_A=rho/nbar

       if (E_over_A.gt.939) then
        print *,'Not a SM model: END !'
        goto 666
       end if

       nbar_p(0)=nbar/MeV**3                            ! fm^-3
       rho_p(0) =rho /MeV**3*(1.d39*1.6d-6/9.d20)       ! g cm^-3
c       pres_p(0)=pres/MeV**3*(1.d39*1.6d-6)            ! erg cm^-3
       pres_p(0)=1.0d0                                  ! erg cm^-3
       ne_p(0)  =ne /MeV**3                             ! fm^-3
       nu_p(0)  =nu /MeV**3                             ! fm^-3
       nd_p(0)  =nd /MeV**3                             ! fm^-3
       ns_p(0)  =ns /MeV**3                             ! fm^-3
       mue_p(0) =mue                                    ! MeV
       muu_p(0) =muu                                    ! MeV
       mud_p(0) =mud                                    ! MeV
       mus_p(0) =mus                                    ! MeV
c       ye=ne/nbar
c       yu=nu/nbar
c       yd=nd/nbar
c       ys=ns/nbar

c     ***************************************************
c     ***************************************************
c     ***************************************************
c     ***************************************************
c     CALCULATE THE EQUATION OF STATE:

       do ia=1,ia_max

       nbar=nbar*(1.0d0+1.d-10*float(ia**4))  ! This starts from the nbar at P=0

 202   mu0=(nbar*pi**2/(1.-2.*alpha_c/pi))**(1./3.)
       muu=mu0
       mud=mu0
       mus=mu0
       mue=1.

 200   continue

       func(1)=eq1(muu,mud,mus,mue)
       func(2)=eq2(muu,mud,mus,mue)
       func(3)=eq3(muu,mud,mus,mue)
       func(4)=eq4(muu,mud,mus,mue)

c     ***************************************************
      
       dfunc(1,1)=-(eq1(     muu      ,mud,mus,mue)-
     1              eq1(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(1,2)=-(eq1(muu,     mud      ,mus,mue)-
     1              eq1(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(1,3)=-(eq1(muu,mud,     mus      ,mue)-
     1              eq1(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(1,4)=-(eq1(muu,mud,mus,     mue      )-
     1              eq1(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(2,1)=-(eq2(     muu      ,mud,mus,mue)-
     1              eq2(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(2,2)=-(eq2(muu,     mud      ,mus,mue)-
     1              eq2(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(2,3)=-(eq2(muu,mud,     mus      ,mue)-
     1              eq2(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(2,4)=-(eq2(muu,mud,mus,     mue      )-
     1              eq2(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(3,1)=-(eq3(     muu      ,mud,mus,mue)-
     1              eq3(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(3,2)=-(eq3(muu,     mud      ,mus,mue)-
     1              eq3(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(3,3)=-(eq3(muu,mud,     mus      ,mue)-
     1              eq3(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(3,4)=-(eq3(muu,mud,mus,     mue      )-
     1              eq3(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

       dfunc(4,1)=-(eq4(     muu      ,mud,mus,mue)-
     1              eq4(muu*(1.+delta),mud,mus,mue))/
     2              (delta*muu)
       dfunc(4,2)=-(eq4(muu,     mud      ,mus,mue)-
     1              eq4(muu,mud*(1.+delta),mus,mue))/
     2              (delta*mud)
       dfunc(4,3)=-(eq4(muu,mud,     mus      ,mue)-
     1              eq4(muu,mud,mus*(1.+delta),mue))/
     2              (delta*mus)
       dfunc(4,4)=-(eq4(muu,mud,mus,     mue      )-
     1              eq4(muu,mud,mus,mue*(1.+delta)))/
     2              (delta*mue)

c     ***************************************************

       do i=1,4
        b(i,1)=-func(i)
        do j=1,4
         a(i,j)=dfunc(i,j)
        end do  
       end do

       call GAUSSJ(a,4,4,b,1,1)

       muu=muu+b(1,1)
       mud=mud+b(2,1)
       mus=mus+b(3,1) 
       mue=mue+b(4,1)

       bb=sqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2+b(4,1)**2)
       mm=sqrt(muu**2+mud**2+mus**2+mue**2)

       if ((bb/mm).ge.epsilon) goto 200
c ***************************************************
       ne=n_el(mue)
       oe=omega_e(mue)
       nu=n_u(muu,alpha_c)
       ou=omega_u(muu,alpha_c)
       nd=n_d(mud,alpha_c)
       od=omega_d(mud,alpha_c)
       ns=n_s(m_s,mus,alpha_c,renorm)
       os=omega_s(m_s,mus,alpha_c,renorm)
       pres=-(oe)-(ou)-(od)-(os)-bbag
       if (pres.le.0.d0) goto 999
       rho=(oe+mue*ne)+(ou+muu*nu)+(od+mud*nd)
     1     +(os+mus*ns)+bbag
c **** Change Units:
       nbar1=nbar/MeV**3
       rho1=rho/MeV**3
       pres1=pres/MeV**3
c       ye=ne/nbar
c       yu=nu/nbar
c       yd=nd/nbar
c       ys=ns/nbar
c **** Store results into arrays:
       nbar_p(ia)=nbar/MeV**3
       rho_p(ia) =rho/MeV**3*(1.d39/9.d20*1.6d-6)
       pres_p(ia)=pres/MeV**3*(1.d39*1.6d-6)
       ne_p(ia)   =ne/MeV**3
       nu_p(ia)   =nu/MeV**3
       nd_p(ia)   =nd/MeV**3
       ns_p(ia)   =ns/MeV**3
       mue_p(ia)  =mue
       muu_p(ia)  =muu
       mud_p(ia)  =mud
       mus_p(ia)  =mus

c       write(6,'(1p11e12.3)')
c     1  nbar_p(ia),rho_p(ia),pres_p(ia),
c     2  ne_p(ia),nu_p(ia),nd_p(ia),ns_p(ia),
c     3  mue_p(ia),muu_p(ia),mud_p(ia),mus_p(ia)

c***** Check if equations are solved:
      x1=(mud-mus)/mud
      x2=(muu+mue-mus)/muu
      x3=(2.*n_u(muu,alpha_c)-n_d(mud,alpha_c)-
     1   n_s(m_s,mus,alpha_c,renorm)-3.*n_el(mue))/nbar
      x4=(n_u(muu,alpha_c)+n_d(mud,alpha_c)+
     1   n_s(m_s,mus,alpha_c,renorm)-3.*nbar)/nbar

      end do

 999  continue

c ****PRINTS OUT RESULTS:

c Write down file name
       ia=int(10.d0*alpha_c)
       ib=int(bbag0)
       ic=int(m_s)
       open(unit=15,file='name.txt',status='new')
        if (ic.ge.100) then
         write(15,'(a13,1i1,a1,1i3,a1,1i3,a4)')
     1     'strange_eos_0',ia,'_',ib,'_',ic,'.dat'
        else if (ic.ge.10) then
         write(15,'(a13,1i1,a1,1i3,a2,1i2,a4)')
     1     'strange_eos_0',ia,'_',ib,'_0',ic,'.dat'
        else
         write(15,'(a13,1i1,a1,1i3,a3,1i1,a4)')
     1     'strange_eos_0',ia,'_',ib,'_00',ic,'.dat'
        end if
       close (unit=15,status='keep')
       open(unit=15,file='name.txt',status='old')
        read(15,'(a26)')filename
       close(unit=15,status='delete')
       print *, filename

      print '(a30,0p3f12.3)','Alpha_c, Bbag^1/4, Ms =',
     1       alpha_c,bbag0,m_s
      print '(a20,0p1f12.3)','E/A =',E_over_A

      open(UNIT=20,FILE=filename,STATUS='NEW')
      write(20,'(0p4f12.3,a60,0p1f12.3)')
     1       alpha_c,bbag0,m_s,renorm,
     2      'Alpha_c    Bbag^1/4     Ms      Pr  ===>  E/A at P=0:',
     3      E_over_A
      write(20,'(0p1f12.3,a80)')
      write(20,*)
      write(20,'(2a18,9a14)')
     1                      'nbar     ','rho      ','presion  ',
     2                      'ne     ','nu     ','nd     ','ns     ',
     3                      'mue    ','muu    ','mud    ','mus    '
      write(20,'(2a18,9a14)')
     1      '[fm^-3]    ','[gm cm^-3]   ','[erg cm^-3]',
     2      '[fm^-3]  ','[fm^-3]  ','[fm^-3]  ','[fm^-3]  ',
     3      ' [MeV]   ',' [MeV]   ',' [MeV]   ',' [MeV]   '
      write(20,*)

      do ia=ia_max,0,-1
       write(20,'(1p2e18.11,1p9e14.5)')
     1  nbar_p(ia),rho_p(ia),pres_p(ia),
     2  ne_p(ia),nu_p(ia),nd_p(ia),ns_p(ia),
     3  mue_p(ia),muu_p(ia),mud_p(ia),mus_p(ia)
c       write(6,'(i5,1p2e18.11,1p19e12.3)')ia,
c     1  nbar_p(ia),rho_p(ia),pres_p(ia),
c     2  ne_p(ia),nu_p(ia),nd_p(ia),ns_p(ia),
c     3  mue_p(ia),muu_p(ia),mud_p(ia),mus_p(ia)
      end do


      close(20)

 666  continue
      end

c     ************************************************
c     ************************************************
c     ************************************************

      function n_el(mu_el)
      implicit none
      real*8 n_el, mu_el, pi
      pi=dacos(-1.d0)
      n_el=mu_el**3/(3.d0*pi**2)
      return
      end
      
      function omega_e(mu_el)
      implicit none
      real*8 omega_e, mu_el, pi
      pi=dacos(-1.d0)
      omega_e=-mu_el**4/(12.d0*pi**2)
      return
      end
      
      function n_u(mu_u,alpha_c)
      implicit none
      real*8 n_u, mu_u, alpha_c, pi, a
      pi=dacos(-1.d0)
      a=1.d0-(2.d0*alpha_c)/(pi)
      n_u=(mu_u**3/(pi**2))*a
      return
      end
      
      function omega_u(mu_u,alpha_c)
      implicit none
      real*8 omega_u, mu_u, alpha_c, pi, a
      pi=dacos(-1.d0)
      a=1.-(2.d0*alpha_c)/(pi)
      omega_u=-(mu_u**4/(4.d0*pi**2))*a
      return
      end
      
      function n_d(mu_d,alpha_c)
      implicit none
      real*8 n_d, mu_d, alpha_c, pi, a
      pi=dacos(-1.d0)
      a=1.d0-(2.d0*alpha_c)/(pi)
      n_d=(mu_d**3/(pi**2))*a
      return
      end

      function omega_d(mu_d,alpha_c)
      implicit none
      real*8 omega_d, mu_d, alpha_c, pi, a
      pi=dacos(-1.d0)
      a=1.d0-(2.d0*alpha_c)/(pi)
      omega_d=-(mu_d**4/(4.d0*pi**2))*a
      return
      end
      
      function n_s(m_s,mu_s,alpha_c,renorm)
      implicit none
      real*8 n_s, m_s, mu_s, alpha_c, pi, b, x, renorm
      pi=dacos(-1.d0)
      x=mu_s**2-m_s**2
      b=mu_s-(3.d0*m_s**2/x**(.5))*(log((mu_s+x**(.5))/renorm))
      n_s=(x/(pi**2))*(x**(.5)-(2.d0*alpha_c/pi)*b)
      return
      end

      function omega_s(m_s,mu_s,alpha_c,renorm)
      implicit none
      real*8 omega_s, m_s, mu_s, alpha_c, renorm
      real*8 pi, as, bs, cs, ds, es, fs, gs, hs, x, y1, y2
      pi=dacos(-1.d0)
      x=mu_s**2-m_s**2
      y1=(mu_s+x**(.5))/m_s
      y2=(mu_s+x**(.5))/mu_s
      as=(1.d0/(4.d0*pi**2))
      bs=mu_s*x**(.5)*(mu_s**2-(2.5d0)*m_s**2)
      cs=(1.5d0)*m_s**4*log(y1)
      ds=3.d0
*(mu_s*(x**(.5))-(m_s**2)*log(y2))**2
      es=6.d0*log(renorm/mu_s)*
     1   (mu_s*(m_s**2)*(x**(.5))-(m_s**4)*log(y1))
      fs=2.d0*alpha_c/pi
      gs=3.d0*(m_s**4)*(log(m_s/mu_s))**2   !ver el log^2
      hs=ds-2.d0*(x**2)-gs+es
      omega_s=-as*(bs+cs-fs*hs)
      return
      end

c      function n_s2(omega_s,epsilon)
c      implicit none
c      real*8 n_s2, omega_s, epsilon, x
c      n_s2=(omega_s*(1.+epsilon)-omega_s)/epsilon
c      return
c      end

c     ************************************************
c     ************************************************
c     ************************************************

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      INTEGER M,MP,N,NP,NMAX
      REAL*8 A(NP,NP),B(NP,NP)
      PARAMETER (NMAX=50)
      INTEGER I,ICOL,IROW,J,K,L,LL,INDXC(NMAX),INDXR(NMAX),IPIV(NMAX)
      REAL*8 BIG,DUM,PIVINV
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

c     ************************************************
c     ************************************************
c     ************************************************
