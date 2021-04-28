c *********************************************************************
c *********************************************************************
c      program check_1
c       implicit real*8(a-h,k-z)
c 100   print *,'P,T,A,Z='
c       read(5,*)P,T,A,Z
c       print *,'Guess for Rho='
c       read(5,*)Rho
c       call density(T,P,A,Z,Rho)
c       print '(a10,1p1e18.6)','Rho =',Rho
c       call pressure(T,Rho,A,Z,Pre)       
c       print '(a10,1p1e18.6)','Pre =',Pre
c       goto 100
c      end
c *********************************************************************
c *********************************************************************
c      program check_2
c       implicit real*8(a-h,k-z)
c 100   print *,'rho,T, A,Z='
c       read(5,*)Rho,T, A,Z
c       nbar=Rho*6.023e23
c       call pressure(T,Rho,A,Z,Pre)       
cc       print '(a15,1p2e18.6)','Nbar, Pres =',nbar*1.d-39,Pre
c       ia=int(A)
c       iz=int(Z)
c       print '(1p1e10.3,1p1e17.3,1p1e16.3,1i9,2i10)',
c     1    nbar*1.d-39,rho,pre,iA,iA,iZ
c       goto 100
c      end
c *********************************************************************
c *********************************************************************
c      program check_3
c       implicit real*8(a-h,k-z)
c 100   print *,'rho, T ='
c       read(5,*)Rho,T
c       A=56
c       Z=26
cc       T=3.e8
c       nbar=Rho*6.023e23
c       call pressure(T,Rho,A,Z,Pre)       
cc       print '(a15,1p2e18.6)','Nbar, Pres =',nbar*1.d-39,Pre
c       ia=int(A)
c       iz=int(Z)
c       print '(1p1e10.3,1p1e17.3,1p1e16.3,1i9,2i10)',
c     1    nbar*1.d-39,rho,pre,iA,iA,iZ
c       goto 100
c      end
c *********************************************************************
c *********************************************************************
c      program plot_pressure
c       implicit real*8 (a-h,k-z)
c       real*4 xplot,yplot
c       dimension xplot(300),yplot(300)
cc ****** draw the box *******************************************    
c123     write(6,*)
c     1    'Screen (1), Laser printer (2) or GIF (3)'
c        read(5,*)output
c        if (output.eq.1.) then
c         call pgbegin (0,'/XSERV',1,1)
c        else if (output.eq.2.) then
c         call pgbegin (0,'EOS.ps/ps',1,1)
c        else if (output.eq.3.) then
c         call pgbegin (0,'EOS.gif/gif',1,1)
c        else
c         goto 123
c        end if      
cc***********************************************
c        call pgslw ( 3 )
c        call pgvsize (1.0,6.0,1.0,6.0)
c        call pgwindow (4.0,10.0,10.0,30.0)
c        call pgbox ('bclnts',1.0,5,'bcnltsv',1.0,5)
c        call pglabel ('\gr (g cm\u-3\d)',
c     1                'Pressure (erg cm\u-3\d)',' ')
c ****** plot the curves ********************************
c 100    print *,'T,A,Z='
c        read(5,*)T,A,Z
c        do ir=40,100
c         lrho=float(ir)/10.
c         Rho=10.**lrho
c         call pressure(T,Rho,A,Z,Pres)
c         xplot(ir)=log10(Rho)
c         yplot(ir)=log10(Pres)
c        end do
c        call pgsls (1)
c        call pgline (61,xplot(40),yplot(40))
c        print *,'Another one ? (0/1)'
c        read(5,*)inext
c        if (inext.eq.1) goto 100
c **********************************************************************
c        call pgend
c      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine density(T,P,A,Z,Rho)
c************************************************************
c given T,P,A,Z, calculates Rho: uses Newton's method
c and Rho as input is taken as a first guess
c************************************************************
c ***** checked
       implicit real*8(a-h,k-z)
        eps=1.d-3
c        print *,'--------------------------------'
c        print '(1a20,0p1f20.8)','--------->',Rho
 100    Rho0=Rho
        call pressure(T,Rho0,A,Z,Pre0)
        Rho1=(1.d0+eps)*Rho0
        call pressure(T,Rho1,A,Z,Pre1)
        f=Pre0-P
        f1=(Pre1-Pre0)/(Rho1-Rho0)
        dRho=-f/f1
        Rho=Rho0+dRho
c        print '(1a20,0p2f20.8,1p3e20.8)','--------->',Rho,Rho1,Pre0,Pre1
        if (abs(dRho/Rho).lt.1.d-5) return
        goto 100
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine pressure(T,Rho,A,Z,Pres)
c **** sept 1991 version, checked on sept 2
       implicit real*8(a-h,k-z)
        hb=1.054588d-27
        kb=1.380662d-16
        c=2.997924d10
        NA=6.022045d23
        me=9.109d-28
        pi=3.141592653d0
C **** Calculate the electron density and ionization
c        call n_electron(T,rho,A,Z,ne)
        ne=Rho*NA*Z/A
        nion=Rho*NA/A
        Zeff=ne/nion
C **** Calculate the ionic pressure
        gamma=2.273E5*Zeff**2*(Rho/A)**(1./3.)/T
        if (gamma.ge.210.) then
         Uion=1.5-0.895929*gamma+3225./gamma**2
        else
         Uion=-0.897744*gamma+0.95043*gamma**.25
     1        +0.18956/gamma**.25-0.81487
        end if
        Pion=nion*kb*T*(1.0d0+Uion/3.d0)      
C **** Get the electron pressure
        call P_electron(T,ne,Pel)
C *********************************
        Pres=Pel+Pion
       return
      end
C ********************************************************************
C ********************************************************************
      subroutine P_electron(T,ne,Pres) 
C*********************************************************************C
C  Calculates the pressure of a perfect electron gas (Fermi-Dirac)    C
C  From Eggleton, Faulkner & Flannery, A&A23 (1973), p. 325-330       C
C  for the approx. which gives thermodynamically consistent D,U & P   C
C                                                                     C
C      CHECKED on nov. 24 1990 and again on July 15 2005              C
C         that it gives the correct perfect gas and                   C
C    degenerate relativistic and non-relativistic limits :            C
C                                                                     C
C      PDEGNR=(3.*PI**2)**(2./3.)/5.*(HB*NE)*(HB/ME)*NE**(2./3.)      C
C      PDEGR=(3.*PI**2)**(1./3.)/4.*(HB*NE)*(  C  )*NE**(1./3.)       C
C      PPERF=NE*KB*T                                                  C
C*********************************************************************C
       implicit real*8(a-h,k-z)
       parameter (pi = 3.14159652d0,NA=6.022045d23)
       parameter(epsilon=1.d-12)
       dimension CD(0:3,0:3),CP(0:3,0:3),CU(0:3,0:3)
       dimension DS(0:3)
       save CP,CD,CU
       save OLDNE,OLDT,OLTF1
C *****
       DATA ((CD(I,J),I=0,3),J=0,3)
     1         / 2.315472, 7.128660, 7.504998, 2.665350,
     2           7.837752,23.507934,23.311317, 7.987465,
     3           9.215560,26.834068,25.082745, 8.020509,
     4           3.693280,10.333176, 9.168960, 2.668248/
       DATA ((CP(I,J),I=0,3),J=0,3)
     1         / 2.315472, 6.748104, 6.564912, 2.132280,
     2           7.837752,21.439740,19.080088, 5.478100,
     3           9.215560,23.551504,19.015888, 4.679944,
     4           3.693280, 8.859868, 6.500712, 1.334124/   
       DATA ((CU(I,J),I=0,3),J=0,3)
     1         / 3.473208,10.122156, 9.847368, 3.198420,
     2          16.121172,43.477194,37.852852,10.496830,
     3          23.971040,60.392810,47.782844,11.361074,
     4          11.079840,26.579604,19.502136, 4.002372/

C *****     
      NEHAT=NE/1.7595d30
C ***** CALCULATE F

      T1=T/5.93d9
      F1=1.D-3*(NE/NA)*(1.d7/T)**3
      IF ( (ABS(OLDNE-NE)/NE.LE.(5.D-1)).AND.
     1        (ABS(OLDT-T)/T.LE.(5.D-1)) )  F1=OLDF1
     
1111  F=ABS(F1)
      G=T1*SQRT(1.D0+F)
      PF=1.D0+F
      PG=1.D0+G

      DO J1=0,3
       DS(J1)=CD(J1,0)+(CD(J1,1)+CD(J1,2)*G+CD(J1,3)*G*G)*G
       DS(J1)=DS(J1)/PF**3
      END DO
      SUM1=DS(0)+(DS(1)+DS(2)*F+DS(3)*F*F)*F

      DO J1=1,3
       DS(J1)=CD(J1,0)+(CD(J1,1)+CD(J1,2)*G+CD(J1,3)*G*G)*G
       DS(J1)=DS(J1)/PF**3
      END DO
      SUM2=DS(1)+(2.*DS(2)+3.*DS(3)*F)*F

      DO J1=0,3
       DS(J1)=(CD(J1,1)+2.*CD(J1,2)*G+3.*CD(J1,3)*G*G)*G
       DS(J1)=DS(J1)/PF**3
      END DO
      SUM3=DS(0)+(DS(1)+DS(2)*F+DS(3)*F*F)*F

      IF(F.LT.1.)THEN
       COEF=F*G**1.5/PF/PG**1.5
       NEF=COEF*SUM1
       NEF1=(COEF/F-
     1       3.25*COEF/PF-
     2       .75*COEF*G/PF/PG) * SUM1 +
     3       COEF*(SUM2+0.5/PF*SUM3)
      ELSE
       COEF=(F/PF)*(G/PG)**1.5
       NEF=COEF*SUM1
       NEF1=COEF*
     1      ( (1./F-3.25/PF-0.75*G/(PF*PG))*SUM1+SUM2+0.5/PF*SUM3 )
      END IF

      F1=ABS(F+(NEHAT-NEF)/NEF1)
     
      IF ((ABS(NEF-NEHAT)/ABS(NEHAT).GT.EPSILON).OR.
     1    (ABS(F1-F)/ABS(F1).GT.EPSILON)) GOTO 1111

C *****  Calculate pressure
     
      G1=T1*SQRT(1.+F1)

      DO J1=0,3
       DS(J1)=CP(J1,0)+(CP(J1,1)+CP(J1,2)*G1+CP(J1,3)*G1*G1)*G1
      END DO
      SUMP1=DS(0)+(DS(1)+DS(2)*F1+DS(3)*F1*F1)*F1

      IF(F1.LT.1)THEN
       P1=1.44E24*F1*G1**2.5/(1.+F1)**4/(1.+G1)**1.5*SUMP1
      ELSE
       P1=1.44E24*(F1**.25/(1.+F1))**4*G1*(G1/(1.+G1))**1.5*SUMP1
      END IF
     
      PRES=P1
     
      OLDNE=NE
      OLDT=T
      OLDF1=F1
     
      RETURN
     
      END

C*********************************************************************
C*********************************************************************
C*********************************************************************

      SUBROUTINE N_ELECTRON(TEMPERATURE,DENSITY,A,Z,NE)

C ******************************************************************C
C     TEMPERATURE in K  DENSITY in gm/cm3  NE in #/cm3              C
C*******************************************************************C
C              CHECKED ON AUGUST 30 1991                            C
C ******************************************************************C

      IMPLICIT REAL*8(A-H,K-Z)
      INTEGER A,Z
      CHARACTER*80 FILENAME

      PARAMETER (ITEMP=50,JRHO=50,EV=11604.)

      DIMENSION TEMP(0:ITEMP),RHO(0:JRHO),NELECT(0:ITEMP,0:JRHO)
      SAVE READ,TEMP,RHO,NELECT

C ****** Read Electron Density Table *********************************

      IF (READ.EQ.1.) GOTO 1234
      FILENAME='[USER.DPAGE.NSTAR.LOSALAMOS]NELECT.DAT'
      OPEN(UNIT=40,FILE=FILENAME,STATUS='OLD')
       DO JUNK=1,6
        READ(40,*)
       END DO
       DO I=0,ITEMP
        READ(40,*)
        DO J=0,JRHO
         READ(40,*)TEMP(I),RHO(J),NELECT(I,J)
        END DO
       END DO
      CLOSE(UNIT=40,STATUS='KEEP')

      READ=1.

1234  CONTINUE

C *******************************************************************

      T=TEMPERATURE/EV
      LT=LOG10(T)
      LR=LOG10(DENSITY)

      IF ((LT.GT.5.).OR.(LR.GT.5.)) THEN
       REAL_A=A
       REAL_Z=Z
       NE=REAL_Z/REAL_A*6.022E23*DENSITY
       RETURN
      END IF

      IF (LT.LE.0.01) THEN
       IT=0
      ELSE IF (LT.GE.4.99) THEN
       IT=ITEMP-1
      ELSE
       IT=10*LT
      END IF

      IF (LR.LE.-4.99) THEN 
       JR=0
      ELSE IF (LR.GE.4.99) THEN
       JR=JRHO-1
      ELSE
       JR=5*(5+LR)
      END IF

      DELTEMP=LOG10(TEMP(IT+1)/TEMP(IT))
      WI1=LOG10(TEMP(IT+1)/T)/DELTEMP
      WI2=1.-WI1
      DELRHO=LOG10(RHO(JR+1)/RHO(JR))
      WJ1=LOG10(RHO(JR+1)/DENSITY)/DELRHO
      WJ2=1.-WJ1

      LNE=WI1*WJ1*LOG10(NELECT(IT,JR))+
     1    WI1*WJ2*LOG10(NELECT(IT,JR+1))+
     2    WI2*WJ1*LOG10(NELECT(IT+1,JR))+
     3    WI2*WJ2*LOG10(NELECT(IT+1,JR+1))

      NE=10.**LNE

      RETURN

      END



c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine density_old(t,p,a,z,dens)

c ***** checked sept 2 1991

      implicit real*8(a-h,k-z)
c      INCLUDE 'files_phys.inc.f'
      INCLUDE 'files.inc.f'

      dimension temp(0:60),pres(0:180),rho(0:60,0:180)
      save read,temp,pres,rho

      if (read.eq.1.) goto 100

c **** read precalculated densities ******************************

      open(unit=33,file=f_density,status='old')
       read(33,*)itmax,ipmax
       do i=1,4
        read(33,*)
       end do
       do it=itmax,0,-1
        read(33,*)
        do ip=ipmax,0,-1
         read(33,*)temp(it),pres(ip),rho(it,ip)
        end do
       end do
      close(unit=33,status='keep')
  
      read=1.

c ****************************************************************

100   continue

c      print *,'t,p=',t,p

      lt=log10(t)
      lp=log10(p)
      
      if(lt.le.log10(temp(0)))then
       it=0
      else if(lt.ge.log10(temp(itmax-1)))then
       it=itmax-1
      else
       it=10*(lt-log10(temp(0)))
      end if

      if(lp.lt.log10(pres(0)))then
       ip=0
      else if(lp.gt.log10(pres(ipmax-1)))then
       ip=ipmax-1
      else
       ip=10*(lp-log10(pres(0)))
      end if

      wt1=(log10(temp(it+1))-lt)/(log10(temp(it+1))-log10(temp(it)))
      wt2=1.-wt1
      wp1=(log10(pres(ip+1))-lp)/(log10(pres(ip+1))-log10(pres(ip)))
      wp2=1.-wp1

      ld=wt1*wp1*log10(rho(it,ip))+
     1   wt1*wp2*log10(rho(it,ip+1))+
     2   wt2*wp1*log10(rho(it+1,ip))+
     3   wt2*wp2*log10(rho(it+1,ip+1))

      dens=10.**ld

      end

