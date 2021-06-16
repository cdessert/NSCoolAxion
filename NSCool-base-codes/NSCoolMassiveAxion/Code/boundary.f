c
c      File cleaned up on Oct. 10, 2007
c

       function fteff(Tb,ifteff,eta,bfield,istep,time,Ts1,Ts2,
     1          Z,A,Rho,logdeltaM,MNS,debug)
c WARNING : Tb must be the local value, not the redshifted one !
c           fteff is also non red-shifted !
c           Ts1 and Ts2 ARE red-shifted !
       implicit real*8 (a-h,k-z)
        if (debug.ge.1.) print *,'Entering fteff'
        if (ifteff.eq.0) then
         fteff=fteff_table(Tb,Ts1,Ts2)
        else if (ifteff.eq.1) then
         fteff=fteff_GPE(Tb)
        else if (ifteff.eq.2) then
         fteff=fteff_NT(Tb)
        else if (ifteff.eq.3) then
         fteff=fteff_acc(Tb,eta)
        else if (ifteff.eq.4) then
         fteff=fteff_field_iron(Tb,bfield)
        else if (ifteff.eq.5) then
         fteff=fteff_ZARho(Tb,Z,A,Rho)
        else if (ifteff.eq.6) then
         fteff=fteff_Beznogov(Tb)
        else if (ifteff.eq.7) then
         fteff=fteff_PCY(Tb,logdeltaM,MNS)
        end if
        if (debug.ge.1.) print *,'Done'
      return
       end
c *********************************************************************
c *********************************************************************

      function fteff_PCY(Tb,logdeltaM,MNS)

c This is the boundary condition for Fe surface with H/He/C
c From https://arxiv.org/pdf/astro-ph/9706148.pdf
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using PCY boundary condition'
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****

      etaPCY=gs14**2.0d0*10.d0**logdeltaM/MNS

      Tb9=Tb/1.d9
      Ts=sqrt(7.0d0*Tb9*sqrt(gs14))
      z=Tb9-Ts/1.d3
      t4_iron=gs14*((7.0d0*z)**2.25+(z/3.0d0)**1.25)
      t4_wacc=gs14*(18.1d0*Tb9)**2.42
      if (etaPCY.gt.1.d-30) then
       a=(1.2d0+(5.3d-6/etaPCY)**0.38)*Tb9**(5./3.)
       t4_acc=(a*t4_iron+t4_wacc)/(a+1.0d0)
      else
       t4_acc=t4_iron
      end if

      fteff_PCY=t4_acc**0.25*1.d6
     
      end
      
c *********************************************************************
c *********************************************************************

      function fteff_Beznogov(Tb)

c This is the Beznogov, Potekhin & Yakovlev 1604.00538
c boundary condition for a H-He, He-C, or C-Fe combined envelope.
c They give Tb(Ts) analytically, so we have generated curves
c and inverted them. This code reads them from the "boundary" folder
c and reads off the numerical relation. 
c Todo: update to interpolation instead of nearest neighbor.
c Note: The C-Fe curves have a typo in the paper, so it is not
c possible to use at this time.

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      integer indlog10tb_lo,indlog10tb_hi
      real indlog10tb
      integer,parameter :: t_dimt = 10000
      real*8,dimension(1:t_dimt) :: t_log10tb
      real*8,dimension(1:t_dimt) :: t_log10ts
      real*8 :: t_log10tb0,t_log10tb1
c      common/gravity/gs14,compactness
      common/tbts_tables/t_log10tb,t_log10ts,t_log10tb0,t_log10tb1
      save print

      if (print.eq.1.) goto 10
       print *,'-----------------------------------------------------'
       print *,'Using BPY boundary condition'
       print *,'-----------------------------------------------------'
       print=1.
 10   continue

      log10tb = log10(tb)
      indlog10tba = int((log10tb-t_log10tb0)/
     1                 (t_log10tb1-t_log10tb0)*t_dimt + 1.5)
      log10tba = t_log10tb(indlog10tba)
      log10tbb = t_log10tb(indlog10tba + 1)
      log10tsa = t_log10ts(indlog10tba)
      log10tsb = t_log10ts(indlog10tba + 1)
      slope = (log10tsb-log10tsa)/(log10tbb-log10tba)
      log10ts = log10tsa + (log10tb-log10tba)*slope
      
c      write(6,*)indlog10tb,t_log10ts(indlog10tb)
      fteff_Beznogov = 10**log10ts

      end

c *********************************************************************
c *********************************************************************

      function fteff_ZARho(Tb,Z,A,Rho)

c This is the Hernquist & Applegate: ApJ 287, 244 (1984), analytical
c boundary condition for an arbitrary density (their Equ (3.12).

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print

      if (print.eq.1.) goto 10
       print *,'-----------------------------------------------------'
       print *,'Using Hernquist & Applegate boundary condition at  '
       print *,'  Rho =',Rho
       print *,'-----------------------------------------------------'
       print=1.
 10   continue

      fteff_ZARho = gs14**0.25 * 
     1              (A**3/Z**4 / 4.25d30 *Tb**6.5/Rho**2)**0.25

      end

c**********************************************************************
c**********************************************************************
       function fteff_table(Tb,Ts1,Ts2)
c checked on nov. 6 1990 and July 2006
c WARNING : Tb must be the local value, not the redshifted one !
c
c Ts1 & Ts2 are two auxiliary temperatures from the envelope calculations,
c They could be, e.g., the Max and min surface T, or anything wanted !
c Ts1 and Ts2 are NOT used for the cooling, they are only informative.
c
       implicit real*8 (a-h,k-z)
       character*100 f_TbTs
       common/bound/Tb_acc0,f_TbTs
       parameter(jtpmax=100)
       dimension tempb(0:jtpmax),temps(0:jtpmax),t2(0:jtpmax),
     1           temps1(0:jtpmax),t21(0:jtpmax),
     2           temps2(0:jtpmax),t22(0:jtpmax)
       common/gravity/gs14,compactness
       save read,tempb,temps,t2
        if (read.eq.1.) goto 999
c ************ read temps and tempb ************************************
         print *,'-------------------------------------------------'
         print *,'Using envelope boundary condition from table'
         print *,'     ',f_TbTs
         print *,'WARNING:   No gs14 correction applied ! '
         print *,'-------------------------------------------------'
        open(unit=15,file=f_TbTs,status='old')
         read(15,*)jmax
         jmax=jmax+1
         if (jmax.gt.jtpmax) then
          print '(a20,i5)',' jtpmax =',jtpmax
          pause 'function_fteff_table: jmax > jtpmax '
         end if
         do j=1,3
          read(15,*)
         end do
         do j=1,jmax
c---------------------
c This assumes Ts1 & Ts2 are defined in the table:
c          read(15,*)temps(j),tempb(j),temps1(j),temps2(j)
c---------------------
c If not, then use this:
          read(15,*)temps(j),tempb(j)
          temps1(j)=tempb(j)
          temps2(j)=tempb(j)
c---------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c HERE DANY: this version assumes the gs14^(1/4) effect is already 
c            included in the table !
c          temps(j)=temps(j)*gs14**0.25
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do
        close(unit=15,status='keep')
        do j=jmax+1,jtpmax
         tempb(j)=tempb(jmax)+
     1            float(j-jmax)*(tempb(jmax)-tempb(jmax-1))
         temps(j)=temps(jmax)+
     1            float(j-jmax)*(temps(jmax)-temps(jmax-1))
         temps1(j)=temps1(jmax)+
     1            float(j-jmax)*(temps1(jmax)-temps1(jmax-1))
         temps2(j)=temps2(jmax)+
     1            float(j-jmax)*(temps2(jmax)-temps2(jmax-1))
        end do
        dt0=(temps(2)-temps(1))/(tempb(2)-tempb(1))
        dt1=(temps(jtpmax)-temps(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps,jtpmax,dt0,dt1,t2)
        dt0=(temps1(2)-temps1(1))/(tempb(2)-tempb(1))
        dt1=(temps1(jtpmax)-temps1(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps1,jtpmax,dt0,dt1,t21)
        dt0=(temps2(2)-temps2(1))/(tempb(2)-tempb(1))
        dt1=(temps2(jtpmax)-temps2(jtpmax-1))/
     1      (tempb(jtpmax)-tempb(jtpmax-1))
        call spline1(tempb,temps2,jtpmax,dt0,dt1,t22)
        read=1.
c ***************************************************************************
999     continue

        call splint1(tempb,temps ,t2 ,jtpmax,Tb,Ts )
        call splint1(tempb,temps1,t21,jtpmax,Tb,Ts1)
        call splint1(tempb,temps2,t22,jtpmax,Tb,Ts2)
       
        fteff_table=Ts

       return
      end
c *********************************************************************
c *********************************************************************

      function fteff_NT(Tb)

c this is the Nomoto & Tsuruta boundary condition at rho=1.e10 g/cm3

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      dimension ltb(11),lts(11)

      data lts/5.698,5.764,5.943,6.103,6.228,
     1         6.363,6.493,6.633,6.773,6.918,7.068/
      data ltb/7.67,7.75,8.00,8.25,8.50,
     2         8.75,9.00,9.25,9.50,9.75,10.00/

      save lts,ltb
      save print
c*********************
      if (print.eq.1.) goto 10
       print *,'-------------------------------------------------'
       print *,'Using Nomoto & Tsuruta boundary condition'
       print *,'-------------------------------------------------'
       print=1.
 10   continue
c*********************

      logtb=dlog10(Tb)

      if (logtb.le.ltb(1)) then
       i1=1
      else if (logtb.ge.ltb(11)) then
       i1=10
      else
       do i=1,10
        if ((logtb.ge.ltb(i)).and.(logtb.le.ltb(i+1))) then
         i1=i
         goto 100
        end if
       end do
100    continue
      end if

      i2=i1+1
      deltb=ltb(i2)-ltb(i1)
      w1=(ltb(i2)-logtb)/deltb
      w2=1.-w1
      logts=lts(i1)*w1+lts(i2)*w2
      logts=logts+0.25*dlog10(gs14)

      fteff_NT=10.**logts      
 
      end

c *********************************************************************
c *********************************************************************

      function fteff_GPE(Tb)

c this is the Gudmundsson, Pethick & Epstein, ApJ 259, L19 (1982)
c boundary condition at rho=1.e10 g/cm3

c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print

      if (print.eq.1.) goto 10
       print *,'-------------------------------------------'
       print *,'Using Gundmundsson et al boundary condition'
       print *,'-------------------------------------------'
       print=1.
 10   continue

      Tb8=Tb/1.e8
      fteff_GPE=0.87e6*gs14**.25*(Tb8)**0.55
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_acc(Tb,eta)

c this is the boundary condition for accreted envelope
c From Potekhin, Chabrier & Yakovlev, A&A 323, 415 (1999)
c CHECKED on MARCH 13, 2001
c WARNING : Tb must be the local value, not the redshifted one !
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using accreted envelope boundary condition'
       print '(a15,1p1e10.2)','with   eta =',eta
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****

      Tb9=Tb/1.d9
      Ts=sqrt(7.0d0*Tb9*sqrt(gs14))
      z=Tb9-Ts/1.d3
      t4_iron=gs14*((7.0d0*z)**2.25+(z/3.0d0)**1.25)
      t4_wacc=gs14*(18.1d0*Tb9)**2.42
      if (eta.gt.1.d-30) then
       a=(1.2d0+(5.3d-6/eta)**0.38)*Tb9**(5./3.)
       t4_acc=(a*t4_iron+t4_wacc)/(a+1.0d0)
      else
       t4_acc=t4_iron
      end if

      fteff_acc=t4_acc**0.25*1.d6
     
      end

c *********************************************************************
c *********************************************************************

      function fteff_field_iron(Tb,bfield)

c this is the boundary condition for iron envelope with B field
c From Potekhin & Yakovlev, A&A 374 (2001), p. 213
c CHECKED on June 5, 2002 
c WARNING : Tb must be the local value, not the redshifted one !
c WARNING: bfield is B at the magnetic pole
c WARNING: good only for Tb > 1e7 K
     
      implicit real*8 (a-h,k-z)
      common/gravity/gs14,compactness
      save print
c*****
      if (print.eq.1.) goto 10
       print *,'------------------------------------------'
       print *,'Using iron envelope with magnetic field:'
       print *,'(a15,1p1e10.2)','B_pole =',bfield
       print *,'------------------------------------------'
       print=1.
 10   continue
c*****
      eta=0.d0
      f_zero=fteff_acc(Tb,eta)

      B12=bfield/1.d12
      T9=Tb/1.d9
      beta=0.074d0*dsqrt(B12)/T9**0.45
      a1=5059.d0*T9**0.75/
     1   (1.d0+20.4d0*T9**0.5+138.d0*T9**1.5+1102.d0*T9**2)**0.5
      a2=1484.d0*T9**0.75/
     1   (1.d0+90.d0*T9**1.5+125.d0*T9**2)**0.5
      a3=5530.d0/dsqrt(1.d0-0.4d0*compactness)*T9**0.75/
     1   (1.d0+8.16d0*T9**0.5+107.8d0*T9**1.5+560.d0*T9**2)**0.5

      ratio=(1.d0+a1*beta**2+a2*beta**3+0.007d0*a3*beta**4)/
     1      (1.d0+a3*beta**2)

      fteff_field_iron=f_zero*ratio**0.25
     
      end

C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINT1(XA,YA,Y2A,IN,X,Y)

C************************************************************************
C Given the arrays XA and YA of length IN, which tabulate a function    *
C (with the XA(i)'s in order), and given the array Y2A, which is the    *
C output from SPLINE, calculate the cubic-spline interpolated value     *
C Y for a given value X.                                                *
C                                                                       *
C From NUMERICAL RECIPES, p.89.                                         *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      DIMENSION XA(IN),YA(IN),Y2A(IN)

      KLO=1
      KHI=IN
1     IF (KHI-KLO.GT.1) THEN
       K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
         KHI=K
        ELSE
         KLO=K
        END IF
       GOTO 1
      END IF
      
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO) + B*YA(KHI)
     1  + ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      RETURN

      END

C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINE1(X,Y,IN,YP1,YPN,Y2)

C************************************************************************
C Given the arrays X and Y of length IN, which tabulate a function      *
C (with the XA(i)'s in order), and given values YP1 and YPN of the      *
C first derivative of the interpolating function at points 1 and N,     *
C this subroutine returns an array Y2 of length IN which contains the   *
C second derivatives of the interpolating function at the tabulated     *
C points X(i)'s.                                                        *
C                                                                       *
C From NUMERICAL RECIPES, p.88                                          *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      PARAMETER (JMAX=100)
      DIMENSION X(IN),Y(IN),Y2(IN),U(JMAX)

      IF (YP1.GT.1.E30) THEN
       Y2(1)=0.
       U(1)=0.
      ELSE
       Y2(1)=-0.5
       U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,IN-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.
       Y2(I)=(SIG-1.)/P
       U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO

      IF (YPN.GT.1.E30) THEN
       QN=0.
       UN=0.
      ELSE
       QN=0.5
       UN=(3./(X(IN)-X(IN-1)))*(YPN-(Y(IN)-Y(IN-1))/(X(IN)-X(IN-1)))
      END IF
      Y2(IN)=(UN-QN*U(IN-1))/(QN*Y2(IN-1)+1.)
      
      DO K=IN-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
      END DO

      RETURN
      END   

