C************************************************************************
C************************************************************************
      function istrlen(st)
      integer        i,istrlen
      character        st*(*)
      i = len(st)
      do while (st(i:i) .eq. ' ')
        i = i - 1
      enddo
      istrlen = i
      return
      end
C************************************************************************
C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINT2(X1A,X2A,YA,Y2A,IM,IN,X1,X2,Y)

C************************************************************************
C Given the arrays X1A of length IM, X2A of length IN, and YA (IMxIN)   *
C which tabulate a two variables function, and given Y2A (IMxIN), which *
C is the output from SPLINE2, calculates the bicubic-spline interpolated*
C value Y for given values X1 and X2.                                   *
C                                                                       *
C From NUMERICAL RECIPES, p.101.                                        *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      PARAMETER (IP=100)
      DIMENSION X1A(IM),X2A(IN),YA(IM,IN),Y2A(IM,IN)
      DIMENSION YTMP(IP),Y2TMP(IP),YYTMP(IP)

      DO IJ=1,IM
       DO IK=1,IN
        YTMP(IK)=YA(IJ,IK)
        Y2TMP(IK)=Y2A(IJ,IK)
       END DO
       CALL SPLINT(X2A,YTMP,Y2TMP,IN,X2,YYTMP(IJ))
      END DO
      CALL SPLINE(X1A,YYTMP,IM,1.D30,1.D30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,IM,X1,Y)

      RETURN

      END

C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINE2(X1A,X2A,YA,IM,IN,Y2A)

C************************************************************************
C Given the arrays X1A of length IM, X2A of length IN, and YA (IMxIN)   *
C which tabulate a two variables function, this subroutine constructs   *
C one-dimensional natural cubic splines of the rows of YA and returns   *
C the second-derivative in the array Y2A.                               *
C                                                                       *
C From NUMERICAL RECIPES, p.100.                                        *
C************************************************************************

      IMPLICIT REAL*8 (A-H,L-Z)
      PARAMETER (IP=100)
      DIMENSION X1A(IM),X2A(IM),YA(IM,IN),Y2A(IM,IN)
      DIMENSION YTMP(IP),Y2TMP(IP)

      DO IJ=1,IM
       DO IK=1,IN
        YTMP(IK)=YA(IJ,IK)
       END DO
       CALL SPLINE (X2A,YTMP,IN,1.D30,1.D30,Y2TMP)
       DO IK=1,IN
        Y2A(IJ,IK)=Y2TMP(IK)
       END DO
      END DO  

      RETURN
      END   

C************************************************************************
C************************************************************************
C************************************************************************

      SUBROUTINE SPLINT(XA,YA,Y2A,IN,X,Y)

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

      SUBROUTINE SPLINE(X,Y,IN,YP1,YPN,Y2)

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

      IF (YP1.GE.1.D30) THEN
       Y2(1)=0.d0
       U(1)=0.d0
      ELSE
       Y2(1)=-0.5d0
       U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,IN-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.
       Y2(I)=(SIG-1.)/P
       U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     1      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      END DO

      IF (YPN.GE.1.D30) THEN
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
C************************************************************************
C************************************************************************
C************************************************************************
