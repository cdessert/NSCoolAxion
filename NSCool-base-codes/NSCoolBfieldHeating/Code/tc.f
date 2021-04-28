cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                   neutron 1s0 superfluidity                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tcn1_sfb(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  SFB  nucl-th/0207004                                                *
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on July 23, 2002

      implicit real*8 (a-h,k-z)
      parameter (imax=18)
      dimension k0(imax),d0(imax),d2(imax)
      save k0,d0,d2

      data k0/0.000d0,
     1        0.100d0,0.200d0,0.300d0,0.400d0,0.500d0,
     1        0.600d0,0.700d0,0.800d0,0.900d0,1.000d0,
     2        1.100d0,1.175d0,1.250d0,1.300d0,1.350d0,
     3        1.400d0,1.450d0/
      data d0/0.000d0,
     1        0.000d0,0.090d0,0.210d0,0.350d0,0.490d0,
     1        0.610d0,0.720d0,0.790d0,0.780d0,0.700d0,
     2        0.570d0,0.440d0,0.280d0,0.190d0,0.100d0,
     3        0.030d0,0.000d0/
c      data d0/0.000d0,0.090d0,0.210d0,0.360d0,0.500d0,
c     1        0.610d0,0.720d0,0.790d0,0.780d0,0.700d0,
c     2        0.580d0,0.450d0,0.280d0,0.190d0,0.100d0,
c     3        0.030d0,0.000d0/
c      data d2/+2.915d1,-4.297d0,+6.040d0,-1.863d0,-4.590d0,
c     1        +2.221d0,-4.296d0,-9.037d0,-7.555d0,-2.741d0,
c     2        -5.480d0,-1.344d1,+1.656d1,-6.667d0,+1.010d1,
c     3        +1.426d1,+2.887d1/
c ************************************
c Generate d2 if first access:
       if (done.ne.1.1111) then 
        call spline_here(k0,d0,imax,0.d0,0.d0,d2)
        done=1.1111
       end if
c ************************************
      if (k.le.k0(1)) then
       tcn1_sfb=0.d0
       return
      else if (k.ge.k0(imax)) then
       tcn1_sfb=0.d0
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*d0(i1)+b*d0(i2)+
     1   ((a**3-a)*d2(i1)+(b**3-b)*d2(i2))*(delk**2)/6.

       tcn1_sfb=t/1.76d0*1.1604d10
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_ccdk(k)

c***********************************************************************
c  calculate the neutron 1S0 superfluidity critical temperature from   *
c  Chen, Clark, Dave & Khodel, Nucl. Phys. A555 (1993), p. 59
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)
      parameter (imax=11)
      dimension k0(imax),d0(imax),t0(imax),t2(imax)
      save k0,t0,t2,done

      data k0/0.10d0,0.20d0,0.30d0,0.40d0,0.50d0,
     1        0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,
     2        1.10d0/

      data d0/0.00d0,0.02d0,0.14d0,0.36d0,0.60d0,
     1        0.83d0,0.86d0,0.67d0,0.35d0,0.07d0,
     2        0.00d0/
c **** Generate t2 if first access: ***************************
      if (done.eq.1.0) goto 100
      do i=1,imax
       t0(i)=d0(i)/1.76d0*1.1604d10
      end do
      call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      done=1.0
c ***************************************************************
 100  continue
      if (k.le.k0(1)) then
       tcn1=0.
      else if (k.ge.k0(imax)) then
       tcn1=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tcn1=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if
      tcn1_ccdk=tcn1

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_wap(k)

c***********************************************************************
c  calculate the neutron 1S0 superfluidity critical temperature from   *
c  Wambach, Ainsworth & Pines Nucl. Phys. A555 (1993), p. 128
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),d0(imax),t0(imax),t2(imax)
      save k0,t0,t2,done

      data k0/0.10d0,0.20d0,0.30d0,0.40d0,0.50d0,
     1        0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,
     2        1.10d0,1.20d0,1.30d0,1.40d0/

      data d0/0.00d0,0.03d0,0.13d0,0.30d0,0.54d0,
     1        0.74d0,0.86d0,0.90d0,0.79d0,0.59d0,
     2        0.35d0,0.14d0,0.03d0,0.00d0/
c **** Generate t2 if first access: ***************************
      if (done.eq.1.0) goto 100
      do i=1,imax
       t0(i)=d0(i)/1.76d0*1.1604d10
      end do
      call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      done=1.0
c ***************************************************************
 100  continue
      if (k.le.k0(1)) then
       tcn1=0.
      else if (k.ge.k0(imax)) then
       tcn1=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tcn1=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if
      tcn1_wap=tcn1

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_GC(k)

c***********************************************************************
c  This attempts to mock-up results of Gezerlis & Carlson, by using the
c  Takatsuka (1972) values scaled down !
c
c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  t72, prog.theor.phys. 48(1972): p.1517 with values as reported in   *
c       table i of t84, prog.theor.phys 71(1984): p.1432.              *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.14,0.20,0.24,0.33,0.45,0.52,
     1        0.61,0.81,0.97,1.09,1.19,1.30,1.50/
      data t0/0.00e00,1.12e09,2.37e09,3.69e09,6.53e09,
     1        1.07e10,1.29e10,1.52e10,1.71e10,1.50e10,
     2        1.15e10,7.51e09,3.30e09,0.00e00/
      data t2/+1.24e+11,+9.50e+10,+3.61e+11,-1.20e+11,+9.13e+10,
     1        -6.95e+10,-6.40e+10,-1.10e+11,-1.34e+11,-1.08e+11,
     2        -1.20e+11,+1.15e+10,+1.54e+11,+1.71e+11/

      if (k.le.k0(1)) then
       tcn1_GC=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_GC=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_GC= 0.7 * t
      end if

      return
      end

c***********************************************************************
c***********************************************************************
      function tcn1_gipsf(k)

c***********************************************************************
c  Calculates the neutron 1S0 superfluidity critical temperature from   *
c  Gondolfi, Illarionov, Pederiva, Schmidt & Fantoni,                  *
c  PRL 111 (2008), 132501                                              *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on February 5, 2009

      implicit real*8 (a-h,k-z)
      parameter (imax=9)
      dimension k0(imax),d0(imax),d2(imax)
      save k0,d0,d2

      data k0/0.000d0,
     1        0.200d0,0.300d0,0.400d0,0.600d0,0.700d0,
     1        0.800d0,1.000d0,1.200d0/
      data d0/0.000d0,
     1        0.300d0,0.900d0,1.500d0,2.100d0,1.900d0,
     1        1.500d0,0.500d0,0.000d0/
c ************************************
c Generate d2 if first access:
       if (done.ne.1.1111) then 
        call spline_here(k0,d0,imax,0.d0,0.d0,d2)
        done=1.1111
       end if
c ************************************
      if (k.le.k0(1)) then
       tcn1_gipsf=0.d0
       return
      else if (k.ge.k0(imax)) then
       tcn1_gipsf=0.d0
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*d0(i1)+b*d0(i2)+
     1   ((a**3-a)*d2(i1)+(b**3-b)*d2(i2))*(delk**2)/6.

       tcn1_gipsf=t/1.76d0*1.1604d10
      end if

      return
      end

c***********************************************************************
c***********************************************************************
c
c  Following functions for neutron 1S0 are old and should NOT be used !
c
c***********************************************************************
c***********************************************************************

      function tcn1_t72(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  t72, prog.theor.phys. 48(1972): p.1517 with values as reported in   *
c       table i of t84, prog.theor.phys 71(1984): p.1432.              *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.14,0.20,0.24,0.33,0.45,0.52,
     1        0.61,0.81,0.97,1.09,1.19,1.30,1.50/
      data t0/0.00e00,1.12e09,2.37e09,3.69e09,6.53e09,
     1        1.07e10,1.29e10,1.52e10,1.71e10,1.50e10,
     2        1.15e10,7.51e09,3.30e09,0.00e00/
      data t2/+1.24e+11,+9.50e+10,+3.61e+11,-1.20e+11,+9.13e+10,
     1        -6.95e+10,-6.40e+10,-1.10e+11,-1.34e+11,-1.08e+11,
     2        -1.20e+11,+1.15e+10,+1.54e+11,+1.71e+11/

      if (k.le.k0(1)) then
       tcn1_t72=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_t72=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_t72=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_ns(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  ns, preprint 1981 as cited in                                       *
c  ao nucl. phys. a437(1986): p. 487, fig 14.                          *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 17 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=17)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.30,0.40,0.50,0.60,0.65,
     1        0.70,0.75,0.80,0.85,0.90,0.95,
     2        1.00,1.05,1.10,1.15,1.20,1.25/
      data t0/0.00e00,1.38e09,3.16e09,5.27e09,6.59e09,
     1        8.50e09,1.07e10,1.42e10,1.61e10,1.74e10,1.81e10,
     2        1.78e10,1.52e10,1.22e10,8.31e09,2.31e09,0.00e00/
      data t2/   +4.48e+11,-6.74e+10,+6.18e+10,+1.83e+10,+4.03e+11,
     1 -2.14e+11,+1.15e+12,-1.26e+12,+3.50e+10,-3.25e+11,-1.77e+11,
     2 -1.37e+12,+1.31e+11,-1.17e+11,-1.80e+12,+2.25e+12,+1.65e+12/

      if (k.le.k0(1)) then
       tcn1_ns=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_ns=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_ns=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_t84(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  t84, prog.theor.phys. 71(1984): p.1432                              *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 17, 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.14,0.20,0.24,0.33,0.45,0.52,
     1        0.61,0.81,0.97,1.09,1.19,1.30,1.50/
      data t0/0.00e00,1.12e09,2.44e09,3.76e09,6.72e09,
     1        1.11e10,1.32e10,1.57e10,1.75e10,1.54e10,
     2        1.18e10,7.25e09,2.97e09,0.00e00/
      data t2/+1.09e+11,+1.25e+11,+3.10e+11,-8.90e+10,+1.12e+11,
     1        -1.44e+11,+3.47e+10,-1.59e+11,-1.17e+11,-1.03e+11,
     2        -2.06e+11,+9.83e+10,+1.71e+11,+1.37e+11/

      if (k.le.k0(1)) then
       tcn1_t84=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_t84=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_t84=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_ao(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  ao, nucl.phys. a437(1985): p.487 from case v(rsc)eff of fig. 11.    *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=21)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.1,0.2,0.3,0.4,0.5,0.6,0.7,
     1        0.8,0.9,1.0,1.1,1.2,1.3,1.4,
     2        1.5,1.6,1.7,1.8,1.9,2.0,2.1/
      data t0/0.00e00,8.89e08,2.64e09,5.93e09,1.12e10,
     1        1.68e10,2.31e10,2.80e10,3.10e10,3.23e10,
     2        3.23e10,3.06e10,2.77e10,2.34e10,1.85e10,
     3        1.22e10,7.58e09,3.96e09,1.65e09,5.25e08,0.00e00/
      data t2/+3.09e+11,-2.43e+10,+1.86e+11,+2.66e+11,-5.97e+10,
     1        +1.71e+11,-2.06e+11,-1.87e+11,-1.84e+11,-9.61e+10,
     2        -2.12e+11,-7.73e+10,-1.99e+11,+3.43e+10,-2.98e+11,
     3   +3.18e+11,+3.58e+10,+1.39e+11,+1.94e+11,-8.32e+10,+2.59e+11/

      if (k.le.k0(1)) then
       tcn1_ao=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_ao=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_ao=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_ccks_var(k)

c************************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from    *
c  ccks nucl.phys. a451(1986): p.509. reid pot. with variational method *
c  uses a cubic spline interpolation.                                   *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c************************************************************************

c checked on feb 17, 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=10)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.5,0.6,0.7,0.8,0.9,
     1        1.0,1.1,1.2,1.3,1.4/
      data t0/1.32e10,1.36e10,1.39e10,1.40e10,1.30e10,
     1        1.02e10,6.06e09,1.98e09,1.32e08,0.00e00/
      data t2/+1.48e+11,-5.52e+10,+1.30e+10,-1.17e+11,-2.05e+11,
     1        -1.43e+11,-2.83e+10,+2.92e+11,+1.99e+11,-6.01e+10/

      if (k.le.0.3)then 
       tcn1_ccks_var=0.
       return
      else if (k.le.k0(1)) then
       w=(k-0.3)/0.2
       tcn1_ccks_var=w*t0(1)
       return
      else if (k.ge.k0(imax)) then
       tcn1_ccks_var=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_ccks_var=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_ccks_cbf(k)

c***********************************************************************
c  calculates the neutron 1s0 superfluidity critical temperature from  *
c  ccks nucl.phys. a451(1986): p.509. reid pot. with cbf method        *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 17, 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=9)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.3,0.4,0.5,0.6,0.7,
     1        0.8,0.9,1.0,1.1/
      data t0/0.00e0,1.25e9,3.43e9,3.49e9,3.03e9,
     1        2.11e9,9.23e8,6.59e7,0.00e0/
      data t2/+2.97e+11,+1.57e+11,-3.65e+11,+3.05e+10,-6.93e+10,
     1        -2.93e+10,+2.62e+10,+1.23e+11,-4.15e+10/

      if (k.le.k0(1)) then
       tcn1_ccks_cbf=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_ccks_cbf=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_ccks_cbf=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_awp_2(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  awp  phys.lett. 222(1989): p.173. case ii, from fig. 3.             *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on feb 17, 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=17)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.1,0.2,0.3,0.4,0.5,0.6,
     1        0.7,0.8,0.9,1.0,1.1,1.2,
     2        1.3,1.4,1.5,1.6,1.7/
      data t0/0.00e0,3.30e8,1.18e9,2.44e9,4.20e9,6.13e9,
     1        7.91e9,9.10e9,9.56e9,9.03e9,7.71e9,5.93e9,
     2        4.15e9,2.50e9,1.12e9,3.61e8,0.00e0/
      data t2/   +7.58e+10,+4.64e+10,+5.06e+10,-2.87e+09,+1.17e+11,
     1 -7.45e+10,-5.28e+10,-6.83e+10,-1.12e+11,-7.75e+10,-5.19e+10,
     2 +9.16e+09,+1.53e+10,+7.71e+09,+1.16e+11,-3.87e+10,+1.58e+11/

      if (k.le.k0(1)) then
       tcn1_awp_2=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_awp_2=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_awp_2=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_awp_3(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  awp  phys.lett. 222(1989): p.173. case iii, from fig. 3.            *
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on feb 17 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=15)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.1,0.2,0.3,0.4,0.5,
     1        0.6,0.7,0.8,0.9,1.0,
     2        1.1,1.2,1.3,1.4,1.5/
      data t0/0.00e0,2.64e8,7.91e8,1.78e9,3.36e9,
     1        5.27e9,6.59e9,7.25e9,7.05e9,5.74e9,
     2        3.96e9,1.98e9,7.91e8,1.32e8,0.00e0/
      data t2/+7.54e+10,+7.60e+09,+5.20e+10,+6.16e+10,+5.62e+10,
     1        -8.83e+10,-5.71e+10,-7.94e+10,-1.41e+11,-2.12e+10,
     2        -5.59e+10,+1.25e+11,+3.07e+10,+7.03e+10,+4.47e+09/

      if (k.le.k0(1)) then
       tcn1_awp_3=0.
       return
      else if (k.ge.k0(imax)) then
       tcn1_awp_3=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn1_awp_3=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn1_bbllp(k)

c***********************************************************************
c  calculate the neutron 1s0 superfluidity critical temperature from   *
c  Broglia et al, Phys.Rev. D***, 1994, p. ****                        *
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)
      parameter (imax=36)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2
c **** Read the data if first access: ***************************
      if (read.eq.1.0) goto 100
      open(unit=35,file='/home/page/nstar/pierre/tc.dat',status='old')
       do i=1,5
        read(35,*)
       end do
       do i=2,imax-1
        read(35,*)k0(i),nothing,t0(i)
       end do
      close(unit=35,status='keep')
      k0(1)=0.5*k0(2)
      t0(1)=0.0
      k0(imax)=2.*k0(imax-1)
      t0(imax)=0.0
      call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      read=1.0
c ***************************************************************
 100  continue
      if (k.le.k0(1)) then
       tcn1=0.
      else if (k.ge.k0(imax)) then
       tcn1=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tcn1=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if
      tcn1_bbllp=tcn1

      pause 'Tc for 1S0 neutron pairing from bbllp badly defined !'

      return
      end
c************************************************************************
c************************************************************************

      function tcn1_sclbl96(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  Schulze, Cugnon, Lejeune, Baldo & Lombardo, preprint                *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sept 26, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      data k0/+0.00E+00,+2.50E-01,+4.00E-01,+6.00E-01,
     1        +7.00E-01,+8.20E-01,+9.00E-01,+1.00E+00,
     2        +1.10E+00,+1.20E+00,+1.30E+00,+1.50E+00/
     1        
      data t0/+0.00E+00,+4.17E+09,+9.06E+09,+1.59E+10,
     1        +1.81E+10,+1.90E+10,+1.86E+10,+1.66E+10,
     2        +1.34E+10,+9.19E+09,+4.50E+09,+0.00E+00/
      data t2/+1.70E+11,+6.08E+10,+3.14E+10,-1.13E+11,
     1        -1.21E+11,-1.40E+11,-1.73E+11,-1.08E+11,
     2        -1.10E+11,-8.58E+10,+1.76E+11,+2.49E+11/
      save k0,t0,t2

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcn1_sclbl96=tc

      return
      end
c************************************************************************
c************************************************************************

      function tcn1_sclbl96_pol(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  Schulze, Cugnon, Lejeune, Baldo & Lombardo, preprint                *
c  with medium polarization effects                                    *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sept 26, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      data k0/+3.00E-01,+5.00E-01,+6.00E-01,+7.00E-01,+9.00E-01,
     1        +1.10E+00,+1.20E+00,+1.25E+00,+1.30E+00,+1.40E+00,
     2        +1.50E+00,+1.55E+00,+1.60E+00,+1.70E+00/
      data t0/+0.00E+00,+3.31E+08,+5.95E+08,+1.46E+09,+3.57E+09,
     1        +5.82E+09,+6.75E+09,+6.81E+09,+6.75E+09,+5.82E+09,
     2        +3.77E+09,+1.92E+09,+7.94E+08,+0.00E+00/
      data t2/+3.31E+10,-1.65E+10,+9.25E+10,+3.57E+09,+2.56E+09,
     1        +6.03E+09,-1.60E+11,-2.24E+09,-1.48E+11,-3.07E+10,
     2        -4.04E+11,+4.99E+11,+1.54E+11,+1.61E+11/
      save k0,t0,t2

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcn1_sclbl96_pol=tc

      return
      end

c**************************************************************************

c**************************************************************************

c**************************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                   neutron 3p2 superfluidity                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tcn3_hgrr(k)

c***********************************************************************
c  calculates the neutron 3p2 superfluidity critical temperature from  *
c  hgrr, phys.rev.lett. 24(1970):p.775.    m*=1                        *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on jan.24 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=21)     
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/1.0,1.2,1.4,1.6,1.8,1.9,2.0,
     1        2.1,2.2,2.3,2.4,2.5,2.6,2.7,
     2        2.8,2.9,3.0,3.2,3.4,3.6,3.8/
      data t0/0.00e0,1.75e8,8.12e8,2.20e9,3.83e9,4.58e9,5.22e9,
     1        5.74e9,6.09e9,6.18e9,6.06e9,5.80e9,5.41e9,4.97e9,
     2        4.41e9,3.83e9,3.19e9,1.97e9,8.12e8,1.48e8,0.00e0/
      data t2/+8.65e09,+8.95e09,+2.48e10,+4.34e09,-5.91e09,
     1        -1.22e10,-1.11e10,-1.52e10,-2.99e10,-2.12e10,
     2        -1.14e10,-1.71e10,+1.67e09,-1.96e10,+4.83e09,
     3        -1.17e10,+5.98e09,-3.09e09,+1.57e10,+1.45e10,+3.87e09/

      if (k.le.k0(1)) then
       tcn3_hgrr=0.
       return
      else if (k.ge.k0(imax)) then
       tcn3_hgrr=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn3_hgrr=t
      end if

      return
      end

c **********************************************************************
c **********************************************************************

      function tcn3_ao(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  a0, nucl.phys. a442(1985):p.163.                                    *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on jan 24 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=26)     
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
     1        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
     2        3.1,3.2,3.3,3.4,3.5,3.6/
      data t0/0.00e0,3.50e7,1.25e8,2.42e8,3.67e8,
     1        4.85e8,5.95e8,7.06e8,8.01e8,8.72e8,
     2        8.80e8,8.31e8,6.79e8,5.12e8,3.95e8,
     3        3.12e8,2.51e8,1.99e8,1.65e8,1.37e8,
     4        1.12e8,8.78e7,6.45e7,4.15e7,1.88e7,0.00e0/
      data t2/+7.63e+09,+5.74e+09,+2.41e+09,+8.37e+08,-9.52e+08,
     1        -1.23e+09,+1.07e+09,-2.46e+09,-8.48e+08,-8.55e+09,
     2        -2.74e+09,-1.47e+10,-3.01e+08,+6.89e+09,+2.72e+09,
     3        +2.61e+09,+3.23e+07,+2.66e+09,+1.29e+08,+4.24e+08,
     4        -2.62e+07,+1.60e+08,-7.58e+07,+3.23e+08,-1.04e+09,
     5        +6.16e+09/

      if (k.le.k0(1)) then
       tcn3_ao=0.
       return
      else if (k.ge.k0(imax)) then
       tcn3_ao=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn3_ao=t
      end if

      return
      end

c **********************************************************************
c **********************************************************************

      function tcn3_t72(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  t72, rep.prog.theor.phys. 48(1972):p.1517.                          *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb 19 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/1.38,1.42,1.46,1.53,1.64,1.73,
     1        1.82,1.89,1.96,1.98,2.02,2.06/
      data t0/0.00e0,1.00e7,1.00e8,3.02e8,6.02e8,7.76e8,
     1        6.02e8,3.49e8,1.52e8,1.00e8,1.00e7,0.00e0/
      data t2/-2.08e+10,+7.91e+10,+4.34e+09,-4.36e+09,
     1        +2.86e+09,-6.03e+10,-1.94e+10,+2.21e+10,
     2        -3.03e+08,-1.02e+10,+8.33e+10,-2.29e+10/

c***** Rescale k because values used to build this function 
c***** were apparently not accurate:
      k1_8=1.68      ! lower value of k where T_c reaches 10^8 K
c                    ! which must be mapped onto 1.46
      k2_8=2.38      ! upper value of k where T_c decreases back to 10^8 K
c                    ! which must be mapped onto 1.98
      slope=(1.98-1.46)/(k2_8-k1_8)
      kk=1.46+slope*(k-k1_8)

      if (kk.le.k0(1)) then
       tcn3_t72=0.
       return
      else if (kk.ge.k0(imax)) then
       tcn3_t72=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.kk) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-kk)/delk
       b=(kk-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn3_t72=t
      end if

      return
      end

c************************************************************************
c************************************************************************

      function tcn3_t72_m1(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  t72, rep.prog.theor.phys. 48(1972):p.1517. from fig.2   m*=1        *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb 20 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=25)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/1.20,1.39,1.56,1.70,1.84,
     1        1.97,2.09,2.20,2.31,2.41,
     2        2.51,2.60,2.69,2.78,2.87,
     3        2.95,3.03,3.11,3.19,3.26,
     4        3.34,3.41,3.48,3.55,3.61/
      data t0/0.00e0,1.66e8,4.85e8,1.05e9,1.73e9,
     1        2.27e9,2.70e9,3.02e9,3.20e9,3.20e9,
     2        3.10e9,2.96e9,2.79e9,2.60e9,2.38e9,
     3        2.15e9,1.91e9,1.65e9,1.38e9,1.14e9,
     4        8.70e8,6.23e8,3.74e8,1.38e8,0.00e0/
      data t2/+1.37e+10,+1.32e+08,+1.95e+10,+6.07e+09,-8.58e+09,
     1        -3.38e+09,-5.16e+09,-1.15e+10,-1.82e+10,-9.07e+09,
     2        -5.51e+09,-3.69e+09,-1.94e+09,-3.36e+09,-6.83e+09,
     3        +5.15e+08,-4.60e+09,-8.46e+08,-1.39e+09,+2.31e+09,
     4        -3.45e+09,-1.04e+09,+5.14e+09,-3.61e+09,+1.17e+11/

      if (k.le.k0(1)) then
       tcn3_t72_m1=0.
       return
      else if (k.ge.k0(imax)) then
       tcn3_t72_m1=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn3_t72_m1=t
      end if

      return
      end

c************************************************************************
c************************************************************************

      function tcn3_ao_m1(k)

c***********************************************************************
c  calculate the neutron 3p2 superfluidity critical temperature from   *
c  a0, nucl.phys. a442(1985):p.163. m*=1 gap from fig.2                *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb 19 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=23)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/1.35,1.4,1.5,1.6,1.7,1.8,
     1         1.9,2.0,2.1,2.2,2.3,2.4,
     2         2.5,2.6,2.7,2.8,2.9,3.0,
     3         3.1,3.2,3.3,3.4,3.5/
      data t0/1.00e0,6.92e7,2.77e8,5.26e8,8.45e8,1.14e9,
     1        1.51e9,1.87e9,2.16e9,2.34e9,2.44e9,2.47e9,
     2        2.40e9,2.23e9,1.95e9,1.56e9,1.16e9,8.10e8,
     3        5.26e8,3.05e8,1.52e8,5.54e7,1.00e0/

      data t2/+8.36e+10,-1.20e+09,+3.42e+09,+1.22e+10,-1.03e+10,
     1        +1.47e+10,-3.60e+09,-6.32e+09,-1.31e+10,-7.25e+09,
     2        -5.91e+09,-1.11e+10,-9.56e+09,-1.06e+10,-1.40e+10,
     3        +6.01e+08,+5.59e+09,+7.02e+09,+5.93e+09,+7.07e+09,
     4        +6.58e+09,+4.33e+08,+1.64e+10/

      if (k.le.k0(1)) then
       tcn3_ao_m1=0.
       return
      else if (k.ge.k0(imax)) then
       tcn3_ao_m1=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcn3_ao_m1=t
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn3_bcll92(k)

c***********************************************************************
c  Calculates the neutron 3p2 superfluidity critical temperature from  *
c  Baldo, Cugnon, Leujeune & Lombardo, N.Ph. A536 (1992), p. 349       *
c  Uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sept. 25, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      data k0/+1.50E00,+1.80E00,+2.00E00,+2.25E00,+2.50E00,+2.60E00,
     1        +2.70E00,+2.80E00,+2.90E00,+3.00E00,+3.50E00,+4.50E00/
      data t0/+0.00E00,+6.89E08,+1.72E09,+4.48E09,+7.64E09,+9.16E09,
     1        +1.03E10,+1.07E10,+1.03E10,+8.95E09,+4.13E09,+0.00E00/
      data t2/+2.17E10,+2.56E09,+4.08E10,-8.32E09,+3.22E10,-5.57E10,
     1        -5.72E10,-8.74E10,-1.30E11,+3.06E10,+2.14E09,+1.13E10/
      save k0,t0,t2

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcn3_bcll92=tc

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn3_eehjo96_nr(k)

c***********************************************************************
c  Calculates the neutron 3p2 superfluidity critical temperature from  *
c  Elgaroy, Engvik, Hjorth-Jensen & Osnes, preprint 1996               *
c  non-relativistic case                                               *
c  Uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sept. 25, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=11)
      dimension k0(imax),t0(imax),t2(imax)
      data k0/+1.06E+00,+1.33E+00,+1.53E+00,+1.68E+00,
     1        +1.81E+00,+1.92E+00,+2.02E+00,+2.12E+00,
     2        +2.20E+00,+2.28E+00,+2.35E+00/
      data t0/+0.00E+00,+3.30E+07,+1.43E+08,+2.59E+08,
     1        +3.09E+08,+3.23E+08,+3.27E+08,+3.14E+08,
     2        +2.51E+08,+1.19E+08,+0.00E+00/
      data t2/+5.69E+07,+2.50E+09,+1.77E+09,-3.92E+09,
     1        -2.10E+09,-3.77E+08,-1.03E+09,-7.04E+09,
     2        -1.22E+10,-1.32E+10,+7.25E+10/
      save k0,t0,t2

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcn3_eehjo96_nr=tc

      return
      end

c***********************************************************************
c***********************************************************************

      function tcn3_eehjo96_r(k)

c***********************************************************************
c  Calculates the neutron 3p2 superfluidity critical temperature from  *
c  Elgaroy, Engvik, Hjorth-Jensen & Osnes, preprint 1996               *
c  non-relativistic case                                               *
c  Uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sept. 25, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=4)
      dimension k0(imax),t0(imax),t2(imax)
      data k0/+1.06E+00,+1.33E+00,+1.53E+00,+1.68E+00/
      data t0/+0.00E+00,+8.30E+07,+7.00E+07,+0.00E+00/
      data t2/+4.62E+09,-2.66E+09,-5.17E+09,+1.15E+10/
      save k0,t0,t2

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcn3_eehjo96_r=tc

      return
      end

c***********************************************************************

c***********************************************************************

c***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                     proton 1s0 superfluidity                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tcp1_ccy_ms(k)

c***********************************************************************
c  calculate the proton 1s0 superfluidity critical temperature from    *
c  ccy, nucl.phys. a179(1972): p.320. from curve msii) fig. 4.         *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb. 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.10,0.20,0.30,0.40,0.50,0.60,
     1        0.65,0.70,0.80,0.90,1.00,1.05,1.10/
      data t0/0.00e0,2.97e8,1.12e9,2.24e9,3.66e9,4.78e9,5.28e9,
     1        5.28e9,5.14e9,4.48e9,3.23e9,1.32e9,1.32e8,0.00e0/
      data t2/+5.91e+10,+5.99e+10,+1.68e+10,+5.12e+10,-4.14e+10,
     1        -6.54e+10,-6.89e+10,-5.60e+10,-4.33e+10,-7.02e+10,
     2        -2.99e+10,-2.06e+11,+7.38e+11,-2.10e+11/

      if (k.le.k0(1)) then
       tcp1_ccy_ms=0.
       return
      else if (k.ge.k0(imax)) then
       tcp1_ccy_ms=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcp1_ccy_ms=t
      end if

      return
      end

c*************************************************************************
c*************************************************************************

      function tcp1_ccy_ps(k)

c***********************************************************************
c  calculate the proton 1s0 superfluidity critical temperature from    *
c  ccy, nucl.phys. a179(1972): p.320. from curve psi) fig.4.           *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb. 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.10,0.20,0.30,0.40,0.50,
     1        0.55,0.60,0.70,0.80,0.90,0.95/
      data t0/0.00e0,2.97e8,1.12e9,2.24e9,3.53e9,4.18e9,
     1        4.28e9,4.20e9,3.59e9,2.24e9,4.94e8,0.00e0/
      data t2/+5.94e+10,+5.94e+10,+1.84e+10,+4.50e+10,
     1        -9.63e+10,-4.37e+10,-8.52e+10,-4.76e+10,
     2        -8.46e+10,-5.79e+10,+7.87e+10,+5.53e+11/

      if (k.le.k0(1)) then
       tcp1_ccy_ps=0.
       return
      else if (k.ge.k0(imax)) then
       tcp1_ccy_ps=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcp1_ccy_ps=t
      end if

      return
      end
c***********************************************************************
c***********************************************************************

      function tcp1_ccdk(k)

c***********************************************************************
c  calculate the neutron 1S0 superfluidity critical temperature from   *
c  Chen, Clark, Dave & Khodel, Nucl. Phys. A555 (1993), p. 59
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),d0(imax),t0(imax),t2(imax)
      save k0,t0,t2,done

      data k0/0.00d0,0.10d0,0.20d0,0.30d0,0.40d0,0.50d0,
     1        0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,1.10d0,
     2        1.20d0,1.30d0/

      data d0/0.00d0,0.05d0,0.19d0,0.41d0,0.66d0,0.84d0,
     1        0.97d0,1.00d0,0.90d0,0.72d0,0.49d0,0.22d0,
     2        0.07d0,0.00d0/
c **** Generate t2 if first access: ***************************
      if (done.eq.1.0) goto 100
      do i=1,imax
       t0(i)=d0(i)/1.76d0*1.1604d10
      end do
      call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      done=1.0
c ***************************************************************
 100  continue
      if (k.le.k0(1)) then
       tcn1=0.
      else if (k.ge.k0(imax)) then
       tcn1=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tcn1=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if
      tcp1_ccdk=tcn1

      return
      end
c*************************************************************************
c*************************************************************************

      function tcp1_t73(k)

c***********************************************************************
c  calculate the proton 1s0 superfluidity critical temperature from    *
c  t73, prog.th.phys. 50(1973):p.1754                                  *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb. 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.30,0.41,0.48,0.53,0.58,
     1        0.63,0.67,0.72,0.76,0.80,0.84/
      data t0/0.00e0,7.58e8,2.27e9,3.03e9,3.29e9,3.00e9,
     1        2.55e9,2.11e9,1.50e9,8.45e8,3.29e7,0.00e0/
      data t2/-2.49e+10,+1.00e+11,-6.79e+10,-5.58e+10,
     1        -3.16e+11,-1.03e+09,-6.41e+10,-1.04e+10,
     2        -5.55e+10,-3.64e+11,+9.21e+11,-3.99e+11/

      if (k.le.k0(1)) then
       tcp1_t73=0.
       return
      else if (k.ge.k0(imax)) then
       tcp1_t73=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcp1_t73=t
      end if

      return
      end

c*************************************************************************
c*************************************************************************

      function tcp1_ns(k)

c***********************************************************************
c  calculate the proton 1s0 superfluidity critical temperature from    *
c  ns as cited by ao, nucl.phys. a437(1985): p.487 from fig. 15.       *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb 16, 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.10,0.20,0.25,0.30,0.35,0.40,
     1        0.45,0.50,0.60,0.70,0.80,0.90/
      data t0/0.00e0,1.09e9,2.54e9,4.30e9,5.52e9,6.20e9,
     1        5.40e9,3.76e9,1.45e9,5.60e8,1.32e8,0.00e0/
      data t2/+1.96e+11,+2.61e+11,+2.13e+11,-3.69e+11,
     1        -3.43e+10,-7.90e+11,-3.57e+11,+2.01e+11,
     2        +1.57e+11,+2.17e+10,+3.32e+10,+2.30e+10/

      if (k.le.k0(1)) then
       tcp1_ns=0.
       return
      else if (k.ge.k0(imax)) then
       tcp1_ns=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcp1_ns=t
      end if

      return
      end

c*************************************************************************
c*************************************************************************

      function tcp1_ao(k)

c***********************************************************************
c  calculate the proton 1s0 superfluidity critical temperature from    *
c  ao, nucl.phys. a437(1985): p.487 from dotted curve fig. 15.         *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in k               *
c***********************************************************************

c checked on feb. 16 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=12)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.00,0.20,0.30,0.40,0.50,0.60,
     1        0.70,0.80,0.90,1.00,1.10,1.20/
      data t0/0.00e0,7.91e8,1.51e9,2.14e9,2.34e9,2.04e9,
     1        1.42e9,8.57e8,4.61e8,1.65e8,6.59e7,0.00e0/
      data t2/+5.13e+10,+1.61e+10,-4.86e+09,-5.00e+10,
     1        -5.30e+10,-3.79e+10,+1.27e+10,+2.14e+10,
     2        +1.82e+09,+3.13e+10,-8.90e+09,+2.42e+10/

      if (k.le.k0(1)) then
       tcp1_ao=0.
       return
      else if (k.ge.k0(imax)) then
       tcp1_ao=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       t=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.

       tcp1_ao=t
      end if

      return
      end

c************************************************************************
c************************************************************************

      function tcp1_bcll92(k)

c***********************************************************************
c  calculate the proton 1S0 superfluidity critical temperature from    *
c  Baldo, Cugnon, Lejeune & Lombardo, NP A536 (1992), p. 349           *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and Tc is given in K               *
c***********************************************************************

c checked on Sep 26, 1996

      implicit real*8 (a-h,k-z)
      parameter (imax=10)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2
      data k0/+5.00E-02,+2.00E-01,+3.10E-01,+4.00E-01,+4.60E-01,
     1        +5.70E-01,+6.30E-01,+7.60E-01,+8.80E-01,+1.10E+00/
      data t0/+0.00E+00,+9.92E+08,+3.77E+09,+5.36E+09,+5.62E+09,
     1        +5.36E+09,+5.03E+09,+3.77E+09,+1.32E+09,+0.00E+00/
      data t2/+1.25E+10,+2.40E+11,-1.33E+11,-2.10E+11,-7.34E+10,
     1        -3.03E+10,-4.51E+09,-1.65E+11,+1.54E+11,+4.80E+09/

      if (k.le.k0(1)) then
       tc=0.
      else if (k.ge.k0(imax)) then
       tc=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tc=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if

      tcp1_bcll92=tc

      return
      end

c***********************************************************************
c***********************************************************************

      function tcp1_eeho(k)

c***********************************************************************
c  calculate the neutron 1S0 superfluidity critical temperature from   *
c  Elgaroy, Engvik, Hjorth-Jensen & Osnes, Nucl. Phys. A604 (1996), p. 466
c  uses a cubic spline interpolation.                                  *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)
      parameter (imax=13)
      dimension k0(imax),d0(imax),t0(imax),t2(imax)
      save k0,t0,t2,done

      data k0/0.00d0,0.10d0,0.20d0,0.30d0,0.40d0,
     1        0.50d0,0.60d0,0.70d0,0.80d0,0.90d0,
     2        1.00d0,1.10d0,1.20d0/

      data d0/0.00d0,0.02d0,0.10d0,0.37d0,0.62d0,
     1        0.79d0,0.87d0,0.86d0,0.70d0,0.49d0,
     2        0.30d0,0.15d0,0.00d0/

c **** Generate t2 if first access: ***************************
      if (done.eq.1.0) goto 100
      do i=1,imax
       t0(i)=d0(i)/1.76d0*1.1604d10
      end do
      call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      done=1.0
c ***************************************************************
 100  continue
      if (k.le.k0(1)) then
       tcn1=0.
      else if (k.ge.k0(imax)) then
       tcn1=0.
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       tcn1=a*t0(i1)+b*t0(i2)+
     1   ((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
      end if
      tcp1_eeho=tcn1

      return
      end

c************************************************************************
c************************************************************************
c************************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                   lambda 1s0 superfluidity                           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tcla1_bb(k,nbar)

c***********************************************************************
c  calculate the Lambda 1s0 superfluidity critical temperature from    *
c  S. Balberg & N. Barnea, nucl-th/9709013                             *
c  T_c depends on k_F and nbar !                                       *
c  T_c at 5rho_0 is exactly as in the paper, at lower densities it is  *
c   an interpolation which is OK at 10%                                *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and tc is given in kelvins         *
c***********************************************************************

c checked on November4, 1997

      implicit real*8 (a-h,k-z)
c      real*4 k0,t0,t2    ! These must be real*4 for use of `spline'
      parameter (imax=14)
      dimension k0(imax),t0(imax),t2(imax)
      save k0,t0,t2

      data k0/0.1,0.2,0.3,0.4,0.5,0.6,0.7,
     2        0.8,0.9,1.0,1.1,1.2,1.3,1.4/
      data t0/0.000,0.130,0.385,0.720,1.040,1.306,1.468,
     2        1.500,1.393,1.156,0.838,0.535,0.175,0.000/

c **** Get t2 if first acces:
      if (read.eq.0.) then
       call spline_here(k0,t0,imax,0.d0,0.d0,t2)
      end if
      read=1.0
c ***************************************************************
      if (k.le.k0(1)) then
       tcla1_bb=0.
       return
      else if (k.ge.k0(imax)) then
       tcla1_bb=0.
       return
      else     
       i1=1
       i2=imax
1      if ((i2-i1).gt.1) then
        i=(i1+i2)/2
        if (k0(i).gt.k) then
         i2=i
        else
         i1=i
        endif
        goto 1
       end if
       delk=k0(i2)-k0(i1)
       a=(k0(i2)-k)/delk
       b=(k-k0(i1))/delk
       gap=a*t0(i1)+b*t0(i2)
     1     +((a**3-a)*t2(i1)+(b**3-b)*t2(i2))*(delk**2)/6.
       tcla1_bb=gap*(0.220+0.156*(nbar/0.16))
     1          /1.76*1.1604e10
c Adjusts gap to actual background density
c Converts gap (in MeV) to Tc and then to K
      end if

      return
      end


c************************************************************************
c************************************************************************

c************************************************************************

c************************************************************************
c************************************************************************

      subroutine spline_here(x,y,in,yp1,ypn,y2)

c***********************************************************************
c calculate the y" coefficients for a spline interpolation procedure   *
c to be used by splint.                                                *
c yp1 and ypn are the first derivative of y(x) at y(1) and y(in)       *
c                                                                      *
c from numerical recipes, p.88                                         *
c***********************************************************************

      implicit real*8 (a-h,k-z)
      parameter (jmax=100)
      dimension x(in),y(in),y2(in),u(jmax)

      if (yp1.gt.1.e30) then
       y2(1)=0.
       u(1)=0.
      else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,in-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn.gt.1.e30) then
       qn=0.
       un=0.
      else
       qn=0.5
       un=(3./(x(in)-x(in-1)))*(ypn-(y(in)-y(in-1))/(x(in)-x(in-1)))
      end if
      y2(in)=(un-qn*u(in-1))/(qn*y2(in-1)+1.)
      
      do ik=in-1,1,-1
       y2(ik)=y2(ik)*y2(ik+1)+u(ik)
      end do

      return
      end   



