       function mstn_pref(k)

c***********************************************************************
c  calculate the neutron effective mass, prefered values from          *
c   analysis of all data                                               *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)

      mstn_pref=min(1.d0,1.09-0.11*k)
      return
      end

c***********************************************************************
c***********************************************************************
c***********************************************************************

      function mstn_awp_3(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  awp  phys.lett. 222(1989): p.173. case iii                          *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=6)
      dimension k0(imax),m0(imax)
      save k0,m0

      data k0/0.6,0.8,1.0,1.2,1.4,1.6/
      data m0/1.004,0.997,0.982,0.961,0.936,0.910/


      if (k.le.k0(1)) then
       mstn_awp_3=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_awp_3=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_awp_3=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_awp_3_d(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  awp  phys.lett. 222(1989): p.173. case iii, direct contribution only*
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=6)
      dimension k0(imax),m0(imax)
      save k0,m0

      data k0/0.6,0.8,1.0,1.2,1.4,1.6/
      data m0/0.958,0.931,0.915,0.894,0.871,0.851/


      if (k.le.k0(1)) then
       mstn_awp_3_d=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_awp_3_d=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_awp_3_d=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_ccks_r_var(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccks nucl.phys. a451(1986): p.509. reid pot. with variational method*
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=11)
      dimension k0(imax),m0(imax)
      save k0,m0


      data k0/0.3,0.4,0.5,0.6,0.7,
     1        0.8,0.9,1.0,1.1,1.2,1.3/
      data m0/0.99,0.98,0.97,0.96,0.95,
     1        0.94,0.92,0.90,0.88,0.86,0.84/


      if (k.le.k0(1)) then
       mstn_ccks_r_var=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_ccks_r_var=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_ccks_r_var=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_ccks_bj_var(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccks nucl.phys. a451(1986): p.509. bj pot. with variational method  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=11)
      dimension k0(imax),m0(imax)
      save k0,m0

      data k0/0.3,0.4,0.5,0.6,0.7,
     1        0.8,0.9,1.0,1.1,1.2,1.3/
      data m0/0.99,0.98,0.96,0.94,0.92,
     1        0.90,0.88,0.85,0.82,0.78,0.75/


      if (k.le.k0(1)) then
       mstn_ccks_bj_var=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_ccks_bj_var=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_ccks_bj_var=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_ccks_r_cbf(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccks nucl.phys. a451(1986): p.509. reid pot. with cbf method        *
c   completed by jkms nucl.phys. a386(1982): p.125.                    *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=18)
      dimension k0(imax),m0(imax),m2(imax)
      save k0,m0,m2

      data k0/0.30,0.40,0.50,0.60,0.70,0.80,
     1        0.90,1.00,1.10,1.20,1.30,1.50,
     2        1.75,2.00,2.25,2.50,2.75,3.00/
      data m0/1.10,1.10,1.095,1.09,1.08,1.07,
     1        1.06,1.05,1.04,1.02,1.00,0.96,
     1        0.89,0.82,0.75,0.69,0.65,0.64/
      data m2/+4.97e-01,-9.95e-01,+4.81e-01,-9.31e-01,+2.41e-01,
     1        -3.40e-02,-1.05e-01,+4.54e-01,-1.71e+00,+3.93e-01,
     2        +1.39e-01,-6.13e-01,+1.76e-01,-9.23e-02,+1.93e-01,
     3        +2.80e-01,+6.06e-01,+1.77e-01/


      if (k.le.k0(1)) then
       mstn_ccks_r_cbf=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_ccks_r_cbf=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)+
     1   ((a**3-a)*m2(i1)+(b**3-b)*m2(i2))*(delk**2)/6.

       mstn_ccks_r_cbf=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_ccks_bj_cbf(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccks nucl.phys. a451(1986): p.509. bj pot. with cbf method          *
c   completed by jkms nucl.phys. a386(1982): p.125.                    *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=18)
      dimension k0(imax),m0(imax),m2(imax)
      save k0,m0,m2

      data k0/0.30,0.40,0.50,0.60,0.70,0.80,
     1        0.90,1.00,1.10,1.20,1.30,1.50,
     2        1.75,2.00,2.25,2.50,2.75,3.00/
      data m0/0.83,0.91,0.97,1.01,1.03,1.04,
     1        1.04,1.03,1.01,0.98,0.94,0.87,
     2        0.77,0.67,0.59,0.53,0.47,0.43/
      data m2/+2.92e+01,-1.04e+01,+3.37e-01,-2.97e+00,-4.73e-01,
     1        -1.14e+00,-9.52e-01,-1.05e+00,-8.58e-01,-1.52e+00,
     2        +9.38e-01,-5.55e-01,+4.79e-02,+3.64e-01,+4.17e-01,
     3        -1.12e-01,+3.21e-02,+1.90e+00/


      if (k.le.k0(1)) then
       mstn_ccks_bj_cbf=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_ccks_bj_cbf=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)+
     1   ((a**3-a)*m2(i1)+(b**3-b)*m2(i2))*(delk**2)/6.

       mstn_ccks_bj_cbf=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_t(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  t72, prog.theor.phys. 48(1972): p.1517 and                          *
c  t84, prog.theor.phys. 71(1984): p.1432                              *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=18)
      dimension k0(imax),m0(imax)
      save k0,m0

      data k0/0.150,0.197,0.242,0.329,0.446,0.518,
     1        0.613,0.811,0.972,1.090,1.190,1.300,
     2        1.370,1.480,1.620,1.920,2.130,2.400/
      data m0/1.00,1.00,1.00,1.00,1.00,0.99,
     1        0.98,0.96,0.95,0.94,0.93,0.91,
     2        0.89,0.86,0.82,0.78,0.75,0.70/


      if (k.le.k0(1)) then
       mstn_t=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_t=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_t=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstn_bks(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  bks phys.lett. 43b(1973): p.263.                                    * 
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=10)
      dimension k0(imax),m0(imax),m2(imax)
      save k0,m0,m2

      data k0/0.62,0.88,1.05,1.20,1.37,
     1        1.48,1.62,1.92,2.13,2.40/
      data m0/1.02,0.99,0.96,0.94,0.91,
     1        0.89,0.87,0.83,0.81,0.80/
c      data m2/-1.25e+00,-1.62e-01,+5.74e-01,-5.42e-01,+1.15e-02,
c     1        +4.88e-01,-8.04e-02,+1.99e-01,+2.38e-01,+2.93e-01/


      if (k.le.k0(1)) then
       mstn_bks=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_bks=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)+
     1   ((a**3-a)*m2(i1)+(b**3-b)*m2(i2))*(delk**2)/6.

       mstn_bks=m
      end if

      return
      end

c***********************************************************************
c***********************************************************************

      function mstn_bbllp(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  Broglia et al, Phys.Rev. D***, 1994, p. ****                        *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=5)
      dimension k0(imax),m0(imax)
      save k0,m0

      data k0/0.250,0.510,0.810,1.080,1.300/
      data m0/0.995,0.980,0.955,0.925,0.910/


      if (k.le.k0(1)) then
       mstn_bbllp=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstn_bbllp=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)

       mstn_bbllp=m
      end if

      return
      end

c************************************************************************

c************************************************************************

c************************************************************************

      function mstp_ccy_ms2(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccy nucl.phys. a179(1972): p.320. from fig.3 case msii)             *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),m0(imax),m2(imax)
      save k0,m0,m2

      data k0/0.0,0.05,0.1,0.15,0.2,0.3,0.4,
     1        0.6,0.8,1.0,1.2,1.4,1.6,1.8/
      data m0/1.000,0.940,0.870,0.820,0.775,0.710,0.665,
     1        0.610,0.575,0.560,0.550,0.555,0.560,0.565/
      data m2/+5.62e+00,-1.12e+01,+1.53e+01,-2.05e+00,+4.87e+00,
     1        +1.42e+00,+1.45e+00,+2.04e-01,+7.39e-01,-1.61e-01,
     2        +6.54e-01,-2.05e-01,+1.66e-01,-4.58e-01/


      if (k.le.k0(1)) then
       mstp_ccy_ms2=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstp_ccy_ms2=m0(imax)
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
c       m=a*m0(i1)+b*m0(i2)
       m=a*m0(i1)+b*m0(i2)+
     1   ((a**3-a)*m2(i1)+(b**3-b)*m2(i2))*(delk**2)/6.

       mstp_ccy_ms2=m
      end if

      return
      end

c************************************************************************
c************************************************************************

      function mstp_ccy_ps1(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  ccy nucl.phys. a179(1972): p.320. from fig.3 case psi)              *
c  uses a cubic spline interpolation.                                  *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on feb 22 1991

      implicit real*8 (a-h,k-z)
      parameter (imax=14)
      dimension k0(imax),m0(imax),m2(imax)
      save k0,m0,m2

      data k0/0.0,0.05,0.1,0.15,0.2,0.3,0.4,
     1        0.6,0.8,1.0,1.2,1.4,1.6,1.8/
      data m0/1.000,0.940,0.872,0.820,0.780,0.710,0.650,
     1        0.555,0.480,0.425,0.385,0.360,0.340,0.320/
      data m2/+4.27e+00,-8.53e+00,+1.07e+01,+4.30e+00,+9.42e-01,
     1        +1.03e+00,+9.58e-01,+3.65e-01,+5.84e-01,+2.99e-01,
     2        +4.69e-01,+7.57e-02,-2.16e-02,+1.08e-02/


      if (k.le.k0(1)) then
       mstp_ccy_ps1=m0(1)
       return
      else if (k.ge.k0(imax)) then
       mstp_ccy_ps1=m0(imax)
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
       m=a*m0(i1)+b*m0(i2)+
     1   ((a**3-a)*m2(i1)+(b**3-b)*m2(i2))*(delk**2)/6.

       mstp_ccy_ps1=m
      end if

      return

      end

c***********************************************************************

c***********************************************************************

c***********************************************************************

      function mstp_awp(k)

c***********************************************************************
c  calculate the neutron effective mass from :                         *
c  awp in nato asi proceeding of heraklion meeting.                    *
c  k is the fermi momentum in fm^-1 and m*=mst is given as the m*/m    *
c***********************************************************************

c checked on 

      implicit real*8 (a-h,k-z)

      mstp_awp=1.35*mstp_ccy_ms2(k)

      return

      end



