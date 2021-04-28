cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                      Yakovlev et al. functions                       c
c                                                                      c
c                                                                      c
c     Values from:                                                     c
c     Yakovlev, Kaminker, & Gnedin, A&A 379, L5   (2001)               c
c     Kaminker, Yakovlev, & Gnedin, A&A 383, 1076 (2002)               c
c     Yakovlev, Kaminker, Haensel, & Gnedin, A&A 389, L24 (2002)       c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c***********************************************************************
c***********************************************************************
c     Proton 1S0:
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_1p(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_1p=Tc_Ioffe(kf,
     x        20.29d0 , 0.d0 , 1.117d0 , 1.241d0 , 0.1473d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_2p(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_2p=Tc_Ioffe(kf,
     x        17.00d0 , 0.d0 , 1.117d0 , 1.329d0 , 0.1179d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_3p(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_3p=Tc_Ioffe(kf,
     x        14.50d0 , 0.d0 , 1.117d0 , 1.518d0 , 0.1179d0)
       return
      end
c***********************************************************************
c***********************************************************************
c     Neutron 1S0:
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_1ns(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_1ns=Tc_Ioffe(kf,
     x        10.20d0 , 0.d0 , 0.6d0 , 1.45d0 , 0.1d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_2ns(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_2ns=Tc_Ioffe(kf,
     x        07.90d0 , 0.d0 , 0.3d0 , 1.45d0 , 0.01d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_3ns(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_3ns=Tc_Ioffe(kf,
     x        1800.d0 , 0.d0 , 21.d0 , 1.45d0 , 0.4125d0)
       return
      end
c***********************************************************************
c***********************************************************************
c     Neutron 3P2:
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_1nt(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_1nt=Tc_Ioffe(kf,
     x        6.461d0 , 1.d0 , 1.961d0 , 2.755d0 , 1.3d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_2nt(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_2nt=Tc_Ioffe(kf,
     x        2.000d0 , 1.d0 , 1.961d0 , 2.755d0 , 1.3d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe_3nt(kf)
       implicit real*8 (a-h,k-z)
        Tc_Ioffe_3nt=Tc_Ioffe(kf,
     x        15.00d0 , 1.d0 , 1.961d0 , 2.755d0 , 1.3d0)
       return
      end
c***********************************************************************
c***********************************************************************
      function Tc_Ioffe(kf,T0,k0,k1,k2,k3)
       implicit real*8 (a-h,k-z)
         if ((kf.gt.k0).and.(kf.lt.k2)) then
          Tc=T0 * (kf-k0)**2/((kf-k0)**2+k1**2)*
     x            (kf-k2)**2/((kf-k2)**2+k3**2)
          Tc_Ioffe=1.d9*Tc
         else
          Tc_Ioffe=1.d0
         end if
       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
