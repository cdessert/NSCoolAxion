      dimension tcn(0:isize),tcp(0:isize),tcla(0:isize),
     2          tcsm(0:isize),tcs0(0:isize),tcsp(0:isize),
     3          tcuu(0:isize),tcdd(0:isize),tcss(0:isize),
     4          tcud(0:isize),tcus(0:isize),tcds(0:isize),
     5          tcu(0:isize),tcd(0:isize),tcs(0:isize)

c April 2004: seperate Tc for flavor AND color:
      dimension tcu1(0:isize),tcu2(0:isize),tcu3(0:isize),
     2          tcd1(0:isize),tcd2(0:isize),tcd3(0:isize),
     3          tcs1(0:isize),tcs2(0:isize),tcs3(0:isize)

      common/superfluid/tcn,tcp,tcla,tcsm,tcs0,tcsp,
     2	                tcuu,tcdd,tcss,tcud,tcus,tcds,
     3		        tcu,tcd,tcs,
     4                  tcu1,tcu2,tcu3,tcd1,tcd2,tcd3,tcs1,tcs2,tcs3,
     4                  sfn1s0,sfn3p2,sfp1s0,sfl1s0,sfquark,
     5                  fn1s0,fn3p2,fp1s0,fl1s0,
     6	                kfmax_n3p2,delkf_n3p2,tcmax_n3p2,isf

c**** Yakovlev et al pairing:
c      common/yak_pairing/T0_ns,k0_ns,k1_ns,k2_ns,k3_ns,
c     2                   T0_nt,k0_nt,k1_nt,k2_nt,k3_nt,
c     3                   T0_p ,k0_p ,k1_p ,k2_ p,k3_p




