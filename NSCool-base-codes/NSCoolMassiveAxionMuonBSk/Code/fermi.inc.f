      dimension mstp(0:isize),mstn(0:isize),mstla(0:isize),
     2          mstsm(0:isize),msts0(0:isize),mstsp(0:isize)
      dimension kfe(0:isize),kfm(0:isize),kfp(0:isize),kfn(0:isize),
     2          kfla(0:isize),kfsm(0:isize),kfs0(0:isize),kfsp(0:isize),
     3          kfqu(0:isize),kfqd(0:isize),kfqs(0:isize)
      common/fermi/mstp,mstn,mstla,mstsm,msts0,mstsp,
     2             kfe,kfm,kfp,kfn,kfla,kfsm,kfs0,kfsp,
     3             kfqu,kfqd,kfqs
      common/effmass/emnco,emncr,emp

c-------------------------------------------------------------------
c The effective masses ms** are the ratio of the 
c Landau effective mass to the free mass.
c The Fermi momenta are in units of fm^-1
c-------------------------------------------------------------------

