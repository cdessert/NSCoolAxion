      dimension rad(0:isize),rrho(0:isize),pres(0:isize),
     2          debar(0:isize),dvol(0:isize),
     3          emas(0:isize),phi(0:isize),rhod(0:isize)
      common/profile_star/rad,rrho,pres,
     2                    debar,dvol,
     3                    emas,phi,rhod
c*************************************************
c  dvol(i) is the volume of the half-shell between 
c   the spheres of radii rad(i) and rad(i+1)
c      (As calculated in subroutine grid)
c*************************************************
      common/star/imax,icore,idrip,ienv

