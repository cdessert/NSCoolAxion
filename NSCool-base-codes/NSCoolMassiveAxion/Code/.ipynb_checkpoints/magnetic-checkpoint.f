      subroutine get_dipole(imax,rad,f,g,h,beta,
     1                      bfield0,st,bfr,bft,bf2,cur)
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       parameter(pi=3.1415926535)
       parameter(c=2.99792e10)
       dimension rad(0:isize),sigma(0:isize),
     1           st(0:isize),bfr(0:isize),bft(0:isize),bf2(0:isize),
     2           cur(0:isize)
       dimension f(0:isize),g(0:isize),h(0:isize)
        bfr(1)=0.
        bft(1)=0.
        cur(1)=0.
        bf2(1)=0.
        do i=3,imax,2
         bfr(i)=beta*bfield0*rad(imax)**2*
     1          st(i)/rad(i)**2
         dsdr=(st(i+2)-st(i-2))/(rad(i+2)-rad(i-2))
         bft(i)=beta*bfield0*rad(imax)**2*
     1          dsqrt(g(i))*dsdr/rad(i)/2.d0
         bf22=4./3.*bfr(i)**2+2./3.*bft(i)**2
c         bf2(i)=(beta*bfield0)**2*rad(imax)**4*
c     2              (4./3.*st(i)**2/rad(i)**4+
c     3               2./3.*g(i)*dsdr**2/rad(i)**2)
         bf2(i)=sqrt(bf22)
         dsdrp=(st(i+2)-st( i ))/(rad(i+2)-rad( i ))
         dsdrm=(st( i )-st(i-2))/(rad( i )-rad(i-2))
         d2sdr2=(dsdrp-dsdrm)/(rad(i+1)-rad(i-1))
         cur(i)=-beta*bfield0*rad(imax)**2*c/4./pi/rad(i)*
     1           (g(i)*d2sdr2+h(i)*dsdr-2.*st(i)/rad(i)**2)
        end do
       return
      end
c ****************************************************************************
c ****************************************************************************
c ****************************************************************************
      subroutine initialize_dipole(i0,i1)
       implicit real*8 (a-h,k-z)
       parameter(pi=3.1415926535)
       parameter(c=2.99792e10)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
c       INCLUDE 'files_in_out.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'mag_field.inc.f'
       character*80 temp_file

       dimension zd(0:isize)
c **** calculate the initial Stoke profile ********************************
      makestoke=2.
      if (makestoke.eq.1.) then
       goto 150                              ! URME
      else if (makestoke.eq.2.) then
       goto 160                              ! Urpin
      else if (makestoke.eq.3.) then
       goto 170                              ! Accreted
      else if (makestoke.eq.4.) then
       goto 180                              ! Accretion
      else if (makestoke.eq.5.) then
       goto 190                              ! Read from Teff file
      else
       pause 'Stoke definition not precised '
      end if
c ULRICH's form: *******************************************
 150  continue
      s0=1.

      hr=1.e3
      zmin=1.e4/hr
c      zmin=0.
      radius=rad(imax)
      do i=0,imax
       zd(i)=zmin+(radius-rad(i))/hr
      end do
      izm=i0
      zmax=zd(izm)
      f1=(2.*zmin)/zmax**2 + hr/radius - (zmin/zmax)**2*hr/radius
      f2=1.-(zmin/zmax)**2
      do i=imax,izm,-1
       stoke(i)=s0*(1.-(zd(i)/zmax)**2) * exp(f1*zd(i)/f2)
      end do
      do i=izm-1,0,-1
       stoke(i)=0.
      end do
      goto 200
c Mirralles, Urpin & Konenkov's form: ***************************
 160  continue
      r0=rad(i0)
      do i=i0,imax
       stoke(i)=(1.-rad(  i )**2/r0**2)/
     1          (1.-rad(imax)**2/r0**2)
      end do
      do i=0,i0-1
       stoke(i)=0.
      end do
      goto 200
c Accreted matter: **********************************************
 170  continue
      do i=0,imax
       stoke(i)=1.d-200
      end do
      mstoke=0.
      do i=i0,i1
       stoke(i)=(rad(i)-rad(i0))*(rad(i1)-rad(i))
       if (stoke(i).gt.mstoke) then
        mstoke=stoke(i)
        im=i
       end if
      end do
      do i=0,imax
       stoke(i)=stoke(i)/mstoke
      end do
      goto 200
c For accretion (as in ApJLett: ************************************
 180  continue
      do i=0,imax
       stoke(i)=1.d-200
       zd(i)=rad(imax)-rad(i)
      end do
      mstoke=0.
      z0=zd(i0)      
      do i=0,imax
       stoke(i)=zd(i)*(rad(imax)-zd(i))/(1.+exp(z0/zd(i)))
       if (stoke(i).gt.mstoke) then
        mstoke=stoke(i)
        im=i
       end if
      end do
      do i=0,imax
       stoke(i)=stoke(i)/mstoke
      end do
      goto 200
c Read stoke fonction from the Temp file of a previous run:
 190  continue
      write(6,*)'Temp file for reading the stoke function='
      read(5,*)temp_file
      open(unit=25,file=temp_file,status='old')
       do i=1,3
        read(25,*)
       end do
       do iread=1,100
         read(25,*,err=999)
         read(25,'(13x,1p1e9.3)',err=999)time
	 print *,'Time=',time
         read(25,*)
         read(25,*)
         read(25,*)
         read(25,*)
        do i=imax,1,-2
         read(25,*,err=999)
     1       ijunk,rr,dd,tt,ll,
     2       bb,ss
         if (iread.eq.10) stoke(i)=ss
        end do
       end do
 999   continue
      close(unit=25,status='keep')
      scale=stoke(imax)
c      scale=1.d0
      print *,'Scale=',scale
      do i=imax,1,-2
       stoke(i)=stoke(i)/scale
      end do
      goto 200
c ***************************************************************

 200   continue

c ***** Add a sin wave to stoke:
c      slen=rad(imax)/(4.*pi)
c      samp=0.0
c      do i=0,imax
c       stoke(i)=stoke(i)*(1.+samp*sin((rad(i)-rad(imax))/slen))
c      end do
c ******
      return
      end
c ****************************************************************************
c ****************************************************************************
c ****************************************************************************
      subroutine solve_induction_dipole(imax,dtime,idt,
     1                                  rad,sigma,velocity,
     2                                  f,g,h,alpha,
     3                                  st,ost)
C     CHECKED on March 19 (with Ulrich !), 1999
       implicit real*8 (a-h,k-z)
c       real*16 st,ost,a,b,c,d,dd,drp,dr,drm,p,q,rrr
       common/stuff/time,istep
       parameter(pi=3.1415926535,cl=2.99792e10)
       INCLUDE 'size.inc.f'
       dimension rad(0:isize),sigma(0:isize),velocity(0:isize),
     1           st(0:isize),ost(0:isize),ost1(0:isize)
       dimension f(0:isize),g(0:isize),h(0:isize)
       dimension a(0:isize),b(0:isize),c(0:isize),d(0:isize),
     1           dd(0:isize)
       dimension p(0:isize),q(0:isize)
c *****
       dt=dtime/float(idt)
c *****
       do i=3,imax-2,2
        drp=rad(i+2)-rad(i)
        dr=rad(i+1)-rad(i-1)
        drm=rad(i)-rad(i-2)
        a(i)=g(i)/drm/dr + 
     1       (+4.*pi*sigma(i)/cl**2*velocity(i)-h(i))/(drp+drm)
        b(i)=-4.*pi*sigma(i)/cl**2/(f(i)*dt)
     1       -2./rad(i)**2-g(i)/drp/dr-g(i)/drm/dr
        c(i)=g(i)/drp/dr + 
     1       (-4.*pi*sigma(i)/cl**2*velocity(i)+h(i))/(drp+drm)
        dd(i)=-4.*pi*sigma(i)/cl**2/(f(i)*dt)
       end do

       do i=1,imax
        st(i)=ost(i)
       end do

       do it=1,idt
        do i=3,imax-2,2
         d(i)=dd(i)*st(i)
        end do

c -------------------------------------------------
c Inner boundary condition: whole star:
        p(1)=0.
        ggg=1./2.
        gamma=((1.-ggg)/(2.-ggg))*rad(1)/(rad(3)-rad(1))
        q(1)=gamma/(1.+gamma)
        do i=3,imax-2,2
         p(i)=(d(i)-a(i)*p(i-2))/(b(i)+a(i)*q(i-2))
         q(i)=       -c(i)      /(b(i)+a(i)*q(i-2))
        end do
c -------------------------------------------------
cc Inner boundary condition: crustal field:
c        do i=1,111
c         p(i)=0.0
c         q(i)=0.0
c        end do     
c        do i=113,imax-2,2
c         p(i)=(d(i)-a(i)*p(i-2))/(b(i)+a(i)*q(i-2))
c         q(i)=       -c(i)      /(b(i)+a(i)*q(i-2))
c        end do
c -------------------------------------------------
c Outer boundary condition:
c With left derivative at the surface:
        rrr=((1.+alpha)*rad(imax)-alpha*rad(imax-2))/rad(imax)
c With centered derivative at the surface:     
c        rrr=((1.+alpha)*rad(imax)+(1.-alpha)*rad(imax-2))/
c     1       ((1.-alpha)*rad(imax)+(1.+alpha)*rad(imax-2))
c
        st(imax)=p(imax-2)/(rrr-q(imax-2))
c
        do i=imax-2,1,-2
         st(i)=p(i)+q(i)*st(i+2)
        end do
       end do

C CHECK MATRIX INVERSION: CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       do i=3,imax-2,2
c        left=a(i)*st(i-2)+b(i)*st(i)+c(i)*st(i+2)
c        right=d(i)
c        if ((istep.eq.10*(istep/10)).or.(istep.eq.1)) then
c         print '(i5,0p1f10.3,3x,1p1e12.3,3x,1p2e16.6,3x,0p1f16.8)',
c     1          i,rad(i)/1.e5,st(i),left,right,left/right
c        end if
c       end do
C CHECK INDUCTION EQUATION SOLUTION: CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       do i=3,imax-2,2
c        drp=rad(i+2)-rad(i)
c        dr=rad(i+1)-rad(i-1)
c        drm=rad(i)-rad(i-2)
c        dsdr_p=(st(i+2)-st( i ))/drp
c        dsdr_m=(st( i )-st(i-2))/drm
c        d2sdr2=(dsdr_p-dsdr_m)/dr
c        left=d2sdr2-2.*st(i)/rad(i)**2
c        right=4.*pi*sigma(i)/cl**2*(st(i)-ost(i))/dt
c        if ((istep.eq.10*(istep/10)).or.(istep.eq.1)) then
c         print '(i5,0p1f10.3,3x,1p1e12.3,3x,1p2e16.6,3x,0p1f16.8)',
c     1          i,rad(i)/1.e5,st(i),left,right,left/right
c        end if
c       end do
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       return
      end
c ****************************************************************************
c ****************************************************************************
c ****************************************************************************
      subroutine get_dipole_force(imax,rad,ostoke,stoke,sigma,dtime,
     1                            bfield0,force_r,force_t)
c NO GR TERMS INCLUDED YET
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       parameter(pi=3.1415926535)
       parameter(c=2.99792e10)
       dimension rad(0:isize),sigma(0:isize),
     1           stoke(0:isize),ostoke(0:isize),
     2           dsdr(0:isize),dsdt(0:isize),
     3           force_r(0:isize),force_t(0:isize)

        force_r(1)=0.d0
        force_t(1)=0.d0
        do i=3,imax-2,2
         dsdt(i)=(stoke(i)-ostoke(i))/dtime
         dsdr(i)=(stoke(i+2)-stoke(i-2))/(rad(i+2)-rad(i-2))
         coeff=sigma(i)/c**2/rad(i)**2*(bfield0*rad(imax)**2)**2
         force_r(i)=-coeff*     dsdt(i)* dsdr(i)
         force_t(i)=+coeff*2.d0*dsdt(i)* stoke(i)/rad(i)
        end do
        force_r(imax)=force_r(imax-2)
        force_t(imax)=force_t(imax-2)
       return
      end
c ****************************************************************************
c ****************************************************************************
c ****************************************************************************



