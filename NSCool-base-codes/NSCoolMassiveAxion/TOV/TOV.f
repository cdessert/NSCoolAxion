      program ov
     
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4
      character*60 title
      parameter (ipmax=10000)
     
      dimension e(0:50000), den(0:50000), r(0:50000), 
     x          p(0:50000), em(0:50000), baryden(0:50000)
      dimension prho(ipmax), phi(0:50000)
     
      common/interp/ i1,i2,i3,limit
      common/const/ pi,g
      common/ans/ bfunc, pmfunc, emfunc, pfunc, phfunc
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
      common/prof/r,den,e,p,em,phi,rhoc,baryden,propmas,delta,stepi
     
      logical iread, fread, cread, dread, sread,
     C        prod, output, stndrd, qseek, seek,
     C        dum_l     
c      data prho/.1,.11,.12,.135,.15,.165,.18,.2,.22,.246,.27,.3,
c     x .33,.37,.4,.45,.5,.55,.6,.67,.74,.82,.9,
c     x 1.,1.1,1.2,1.35,1.5,1.65,1.8,2.,2.2,2.46,2.7,3./
     
      prho(1)=0.1
c      prho(1)=0.31
c      prho(1)=0.277
      factor=30.**(1./float(ipmax))
      do i=2,ipmax
       prho(i)=prho(1)*factor**i
      end do

c      do i=1,ipmax
c       prho(i)=0.2*(float(i)/float(ipmax))
c      end do

      pi=acos(-1.)
      g=1.484
     
C ***** Defines default values *******************************************
     
      stepi=.01
      delta=.5
      title=' '
      limit=96
      output=.false.
      rhoc=.4
      rhol=.1324
c This correspond to a density of 2.2e14 g/cm3 for the core limit
      rhod=2.573e-4
c This correspond to a density of 4.3e11 g/cm3 for neutron drip
      const = 1.d0
      c2=1.d0
      target=1.44
     
1     prod=.false.
      stndrd=.false.
     
C ***** Initializes result arrays ****************************************
     
      do i=0,50000
      em(i)=0.
      den(i)=0.
      p(i)=0.
      r(i)=0.
      e(i)=0.
      phi(i)=0.
      end do
     
C ***** Waits for instructions *********************************************
     
5     call gread(5,6,.true.,.true.)
      if(cread('go')) go to 10
      if(cread('production')) then
        open (unit=7,file=
     1        'Production/prod.dat',
     2        status='new',err=1090)
        prod=.true.
        goto 1091
 1090   print *,' > `prod.dat'' already exists: delete please it'
 1091   continue
      endif
c      if(cread('standard')) stndrd = .true.
      if(dread('targetmass', target, dum, 1, stndrd)) goto 5
      if(cread('EOF')) call exit
      if(cread('eof')) call exit
      if(cread('dump')) output=.true.
      if(cread('nod'))  output=.false.
      if(cread('SMU')) then
        const = 1.0d0
        c2=1.0d0
      end if
      if(cread('MeV')) then
        const = 1.0d0 / 1.122d6
        c2=1.0d0
      end if
      if(cread('cgs')) then
        const = 1.0d0 / 1.989d18
        c2=9.d20
      endif
      if(cread('units?')) write(6,*) ' const = ', 1./const
      if(sread('eos', title, dum)) goto 5

      if(cread('how')) write(6,*) ' prod?', prod, '  dump?', output,
     c                            '   standard?', stndrd
      if(cread('what?')) write(6,*) ' eos = ',title
      if(cread('geteos'))  call eos(title,c2)
      if(iread('limit',limit, dum, 1, dum)) goto 5
      if(cread('limit?')) write(6,*) limit
      if(dread('omegac', omebc, dum, 1, dum_l)) goto 5
      if(cread('omegac?')) write(6,*) omebc
      if(dread('rhoc', rhoc, dum, 1, dum_l)) goto 5
      if(cread('rhoc?'))write(6,*) rhoc
      if(dread('rhol', rhol, dum, 1, dum_l)) goto 5
      if(cread('rhol?'))write(6,'(0p1f10.4,1p1e16.3)')
     c                  rhol, 2.2e14*rhol/.1324
      if(dread('fs', stepi, dum, 1, dum_l)) goto 5
      if(cread('fs?')) write(6,*)  stepi
      if(dread('ss', delta, dum, 1, dum_l)) goto 5
      if(cread('ss?')) write(6,*)  delta
      goto 5
     
10    iii=1
11    if(prod) rhoc=prho(iii)
12    continue
     
      i1=2
      i2=2
      i3=2
     
C ***** Do the first step ********************************************
     
      den(0)=rhoc
      p(0)=pres(den(0))
      e(0)=ener(p(0))
      phi(0)=0.
      icore=0
      idrip=0
     
      k1=0.
      l1=0.
      m1=0.
      n1=0.
      p1=0.
      step=stepi
     
      call twostep(r(0)+step/2.,p(0)+k1/2.,em(0)+l1/2.,phi(0)+p1/2.)
      k2 = step *  pfunc
      l2 = step * emfunc
      m2 = step *  bfunc
      n2 = step * pmfunc
      p2 = step * phfunc
     
      call twostep(r(0)+step/2.,p(0)+k2/2.,em(0)+l2/2.,phi(0)+p2/2.)
      k3 = step *  pfunc
      l3 = step * emfunc
      m3 = step *  bfunc
      n3 = step * pmfunc
      p3 = step * phfunc
     
      call twostep(r(0)+step,p(0)+k3,em(0)+l3,phi(0)+p3)
      k4 = step *  pfunc
      l4 = step * emfunc
      m4 = step *  bfunc
      n4 = step * pmfunc
      p4 = step * phfunc
     
      p(1)  =  p(0) + (k1+2.*k2+2.*k3+k4)/6.
      em(1) = em(0) + (l1+2.*l2+2.*l3+l4)/6.
      r(1)  =  r(0) + step
      den(1)= rho(p(1))
      baryden(1) = (m1+2.*m2+2.*m3+m4)/6.
      propmas = (n1+2.*n2+2.*n3+n4)/6.
      e(1)  = ener(p(1))
      phi(1)   = phi(0)   + (p1+2.*p2+2.*p3+p4)/6.
     
C ***** Do the next step **********************************************
     
      do i=1,50000-1
     
      imax=i+1
     
      call twostep(r(i),p(i),em(i),phi(i))
     
      step = delta/(emfunc/em(i)-pfunc/p(i))
     
      k1 = step *  pfunc
      l1 = step * emfunc
      m1 = step *  bfunc
      n1 = step * pmfunc
      p1 = step * phfunc
     
      call twostep(r(i)+step/2.,p(i)+k1/2.,em(i)+l1/2.,phi(i)+p1/2.)
      k2 = step *  pfunc
      l2 = step * emfunc
      m2 = step *  bfunc
      n2 = step * pmfunc
      p2 = step * phfunc
     
      call twostep(r(i)+step/2.,p(i)+k2/2.,em(i)+l2/2.,phi(i)+p2/2.)
      k3 = step *  pfunc
      l3 = step * emfunc
      m3 = step *  bfunc
      n3 = step * pmfunc
      p3 = step * phfunc
     
      call twostep(r(i)+step,p(i)+k3,em(i)+l3,phi(i)+p3)
      k4 = step *  pfunc
      l4 = step * emfunc
      m4 = step *  bfunc
      n4 = step * pmfunc
      p4 = step * phfunc
     
      p(i+1)  = p(i) + (k1+2.*k2+2.*k3+k4)/6.
      em(i+1) = em(i) + (l1+2.*l2+2.*l3+l4)/6.
      r(i+1)  = r(i) + step
      den(i+1)= rho(p(i+1))
      baryden(i+1) = baryden(i) + (m1+2.*m2+2.*m3+m4)/6.
      propmas = propmas + (n1+2.*n2+2.*n3+n4)/6.
      e(i+1)  = ener(p(i+1))
      phi(i+1)   = phi(i)   + (p1+2.*p2+2.*p3+p4)/6.
     
      if ((icore .eq. 0).and.(den(i+1) .le. rhol)) icore=i
     
      if ((idrip .eq. 0).and.(den(i+1) .lt. rhod)) idrip=i
     
      if (den(i+1).lt.deos(limit-1)) goto 1100
     
      end do
     
      write(6,*) ' '
      write(6,*) ' The star surface has not been reached! '
      write(6,*) ' '
     
C ***** The calculation is finished ****************************************
     
1100  continue
     
C ***** Rescales phi ********************************************************
     
      const1=.5*dlog(1.-2.*g*em(imax)/r(imax))-phi(imax)
      do i=0,imax
        phi(i)=phi(i)+const1
      end do
     
C ***** Outputs the profile if 'dump' asked *********************************
     
      if(output) call profile(title,imax,icore,idrip)
     
C ***** Outputs the results if 'production' was asked **********************
     
      if(prod) then
     
        if(iii.eq.1) then
          write(7,*) 5, ipmax
          write(7,*)
          write(7,*) '      EOS file : ', title
          write(7,*)
          write(7,*) '   rhoc', '        rad', '    grav. mass',
     c   '   b. mass ', ' red shift ', ' % binding '
          write(7,*)
        endif
     
        write(7,'( 1p6e11.3 )' ) rhoc,r(imax),em(imax),baryden(imax),
     c         sqrt(1.-2.*g*em(imax)/r(imax)),
     c         (baryden(imax)-em(imax))/em(imax)
      endif
     
C **********************************************************************
     
      write(6,*) ' '
      write(6,*) title
      write(6,*) ' '
      write(6,'('' rhoc:'', 1pe12.4, '' radius:'', 1pe12.4,/,
     x '' masses(g,b,p):'', 1p3e17.9)' )
     x    rhoc,r(imax),em(imax),baryden(imax),propmas
      write(6,*) ' '
      write(6, '( '' step size:'', 1pe12.4, '' initial step:'',
     x 1pe12.4)' ) delta, stepi
      write(6,*) ' '
      write(6,*) ' Here comes another one! '
      write(6,*) ' '

C *********************************************************************
     
      iii=iii+1
      if (stndrd) then
ccccccccccccccccccccccccccccccccccccccccccccccc
        emax = em(imax)-target
c        emax = baryden(imax)-target
ccccccccccccccccccccccccccccccccccccccccccccccc
c        if(.not.seek(emax,rhoc,emold,rhoold,1.e-5,iii-1)
c     1     .and. iii.le.ipmax)  go to 12
        accuracy=1.e-7
        if(.not.qseek(emax,rhoc,emold1,rhoold1,emold2,rhoold2,
     1                accuracy,iii-1)
     1     .and. iii.le.30)  go to 12
     
      elseif (prod) then
        if(iii.le.ipmax) go to 11
        close(unit=7,status='keep')
      endif
     
      go to 1
      end
     
c *************************************************************************
c *************************************************************************
     
      function ener(p)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     ener=0.
     
      if ( p.ge.peos(i1) ) then
        ener = eeos(i1) * exp( log(p/peos(i1)) *
     x                         log(eeos(i1-1)/eeos(i1)) /
     x                         log(peos(i1-1)/peos(i1)) )
        return
      else
        i1=i1+1
      endif
     
      if(i1.gt.limit) then
        i1=limit
        return
      endif
     
      goto 1
     
      end
     
c ************************************************************************
c ************************************************************************
     
      function pres(d)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     pres=0.
     
      if ( d.ge.deos(i2) ) then
        pres = peos(i2) * exp( log(d/deos(i2)) *
     x                         log(peos(i2-1)/peos(i2)) /
     x                         log(deos(i2-1)/deos(i2)) )
        return
      else
        i2=i2+1
      endif
     
      if(i2.gt.limit) then
        i2=limit
        return
      endif
     
      goto 1
     
      end
     
c **************************************************************************
c **************************************************************************
     
      function rho(p)
     
      implicit real*8 (a-h,o-z)
     
      common/interp/ i1,i2,i3,limit
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
     
1     rho=0.
     
      if ( p.ge.peos(i3) ) then
        rho = deos(i3) * exp( log(p/peos(i3)) *
     x                         log(deos(i3-1)/deos(i3)) /
     x                         log(peos(i3-1)/peos(i3)) )
        return
      else
        i3=i3+1
      endif
     
      if(i3.gt.limit) then
        i3=limit
        return
      endif
     
      goto 1
     
      end
     
c *************************************************************************
c *************************************************************************
c *************************************************************************
     
      subroutine eos(title,c2)
     
      implicit real*8 (a-h,o-z)
      character*60 title
     
      common/eos_dat/ peos(50000),eeos(50000),deos(50000),const
      common/interp/ i1,i2,i3,limit
c      write(6,*) c2
      title='/Users/dany/Work/NStar/EOS/'//title
      open (unit=15, file=title, status='old')
c       read(15,*) itext,limit
       itext=6
       limit=1000
       do i=1,itext
        read(15,*)
       end do
       do i=1,limit
        read(15,*,err=100,end=100) x1,x2,x3
c Check order of columns in the EOS file:
        if(i.eq.1) then
         if ((x3.le.10.).and.(x2.ge.1e30)) then
          ilist=1
         else if ((x1.le.10.).and.(x3.ge.1.e30))then
          ilist=2
         else
          pause 'Check EOS column ordering !'
         end if
        end if
        if (ilist.eq.1) then
         eeos(i)=x1
         peos(i)=x2
         deos(i)=x3
        else if (ilist.eq.2) then
         eeos(i)=x2
         peos(i)=x3
         deos(i)=x1
        else
         pause 'Check EOS column ordering !'
        end if
        peos(i)=peos(i) * const / c2
        eeos(i)=eeos(i) * const
       end do
 100  continue
      limit=i-1
      close(unit=15, status='keep')

      do i=1,limit
       write(6,10010) i,eeos(i)/const,peos(i)/const*c2,deos(i)
      end do
10010 format(1i5,1p3e16.4)
     
      return
      end
     
c ************************************************************************
c ************************************************************************
c ************************************************************************
     
      subroutine twostep (r, p, em, phi )
     
      implicit real*8 (a-h,o-z)
      real*8 infunc, nodrag,  nonrel
     
      common/ans/ bfunc, pmfunc, emfunc, pfunc, phfunc
      common/const/ pi,g
     
      dens = rho(p)
      energy = ener(p)
      vol = 4.*pi*r*r*r
     
      bfunc = dens *4.*pi*r*r / sqrt( 1.-2.*g*em / r )
      bfunc = bfunc * 8.42e-4
     
c  this converts to solar masses from [#-km**3/fm**3]
c  I hope. this is the baryonic mass and is good for
c  computing the gravitational binding energy of a star.
     
      pmfunc = energy *4.*pi*r*r / sqrt( 1.-2.*g*em / r )
     
      pfunc=-g * (energy+p) * (em+vol*p) / ( r*(r-2.*em*g) )
     
      emfunc= 4.*pi*r*r * energy
     
      phfunc = g * (em+vol*p) / ( r*(r-2.*em*g) )
     
      return
     
      end
     
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
     
      subroutine profile(title,imax,icore,idrip)
     
      implicit real*8 (a-h,o-z)
      character*60 title
      character*2 error
     
      dimension e(0:50000), den(0:50000), r(0:50000), 
     x          p(0:50000), em(0:50000), baryden(0:50000)
      dimension phi(0:50000)
     
      common/prof/r,den,e,p,em,phi,rhoc,baryden,propmas,delta,stepi
     
      write(6,*) ' Profile being output. '
     
 455  open (unit=15,file=
     1      'Profile/prof.dat', 
     2      status='new',err=456)
      write(15,10001) 6, imax, icore, idrip
      write(15,*)
      write(15,*) '    EOS file :  ',title
      write(15,*)
      write(15,10002)'step','radius ','baryon#  ','density   ',
     x               '  pressure  ','encl. mass   ','phi   ',
     x               'encl. bar. mass '
      write(15,10002) '   ','  (m)  ','(#/fm3)  ','(g/cm3)   ',
     x               '  (dyn/cm2) ',' (sol. mass)  ','     ',
     x               ' (sol. mass)   '
      write(15,*)
     
      do i=0,imax
       write(15,10000) i,r(i)*1.e3,den(i),e(i)*1.989e18,
     x                 p(i)*1.989e18*9e20,em(i),phi(i),
     x                 baryden(i)
      end do
     
      close (unit=15, status='keep')
     
10000 format(i6,0pf15.6,1p1e15.6,1e15.6,1e15.5,1e18.9,1e15.6,1e18.9)
10001 format(4i8)
10002 format(a6,a15,a15,a15,a15,a18,a15,a18)

      goto 457

 456  print *,
     1  'prof.dat file already exists: delete it and type `go'' '
      read(5,*)error
      goto 455

 457  return
     
      end
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c          Below is included content of old file "seek.for"
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************    
c **
      logical function seek(f1,x1,f2,x2,tol,i)
c **
      implicit real*8 (a-h,o-z)
c **
      seek=.true.
      if(abs(f1).lt.tol) return
c **
      seek=.false.
      if(i.eq.1) then
        xnew=.9*x1
      else
        xnew=(f2*x1-f1*x2)/(f2-f1)
      endif
c **
      f2=f1
      x2=x1
      x1=xnew
c **
      return
      end
c **
c **
      logical function qseek(f1,x1,f2,x2,f3,x3,tol,i)
c **
      implicit real*8 (a-h,o-z)
      logical seek
c **
      qseek=.true.
      if(abs(f1).lt.tol) return
c **
      qseek=.false.
      if (i.gt.2) then
c **
        a=(x1-x2)*f3+(x2-x3)*f1+(x3-x1)*f2
        a=-a/((x1-x2)*(x2-x3)*(x3-x1))
        b=(f1-f3)/(x1-x3)-(x1+x3)*a
        c=f2-a*x2*x2-b*x2
        d=b*b-4.0*a*c
        xnew=-0.5*b/a
c **
        if (d.lt.0.0) then
          if (abs(xnew-x1).lt.tol.or.abs(f1-f2).lt.tol) qseek=.true.
        else if (xnew.lt.x1) then
          xnew=xnew + 0.5*abs(sqrt(d)/a)
        else
          xnew=xnew - 0.5*abs(sqrt(d)/a)
        endif
c **
        f3=f2
        f2=f1
        x3=x2
        x2=x1
        x1=xnew
c **
      else
c **
        f3=f2
        x3=x2
        qseek=seek(f1,x1,f2,x2,tol,i)
        x1=max( x2*x2/x1, min( x1, x1*x1/x2 ) )
c **
      endif
c **
      return
      end
c **
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c          Below is included content of old file "cname.for"
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************    
      logical function gread(unit1,init2,echo,object)
      integer unit1,unit2,quote,iary(*)
      real rary(*)
      real*8 d,dary(*)
c      logical object,echo,dblank,image/.false./,new
      logical object,echo,dblank,image,new
      logical aread,bread,cread,dread,fread,iread,sread,tread
      character*(*) name,sary
      character*80 card
c      character*4 eof/'eof '/,blanks/'    '/
      character*4 eof,blanks
c      character*1 blank/' '/,comma/','/
      character*1 blank,comma
c
      save unit2,card,length,ieq,quote,image
c
      image=.false.
      eof='eof'
      blanks='    '
      blank=' '
      comma=','
c
      unit2=init2
      if(image.and.object) write(unit2,102)
      write(unit2,'('' >'',$)')
      read(unit1,100,end=2) card
      if(echo) write(unit2,101) card
      image=dblank(card,80,length,quote)
      if(length.eq.0) go to 2
      ieq=index(card,'=')
      return
2     card(1:4)=eof
      write(unit2,103) eof
      length=3
      image=.true.
      quote=0
      ieq=0
      return
c
      entry bread(rary,lf,lmax,new)
      gread=.false.
      if(ieq.ne.0.or.quote.ne.0.or..not.image) return
      int=2
      go to 8
c
      entry iread(name,iary,lf,lmax,new)
      int=1
      go to 7
c
      entry fread(name,rary,lf,lmax,new)
      int=2
      go to 7
c
      entry dread(name,dary,lf,lmax,new)
      int=3
c
7     l=len(name)
      gread=.false.
      if(ieq-1.ne.l.or..not.image) return
      if(name(:l).ne.card(:l)) return
      if(quote.ne.0) then
         write(unit2,104)
         image=.false.
         return
      endif
8     new=.true.
      gread=.true.
      call getdp(d,card,ieq,length,.false.)
      do 6 i=1,lmax
      if(int.eq.2) then
         d=rary(i)
      else if(int.eq.1) then
         d=iary(i)
      else if(int.eq.3) then
         d=dary(i)
      endif
      call getdp(d,card,ieq,length,image)
      if(int.eq.2) then
         rary(i)=d
      else if(int.eq.1) then
         iary(i)=d
      else if(int.eq.3) then
         dary(i)=d
      endif
      lf=i
      if(length.eq.0) write(unit2,104)
      if(.not.image) return
6     continue
      lf=lmax
      image=.false.
      return
c
      entry cread(name)
      l=len(name)
      gread=.false.
      if(ieq.ne.0.or.length.ne.l.or..not.image) return
      image=name(:l).ne.card(:l)
      gread=.not.image
      if(.not.(gread.and.quote.ne.0)) return
      write(unit2,104)
      gread=.false.
      return
c
      entry sread(name,sary,new)
      l=len(name)
      gread=.false.
      if(ieq-1.ne.l.or..not.image) return
      if(name(:l).ne.card(:l)) return
      image=.false.
      if(quote.ne.ieq+1) then
         write(unit2,104)
         return
      endif
3     new=.true.
      gread=.true.
      sary=card(quote:length)
      return
c
      entry tread(sary,new)
      gread=.false.
      if(ieq.ne.0.or.quote.eq.0.or..not.image) return
      image=.false.
      go to 3
c
100   format(a80)
101   format(1x,a80)
102   format(1x,'??????????')
103   format(1x,a4)
104   format(1x,'data format error.')
      end
c
c
      subroutine getdp(d,name,iprev,imax,notmt)
      logical notmt,nextch
      character*80 name
      character*1 char,digit(10)
      real*8 d,dd,pow
      save nn,dd
      data digit/'0','1','2','3','4','5','6','7','8','9'/,
     x pow/10.d0/
      if(.not.notmt) then
         nn=0
         return
      else if(nn.ge.1) then
         nn=nn-1
         d=dd
         notmt=imax.gt.iprev.or.nn.ge.1
         return
      endif
      if(.not.nextch(name,char,iprev,imax)) go to 120
      if(char.eq.',') return
      nn=0
10    sn=1
      se=1
      ne=0
      id=0
      nd=0
      dd=0.d0
      if(char.eq.'-') then
         sn=-1
         go to 40
      else if(char.eq.'+') then
         go to 40
      endif
20    if(char.eq.'d'.or.char.eq.'e'  .or.
     x   char.eq.'d'.or.char.eq.'e') go to 50
      do 30 i=1,10
      if(char.eq.digit(i)) then
         dd=10.d0*dd+i-1.d0
         nd=nd+id
         go to 40
      endif
30    continue
      if(char.eq.'.') then
         if(id.ne.0) go to 110
         id=1
         go to 40
      else if(char.eq.',') then
         go to 90
      else if(char.eq.'*') then
         if(nn.ne.0.or.id.ne.0.or.int(dd).eq.0) go to 110
         if(.not.nextch(name,char,iprev,imax)) go to 110
         if(char.eq.',') go to 110
         nn=dd
         go to 10
      else
         go to 110
      endif
40    if(nextch(name,char,iprev,imax)) go to 20
      go to 90
50    if(.not.nextch(name,char,iprev,imax)) go to 90
      if(char.eq.'-') then
         se=-1
         go to 80
      else if(char.eq.'+') then
         go to 80
      endif
60    do 70 i=1,10
      if(char.eq.digit(i)) then
         ne=10*ne+i-1
         go to 80
      endif
70    continue
      if(char.eq.',') go to 90
      go to 110
80    if(nextch(name,char,iprev,imax)) go to 60
90    ne=se*ne-nd
      dd=sn*dd
      if(ne.lt.0) pow=.1d0
      do 100 i=1,abs(ne)
      dd=dd*pow
100   continue
      nn=nn-1
      d=dd
      notmt=imax.gt.iprev.or.nn.ge.1
      return
110   imax=0
120   notmt=.false.
      return
      end
c
c
      logical function dblank(name,l,length,nquo)
      character*80 name
      nquo=0
      length=0
      do 20 i=1,l
      if(name(i:i).eq.'''') then
         nquo=length+1
         do 10 j=i+1,l-1
         if(name(j:j).eq.'''') then
              if(name(j+1:j+1).eq.'''') then
                   name(j:l-1)=name(j+1:l)
                   name(l:l)=' '
              else
                   go to 30
              endif
         endif
         length=length+1
         name(length:length)=name(j:j)
10       continue
         if(name(l:l).ne.'''') then
              length=length+1
              name(length:length)=name(l:l)
         endif
         go to 30
      else if(name(i:i).ne.' ') then
         length=length+1
         name(length:length)=name(i:i)
      endif
20    continue
30    dblank=length.ne.0
      return
      end
c
cc
c      integer function ident(name,l,iden,ist)
c      character*80 name
c      character*1 iden
c      do 10 i=ist,l
c      ident=i
c      if(name(i:i).eq.iden) return
c10    continue
c      ident=0
c      return
c      end
cc
c
      logical function nextch(name,char,iprev,imax)
      character*80 name
      character*1 char
      nextch=.false.
      iprev=iprev+1
      if(iprev.gt.imax) return
      char=name(iprev:iprev)
      nextch=.true.
      return
      end

