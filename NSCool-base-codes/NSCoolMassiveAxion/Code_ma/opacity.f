      subroutine opacity(T,Rho,A,Z,kappa)

c ******************************************************************c
c     temperature in K  density in gm/cm3  op given in cm2/gm       c
c*******************************************************************c
c     Uses Kramer opacity only, from Shapiro & Teukolsky 
c ******************************************************************c
       implicit real*8(a-h,k-z)
       INCLUDE 'control_con.inc.f'
        if (iopacity.eq.0) then
         kappa=1.d200
        else
         if (Rho.ge.1.e14) then
          A=100.
          Z=32.
         end if
         kappa = 0.645d23*Z**3/A**2 * Rho/T**3.5
        end if
       return
      end
c ******************************************************************c
c ******************************************************************c
      subroutine opacity_old(temperature,density,ia,iz,op)

c ******************************************************************c
c     temperature in k  density in gm/cm3  op given in cm2/gm       c
c*******************************************************************c
c              checked on august 28 1991                            c
c ******************************************************************c

      implicit real*8(a-h,k-z)
      INCLUDE 'files.inc.f'
      INCLUDE 'control_con.inc.f'
      parameter (itemp=60,jrho=75,ev=11604.)
      dimension temp(0:itemp),rho(0:jrho),kappa(0:itemp,0:jrho)
      save read,temp,rho,kappa

      if (iopacity.eq.0) then
       op=1.d200
       return
      end if

c ****** read opacity table *****************************************

      if (read.eq.1.) goto 1234
      open(unit=40,file=f_opacity,status='old')
       do junk=1,6
        read(40,*)
       end do
       do i=0,itemp
        read(40,*)
        do j=0,jrho
         read(40,*)temp(i),rho(j),kappa(i,j)
        end do
       end do
      close(unit=40,status='keep')

      read=1.

1234  continue

c *******************************************************************

      t=temperature/ev
      lt=log10(t)
      lr=log10(density)

      if (lt.le.0.01) then
       it=0
      else if (lt.ge.5.99) then
       it=itemp-1
      else
       it=10*lt
      end if

      if (lr.le.-4.99) then 
       jr=0
      else if (lr.ge.9.99) then
       jr=jrho-1
      else
       jr=5*(5+lr)
      end if

      deltemp=log10(temp(it+1)/temp(it))
      wi1=log10(temp(it+1)/t)/deltemp
      wi2=1.-wi1
      delrho=log10(rho(jr+1)/rho(jr))
      wj1=log10(rho(jr+1)/density)/delrho
      wj2=1.-wj1

      lop=wi1*wj1*log10(kappa(it,jr))+
     1    wi1*wj2*log10(kappa(it,jr+1))+
     2    wi2*wj1*log10(kappa(it+1,jr))+
     3    wi2*wj2*log10(kappa(it+1,jr+1))

      lopmax=20.
      lop=min(lop,lopmax)

      op=10.**lop

      return

      end

