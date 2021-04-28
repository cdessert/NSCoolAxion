c **** calculate the t & l profile ************************************

      if (debug.ge.1.) print *,'Calculating Initial T profile'
       if (ifteff.ne.15) then
       tsurface=ephi(imax) *1.d11
       tdrip   =ephi(idrip)*3.d11
       tcore   =ephi(0)    *1.d12
       if (tempini.ne.0.) then
        tsurface=0.5 * ephi(imax) *tempini
        tdrip   =0.8 * ephi(idrip)*tempini
        tcore   =1.0 * ephi(0)    *tempini
       else
        tsurface=1.e9
        tdrip   =2.e10
        tcore   =1.e11
       end if
      else                     ! this is for fixed T_b = tb_acc0
       tb_acc0=tb_acc0*ephi(imax)
       tsurface=tb_acc0
       tdrip   =tb_acc0
       tcore   =tb_acc0
      end if

      dec=1.1

      do i=0,icore
       temp(i)=tcore
      end do
      do i=icore+1,idrip
       w1=(dlog10(rad(idrip))-dlog10(rad(  i  )))/
     1    (dlog10(rad(idrip))-dlog10(rad(icore)))
       w2=1.d0-w1
       lt=w1*dlog10(tcore)+w2*dlog10(tdrip)
       temp(i)=10.d0**lt
      end do
      do i=idrip+1,imax
       w1=(dlog10(rad(imax))-dlog10(rad(  i  )))/
     1    (dlog10(rad(imax))-dlog10(rad(idrip)))
       w2=1.d0-w1
       lt=w1*dlog10(tdrip)+w2*dlog10(tsurface)
       temp(i)=10.d0**lt
      end do

      dtemp(0)=0.
      do i=2,imax-1,2
       dtemp(i)=(temp(i+1)-temp(i-1))/(debar(i)+debar(i+1))
      end do

c *********************************************************************
c ***** calculate the inner envelope profile **************************

      if (debug.ge.1.) print *,'Calculating envelope profile'
      do i=imax-1,ienv+1,-2
       x=(dlog(rrho(i+1))-dlog(rrho(i)))/
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       y=(dlog(rrho(i))-dlog(rrho(i-1)))/    
     1   (dlog(rrho(i+1))-dlog(rrho(i-1)))
       ltemp=y*dlog(temp(i+1))+x*dlog(temp(i-1))
c       ltemp=y*log(temp(i+1))+x*log(temp(i-1))
       temp(i)=dexp(ltemp)
      end do

      do i=ienv,imax
       if (temp(i).lt.tcon) then      
        call density(temp(i)/ephi(i),pres(i),a_ion(i),z_ion(i),rrho(i))
        rrho(i)=min(rrho(i),rhod(i))
        bar(i)=6.022d-16*rrho(i)
        dr=debar(i)/rrho(i)*factor
        rad(i+1)=rad(i)+dr
       end if
      end do

c *********************************************************************

      if (debug.ge.1.) print *,'Calculating initial L profile'
      do i=imax,1,-2
       call conduct(i,temp(i)/ephi(i),rrho(i),
     1              a_cell(i),a_ion(i),z_ion(i),qimp,
     2              nbfield2(i),
     3              sig,lambda(i),debug,
     4              nu_e_s,nu_e_l)
       call opacity(temp(i)/ephi(i),rrho(i),a_cell(i),z_ion(i),kappa(i))
       acd=7.56d-15*c/(3.d0*kappa(i)*rrho(i))
       fp(i)=(lambda(i)+4.d0*acd*(temp(i)/ephi(i))**3)*bar(i)/lsol
      end do

      lum(0)=0.
      do i=2,imax-1,2
       lum(i)=-(fp(i+1)+fp(i-1))/2.d0*a2ephin(i)*dtemp(i)
       if (lum(i).eq.0.0) lum(i)=1.d-3
      end do

      do i=1,imax-2,2
       dlum(i)=(nlum(i+1)-nlum(i-1))/(debar(i)+debar(i+1))
      end do
c look +++++++++++++++++++++++++++++++++++++++++++++
      dlum(imax)=0.
c      dlum(imax)=dlum(imax-2)
c ++++++++++++++++++++++++++++++++++++++++++++++++++
     
      do i=1,imax
       otemp(i)=temp(i)
       olum(i-1)=lum(i-1)
      end do

c *********************************************************************

      do i=0,imax
       orad(i)=rad(i)
       orrho(i)=rrho(i)
       rrho1(i)=rrho(i)
       obar(i)=bar(i)
       bar1(i)=bar(i)
      end do

c *********************************************************************
