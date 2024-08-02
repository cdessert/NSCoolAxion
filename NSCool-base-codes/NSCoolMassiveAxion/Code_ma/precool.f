c *********************************************************************
c *********************************************************************
      subroutine grid(idec,rhocore,rhodrip,rhoenv,rhosurf,Model)
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
c       INCLUDE 'profile_comp.inc.f'
       character*79 Model
       dimension rad_t(10000),bar_t(10000),rho_t(10000),pres_t(10000),
     1            emas_t(10000),phi_t(10000)
       parameter (pi=3.1415926535d0)
c Define zone indices: icore, idrip & isurf ***
       icore=2*((icore-1)/2)+1           ! Makes sure icore is odd
       idel1=int(Log10(rhocore/rhodrip)*float(idec))
       idel1=2*(idel1/2)                 ! Makes sure idel1 is even
       idrip=icore+idel1
       idel2=int(Log10(rhodrip/rhosurf)*float(idec))
       idel2=2*(idel2/2)                 ! Makes sure idel2 is even
       isurf=idrip+idel2
c       print *,icore,idrip,isurf
c       read(5,*)
c Read Star structure file: ********************************
        open(unit=20,file=f_profile,status='old')
         read(20,*)jtext,jmax
         jmax=jmax+1         ! Line indices in file star at j=0 instead of 1
         read(20,*)
         read(20,*)Model
         do j=3,jtext
          read(20,*)
         end do
         jdrip=0
         jcore=0
         do j=1,jmax
          read(20,*)jj,
     1     rad_t(j),bar_t(j),rho_t(j),pres_t(j),emas_t(j),phi_t(j)
           rad_t(j)=rad_t(j)*100.d0
           if ((rho_t(j).lt.rhocore).and.(jcore.eq.0)) then
            jcore=j-1
           end if
         end do
        close(unit=20,status='keep')
c Get the core radius exactly:
        drho=rho_t(jcore)-rho_t(jcore+1)
        w1=(rhocore-rho_t(jcore+1))/drho
        w2=1.d0-w1
        rad_core=w1*rad_t(jcore)+w2*rad_t(jcore+1)
c
c Define Star grid: ****************************************
c
c         dvol = physical volume between rad(i-1) and rad(i)
c
c CORE: zoning with (approximate) constant volume: *********
        do i=0,icore
         rad(i)=(float(i)/float(icore))**(1./3.) * rad_core
        end do
c        debar(0)=0.d0
c        bar(0)=bar_t(1)                         ! bar is better defined in get_*_chemistry, from the EOS
        rrho(0)=rho_t(1)
        emas(0)=0.d0
        phi(0)=phi_t(1)
        pres(0)=pres_t(1)
        dvol(0)=0.d0
        j=0
        do i=1,icore
 500     j=j+1
         if (rad_t(j).lt.rad(i)) goto 500
         delrad=rad_t(j)-rad_t(j-1)
         w1=(rad_t(j)-rad(i))/delrad
         w2=1.d0-w1
c         bar(i) =w1*bar_t(j-1) +w2*bar_t(j)     ! bar is better defined in get_*_chemistry, from the EOS
         rrho(i)=w1*rho_t(j-1) +w2*rho_t(j)
         emas(i)=w1*emas_t(j-1)+w2*emas_t(j)
         phi(i) =w1*phi_t(j-1) +w2*phi_t(j)
         pres(i)=w1*pres_t(j-1)+w2*pres_t(j)
         if (i.eq.1) then
          dvol(i)=4./3.*pi*rad(i)**3
         else
          dvol(i)=4.*pi*((rad(i-1)+rad(i))/2.)**2*
     x                   (rad(i)-rad(i-1)) /
     x            dsqrt(1.-2.92d5*(emas(i-1)+emas(i))/(rad(i-1)+rad(i)))
         end if
c         debar(i)=(bar(i-1)+bar(i))/2. * dvol(i)
         j=j-1
        end do
c ***   Just to make sure (exact accuracy, not from above interpolations)
        rrho(icore)=rhocore
c Crust: zoning with "idec" zone per decade in density: ****
        dlogrho=log10(rhocore/rhodrip)
        dlrho=dlogrho/float(idrip-icore)
        do i=icore+1,idrip
         lrho=log10(rhocore)-float(i-icore)*dlrho
         rrho(i)=10.d0**lrho
        end do
c ***
        dlogrho=log10(rhodrip/rhosurf)
        dlrho=dlogrho/float(isurf-idrip)
        i=idrip
        do i=idrip+1,isurf
         lrho=log10(rhodrip)-float(i-idrip)*dlrho
         rrho(i)=10.d0**lrho
        end do
        j=0
        do i=icore+1,isurf
 600     j=j+1
         if (rho_t(j).gt.rrho(i)) goto 600
         dellrho=log10(rho_t(j-1))-log10(rho_t(j))
         w1=(log10(rrho(i))-log10(rho_t(j)))/dellrho
         w2=1.d0-w1
c         bar(i) =w1*bar_t(j-1) +w2*bar_t(j)   ! bar is better defined in get_*_chemistry, from the EOS
         rad(i) =w1*rad_t(j-1) +w2*rad_t(j)
         emas(i)=w1*emas_t(j-1)+w2*emas_t(j)
         phi(i) =w1*phi_t(j-1) +w2*phi_t(j)
         pres(i)=w1*pres_t(j-1)+w2*pres_t(j)
         dvol(i)=4.*pi*((rad(i-1)+rad(i))/2.)**2*
     x                  (rad(i)-rad(i-1)) /
     x            dsqrt(1.-2.92d5*(emas(i-1)+emas(i))/(rad(i-1)+rad(i)))
c         debar(i)=(bar(i-1)+bar(i))/2. * dvol(i)
         j=j-1
        end do
        dvol(isurf+1)=dvol(isurf)
c ***   Just to make sure (exact accuracy, not from above interpolations)
        rrho(idrip)=rhodrip
        rrho(isurf)=rhosurf
c Find the envelope boundary: ******************************
        ienv=isurf+2
        do i=isurf,idrip,-2
         if (rrho(i).lt.rhoenv) ienv=i
        end do
c ***
        imax=isurf
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine grid_strange(idec,rhocore,rhodrip,rhoenv,rhosurf,Model)
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
c HERE DANY: looks like it's not used:
c       INCLUDE 'profile_comp.inc.f'
       character*79 Model
       dimension rad_t(10000),bar_t(10000),rho_t(10000),pres_t(10000),
     1            emas_t(10000),phi_t(10000)
       parameter (pi=3.1415926535d0)
c Read Star structure file: ********************************
        open(unit=20,file=f_profile,status='old')
         read(20,*)jtext,jmax
         jmax=jmax+1         ! Line indices in file star at j=0 instead of 1
         read(20,*)
         read(20,*)Model
         do j=3,jtext
          read(20,*)
         end do
         jdrip=0
         jcore=0
         do j=1,jmax
          read(20,*)jj,
     1     rad_t(j),bar_t(j),rho_t(j),pres_t(j),emas_t(j),phi_t(j)
           rad_t(j)=rad_t(j)*100.d0
           if ((rho_t(j).lt.rhocore).and.(jcore.eq.0)) then
            jcore=j-1
           end if
         end do
        close(unit=20,status='keep')
c Get the core radius exactly:
c        drho=rho_t(jcore)-rho_t(jcore+1)
c        w1=(rhocore-rho_t(jcore+1))/drho
c        w2=1.d0-w1
c        rad_core=w1*rad_t(jcore)+w2*rad_t(jcore+1)
       rad_core=rad_t(jcore)
       rhocore=rho_t(jcore-1)
       rhodrip=rho_t(jcore+1)
c HERE DANY ccccccccccccccccccccc
c       print '(a40,1p4e12.3)','rhocore, rhodrip, rhoenv, rhosurf =',
c     x       rhocore,rhodrip,rhoenv,rhosurf
ccccccccccccccccccccccccccccccccc
c Define zone indices: icore, idrip & isurf ***
       icore=2*((icore-1)/2)+1           ! Makes sure icore is odd
c       idel1=int(Log10(rhocore/rhodrip)*float(idec))
c       idel1=2*(idel1/2)                 ! Makes sure idel1 is even
c       idrip=icore+idel1
       idrip=icore+2
       idel2=int(Log10(rhodrip/rhosurf)*float(idec))
       idel2=2*(idel2/2)                 ! Makes sure idel2 is even
       isurf=idrip+idel2
c       print *,icore,idrip,isurf
c       read(5,*)
c
c Define Star grid: ****************************************
c
c         dvol = physical volume between rad(i-1) and rad(i)
c
c CORE: zoning with (approximate) constant volume: *********
        do i=0,icore
         rad(i)=(float(i)/float(icore))**(1./3.) * rad_core
        end do
c        debar(0)=0.d0
c        bar(0)=bar_t(1)                         ! bar is better defined in get_*_chemistry, from the EOS
        rrho(0)=rho_t(1)
        emas(0)=0.d0
        phi(0)=phi_t(1)
        pres(0)=pres_t(1)
        dvol(0)=0.d0
        j=0
        do i=1,icore
 500     j=j+1
         if (rad_t(j).lt.rad(i)) goto 500
         delrad=rad_t(j)-rad_t(j-1)
         w1=(rad_t(j)-rad(i))/delrad
         w2=1.d0-w1
c         bar(i) =w1*bar_t(j-1) +w2*bar_t(j)     ! bar is better defined in get_*_chemistry, from the EOS
         rrho(i)=w1*rho_t(j-1) +w2*rho_t(j)
         emas(i)=w1*emas_t(j-1)+w2*emas_t(j)
         phi(i) =w1*phi_t(j-1) +w2*phi_t(j)
         pres(i)=w1*pres_t(j-1)+w2*pres_t(j)
         if (i.eq.1) then
          dvol(i)=4./3.*pi*rad(i)**3
         else
          dvol(i)=4.*pi*((rad(i-1)+rad(i))/2.)**2*
     x                   (rad(i)-rad(i-1)) /
     x            dsqrt(1.-2.92d5*(emas(i-1)+emas(i))/(rad(i-1)+rad(i)))
         end if
c         debar(i)=(bar(i-1)+bar(i))/2. * dvol(i)
         j=j-1
        end do
c ***   Just to make sure (exact accuracy, not from above interpolations)
        rrho(icore)=rhocore
c Crust: zoning with "idec" zone per decade in density: ****
c        dlogrho=log10(rhocore/rhodrip)
c        dlrho=dlogrho/float(idrip-icore)
c        do i=icore+1,idrip
c         lrho=log10(rhocore)-float(i-icore)*dlrho
c         rrho(i)=10.d0**lrho
c        end do
c ***
c       This is temporary: redfined below, but necesary for the time being:
        rrho(idrip-1)=1.1*rhodrip
c ****
        rrho(idrip)=rhodrip
        dlogrho=log10(rhodrip/rhosurf)
        dlrho=dlogrho/float(isurf-idrip)
        i=idrip
        do i=idrip+1,isurf
         lrho=log10(rhodrip)-float(i-idrip)*dlrho
         rrho(i)=10.d0**lrho
        end do
c        j=jcore-1
        j=0
        do i=icore+1,isurf
 600     j=j+1
         if (rho_t(j).gt.rrho(i)) goto 600
         dellrho=log10(rho_t(j-1))-log10(rho_t(j))
         w1=(log10(rrho(i))-log10(rho_t(j)))/dellrho
         w2=1.d0-w1
c         bar(i) =w1*bar_t(j-1) +w2*bar_t(j)   ! bar is better defined in get_*_chemistry, from the EOS
         rad(i) =w1*rad_t(j-1) +w2*rad_t(j)
         emas(i)=w1*emas_t(j-1)+w2*emas_t(j)
         phi(i) =w1*phi_t(j-1) +w2*phi_t(j)
         pres(i)=w1*pres_t(j-1)+w2*pres_t(j)
         dvol(i)=4.*pi*((rad(i-1)+rad(i))/2.)**2*
     x                  (rad(i)-rad(i-1)) /
     x            dsqrt(1.-2.92d5*(emas(i-1)+emas(i))/(rad(i-1)+rad(i)))
c         debar(i)=(bar(i-1)+bar(i))/2. * dvol(i)
         j=j-1
        end do
        dvol(isurf+1)=dvol(isurf)
c ***   Just to make sure (exact accuracy, not from above interpolations)
        rrho(idrip)=rhodrip
        rrho(isurf)=rhosurf
c ***   Define stuff at icore+1 = idrip-1
        rrho(idrip-1)=rhodrip
        rad(idrip-1) =0.5*(rad(icore)+rad(idrip))
        emas(idrip-1)=0.5*(emas(icore)+emas(idrip))
        phi(idrip-1) =0.5*(phi(icore)+phi(idrip))
        pres(idrip-1)=0.5*(pres(icore)+pres(idrip))
        dvol(idrip-1)=dvol(idrip)
c Find the envelope boundary: ******************************
        ienv=isurf+2
        do i=isurf,idrip,-2
         if (rrho(i).lt.rhoenv) ienv=i
        end do
c ***
        imax=isurf
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        print '(a40,1p4e12.3)','rhocore,rhodrip,rhoenv,rhosurf =',
c     x        rhocore,rhodrip,rhoenv,rhosurf
c        print '(a40,4i12)','icore, idrip, ienv, imax =',
c     x       icore,idrip,ienv,imax
c        do i=0,imax
c         print '(i10,0p1f15.5,1p4e12.3)',
c     x     i,rad(i)/100.,rrho(i),pres(i),emas(i),dvol(i)
c        end do
c        print *,'that''s it !'
c        read(5,*)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine get_core_chemistry
c *********************************************************************
c      This subroutine calculates the concentrations Y's of 
c      all particles in the core using the EOS table.
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
c       INCLUDE 'files_in_out.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'fermi.inc.f'
c       INCLUDE 'control_nu.inc.f'
       INCLUDE 'quark.inc.f'
       dimension rho_t(500),nbar_t(500)
       dimension yneutr_t(500),yprot_t(500),
     2           yelect_t(500),ymuon_t(500),
     3           ylambda_t(500),
     4           ysminus_t(500),yszero_t(500),ysplus_t(500),
     5           yquarku_t(500),yquarkd_t(500),yquarks_t(500),
     6           fhad_t(500)
      dimension mstp_t(0:isize),mstn_t(0:isize),mstla_t(0:isize),
     2          mstsm_t(0:isize),msts0_t(0:isize),mstsp_t(0:isize)
       dimension theta_k_t(500),theta_p_t(500)

c ***** Read EOS file: *************************
        if (version.eq.'old') then
c old:
         open(unit=59,file=f_stareos,status='old')
          read(59,*)it,i0,ix
          do i=1,it
           read(59,*)
          end do
          do i=1,ix
           read(59,*)rho_t(i),x1,x2,yprot_t(i),ymuon_t(i)
           yelect_t(i)=yprot_t(i)-ymuon_t(i)
           yneutr_t(i)=1.d0-yprot_t(i)
           ylambda_t(i)=0.0d0
           ysminus_t(i)=0.0d0
           yszero_t(i)=0.0d0
           ysplus_t(i)=0.0d0
           theta_k_t(i)=0.0d0
           theta_p_t(i)=0.0d0
           yquarku_t(i)=0.0d0
           yquarkd_t(i)=0.0d0
           yquarks_t(i)=0.0d0
          end do
         close(unit=59,status='keep')
        else if (version.eq.'new')then
c new:
         open(unit=59,file=f_stareos,status='old')
          read(59,*)it,i0,ix
          do i=1,it
           read(59,*)
          end do
          do i=1,ix
c This ynut is for formats which list fhad !
c           read(59,*)rho_t(i),xnut,nbar_t(i),ynut,
           read(59,*)rho_t(i),xnut,nbar_t(i),
     2         yelect_t(i),ymuon_t(i),yneutr_t(i),yprot_t(i),
     3         ylambda_t(i),ysminus_t(i),yszero_t(i),ysplus_t(i)
           theta_k_t(i)=0.0d0
           theta_p_t(i)=0.0d0
           yquarku_t(i)=0.0d0
           yquarkd_t(i)=0.0d0
           yquarks_t(i)=0.0d0
          end do
        else if (version.eq.'NEW')then
c NEW:
         open(unit=59,file=f_stareos,status='old')
          read(59,*)it,i0,ix
          do i=1,it
           read(59,*)
          end do
          do i=1,ix
           read(59,*)rho_t(i),xnut,nbar_t(i),
     2         yelect_t(i),ymuon_t(i),yneutr_t(i),yprot_t(i),
     3         ylambda_t(i),ysminus_t(i),yszero_t(i),ysplus_t(i),
     4         mstp_t(i),mstn_t(i),mstla_t(i),
     5         mstsm_t(i),msts0_t(i),mstsp_t(i)
           theta_k_t(i)=0.0d0
           theta_p_t(i)=0.0d0
           yquarku_t(i)=0.0d0
           yquarkd_t(i)=0.0d0
           yquarks_t(i)=0.0d0
          end do
         close(unit=59,status='keep')
        else if (version.eq.'QRK')then
c QRK:
         open(unit=59,file=f_stareos,status='old')
          read(59,*)it,i0,ix,alpha_c,strange_mass
          do i=1,it
           read(59,*)
          end do
          do i=1,ix
           read(59,*)rho_t(i),
     1         xnut,nbar_t(i),fhad_t(i),
     2         yelect_t(i),ymuon_t(i),yneutr_t(i),yprot_t(i),
     3         ylambda_t(i),ysminus_t(i),yszero_t(i),ysplus_t(i),
     4         yquarku_t(i),yquarkd_t(i),yquarks_t(i)
           theta_k_t(i)=0.0d0
           theta_p_t(i)=0.0d0
          end do
         close(unit=59,status='keep')
        else
         pause 'Cannot read the f_stareos file'
        end if
c ***** Interpolate the particle concentrations: ********
        i1=1
        do i0=0,icore
         if(rrho(i0).ge.rho_t(1))then
          i1=1
          i2=2
         else if(rrho(i0).le.rho_t(ix))then
          i1=ix-1
          i2=ix
         else
          i=i1-1
5768      i=i+1
          if((rrho(i0).ge.rho_t(i+1)).and.(rrho(i0).le.rho_t(i)))then
           i1=i
           i2=i+1
          else
           goto 5768
          end if
         end if
         x1=(dlog(rho_t(i2))-dlog(rrho(i0)))/
     1      (dlog(rho_t(i2))-dlog(rho_t(i1)))
         x2=(dlog(rrho(i0))-dlog(rho_t(i1)))/
     1      (dlog(rho_t(i2))-dlog(rho_t(i1)))
c         print '(1p3e12.3,5x,0p2f10.5,2i5)',
c     1       rho_t(i1),rrho(i0),rho_t(i2),x1,x2,i1,i2
         bar(i0)    =x1*nbar_t(i1)   +x2*nbar_t(i2)
         yelect(i0) =x1*yelect_t(i1) +x2*yelect_t(i2)
         ymuon(i0)  =x1*ymuon_t(i1)  +x2*ymuon_t(i2)
         yneutr(i0) =x1*yneutr_t(i1) +x2*yneutr_t(i2)
         yprot(i0)  =x1*yprot_t(i1)  +x2*yprot_t(i2)
         ylambda(i0)=x1*ylambda_t(i1)+x2*ylambda_t(i2)
         ysminus(i0)=x1*ysminus_t(i1)+x2*ysminus_t(i2)
         yszero(i0) =x1*yszero_t(i1) +x2*yszero_t(i2)
         ysplus(i0) =x1*ysplus_t(i1) +x2*ysplus_t(i2)
         yquarku(i0)=x1*yquarku_t(i1)+x2*yquarku_t(i2)
         yquarkd(i0)=x1*yquarkd_t(i1)+x2*yquarkd_t(i2)
         yquarks(i0)=x1*yquarks_t(i1)+x2*yquarks_t(i2)
         theta_k(i0)=x1*theta_k_t(i1)+x2*theta_k_t(i2)
         theta_p(i0)=x1*theta_p_t(i1)+x2*theta_p_t(i2)
         if ((version.eq.'old').or.(version.eq.'new').or.
     1        (version.eq.'NEW')) then
          fhad(i0)=1.d0
         else if (version.eq.'QRK') then
          fhad(i0)=   x1*fhad_t(i1)   +x2*fhad_t(i2)
          if (fhad(i0).gt.1.0) fhad(i0)=1.d0
          if (fhad(i0).lt.0.0) fhad(i0)=0.d0
         end if
c NOTE: the Y's are defined such that 
c       Y_i*bar=number of particle of type i per fm^3
c       but the density of baryons in the baryon phase is then
c       Y_i*bar/fhad 
c       since they occupy only a fraction fhad of the volume
c       and for quarks it is
c       Y_q*bar/(1-fhad)
c *** Check for consistency:
         bnuc=yneutr(i0)+yprot(i0)
         bhyp=ylambda(i0)+ysminus(i0)+yszero(i0)+ysplus(i0) 
         bqua=1./3.*(yquarku(i0)+yquarkd(i0)+yquarks(i0))
         btot=bnuc+bhyp+bqua
         qlep=-yelect(i0)-ymuon(i0)
         qnuc=yprot(i0)
         qhyp=ysplus(i0)-ysminus(i0)
         qqua=1./3.*(2.*yquarku(i0)-yquarkd(i0)-yquarks(i0))
         qtot=qlep+qnuc+qhyp+qqua
         if (abs(btot-1.0).gt.1.e-2) then
          print '(a30,i5,1p1e12.3,5x,0p2f10.5)',
     1          'i, rho, Btot, Qtot = ',
     2           i0,rrho(i0),btot,qtot
          pause 'Btot not equal to 1 !'
         end if
         if (abs(qtot).gt.1.e-2) then
          print '(a30,i5,1p1e12.3,5x,0p2f10.5)',
     1          'i, rho, Btot, Qtot = ',
     2           i0,rrho(i0),btot,qtot
          pause 'Qtot not equal to 0 !'
         end if
c Get the baryon effective masses if NEW:
         if (version.eq.'NEW') then
          mstp(i0) =x1*mstp_t(i1) +x2*mstp_t(i2)
          mstn(i0) =x1*mstn_t(i1) +x2*mstn_t(i2)
          mstla(i0)=x1*mstla_t(i1)+x2*mstla_t(i2)
          mstsm(i0)=x1*mstsm_t(i1)+x2*mstsm_t(i2)
          msts0(i0)=x1*msts0_t(i1)+x2*msts0_t(i2)
          mstsp(i0)=x1*mstsp_t(i1)+x2*mstsp_t(i2)
         end if
        end do
c *** Clean up Y's in the crust, just in case:
c (yelect & yneutr will be calculated in "get_crust_chemistry")
        do i0=icore+1,imax
         yelect(i0) =0.d0
         ymuon(i0)  =0.d0
         yneutr(i0) =0.d0
         yprot(i0)  =0.d0
         ylambda(i0)=0.d0
         ysminus(i0)=0.d0
         yszero(i0) =0.d0
         ysplus(i0) =0.d0
         yquarku(i0)=0.d0
         yquarkd(i0)=0.d0
         yquarks(i0)=0.d0
         theta_k(i0)=0.d0
         theta_p(i0)=0.d0
         fhad(i0)=1.d0              ! Also just in case
        end do
c *****
       return
      end
c *************************************************************************
c *************************************************************************
      subroutine get_core_chemistry_strange(debug)
c
c     Specifically for pure quark matter cores !
c
c *********************************************************************
c      This subroutine calculates the concentrations Y's of 
c      all particles in the core using the EOS table.
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
c       INCLUDE 'files_in_out.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'fermi.inc.f'
c       INCLUDE 'control_nu.inc.f'
       INCLUDE 'quark.inc.f'
       dimension rho_t(500),nbar_t(500),pres_t(500)
       dimension yelect_t(500),
     1           yquarku_t(500),yquarkd_t(500),yquarks_t(500)
c *****
        if (debug.ge.1.) then
         print *,'Entering subroutine get_core_chemistry_strange'
        end if
c ***** Read EOS file: *************************
         open(unit=59,file=f_stareos,status='old')
          read(59,*)alpha_c,b_bag,strange_mass
          it=5
          do i=1,it
           read(59,*)
          end do
          do i=1,1000
           read(59,*,end=999,err=999)
     1         nbar_t(i),rho_t(i),pres_t(i),
     2         yelect_t(i),yquarku_t(i),yquarkd_t(i),yquarks_t(i)
c          Adjust, since file actually gives ne,nu,nd & ns:
           yelect_t(i) =yelect_t(i) /nbar_t(i)
           yquarku_t(i)=yquarku_t(i)/nbar_t(i)
           yquarkd_t(i)=yquarkd_t(i)/nbar_t(i)
           yquarks_t(i)=yquarks_t(i)/nbar_t(i)
          end do
 999     ix=i-1
         close(unit=59,status='keep')
c ***** Interpolate the particle concentrations: ********
        i1=1
        do i0=0,icore
         if(pres(i0).ge.pres_t(1))then
          i1=1
          i2=2
         else if(pres(i0).le.pres_t(ix))then
          i1=ix-1
          i2=ix
         else
          i=i1-1
5768      i=i+1
          if((pres(i0).ge.pres_t(i+1)).and.
     1       (pres(i0).le.pres_t(i)))then
           i1=i
           i2=i+1
          else
           goto 5768
          end if
         end if
         x1=(dlog(pres_t(i2))-dlog(pres(i0)))/
     1      (dlog(pres_t(i2))-dlog(pres_t(i1)))
         x2=(dlog(pres(i0))-dlog(pres_t(i1)))/
     1      (dlog(pres_t(i2))-dlog(pres_t(i1)))
         yelect(i0) =x1*yelect_t(i1) +x2*yelect_t(i2)
         yquarku(i0)=x1*yquarku_t(i1)+x2*yquarku_t(i2)
         yquarkd(i0)=x1*yquarkd_t(i1)+x2*yquarkd_t(i2)
         yquarks(i0)=x1*yquarks_t(i1)+x2*yquarks_t(i2)
         bar(i0)    =x1*    nbar_t(i1)+x2*    nbar_t(i2)
c         print *
c         print '(i5,1p3e28.18)',i0,pres_t(i1),pres(i0),pres_t(i2)
c         print '(i5,0p3f28.18)',i0,nbar_t(i1),nnb,nbar_t(i2)
c *** Check for consistency:
         bqua=1./3.*(yquarku(i0)+yquarkd(i0)+yquarks(i0))
         qlep=-yelect(i0)
         qqua=1./3.*(2.*yquarku(i0)-yquarkd(i0)-yquarks(i0))
         qtot=qlep+qqua
         if (abs(bqua-1.0).gt.1.e-2) then
          print '(a30,i5,1p1e12.3,5x,0p2f10.5)',
     1          'i, rho, Btot, Qtot = ',
     2           i0,rrho(i0),btot,qtot
          pause 'Bqua not equal to 1 !'
         end if
         if (abs(qtot).gt.1.e-2) then
          print '(a30,i5,1p1e12.3,5x,0p2f10.5)',
     1          'i, rho, Btot, Qtot = ',
     2           i0,rrho(i0),btot,qtot
          pause 'Qtot not equal to 0 !'
         end if
        end do
c *************************************
c Clean up, just in case:
        do i=0,icore
         yneutr(i)=0.d0
         yprot(i)=0.d0
         ylambda(i)=0.d0
         ysminus(i)=0.d0
         yszero(i)=0.d0
         ysplus(i)=0.d0
         fhad(i)=0.d0
         ymuon(i)=0.d0
         theta_k(i)=0.d0
         theta_p(i)=0.d0
        end do
        do i=icore+1,imax
         yneutr(i)=0.d0
         yprot(i)=0.d0
         ylambda(i)=0.d0
         ysminus(i)=0.d0
         yszero(i)=0.d0
         ysplus(i)=0.d0
         fhad(i)=1.d0
         ymuon(i)=0.d0
         theta_k(i)=0.d0
         theta_p(i)=0.d0
         yquarku(i)=0.d0
         yquarkd(i)=0.d0
         yquarks(i)=0.d0
        end do
c *****
        if (debug.ge.1.) then
         print *,'Exiting subroutine get_core_chemistry_strange'
        end if
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_crust_chemistry(debug)
c *********************************************************************
c      This subroutine calculates the chemical composition, A, A' & Z
c      and the e & n Y's in the crust using the crust_cc table.
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       dimension rho_t(500),pres_t(500),bar_t(500),
     1           A_cell_t(500),A_ion_t(500),Z_ion_t(500)
c HERE DANY ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this common block is , apparently, not used !
c       common/temporary/rho_t,a_cell_t,a_ion_t,z_ion_t,jmax
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *****
        if (debug.ge.1.) then
         print *,'Entering subroutine get_crust_chemistry'
        end if
c *****
        open (unit=20,file=f_crusteos,status='old')
         read(20,*)jtext,jmax
         do j=1,jtext
          read(20,*)
         end do
         jget_drip=0
         do j=jmax,1,-1
          read(20,*)rho_t(j),pres_t(j),bar_t(j),
     1              A_cell_t(j),A_ion_t(j),Z_ion_t(j)
          if (jget_drip.eq.0) then
           if (A_cell_t(j).ne.A_ion_t(j)) then
            jdrip=j
            jget_drip=1
           end if
          end if
         end do
        close(unit=20,status='keep')
c **************************************************************
c Make sure that rho_t and bar_t at jmax are smaller than in core:
        jjmax=jmax
        do j=jmax,1,-1
         if (rho_t(j).ge.rrho(icore)) then
          jjmax=j-1
         end if
        end do
        jmax=jjmax
c **************************************************************
c This is for interpolation from the last crust_EOS line
c up to the core, in a, hopefully, consistent way:
c ***************
c First version:
c - Take rho & bar from core EOS:
        jmax=jmax+1
        rho_t(jmax) =rrho(icore)
        bar_t(jmax) =bar(icore)
        pres_t(jmax)=pres(icore)
c - Keep Z constant:
        Z_ion_t(jmax)=Z_ion_t(jmax-1)
c - Take A_ion constant: results in that the driped neutron
c   fraction grows and kfn only has a small discontinuity
c   (this is the same recipe as below for lower densities)
        A_ion_t(jmax) =A_ion_t(jmax-1)
c - Take Z/A_cell consistent with Ye in core_EOS:
        A_cell_t(jmax)=Z_ion_t(jmax) / yelect(icore)
c ***************
c Second version:
c
c To be done !
c
c **************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do j=jmax,1,-1
c         print '(i5,1p3e12.3,5x,0p3f12.1)',
c     x          j,rho_t(j),pres_t(j),bar_t(j),
c     x          A_cell_t(j),A_ion_t(j),Z_ion_t(j)
c        end do
c        read(5,*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        j=jmax
        do i=icore+1,imax
 100     j=j-1
         if (rrho(i).le.rho_t(j)) goto 100
         dd=(rho_t(j+1)-rho_t(j))
         w2=(rrho(i)-rho_t(j))/dd
         w1=1.d0-w2
         if (A_ion_t(j+1).eq.A_cell_t(j+1)) then
          A_ion(i) =A_ion_t(j+1)
          A_cell(i)=A_cell_t(j+1)
          Z_ion(i) =Z_ion_t(j+1)
         else
c This is the old version:
c          A_cell(i)=w1*A_cell_t(j)+w2*A_cell_t(j+1)
c          A_ion(i) =w1*A_ion_t(j) +w2*A_ion_t(j+1)
c          Z_ion(i) =Z_ion_t(j+1)
c New version: A & Z constant, only A_cell changes:
          A_cell(i)=w1*A_cell_t(j)+w2*A_cell_t(j+1)
          A_ion(i) =A_ion_t(j+1)
          Z_ion(i) =Z_ion_t(j+1)
         end if
c ****  Get the baryon number density:
         bar(i)=w1*bar_t(j)+w2*bar_t(j+1)
c ****  Calculate the fraction of volume occupied by ions:
         r1=1.1d0                                     ! Scale parameter, in fm
         vion=4.d0/3.d0*3.14159d0 * r1**3 * a_ion(i)  ! vion in fm^3
         vion=1.d-39 * vion                           ! in cm^3
         nion=rrho(i)/(1.66d-24*a_ion(i))             ! ion density per cm^3
         v_ion(i)=nion*vion
         v_ion(i)=min(1.d0,v_ion(i))                  ! Make sure it's < 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       j+1,rho_t(j+1),bar_t(j+1),
c     2       a_cell_t(j+1),a_ion_t(j+1),z_ion_t(j+1)
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       i ,rrho( i ),bar(  i ),
c     2        a_cell(  i ),a_ion(  i ),z_ion(  i )
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       j ,rho_t( j ),bar_t( j ),
c     2        a_cell_t( j ),a_ion_t( j ),z_ion_t( j )
c      read(5,*)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         j=j+1
        end do
c$$$c Find the neutron drip point :
c$$$c HERE DANY: idrip is now defined in subroutine grid
c$$$        do i=imax,icore+2,-2
c$$$         if (A_ion(i).eq.A_cell(i)) idrip=i-2
c$$$        end do
c$$$c       Following is just in case strange stars are considered:
c$$$c       and garantees that neutron drip is in the crust !
c$$$        if (idrip.eq.icore) idrip=icore+2
c$$$        rhodrip=rrho(idrip)
c Calculate the Y's of the e & n:
        do i=icore+1,imax
         yelect(i) =Z_ion(i)/A_cell(i)
         yneutr(i) =(A_cell(i)-A_ion(i))/A_cell(i)
        end do
c Clean up the core, just in case:
        do i=0,icore
         Z_ion(i) =0.d0
         A_ion(i) =0.d0
         A_cell(i)=0.d0
         v_ion(i) =0.d0
        end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do i=icore+1,imax,2
c         if (i.eq.idrip) print*,'+++++++++++++'
c         print '(i5,1p3e12.3,3x,0p3f12.1,3x,1p3e12.3)',
c     x         i,rrho(i),pres(i),bar(i),
c     x         A_cell(i),A_ion(i),Z_ion(i),v_ion(i),yelect(i),yneutr(i)
c         if (i.eq.idrip) print*,'+++++++++++++'
c        end do
c        read(5,*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c *********************************************************************
      debar(0)=0.d0
      do i=1,imax
       debar(i)=(bar(i-1)+bar(i))/2. * dvol(i)
      end do
c *********************************************************************

c *****
        if (debug.ge.1.) then
         print *,'Exiting subroutine get_crust_chemistry'
        end if
c *****
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_crust_chemistry_OLD
c *********************************************************************
c      This subroutine calculates the chemical composition, A, A' & Z
c      and the e & n Y's in the crust using the crust_cc table.
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'rho_limits.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       dimension rho_t(500),pres_t(500),bar_t(500),
     1           a_cell_t(500),a_ion_t(500),z_ion_t(500)
       common/temporary/rho_t,a_cell_t,a_ion_t,z_ion_t,jmax
        open (unit=20,file=f_crusteos,status='old')
         read(20,*)jtext,jmax
         do j=1,jtext
          read(20,*)
         end do
         do j=jmax,1,-1
          read(20,*)rho_t(j),pres_t(j),bar_t(j),
     1              a_cell_t(j),a_ion_t(j),z_ion_t(j)
         end do
         jmax=jmax+1
         rho_t(jmax)=2.8d14
         bar_t(jmax)=0.16d0
         a_cell_t(jmax)=a_cell_t(jmax-1)
         a_ion_t(jmax) =a_ion_t(jmax-1)
         z_ion_t(jmax)=z_ion_t(jmax-1)
        close(unit=20,status='keep')

        j=jmax
        do i=icore+1,imax
 100     j=j-1
         if (bar(i).le.bar_t(j)) goto 100
         if (a_ion_t(j+1).eq.a_cell_t(j+1)) then
          a_ion(i) =a_ion_t(j+1)
          a_cell(i)=a_cell_t(j+1)
          z_ion(i) =z_ion_t(j+1)
         else
          db=(bar_t(j+1)-bar_t(j))
          w2=(bar(i)-bar_t(j))/db
          w1=1.d0-w2
          a_ion(i) =w1*a_ion_t(j) +w2*a_ion_t(j+1)
          a_cell(i)=w1*a_cell_t(j)+w2*a_cell_t(j+1)
          z_ion(i) =z_ion_t(j+1)
         end if
c ****  Calculate the fraction of volume occupied by ions:
         r1=1.1        ! Scale parameter, in fm
         vion=4.d0/3.d0*3.14159 * r1**3 * a_ion(i)    ! vion in fm^3
         vion=1.d-39 * vion                           ! in cm^3
         nion=rrho(i)/(1.66d-24*a_ion(i))
         v_ion(i)=nion*vion
         v_ion(i)=min(1.d0,v_ion(i))    ! Make sure it's < 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       j+1,rho_t(j+1),bar_t(j+1),
c     2       a_cell_t(j+1),a_ion_t(j+1),z_ion_t(j+1)
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       i ,rrho( i ),bar(  i ),
c     2        a_cell(  i ),a_ion(  i ),z_ion(  i )
c      print '(i5,1p1e12.3,0p1f17.12,0p3f10.1)',
c     1       j ,rho_t( j ),bar_t( j ),
c     2        a_cell_t( j ),a_ion_t( j ),z_ion_t( j )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         j=j+1
        end do
c Find the neutron drip point :
        do i=imax,icore+1,2
         if (A_ion(i).eq.A_cell(i)) idrip=i-2
        end do
c       Following is just in case strange stars are considered:
c       and garantees that neutron drip is in the crust !
        if (idrip.eq.icore) idrip=icore+2
        rhodrip=rrho(idrip)
c Calculate the Y's of the e & n:
        do i=icore+1,imax
         yelect(i) =z_ion(i)/a_cell(i)
         yneutr(i) =(a_cell(i)-a_ion(i))/a_cell(i)
c         print '(i5,1p1e12.3,0p1f10.5,5x,0p2f12.5)',
c     1      i,rrho(i),bar(i),yelect(i),yneutr(i)
        end do
c Clean up the core, just in case:
        do i=0,icore
         z_ion(i) =0.d0
         a_ion(i) =0.d0
         a_cell(i)=0.d0
         v_ion(i) =0.d0
        end do
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_fermi_momenta
c *********************************************************************
c      This subroutine calculates the Fermi momentum of all fermions
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'control_nu.inc.f'
       parameter(pi=3.14159265d0)
c ****  Calculate the fermi momenta in the core:
        do j=0,icore
c----------------------------------------------
         if (fhad(j).ne.0.d0) then
          nn =max(0.e0,yneutr(j) *bar(j)) / fhad(j)
          np =max(0.e0,yprot(j)  *bar(j)) / fhad(j)
          nla=max(0.e0,ylambda(j)*bar(j)) / fhad(j)
          nsm=max(0.e0,ysminus(j)*bar(j)) / fhad(j)
          ns0=max(0.e0,yszero(j) *bar(j)) / fhad(j)
          nsp=max(0.e0,ysplus(j) *bar(j)) / fhad(j)
         else
          nn =0.d0
          np =0.d0
          nla=0.d0
          nsm=0.d0
          ns0=0.d0
          nsp=0.d0
         end if
         if (fhad(j).ne.1.d0) then
          nqu=max(0.e0,yquarku(j)*bar(j)) / (1.d0-fhad(j))
          nqd=max(0.e0,yquarkd(j)*bar(j)) / (1.d0-fhad(j))
          nqs=max(0.e0,yquarks(j)*bar(j)) / (1.d0-fhad(j))
         else
          nqu=0.d0
          nqd=0.d0
          nqs=0.d0
         end if
         kfn(j) =(3.d0*pi**2*nn )**(1.d0/3.d0)
         kfp(j) =(3.d0*pi**2*np )**(1.d0/3.d0)
         kfla(j)=(3.d0*pi**2*nla)**(1.d0/3.d0)
         kfsm(j)=(3.d0*pi**2*nsm)**(1.d0/3.d0)
         kfs0(j)=(3.d0*pi**2*ns0)**(1.d0/3.d0)
         kfsp(j)=(3.d0*pi**2*nsp)**(1.d0/3.d0)
         kfqu(j)=(3.d0*pi**2*nqu)**(1.d0/3.d0)
         kfqd(j)=(3.d0*pi**2*nqd)**(1.d0/3.d0)
         kfqs(j)=(3.d0*pi**2*nqs)**(1.d0/3.d0)
c----------------------------------------------
         ne=abs(yelect(j)*bar(j))
         nm=abs( ymuon(j)*bar(j))
         kfe(j)=(3.d0*pi**2*ne)**(1.d0/3.d0)
         if (yelect(j).le.0.d0) kfe(j)=-kfe(j)
         kfm(j)=(3.d0*pi**2*nm)**(1.d0/3.d0)
         if ( ymuon(j).le.0.d0) kfm(j)=-kfm(j)
c----------------------------------------------
c Check all this stuff:
         coeff=3.d0*pi**2
         nn = kfn(j)**3/coeff * fhad(j)
         np = kfp(j)**3/coeff * fhad(j)
         nla=kfla(j)**3/coeff * fhad(j)
         nsm=kfsm(j)**3/coeff * fhad(j)
         ns0=kfs0(j)**3/coeff * fhad(j)
         nsp=kfsp(j)**3/coeff * fhad(j)
         nqu=kfqu(j)**3/coeff * (1.d0-fhad(j))
         nqd=kfqd(j)**3/coeff * (1.d0-fhad(j))
         nqs=kfqs(j)**3/coeff * (1.d0-fhad(j))
         ne = kfe(j)**3/coeff * 1.d0
         nm = kfm(j)**3/coeff * 1.d0
         charge_l=-ne-nm
         charge_h=(np+nsp-nsm)
         charge_q=(2./3.*nqu-1./3.*nqd-1./3.*nqs)
         charge=(charge_l+charge_h+charge_q)
         baryon_h=(nn+np+nla+nsm+ns0+nsp)
         baryon_q=1./3.*(nqu+nqd+nqs)
         baryon=(baryon_h+baryon_q)
         if (abs(charge).ge.1.d-2) then
          print *,'Charge neutrality violated at:'
          print '(i5,a6,1p1e12.3,a18,0p1f12.5)',
     1          j,'Rho= ',rrho(j),': charge/fm3= ',charge
          pause
         end if
         barrel=baryon/bar(j)
         if ((abs(barrel)-1.d0).ge.1.d-2) then
          print *,'Baryons do not sum up to baryon density at:'
          print '(i5,a6,1p1e12.3,a26,0p1f12.5)',
     1        j,'Rho= ',rrho(j),': sum(baryons)/baryon#= ',baryon
          pause
         end if
        end do
c ****  Calculate the fermi momenta in the crust:
        do j=icore+1,imax
         ne=yelect(j)*bar(j)
         nn=yneutr(j)*bar(j)
         kfe(j) =(3.d0*pi**2*ne)**(1.d0/3.d0)
         kfm(j) =0.d0
         kfn(j) =(3.d0*pi**2*nn)**(1.d0/3.d0)
         kfp(j) =0.d0
         kfla(j)=0.d0
         kfsm(j)=0.d0
         kfs0(j)=0.d0
         kfsp(j)=0.d0
         kfqu(j)=0.d0
         kfqd(j)=0.d0
         kfqs(j)=0.d0
        end do
c ******************************************************
c CHECK FOR DIRECT URCAS *******************************
c ******************************************************
        do j=0,icore
c n-p:
         if    ( (kfp(j) .lt.kfn(j) +kfe(j) ).and.
     1           (kfn(j) .lt.kfp(j) +kfe(j) ).and.
     2           (kfe(j) .lt.kfp(j) +kfn(j) )     ) then
          idurca_np(j)=1
          if    ( (kfp(j) .lt.kfn(j) +kfm(j) ).and.
     1            (kfn(j) .lt.kfp(j) +kfm(j) ).and.
     2            (kfm(j) .lt.kfp(j) +kfn(j) )     ) then
           idurca_np(j)=2
          end if
         else
          idurca_np(j)=0
         end if
c la-p:
         if     ( (kfp(j) .lt.kfla(j)+kfe(j) ).and.
     1            (kfla(j).lt. kfp(j)+kfe(j) ).and.
     2            (kfe(j) .lt. kfp(j)+kfla(j))     ) then
          idurca_lap(j)=1
          if    ( (kfp(j) .lt.kfla(j)+kfm(j) ).and.
     1          (kfla(j).lt. kfp(j)+kfm(j) ).and.
     2          (kfm(j) .lt. kfp(j)+kfla(j))     ) then
           idurca_lap(j)=2
          end if
         else
          idurca_lap(j)=0
         end if
c sm-n:
         if     ( (kfsm(j).lt.kfn(j) +kfe(j) ).and.
     1            (kfn(j) .lt.kfsm(j)+kfe(j) ).and.
     2            (kfe(j) .lt.kfsm(j)+kfn(j) )     ) then
          idurca_smn(j)=1
          if    ( (kfsm(j).lt.kfn(j) +kfm(j) ).and.
     1          (kfn(j) .lt.kfsm(j)+kfm(j) ).and.
     2          (kfm(j) .lt.kfsm(j)+kfn(j) )     ) then
           idurca_smn(j)=2
          end if
         else
          idurca_smn(j)=0
         end if
c sm-la:
         if     ( (kfsm(j).lt.kfla(j)+kfe(j) ).and.
     1            (kfla(j).lt.kfsm(j)+kfe(j) ).and.
     2            (kfe(j) .lt.kfsm(j)+kfla(j))     ) then
          idurca_smla(j)=1
          if    ( (kfsm(j).lt.kfla(j)+kfm(j) ).and.
     1            (kfla(j).lt.kfsm(j)+kfm(j) ).and.
     2            (kfm(j) .lt.kfsm(j)+kfla(j))    ) then
           idurca_smla(j)=2
          end if
         else
          idurca_smla(j)=0
         end if
c sm-s0:
         if     ( (kfsm(j).lt.kfs0(j)+kfe(j) ).and.
     1            (kfs0(j).lt.kfsm(j)+kfe(j) ).and.
     2            (kfe(j) .lt.kfsm(j)+kfs0(j))     ) then
          idurca_sms0(j)=1
          if    ( (kfsm(j).lt.kfs0(j)+kfm(j) ).and.
     1           (kfs0(j).lt.kfsm(j)+kfm(j) ).and.
     2           (kfm(j) .lt.kfsm(j)+kfs0(j))     ) then
           idurca_sms0(j)=2
          end if
         else
          idurca_sms0(j)=0
         end if
c qu-qd:
         if     ( (kfqu(j).lt.kfqd(j)+kfe(j) ).and.
     1            (kfqd(j).lt.kfqu(j)+kfe(j) ).and.
     2            (kfe(j) .lt.kfqu(j)+kfqd(j))     ) then
          idurca_quqd(j)=1
          if    ( (kfqu(j).lt.kfqd(j)+kfm(j) ).and.
     1            (kfqd(j).lt.kfqu(j)+kfm(j) ).and.
     2            (kfm(j) .lt.kfqu(j)+kfqd(j))     ) then
           idurca_quqd(j)=2
          end if
         else
          idurca_quqd(j)=0
         end if
c qu-qs:
         if     ( (kfqu(j).lt.kfqs(j)+kfe(j) ).and.
     1            (kfqs(j).lt.kfqu(j)+kfe(j) ).and.
     2            (kfe(j) .lt.kfqu(j)+kfqs(j))     ) then
          idurca_quqs(j)=1
          if    ( (kfqu(j).lt.kfqs(j)+kfm(j) ).and.
     1            (kfqs(j).lt.kfqu(j)+kfm(j) ).and.
     2            (kfm(j) .lt.kfqu(j)+kfqs(j))     ) then
           idurca_quqs(j)=2
          end if
         else
          idurca_quqs(j)=0
         end if
        end do     
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_effective_masses
c *********************************************************************
c      This subroutine calculates the baryon effective masses
c      (in a lousy way) incase thery are not given by the EOS table !
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'fermi.inc.f'
c       INCLUDE 'files_in_out.inc.f'
       INCLUDE 'files.inc.f'
c Crust neutrons:
        if (emncr.lt.1.) then
         do i=icore+1,idrip
          mstn(i)=emncr
         end do
        else if (emncr.eq.1.) then
         do i=icore+1,idrip
          mstn(i)=mstn_ccks_r_cbf(kfn(i))
         end do
        else if (emncr.eq.2.) then
         do i=icore+1,idrip
          mstn(i)=mstn_ccks_bj_cbf(kfn(i))
         end do
        else if (emncr.eq.3.) then
         do i=icore+1,idrip
          mstn(i)=mstn_bks(kfn(i))
         end do
        else if (emncr.eq.4.) then
         do i=icore+1,idrip
          mstn(i)=mstn_t(kfn(i))
         end do
        else if (emncr.eq.5.) then
         do i=icore+1,idrip
          mstn(i)=mstn_pref(kfn(i))
         end do
        else if (emncr.eq.6.) then
         do i=icore+1,idrip
          mstn(i)=mstn_awp_3(kfn(i))
         end do
        else if (emncr.eq.7.) then
         do i=icore+1,idrip
          mstn(i)=mstn_awp_3_d(kfn(i))
         end do
        else if (emncr.eq.8.) then
         do i=icore+1,idrip
          mstn(i)=mstn_ccks_r_var(kfn(i))
         end do
        else if (emncr.eq.9.) then
         do i=icore+1,idrip
          mstn(i)=mstn_ccks_bj_var(kfn(i))
         end do
        end if
c If version = NEW exit because core mst's already gotten from EOS table
        if ((version.eq.'NEW').or.(version.eq.'QRK')) return
c Core neutrons:
        if (emnco.lt.1.) then
         do i=0,icore
          mstn(i)=emnco
         end do
        else if (emnco.eq.1.) then
         do i=0,icore
          mstn(i)=mstn_ccks_r_cbf(kfn(i))
         end do
        else if (emnco.eq.2.) then
         do i=0,icore
          mstn(i)=mstn_ccks_bj_cbf(kfn(i))
         end do
        else if (emnco.eq.3.) then
         do i=0,icore
          mstn(i)=mstn_bks(kfn(i))
         end do
        else if (emnco.eq.4.) then
         do i=0,icore
          mstn(i)=mstn_t(kfn(i))
         end do
        else if (emnco.eq.5.) then
         do i=0,icore
          mstn(i)=mstn_pref(kfn(i))
         end do
        end if
c Protons:
        if (emp.lt.1.) then
         do i=0,icore
          mstp(i)=emp
         end do
        else if (emp.eq.1.) then
         do i=0,icore
          mstp(i)=mstp_ccy_ms2(kfp(i))
         end do
        else if (emp.eq.2.) then
         do i=0,icore
          mstp(i)=mstp_ccy_ps1(kfp(i))
         end do
        else if (emp.eq.3.) then
         do i=0,icore
          mstp(i)=mstp_awp(kfp(i))
         end do
        end if
c Hyperons:
        do i=0,icore
         mstla(i)=0.8d0
         mstsm(i)=0.8d0
         msts0(i)=0.8d0
         mstsp(i)=0.8d0
        end do
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_spec_heat_degenerate
c *********************************************************************
c      This subroutine calculates Cv/T for degenerate particles
c *********************************************************************
       implicit real*8 (a-h,k-z)
       parameter(pi=3.14159265)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'spec_heat.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'quark.inc.f'
        do j=0,imax
         pfe=kfe(j)*197.d0
         me =sqrt(0.511d0**2+pfe**2)   ! electron effective mass
         pfm=kfm(j)*197.d0
         mm =sqrt(105.d0**2+pfm**2)    ! muon effective mass
         pfn=kfn(j)*197.d0
         mn =939.56d0*mstn(j)
         pfp=kfp(j)*197.d0
         mp =938.27d0*mstp(j)
         pfla=kfla(j)*197.d0
         mla=1116.0d0*mstla(j)
         pfsm=kfsm(j)*197.d0
         msm=1193.0d0*mstsm(j)
         pfs0=kfs0(j)*197.d0
         ms0=1193.0d0*msts0(j)
         pfsp=kfsp(j)*197.d0
         msp=1193.0d0*mstsp(j)
         pfqu=kfqu(j)*197.d0
         mqu=sqrt(5.d0**2+pfqu**2)            ! up effective mass
         pfqd=kfqd(j)*197.d0
         mqd=sqrt(8.d0**2+pfqd**2)            ! down effective mass
         pfqs=kfqs(j)*197.d0
         mqs=sqrt(strange_mass**2+pfqs**2)    ! strange effective mass
         cve(j) =cvt_deg(pfe,me)
         cvm(j) =cvt_deg(pfm,mm)
         cvn(j) =cvt_deg(pfn,mn)   * fhad(j)
         cvp(j) =cvt_deg(pfp,mp)   * fhad(j)
         cvla(j)=cvt_deg(pfla,mla) * fhad(j)
         cvsm(j)=cvt_deg(pfsm,msm) * fhad(j)
         cvs0(j)=cvt_deg(pfs0,ms0) * fhad(j)
         cvsp(j)=cvt_deg(pfsp,msp) * fhad(j)
         cvqu(j)=cvt_deg(pfqu,mqu) * (1.d0-fhad(j))
         cvqd(j)=cvt_deg(pfqd,mqd) * (1.d0-fhad(j))
         cvqs(j)=cvt_deg(pfqs,mqs) * (1.d0-fhad(j))
        end do   
       return
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_Tc
c *********************************************************************
c      This subroutine calculate Tc for all baryons and quarks
c *********************************************************************
       implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_star.inc.f'
       INCLUDE 'profile_comp.inc.f'
       INCLUDE 'fermi.inc.f'
       INCLUDE 'pairing.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       dimension tcrit(10),rho_lo(10),rho_hi(10)
       integer,parameter :: p_kfd = 1000
       real*8,dimension(1:p_kfd) :: kfInt
       real*8,dimension(1:p_kfd) :: DeltaInt 
       integer kfind,p_I

c Just to be safe:
       do i=0,imax
        tcn(i) =1.0d0
        tcp(i) =1.0d0
        tcla(i)=1.0d0
        tcuu(i)=1.0d0
        tcdd(i)=1.0d0
        tcss(i)=1.0d0
        tcud(i)=1.0d0
        tcus(i)=1.0d0
        tcds(i)=1.0d0
        tcu(i) =1.0d0
        tcd(i) =1.0d0
        tcs(i) =1.0d0
       end do
c ***** 1s0 neutron superfluidity *************************************
      if (sfn1s0.eq.1.) then
       do i=0,idrip
        tcn(i)=max(1.d0,tcn1_sfb(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.2.) then
       do i=0,idrip
        tcn(i)=max(1.d0,tcn1_ccdk(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.3.) then
       do i=0,idrip
        tcn(i)=max(1.d0,tcn1_wap(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.4.) then
       do i=0,idrip
        tcn(i)=max(1.d0,tcn1_gc(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.5.) then
       do i=0,idrip
        tcn(i)=max(1.d0,tcn1_gipsf(kfn(i)))*fn1s0
       end do
c     Ioffe gaps:
      else if (sfn1s0.eq.201.) then
       do i=0,idrip
        tcn(i)=max(1.d0,Tc_Ioffe_1ns(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.202.) then
       do i=0,idrip
        tcn(i)=max(1.d0,Tc_Ioffe_2ns(kfn(i)))*fn1s0
       end do
      else if (sfn1s0.eq.203.) then
       do i=0,idrip
        tcn(i)=max(1.d0,Tc_Ioffe_3ns(kfn(i)))*fn1s0
       end do
      end if
c ***** 3p2 neutron superfluidity *************************************
      isf=-1
c isf will be the largest zone number at which triplet pairing is present
      if(sfn3p2.eq.1.)then
       do i=0,idrip
        temp=tcn3_hgrr(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.2.)then
       do i=0,idrip
        temp=tcn3_ao(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.3.)then
       do i=0,idrip
        temp=tcn3_ao_m1(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.4.)then
       do i=0,idrip
        temp=tcn3_t72(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.5.)then
       temp=tcn3_t72_m1(kfn(i))*fn3p2
       do i=0,idrip
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.6.)then
       do i=0,idrip
        temp=tcn3_bcll92(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.7.)then
       do i=0,idrip
        temp=tcn3_eehjo96_nr(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if(sfn3p2.eq.8.)then
       do i=0,idrip
        temp=tcn3_eehjo96_r(kfn(i))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
c     Minimal Cooling paper gaps:
      else if((sfn3p2.ge.100.).and.(sfn3p2.lt.200.)) then
       if (sfn3p2.eq.101.) then          ! Gap "a"
        kfmax_n3p2=1.8d0
        delkf_n3p2=0.5d0
        tcmax_n3p2=1.0d9
       else if(sfn3p2.eq.102.) then      ! Gap "b"
        kfmax_n3p2=2.0d0
        delkf_n3p2=0.5d0
        tcmax_n3p2=3.0d9
       else if(sfn3p2.eq.103.) then      ! Gap "c"
        kfmax_n3p2=2.5d0
        delkf_n3p2=0.7d0
        tcmax_n3p2=1.0d10
       else if(sfn3p2.eq.150.) then
        ! Nothing to do: parameters were read from I_Pairing*.dat file
       end if

       open(unit=20,file='SCGF.txt')
       do p_I = 1,p_kfd
        read(20,*) kfInt(p_I)
        read(20,*) DeltaInt(p_I)
       end do
       close(20)

       do i=0,idrip
        if(sfn3p2.eq.104.) then      ! new Gap model SCGF
         kfind=int(kfn(i)/10d0*p_kfd)
         temp=DeltaInt(kfind) * 1.16d10 / 0.8416d0 
c         write(*,*) i,kfind,kfn(i),temp,DeltaInt(kfind)
        else
         temp=tcmax_n3p2*
     1       exp(-(kfn(i)-kfmax_n3p2)**2/delkf_n3p2**2)*fn3p2
        end if
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         if (isf.eq.i-1) isf=i
        end if
       end do
c     Ioffe gaps:
      else if (sfn3p2.eq.201.) then
       do i=0,idrip
        temp=max(1.d0,Tc_Ioffe_1nt(kfn(i)))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if (sfn3p2.eq.202.) then
       do i=0,idrip
        temp=max(1.d0,Tc_Ioffe_2nt(kfn(i)))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
      else if (sfn3p2.eq.203.) then
       do i=0,idrip
        temp=max(1.d0,Tc_Ioffe_3nt(kfn(i)))*fn3p2
        if (temp.ge.tcn(i)) then
         tcn(i)=temp
         isf=i
        end if
       end do
c     Uniform Tc gap:
      else if(sfn3p2.ge.1.e3)then
       do i=0,icore
        tcn(i)=sfn3p2
       end do
       isf=icore
      end if
c ***** 1s0 proton superconductivity **********************************
      if(sfp1s0.eq.1.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_ccy_ms(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.2.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_ccy_ps(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.3.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_t73(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.4.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_ns(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.5.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_ao(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.6.)then  
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_bcll92(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.7.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcp1_ccdk(kfp(i)))*fp1s0
       end do
c     Using Neutron 1S0 gaps for Protons:
      else if(sfp1s0.eq.21.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcn1_t72(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.22.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcn1_awp_2(kfp(i)))*fp1s0
       end do
      else if(sfp1s0.eq.23.)then
       do i=0,icore
        tcp(i)=max(1.d0,tcn1_awp_3(kfp(i)))*fp1s0
       end do
c     Ioffe gaps:
      else if (sfp1s0.eq.201.) then
       do i=0,idrip
        tcp(i)=max(1.d0,Tc_Ioffe_1p(kfn(i)))*fp1s0
       end do
      else if (sfp1s0.eq.202.) then
       do i=0,idrip
        tcp(i)=max(1.d0,Tc_Ioffe_2p(kfn(i)))*fp1s0
       end do
      else if (sfp1s0.eq.203.) then
       do i=0,idrip
        tcp(i)=max(1.d0,Tc_Ioffe_3p(kfn(i)))*fp1s0
       end do
c     Uniform Tc gap:
      else if(sfp1s0.ge.1.e3)then
       do i=0,icore
        tcp(i)=sfp1s0
       end do
      end if
c ***** 1s0 Lambda superfluidity **************************************
      if (sfl1s0.eq.1.) then
       do i=0,icore
        tcla(i)=max(1.d0,tcla1_bb(kfla(i),bar(i)))*fl1s0
       end do
      end if 
c ***** Quark pairing *************************************************
      open(unit=32,file=f_quark_tc,status='old',err=6616)
       read(32,*)
       read(32,*)
       read(32,*)
       read(32,*)pf0_uu,dpf_uu,gap_uu
       read(32,*)pf0_dd,dpf_dd,gap_dd
       read(32,*)pf0_ss,dpf_ss,gap_ss
       read(32,*)pf0_ud,dpf_ud,gap_ud
       read(32,*)pf0_us,dpf_us,gap_us
       read(32,*)pf0_ds,dpf_ds,gap_ds
      close(unit=32,status='keep')
      do i=0,icore
c **** uu:
       pf_uu=kfqu(i)*.197
       gap=dexp(-(pf_uu-pf0_uu)**2/dpf_uu**2)
       tcuu(i)=1.1604e13*gap_uu*gap
       if (yquarku(i).eq.0.) tcuu(i)=0.0
c **** dd:
       pf_dd=kfqd(i)*.197
       gap=dexp(-(pf_dd-pf0_dd)**2/dpf_dd**2)
       tcdd(i)=1.1604e13*gap_dd*gap
       if (yquarkd(i).eq.0.) tcdd(i)=0.0
c **** ss:
       pf_ss=kfqs(i)*.197
       gap=dexp(-(pf_ss-pf0_ss)**2/dpf_ss**2)
       tcss(i)=1.1604e13*gap_ss*gap
       if (yquarks(i).eq.0.) tcss(i)=0.0
c **** ud:
       pf_ud=(kfqu(i)+kfqd(i))/2.*.197
       gap=dexp(-(pf_ud-pf0_ud)**2/dpf_ud**2)
       tcud(i)=1.1604e13*gap_ud*gap
       if ((yquarku(i).eq.0.).or.(yquarkd(i).eq.0.)) tcud(i)=0.0
c *** us:
       pf_us=(kfqu(i)+kfqs(i))/2.*.197
       gap=dexp(-(pf_us-pf0_us)**2/dpf_us**2)
       tcus(i)=1.1604e13*gap_us*gap
       if ((yquarku(i).eq.0.).or.(yquarks(i).eq.0.)) tcus(i)=0.0
c *** ds:
       pf_ds=(kfqd(i)+kfqs(i))/2.*.197
       gap=dexp(-(pf_ds-pf0_ds)**2/dpf_ds**2)
       tcds(i)=1.1604e13*gap_ds*gap
       if ((yquarkd(i).eq.0.).or.(yquarks(i).eq.0.)) tcds(i)=0.0
c *** Take only the maximum Tc:
       tcu(i)=max(tcuu(i),tcud(i),tcus(i))
       tcd(i)=max(tcdd(i),tcud(i),tcds(i))
       tcs(i)=max(tcss(i),tcus(i),tcds(i))
       if (fhad(i).eq.1.0) then
        tcu(i)=1.d0
        tcd(i)=1.d0
        tcs(i)=1.d0
       end if
      end do
      goto 6617
c Following is in case the file "f_quark_tc" is not found:
 6616 continue
c      print *,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
c      print *,'WARNING:'
c      print *,'  File ',f_quark_tc(1:istrlen(f_quark_tc)),
c     x        'not found ==> Quark Tc=0 !'
c      print *,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      do i=0,icore
       tcu(i)=1.d0
       tcd(i)=1.d0
       tcs(i)=1.d0
      end do
 6617 continue
c *********************************************
c Just in case the above formulas give Tc<0 instead of Tc=0
c at the edges where Tc is almost 0
        do i=0,idrip
         tcn(i) =abs( tcn(i))
         tcp(i) =abs( tcp(i))
         tcla(i)=abs(tcla(i))
         tcuu(i)=abs(tcuu(i))
         tcdd(i)=abs(tcdd(i))
         tcss(i)=abs(tcss(i))
         tcud(i)=abs(tcud(i))
         tcus(i)=abs(tcus(i))
         tcds(i)=abs(tcds(i))
         tcu(i) =abs( tcu(i))
         tcd(i) =abs( tcd(i))
         tcs(i) =abs( tcs(i))
        end do
       return
c *********************************************************************
      end
c *********************************************************************
c *********************************************************************
c *********************************************************************
      subroutine get_degenerate_density
       Implicit real*8 (a-h,k-z)
       INCLUDE 'size.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'profile_star.inc.f'
       dimension rho2(0:isize),pres2(0:isize)
       dimension theta_vib(4,0:isize)
c *** Read crust chemical composition:
        open(unit=16,file=f_crusteos,status='old')
         read(16,*)itext,idata
         do i=1,itext
          read(16,*)
         end do
         do j=idata,1,-1
          read(16,*)rho2(j),pres2(j)
         end do
        close(unit=16,status='keep')
c ***
        do i=0,ienv-1
         rhod(i)=rrho(i)
        end do
        do i=imax,ienv,-1
         j=0
765      j=j+1
         if ((pres(i).ge.pres2( j )).and.
     1       (pres(i).le.pres2(j+1))     ) then
          x=(dlog(pres(i))-dlog(pres2(j)))/
     1      (dlog(pres2(j+1))-dlog(pres2(j)))
          y=(dlog(pres2(j+1))-dlog(pres(i)))/
     1      (dlog(pres2(j+1))-dlog(pres2(j)))
          lrhod=y*dlog(rho2(j))+x*dlog(rho2(j+1))
          rhod(i)=max(dexp(lrhod),rrho(i))
          if ((i.lt.imax).and.(rrho(i+1).eq.rhod(i+1))) then
           rhod(i)=rrho(i)
          end if
         else
          goto 765
         end if
        end do
       return
      end
c *********************************************************************
c *********************************************************************
      subroutine get_ion_vibration_frequency
       Implicit real*8 (a-h,k-z)
c       INCLUDE 'files_in_out.inc.f'
c       INCLUDE 'files_phys.inc.f'
       INCLUDE 'files.inc.f'
       INCLUDE 'size.inc.f'
       INCLUDE 'profile_comp.inc.f'
       dimension theta_vib(4,0:isize)
       common/cvn_vib/theta_vib
c ****** to be checked !!!!!!!!!!!!!!!
        pi=dacos(-1.d0)
        hb=1.054588d-27
        kb=1.380662d-16
        c=2.997924d10
        mev=1.602d-6
        do i=icore+2,imax,2
         y=1.d0*z_ion(i)/a_ion(i)
         les=3.251d0-(36.810d0+y**3.10d0)**0.3226d0
         es=10.d0**les*mev
         r0=1.2d0*a_ion(i)**(1.d0/3.d0)
         x1=r0**2*es
         x2=3.d0/2.d0/pi*(4.803d-10*z_ion(i))**2/r0*1.d13
         rbar=(4.d0/3.d0*pi*bar(i))**(-1.d0/3.d0)
         x3=(r0/rbar)**3
         c1=0.0d0*x1-x2*(0.d0/3.d0-.5d0*x3)
         c2=4.0d0*x1-x2*(1.d0/5.d0-.5d0*x3)
         c3=10.d0*x1-x2*(2.d0/7.d0-.5d0*x3)
         mu=1.6d-24
         b0=3.d0*mu*a_ion(i)*r0**2/4.d0/pi
         b1=b0/1.d0
         b2=b0/2.d0
         b3=b0/3.d0
         om1=sqrt(c1/b1)*1.d13
         om2=sqrt(c2/b2)*1.d13
         om3=sqrt(c3/b3)*1.d13
         theta_vib(1,i)=     hb*om1/kb
         theta_vib(2,i)=     hb*om2/kb
         theta_vib(3,i)=3.d0*hb*om2/kb
         theta_vib(4,i)=     hb*om3/kb
        end do
       return
      end
c***************************************************************************
c***************************************************************************
      function cvt_deg(pf,m)
c**************************************************************c
c  Calculates the degenerate specific heat per cm^3 over T:    c
c                       Cv/T                                   c
c  For spin 1/2 fermions !                                     c
c  pf and m must be in MeV, but cvt is returned in cgs units.  c
c  m must be the Landau effective mass, i.e., m* for baryons   c
c  and sqrt(m**2+pf**2) for leptons !                          c
c              CHECKED ON MARCH 913, 2001                      c
c**************************************************************c
       implicit real*8(a-h,k-z)
       parameter(pi=3.14159265d0)
       parameter (kb=1.38d-16, MeV=1.602d-6)
        if (pf.eq.0.0d0) then
         cvt_deg=0.0d0
        else
         N0 = 2.d0 * m*pf/2.d0/pi**2
         cvt_deg=pi**2/3.d0 * N0
        end if
        cvt_deg=cvt_deg * kb**2/MeV/197.d0**3*1.d39 ! converts to cgs units
       return
      end
c***************************************************************************
c***************************************************************************
      function cvt_deg_old1(n,m)
c**************************************************************c
c  Calculates the degenerate specific heat over T: Cv/T        c
c  For spin 1/2 fermions !                                     c
c  n must be in #/fm3 and m in MeV. cvt is given in cgs units. c
c  m must be the Landau effective mass, i.e., m* for baryons   c
c  and sqrt(m**2+kf**2) for leptons !                          c
c              CHECKED ON MARCH 9 2001                         c
c**************************************************************c
       implicit real*8(a-h,k-z)
       parameter(pi=3.14159265d0)
       parameter (kb=1.38d-16,MeV=1.602d-6,fm3=1.d39)
        if (n.eq.0.0d0) then
         cvt_deg_old1=0.0d0
        else
         pf=(3.d0*pi**2*n)**(1.d0/3.d0)*197.d0
         coeff=m/pf**2/MeV
         cvt_deg_old1=pi**2 * n*fm3 * kb**2 * coeff
        end if
       return
      end
c***************************************************************************
c***************************************************************************
      function cvt_deg_old2(n,m)
c**************************************************************c
c  Calculates the degenerate specific heat over T: Cv/T        c
c  For spin 1/2 fermions !                                     c
c  n must be in #/fm3 and m in MeV. cvt is given in cgs units. c
c  The function works for massless fermions also.              c
c  WRONG: for interacting particle, i.e., baryons: double      c
c         counting of the effective mass !                     c
c**************************************************************c
       implicit real*8(a-h,k-z)
       parameter(pi=3.14159265d0)
       parameter (kb=1.38d-16,MeV=1.602d-6,fm3=1.d39)
        if (n.eq.0.0d0) then
         cvt_deg_old2=0.0d0
        else
         pf=(3.d0*pi**2*n)**(1.d0/3.d0)*197.d0
         coeff=sqrt(m**2+pf**2)/pf**2/MeV
         cvt_deg_old2=pi**2 * n*fm3 * kb**2 * coeff
        end if
       return
      end
c***************************************************************************
c***************************************************************************
      function cvt_deg_old3(pf,m)
c**************************************************************c
c  Calculates the degenerate specific heat over T: Cv/T        c
c  For spin 1/2 fermions !                                     c
c  pf and m in MeV. cvt is given in cgs units.                 c
c  The function works for massless fermions also.              c
c  WRONG: for interacting particle, i.e., baryons: double      c
c         counting of the effective mass !                     c
c**************************************************************c
       implicit real*8(a-h,k-z)
       parameter(pi=3.14159265d0)
       parameter (kb=1.38d-16,MeV=1.602d-6,fm3=1.d39)
        if (pf.eq.0.0d0) then
         cvt_deg_old3=0.0d0
        else
         n=(pf/197.)**3/3./pi*2
         coeff=sqrt(m**2+pf**2)/pf**2/MeV
         cvt_deg_old3=pi**2 * n*fm3 * kb**2 * coeff
        end if
       return
      end
c***************************************************************************
c***************************************************************************

