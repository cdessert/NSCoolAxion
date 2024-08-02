c *********************************************************************
c *********************************************************************
c *********************************************************************
c                                  TEMP
c *********************************************************************
c *********************************************************************
c *********************************************************************
      if(ptemp.eq.1.)then
       if (debug.ge.1.) print *,'Printing in Temp file'
       if ( ( ((time+dtime)/year) .ge. tprint(itprint)) .and.
     1     (itprint.le.itpmax)) then
  
        itprint=itprint+1
  
        w1=((time+dtime)-tprint(itprint-1)*year)/dtime
        w2=1.d0-w1
        logtemp=w1*dlog(sign_l*oteffective)+
     1          w2*dlog(sign_l* teffective)
        t_effective=sign_l*dexp(logtemp)

        write(10,514)' Time=',tprint(itprint-1),
     1               ' Te_inf=',t_effective
514     format(2(5x,1a7,1p1e10.3))
        write(10,*)
        write(10,516)
     1       ' zone ',
     2       'Rad [m]  ','Rho [gm/cc]','exp(phi)','dvol ',
     3       'Temp [K]  ','Lum [Lsun]  '
516     format(
     1       a6,
     2       4a15,2x,
     3       2a15,2x)
        write(10,'(68x,2a15)')'no exp(phi) ','no exp(2phi) '
        write(10,*)
        do i=imax,1,-2
         logtemp=w1*log(otemp(i))+w2*log(temp(i))
         temperature=dexp(logtemp)
         lumino=0.0d0
         if (i-1.ne.0) then
         loglum=w1*log(abs(olum(i-1)))+w2*log(abs(lum(i-1)))
         lumino=dexp(loglum)
         if (lum(i-1).le.0.0) lumino=-lumino
         end if
         write(10,517)
     1        i,
     2        rad(i)/1.d2,rrho(i),ephi(i),dvol(i)+dvol(i+1),
     3        temperature/ephi(i),lumino/e2phi(i-1)
        end do
        write(10,*)
517     format(
     1        i5,1x,
     2        1p4e15.7,2x,
     3        1p2e15.7,2x)
       end if
      end if
c *********************************************************************

c *********************************************************************
c *********************************************************************
c *********************************************************************
c                               TEFF
c *********************************************************************
c *********************************************************************
c *********************************************************************

      if (pteff.ge.1.0) then
       if (debug.ge.1.) print *,'Printing in Teff file'
       if ((float(idump1)/2.).ne.float(idump1/2)) then
c        temp1=temp(idump1)/ephi(idump1)
        temp1=temp(idump1)
       else
c        temp1=(temp(idump1-1)+temp(idump1+1))/2./ephi(idump1)
        temp1=(temp(idump1-1)+temp(idump1+1))/2.
       end if
       if ((float(idump2)/2.).ne.float(idump2/2)) then
c        temp2=temp(idump2)/ephi(idump2)
        temp2=temp(idump2)
       else
c        temp2=(temp(idump2-1)+temp(idump2+1))/2./ephi(idump2)
        temp2=(temp(idump2-1)+temp(idump2+1))/2.
       end if
       if ((float(idump3)/2.).ne.float(idump3/2)) then
c        temp3=temp(idump3)/ephi(idump3)
        temp3=temp(idump3)
       else
c        temp3=(temp(idump3-1)+temp(idump3+1))/2./ephi(idump3)
        temp3=(temp(idump3-1)+temp(idump3+1))/2.
       end if
       cv_core =0.d0
       cv_crust=0.d0
       do i=1,icore,2
        cv_core=cv_core+(dvol(i)+dvol(i+1))*cv(i)
       end do
       do i=icore+2,imax,2
        cv_crust=cv_crust+(dvol(i)+dvol(i+1))*cv(i)
       end do
c*********************************************************************
       if (pteff.eq.1.0) then
        write(19,601) istep,(time+dtime)/year,teffective,
     1               lum(imax-1)*lsol,lnu*lsol,lh*lsol,Bcurr,tcenter
        write(*,*) 'res: ',(time+dtime)/year,lum(imax-1),lnu
 
 601    format(i8,1p1e12.3,5x,1p4e12.3,1p4e12.3,1p4e12.3)
c*********************************************************************
       else
        pause 'WARNING: Not Teff print out defined !'
       end if
c*********************************************************************
c       write(19,615) istep,(time+dtime)/year,teffective,
c     1               bf_r(imax),period,
c     2               dtime/year,dtime/odtime,max_dtemp,itrial,
c     3               omegab_tau1,omegab_tau2,omegab_tau3,
c     4               icycle,m_dot*3.15e7/2.e33,
c     5               htot*lsol,lh*lsol,
c     6               lcore+lexo,lcrust
c*********************************************************************
c       write(19,615) istep,(time+dtime)/year,teffective,
c     1               bf_r(imax),period,
c     2               dtime/year,dtime/odtime,max_dtemp,itrial,
c     3               temp1,temp2,temp3,
c     4               icycle,m_dot*3.15e7/2.e33,
c     5               htot*lsol,lh*lsol,(lcore+lexo)*lsol,lcrust*lsol,
c     6               cv_core,cv_crust
c615    format(i8 , 1p1e22.15 , 1p1e16.8 ,
c     1        1p1e12.3 , 1p1e12.3 ,        
c     2        1p1e10.2 , 0p2f8.3, i8 ,
c     3        1p3e16.8 ,
c     4        i8 , 1p1e14.5,
c     5        1p4e12.3,
c     6        1p2e12.3)       
c*********************************************************************
c       write(19,615) istep,(time+dtime)/year,teffective,
c     1               bf_r(imax),
c     2               lum(imax-1)*lsol,lnu_tot,
c     3               lmurca_nucl,lbrem_nucl,lplasma,
c     4               lpbf_n1S0,lpbf_n3P2,lpbf_p1S0
c615    format(i8 , 1p1e12.3 , 1p1e12.3 , 3x ,
c     1        1p1e12.3, 5x ,
c     2        1p1e12.3 , 2x , 1p1e12.3 , 2x ,
c     3        1p3e12.3 ,
c     4        1p3e12.3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       menv= (emas(imax)/gs14**2) * eta
c       write(19,615) istep,(time+dtime)/year,teffective,ts1,ts2,
c     1               temp1,temp2,temp3,
c     2               bf_r(imax),
c     3               lum(imax-1)*lsol,lnu_tot,
c     4               lmurca_nucl,lbrem_nucl,lplasma,lnpb,
c     5               lpbf_n1S0,lpbf_n3P2,lpbf_p1S0
c615    format(i8 , 1p1e12.3 , 1p3e12.3 , 3x ,
c     1        1p3e12.3,        
c     2        1p1e12.3, 5x ,
c     3        1p1e12.3 , 2x , 1p1e12.3 , 2x ,
c     4        1p4e12.3 ,
c     5        1p3e12.3)
c*********************************************************************
      end if

