      program Plot_Teff

      implicit real*8 (a-h,k-z)
      real*4 xplot,yplot
      parameter (year=3.15576e7)
      character*60 filename
      character*2 next
      dimension xplot(100000),yplot(100000)
c ****** draw the box *******************************************    
123   write(6,*)
     1    'Screen (1) or Laser printer (2)'
      read(5,*)output
      if (output.eq.1.) then
       call pgbegin (0,'/XSERV',1,1)
      else if (output.eq.2.) then
       call pgbegin (0,'Plot_Teff.ps/ps',1,1)
      else
       goto 123
      end if      
c***********************************************
      call pgslw(3)
      call pgsls(1)
      call pgvsize (2.0,6.0,2.0,6.0)
      call pgwindow (0.0,7.0,5.0,6.7)
      call pgsch(1.5)
      call pgbox ('bcnts',1.0,5,'bcntsv',0.5,5)
      call pglabel ('Log t [yrs]',' ',' ')
      call pgmtext ('l',3.,0.5,0.5,'Log T\de\u\u\(2270)\d [K]')
      call pgsch(1.0)
c ****** plot the curves ********************************
       goto 110
 100   print *,'I can''t open the file ! Try again.'
 110   write(6,*)'filename ='
       read(5,*,err=110)filename
 111   write(6,*)'Line style ='
       read(5,*,err=111)isls

       open(unit=20,file=filename,status='old',err=100)
        do i=1,30
         read(20,*)
        end do
        do i=1,100000
         read(20,*,err=90,end=90)ijunk,time,teff,l_phot,l_nu
         xplot(i)=log10(time)
         yplot(i)=log10(teff)
        end do
 90     continue
       close(unit=20,status='keep')
       imax=i-1
       call pgsls(isls)
       call pgline (imax,xplot,yplot)
c **********************************************************************
 900   write(6,*)'Another one ? [y/n]'
      read(5,*)next
      if (next.eq.'y') then
       goto 110
      else if (next.eq.'n') then
       goto 999
      else
       print *,'I didn''t get it ! Tell me again.'
       goto 900
      end if
c **********************************************************************
 999  continue
      call pgend

      end

c **********************************************************************
c **********************************************************************
