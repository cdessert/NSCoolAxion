      if (ifield.ne.0) then
      if (debug.ge.1.) print *,'Initializing B field'
       call get_dipole(imax,rad,
     1          f_gr_field,g_gr_field,h_gr_field,beta_gr_field,
     2          bfield0,stoke,bf_r,bf_t,bfield2,current)
      else
       do i=1,imax,2
        stoke(i)=1.d0
c        bf_r(i)=1./dsqrt(2.d0)*bfield0
c        bf_t(i)=1./dsqrt(2.d0)*bfield0
        bf_r(i)=bfield0
        bf_t(i)=bfield0
        bfield2(i)=bfield0
        current(i)=1.d0
       end do
      end if

      do i=1,imax,2
       ostoke(i)=stoke(i)
       obf_r(i)=bf_r(i)
       obf_t(i)=bf_t(i)
       obfield2(i)=bfield2(i)
       ocurrent(i)=current(i)
      end do

      do i=1,imax,2
       nstoke(i)=stoke(i)
       onstoke(i)=stoke(i)
       nbf_r(i)=bf_r(i)
       nbf_t(i)=bf_t(i)
       nbfield2(i)=bfield2(i)
       ncurrent(i)=current(i)
      end do

      bfield_energy=0.d0
      do i=1,imax,2
       bfield_energy=bfield_energy+
     1      bfield2(i)**2/8.d0/pi*(dvol(i-1)+dvol(i))
      end do
