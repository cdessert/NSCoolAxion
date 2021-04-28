      if (ifield.ne.0) then
       if (ifield.eq.1.) extra=0.0d0
       if (ifield.eq.2.) extra=1.0d0
       do i=1,imax,2
         nstoke(i)  =stoke(i) +
     1               extra*(stoke(i)  -ostoke(i)  )*dtime/odtime
         onstoke(i)=nstoke(i)
       end do
       call get_dipole(imax,rad,
     1          f_gr_field,g_gr_field,h_gr_field,beta_gr_field,
     2          bfield0,nstoke,nbf_r,nbf_t,nbfield2,ncurrent)
      end if
