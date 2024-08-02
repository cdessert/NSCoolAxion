c GR coefficients for the dipole evolution:

      if (i_gr_field.eq.1) then
       f_gr_field(0)=ephi(0)
       g_gr_field(0)=1.
       h_gr_field(0)=0.
       do i=1,imax
        f_gr_field(i)=ephi(i)
        g_gr_field(i)=1.-2.d0*g*msol*emas(i)/rad(i)/c**2
        h_gr_field(i)=(2.d0*g*msol*emas(i)/rad(i)**2/c**2-
     1                4.d0*pi*g*rad(i)/c**2*(rrho(i)-pres(i)/c**2))
       end do
       x=rad(imax)*c**2/2./g/msol/emas(imax)
c For the outer boundary condition:
       alpha_gr_field=-x*(2.d0*x*log(1.d0-1.d0/x)+
     1                    (2.d0*x-1.d0)/(x-1.d0)  )/
     2                (x**2*log(1.-1./x)+x+0.5)
c For obtaining the surface field by extrapolation form the light cylinder,
c i.e., included GR effects outside the neutron star:
c       beta_gr_field=-3.d0*x**3*
c     1                 (log(1.d0-1.d0/x)+1.d0/2.d0/x*(2.d0+1.d0/x))
       beta_gr_field=1.d0
      else if (i_gr_field.eq.0) then               ! If no GR effects:
       do i=0,imax
        f_gr_field(i)=1.d0
        g_gr_field(i)=1.d0
        h_gr_field(i)=0.d0
       end do
       alpha_gr_field=1.d0
       beta_gr_field=1.d0
      else
       pause 'I_GR_FIELD improperly defined !'
      end if
