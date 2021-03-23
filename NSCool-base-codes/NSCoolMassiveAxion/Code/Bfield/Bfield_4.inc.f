c *********************************************************************
c ****** get nstoke in case of exact Field evolution ******************
c *********************************************************************

      if (ifield.eq.2) then
      if (debug.ge.1.) print *,'Getting NStoke (exact evolution)'
       if ((time+dtime)/year.ge.start_b_diffusion) then
        do i=1,imax,2
         t=ntemp(i)/ephi(i)
         d=rrho(i)
         a=a_cell(i)
         a1=a_ion(i)
         z=z_ion(i)
         q=qimp
         call conduct(i,t,d,a,a1,z,q,nbfield2(i),i_eleconb,
     1                sigma(i),lll,debug,
     2                nu_e_s,nu_e_l)
c         call elecon(i,t,d,a,a1,z,q,tcp(i),nbfield2(i),i_eleconb,
c     1               sigma(i),print)
        end do
        call solve_induction_dipole(imax,dtime,1,rad,sigma,v_acc,
     1             f_gr_field,g_gr_field,h_gr_field,alpha_gr_field,
     2                              nstoke,stoke)

        call get_dipole(imax,rad,
     1           f_gr_field,g_gr_field,h_gr_field,beta_gr_field,
     2           bfield0,nstoke,nbf_r,nbf_t,nbfield2,ncurrent)
       end if
      end if
