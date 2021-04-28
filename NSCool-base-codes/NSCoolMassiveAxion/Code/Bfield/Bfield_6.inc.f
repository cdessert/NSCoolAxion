c *********************************************************************
c ****** get nstoke in case of retarded Field evolution ***************
c *********************************************************************

      if (ifield.eq.1) then
      if (debug.ge.1.) print *,'Getting Nstoke(retarded evolution)'
       if ((time+dtime)/year.ge.start_b_diffusion) then
        do i=1,imax,2
         t=ntemp(i)/ephi(i)
         d=rrho(i)
         a=a_cell(i)
         a1=a_ion(i)
         z=z_ion(i)
         q=qimp
         print=0.d0
         if (istep.eq.(istep/50)*50) print=1.d0
         call conduct(i,t,d,a,a1,z,q,nbfield2(i),i_eleconb,
     1                sigma(i),lll,debug,
     2                nu_e_s,nu_e_l)
c         call elecon(i,t,d,a,a1,z,q,tcp(i),nbfield2(i),i_eleconb,
c     1               sigma(i),print)
        end do

        call solve_induction_dipole(imax,dtime,idt,rad,sigma,v_acc,
     1             f_gr_field,g_gr_field,h_gr_field,alpha_gr_field,
     2                              nstoke,stoke)
        call get_dipole(imax,rad,
     1           f_gr_field,g_gr_field,h_gr_field,beta_gr_field,
     2           bfield0,nstoke,nbf_r,nbf_t,nbfield2,ncurrent)
        do i=1,imax,2
         ostoke(i)  =stoke(i)
         obf_r(i) =bf_r(i)
         obf_t(i) =bf_t(i)
         obfield2(i)=bfield2(i)
         ocurrent(i)=current(i)
        end do
       end if
      end if

c Control the solve_induction DTIME cut `idt':
c With accretion:
c        if (m_dot.ne.0.) then
c         tau0=abs((rad(imax)-rad(imax-2))/v_acc(imax))
c         dtaa=tau0/1.d1
c         if (dtaa.ge.dtime) then
c          idt=1
c         else
c          idt=int(dtime/dtaa)
c         end if
c        end if
cc Without acretion:
c      if (m_dot.eq.0.d0) then
c       ichs=0
c       mdstoke=0.d0
c       do i=1,imax-2,2
c        mds=abs(stoke(i)-nstoke(i))/(abs(nstoke(i))+1.d-5)
c        if (mds.gt.mdstoke)then
c         mdstoke=mds
c         ichs=i
c        end if
c       end do
c       dsvar=svar-1.d0
c       dsto=mdstoke-1.d0
c       if(dsvar.lt.dsto)then
c        idt=int(dsto/dsvar)+1
c       else
c        idt=1
c       end if
c      end if
