c *********************************************************************
c ************     Read Filenames for Physics    **********************
c *********************************************************************

      open(unit=21,file='Code/Data_Files/FILES_PHYS.dat',
     x     status='old')
       read(21,*)f_concryst
       read(21,*)f_opacity
       read(21,*)f_npbl
       read(21,*)f_npbc
       read(21,*)f_density
      close(unit=21,status='keep')

c *********************************************************************
c **************     Read Numerical Parameters    *********************
c *********************************************************************

      open(unit=22,file='Code/Data_Files/NUM_PARAM.dat',
     y     status='old')
       read(22,*)
       read(22,*)time0,timemax,istepmax,itrial_max,itrial_opt,tcut
       read(22,*)
       read(22,*)
       read(22,*)dtime,dtlimit,scale_dt0,scale_dt1,repeat
       read(22,*)
       read(22,*)
       read(22,*)istart,mratt,mratl,mrats,tvar,svar
       read(22,*)
      close(unit=22,status='keep')

      dtime=dtime*year
      time0=time0*year
      odtime=dtime
      scale_dt=scale_dt1

      tcon=1.e12

c *********************************************************************
c ********     Initialize the cooling calculation     *****************
c *********************************************************************
      open(unit=20,file=f_i,status='old')
       read(20,*)
       read(20,*)pscreen
       read(20,*)debug_keep
       debug=debug_keep
       read(20,*)istep_debug
       read(20,*)pteff
       read(20,*)ptemp
       read(20,*)pstar
       read(20,*)idump1
       read(20,*)idump2
       read(20,*)idump3
       read(20,*)tempmin
       read(20,*)tempini
       read(20,*)
       read(20,*)icvel_nodeg
       read(20,*)
       read(20,*)emnco
       read(20,*)emncr
       read(20,*)emp
       read(20,*)            
       read(20,*)p0
       period=p0
       read(20,*)  
       tprint(0)=0.0
       do i=1,50
        read(20,*,err=66,end=66)tprint(i)
       end do
 66    itpmax=i-1
      close(unit=20,status='keep')

c READ OUTER BOUNDARY PARAMETERS: ************************************
      if (debug.ge.1.) print *,'Opening I_Bound*.dat'
      open(unit=20,file=f_Bound,status='old')
       read(20,*)
       read(20,*)ifteff
       read(20,*)eta
       eta=eta/0.443d0
       eta_0=eta
       read(20,*)mag_coeff
       read(20,*)tb_acc0
       read(20,*)f_TbTs
      close(unit=20,status='keep')
      if (ifteff.eq.15) then
       print '(a40,1p1e10.1)',
     2    'Fixed outer boundary temperature T_b =',tb_acc0     
      end if

c READ PAIRING PARAMETERS: ********************************************
      if (debug.ge.1.) print *,'Opening I_Pairing*.dat'
      open(unit=20,file=f_Pairing,status='old')
       read(20,*)
       read(20,*)sfn1s0
       read(20,*)sfn3p2
       read(20,*)sfp1s0
       read(20,*)sfl1s0
       read(20,*)fn1s0
       read(20,*)fn3p2
       read(20,*)fp1s0
       read(20,*)fl1s0
       read(20,*)sfquark
       if (sfquark.eq.1.) read(20,*)f_quark_tc       ! File for Quark pairing
       if (sfn3p2.eq.150.) then
        read(20,*)
        read(20,*)kfmax_n3p2
        read(20,*)delkf_n3p2
        read(20,*)tcmax_n3p2
       end if
      close(unit=20,status='keep')

c READ NEUTRINO PARAMETERS: ********************************************
      if (debug.ge.1.) print *,'Opening I_Neutrino*.dat'
      open(unit=20,file=f_Neutrino,status='old')
       read(20,*)
       read(20,*)murca_increase
       read(20,*)inu_durca
       read(20,*)inu_eion
       read(20,*)inu_plasma
       read(20,*)inu_synch
       read(20,*)inu_n1s0_pbf
       read(20,*)inu_n3p2_pbf
       read(20,*)inu_p_pbf
       read(20,*)inu_bubble
       read(20,*)inu_photo
       read(20,*)inu_pair
       read(20,*)inu_nuts1       ! Not used so far
       read(20,*)inu_nuts2       ! Not used so far
       read(20,*)inu_nuts3       ! Not used so far
       read(20,*)inu_nuts4       ! Not used so far
       read(20,*)inu_nuts5       ! Not used so far
       read(20,*)
       read(20,*)rhoexo
       read(20,*)cexo
       read(20,*)pexo
       read(20,*)pexosn
       read(20,*)pexosp
       read(20,*)nonothing1      ! Not used so far
       read(20,*)nonothing2      ! Not used so far
       read(20,*)nonothing3      ! Not used so far
       read(20,*)nonothing4      ! Not used so far
      close(unit=20,status='keep')

c READ CONDUCTIVITY PARAMETERS: ********************************************
      if (debug.ge.1.) print *,'Opening I_Conduct*.dat'
      open(unit=20,file=f_Conduct,status='old')
       read(20,*)
       read(20,*)iopacity
       read(20,*)icon_crust
       read(20,*)icon_core
       read(20,*)iconnothing2      ! Not used so far
       read(20,*)iconnothing3      ! Not used so far
       read(20,*)iconnothing4      ! Not used so far
       read(20,*)iconnothing5      ! Not used so far
       read(20,*)connothing1       ! Not used so far
       read(20,*)connothing2       ! Not used so far
       read(20,*)qimp
       read(20,*)connothing4       ! Not used so far
       read(20,*)connothing5       ! Not used so far
      close(unit=20,status='keep')

c READ HEATING PARAMETERS: ******************************************
      if (debug.ge.1.) print *,'Opening I_Heat*.dat'
      open(unit=20,file=f_Heat,status='old')
       read(20,*)
c      Deep crustal heating:
       read(20,*)i_heat_deep_crust
c      Sudden heating:
       read(20,*)i_heat_deposit
       read(20,*)t_dep
       read(20,*)del_t_dep
       read(20,*)total_heat
       read(20,*)i_dep
c      Baryon --> SQM conersion heating:
       read(20,*)i_heat_convert           ! these values are overwritten by the ones
       read(20,*)MeV_neutron              ! read from the file I_Strange
c      Vortex creep heating:
       read(20,*)i_heat_vortex_creep
       read(20,*)j_44
c      Joule heating:
       read(20,*)i_heat_joule
c      Field decay heating:
       read(20,*)i_heat_field_decay
c      Heating parameters for future use:
       read(20,*)i_heat_cold1
       read(20,*)i_heat_cold2
       read(20,*)i_heat_cold3
       read(20,*)i_heat_cold4
       read(20,*)heat_cold1
       read(20,*)heat_cold2
       read(20,*)heat_cold3
       read(20,*)heat_cold4
      close(unit=20,status='keep')

c READ MAGNETIC FIELD PARAMETERS: **********************************
      if (debug.ge.1.) print *,'Opening I_Bfield*.dat'
      open(unit=20,file=f_Bfield,status='old')
       read(20,*)
       read(20,*)ifield
       read(20,*)i_gr_field
       read(20,*)bfield0
       read(20,*)i0
       read(20,*)i1
       read(20,*)nothing
       read(20,*)start_b_diffusion
       read(20,*)i_joule_heat
       read(20,*)i_conb
       read(20,*)i_eleconb
      close(unit=20,status='keep')

c READ ACCRETION PARAMETERS: **************************************
      if (debug.ge.1.) print *,'Opening I_Accretion*.dat'
      open(unit=20,file=f_Accretion,status='old')
       read(20,*)        
       read(20,*)i_acc       
       read(20,*)m_dot0
       read(20,*)t_acc0
       t_acc0=t_acc0*3.15576d7  ! convert year=365.25 days to seconds
       read(20,*)t_acc1
       t_acc1=t_acc1*3.15576d7 
       read(20,*)t_acc2
       t_acc2=t_acc2*3.15576d7 
       read(20,*)alpha_acc
       read(20,*)time_step_min
       time_step_min=time_step_min*3.15576d7 
       read(20,*)eta_Edd
       read(20,*)X_Edd
      close(unit=20,status='keep')
      eta=max(1.d-50,eta)
      eta0=max(1.d-50,eta0)

c READ STRANGE QUARK MATTER PARAMETERS: ***************************
      if (istrange.eq.1) then
       if (debug.ge.1.) print *,'Opening I_Stange*.dat'
       open(unit=20,file=f_Strange,status='old')
        read(20,*)c_nu_str          ! SQM neutrino emission
        read(20,*)p_nu_str          ! SQM neutrino emission
        read(20,*)c_con_str         ! SQM thermal conductivity
        read(20,*)p_con_str         ! SQM thermal conductivity
        read(20,*)c_cv_str          ! SQM specific heat
        read(20,*)i_heat_convert    ! Heating from baryon -> SQM conversion
        read(20,*)MeV_neutron       ! Heating from baryon -> SQM conversion
       close(unit=20,status='keep')
      else if (istrange.eq.0) then
        c_nu_str      =0.d0
        p_nu_str      =0.d0
        c_con_str     =0.d0
        p_con_str     =0.d0
        c_cv_str      =0.d0
        i_heat_convert=0.d0
        MeV_neutron   =0.d0
      else
       print *,'NSCool_READ: istrange not defined !'
       pause
      end if

c READ STAR STRUCTURE LAYOUT FILE: ********************************
      if (debug.ge.1.) print *,'Opening I_Structure*.dat'
      open(unit=20,file=f_Structure,status='old')
       read(20,*)rhocore,rhodrip,rhoenv,rhosurf
       read(20,*)icore,idec
      close(unit=20,status='keep')

