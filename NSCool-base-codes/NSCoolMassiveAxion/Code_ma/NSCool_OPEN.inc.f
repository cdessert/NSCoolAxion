c *********************************************************************
c PRINT OUT STAR PROPERTIES *******************************************
c *********************************************************************
      if (pstar.eq.1.) then
       open(unit=49,file=f_Star,status='new')
        itext=5
        write(49,'(4i10)')itext,imax,icore,idrip
        write(49,*)
        write(49,'(a79)')model
        write(49,*)
        write(49,
     1   '(a5,2a15,3a12,3x,
     2     8a12,3x,
     3     3a12,3x,
     4     a5,
     5     a5,3a7)')
     1     'i  ','rad  ','emas  ','rho  ','pres ','nb ',
     2     ' kf(e)','kf(mu)',' kf(p)',' kf(n)',
     2        'kf(la)','kf(S-)','kf(S0)','kf(S+)',
     3     'Tc(n) ','Tc(p) ','Tc(la) ',
     4     'Durca',
     5     'nSF',' Acel',' Aion',' Zion'
        write(49,*)
c       Print out core:
        do i=0,icore
         if (i.gt.isf)       what='  1s0'
         if (i.le.isf)       what='  3p2'
         if (tcn(i).eq.1.d0) what='  no '
         write(49,
     1    '(i5,1p2e15.6,1p3e12.3,3x,
     2      1p8e12.3,3x,
     3      1p3e12.3,3x,
     4      5i1,
     5      a5,3a7)')
     1      i,rad(i)/100.,emas(i),rrho(i),pres(i),bar(i),
     2      kfe(i),kfm(i),kfp(i),kfn(i),kfla(i),kfsm(i),kfs0(i),kfsp(i),
     3      tcn(i),tcp(i),tcla(i),
     4      idurca_np(i),
     4       idurca_lap(i),idurca_smn(i),idurca_smla(i),idurca_sms0(i),
     5      what,' N/A ',' N/A ',' N/A '
        end do
c       Print out inner crust:
        do i=icore+1,idrip
         if (i.gt.isf)       what='  1s0'
         if (i.le.isf)       what='  3p2'
         if (tcn(i).eq.1.d0) what='  no '
         write(49,
     1     '(i5,1p2e15.6,1p3e12.3,3x,
     2       1p1e12.3,24x,1p1e12.3,48x,3x,
     3       1p1e12.3,24x,3x,
     4       5x,       
     5       a5,0p3f7.1)')       
     1       i,rad(i)/100.,emas(i),rrho(i),pres(i),bar(i),
     2       kfe(i),kfn(i),
     3       tcn(i),
     5       what,a_cell(i),a_ion(i),z_ion(i)
        end do
c       Print out outer crust:
        do i=idrip+1,imax
         what=' '
         write(49,
     1     '(i5,1p2e15.6,1p3e12.3,3x,
     2       1p1e12.3,84x,3x,
     3       36x,3x,
     4       5x,
     5       a5,0p3f7.1)')
     1       i,rad(i)/100.,emas(i),rrho(i),pres(i),bar(i),
     2       kfe(i),
     5       what,a_cell(i),a_ion(i),z_ion(i)
        end do
       close(unit=49,status='keep')
      end if
c ***********************************************************************
c **** OPEN TEMP file: **************************************************
c ***********************************************************************
c     ptemp=1 : print out at times defined by TPRINT read from I_*.dat
c ***********************************************************************
      if (ptemp.eq.1) then
       if (debug.ge.1.) print *,'Opening Temp file'
       open(unit=10,file=f_Temp,status='new')
       write(10,*)'   Time in years,  Te_inf includes ephi'
       write(10,*)'   Temp and Lum are local values, i.e. no ephi'
       write(10,*)
       write(10,*)
      end if
c ***********************************************************************
c **** OPEN TEFF file: **************************************************
c ***********************************************************************
      if (pteff.ge.1.) then
      if (debug.ge.1.) print *,'Opening Teff file'
       rstar=rad(imax)
       mstar=emas(imax)
       rdurca=-1.d0
       mdurca=-1.d0
       rhodurca=-1.d0
       do i=0,icore
        if(idurca_np(i).eq.1)then
         rdurca=rad(i)
         mdurca=emas(i)
         rhodurca=rrho(i)
        end if
       end do
       if(jexo.gt.-1) then
        rexo=rad(jexo)
        mexo=emas(jexo)
       end if
       rcore=rad(icore)
       mcore=emas(icore)
       rincrust=rad(idrip)-rcore
       mincrust=emas(idrip)-mcore
       routcrust=rad(imax)-rad(idrip)
       moutcrust=pres(idrip)*4.d0*pi*(rad(idrip)/msol**.25)**4*
     1           root/g/emas(idrip)/msol

       open(unit=19,file=f_Teff,status='new')

       write(19,*)
       write(19,*)model
       write(19,*)
       write(19,1230)'Mstar     =',    mstar,'Rstar     =',    rstar
       if(rdurca.gt.0.)then
        write(19,1230)'Mdurca    =',   mdurca,'Rdurca    =',   rdurca
       end if
       if(jexo.gt.-1)then
        write(19,1230)'Mexo      =',     mexo,'Rexo      =',     rexo
       end if
       write(19,1230)'Mcore     =',    mcore,'Rcore     =',    rcore
       write(19,1230)'Mincrust  =', mincrust,'Rincrust  =', rincrust
       write(19,1230)'Moutcrust =',moutcrust,'Routcrust =',routcrust
       write(19,*)
1230   format(10x,2(1a15,1p1e10.3))

       write(19,9750)
     1   'icore=',icore,'sfn1s0=',sfn1s0,               'emnco=',emnco,
     2     'rhoexo=',rhoexo
       write(19,9751)
     1   'idrip=',idrip,'sfn3p2=',sfn3p2,'fn3p2=',fn3p2,'emncr=',emncr,
     2     '  cexo=',  cexo
       write(19,9751)
     1   ' imax=', imax,'sfp1s0=',sfp1s0,'fp1s0=',fp1s0,'  emp=',  emp,
     2     '  pexo=',  pexo
       write(19,*)
 9750  format(a8,i4,a10,0p1f3.0,8x,  5x    ,a10,0p1f3.0,a10,1p1e12.3)
 9751  format(a8,i4,a10,0p1f3.0,a8,0p1f05.2,a10,0p1f3.0,a10,1p1e12.3)
 9752  format(a8,i4,a10,0p1f3.0,8x,  5x    ,a10,0p1f3.0,a10,0p1f05.0)

       if (rdurca.le.0.)then
          write(19,9791)'DUrca not possible'
       else
          if (inu_durca.eq.0) then
             write(19,9791)'DUrca possible but TURNED OFF'
          else 
             write(19,9792)'DUrca active above Rho =',rhodurca,
     2                     'gm/cm^3'
          end if
       end if
       write(19,9790)'e-ion bremstrahlung:',inu_eion,
     2                  'bubble neutrinos:',inu_bubble
       write(19,9790)'plasma neutrinos:',inu_plasma,
     2               ' neutron 1S0 pbf:',inu_n1s0_pbf
       write(19,9790)'  pair neutrinos:',inu_pair,
     2               ' neutron 3P2 pbf:',inu_n3p2_pbf
       write(19,9790)' photo neutrinos:',inu_photo,
     2               '  proton 1S0 pbf:',inu_p_pbf
       write(19,*)
 9790  format (2(1a30,1i3))
 9791  format(a30)
 9792  format(a30,1p1e12.3,a8)

       if (ifteff.eq.0) then
        write(19,'(a30,a40)')'Envelope from file: ',f_TbTs
       else if (ifteff.eq.1) then
        write(19,'(a40)')' Iron envelope from Gundmundsson et al'
       else if (ifteff.eq.2) then
        write(19,'(a40)')' Iron envelope from Nomoto & Tsuruta'
       else if (ifteff.eq.3) then
        etap=eta
        if (eta.le.1.d-40) etap=0.0
        write(19,'(a58,1p1e8.2)')
     2    'Accreted envelope from Potekhin et al. with Eta = ',etap
       else if (ifteff.eq.10) then
        write(19,'(a40,1p1e10.1,a3)')
     2    ' Magnetized Iron envelope with',bf_r(imax),' G'
       else if (ifteff.eq.11) then
        write(19,'(a40,1p1e10.1,a3)')
     2    ' Magnetized Accreted envelope with',bf_r(imax),' G'
       else if (ifteff.eq.15) then
        write(19,'(a40,1p1e10.1)')
     2    'Fixed outer boundary temperature T_b =',tb_acc0     
       end if
       write(19,*)
c ***********************************************************************
       write(19,*)
       if (pteff.eq.1.0) then
        write(19,1231)
     1   'Step',' Time  ',' Teff   ',
     2   ' L_phot  ',' L_nu   ',' L_heat '
        write(19,1231)
     1    '   ','[years]',' at inf [K] ',
     2    '[erg/sec]','[erg/sec]','[erg/sec]'
        write(19,*)
 1231   format(a8,a12,5x,4a12)
       else
        write(19,1235)
     1  'step','Time      ','Teff   ',
     2  'B field ','Period ',
     3  'dt  ','dt/odt','dtemp','itrial',
     4  'Temp(i=',idump1,')','Temp(i=',idump2,')','Temp(i=',idump3,')',
     5  'Cycle','M_dot  ',
     6  'Heat  ','L_H   ','L_nu(core)','L_nu(crust',
     7  'Cv_core ','Cv_crust '
1235    format(
     1   1a8 , 1a22 , 1a16 ,
     2   1a12 , 1a12 ,
     3   1a10 , 2a8 , 1a8 ,
     4   3(a12,i3,a1) ,
     5   1a8 , 1a14 ,
     6   4a12,
     7   2a12    )
        write(19,1236)
     1   ' ','[years]     ','at inf [K] ',
     2   '[G]   ','[sec.] ',
     3   ' ',' ',' ',' ',   
     4   rrho(idump1),rrho(idump2),rrho(idump3),
     5   ' ','[Msun/yr]',
     6   '[erg] ','[erg/sec]','[erg/sec]','[erg/sec]'
 1236   format(
     1   1a8 , 1a22 , 1a16 ,
     2   1a12 , 1a12 ,
     3   1a10 , 2a8 , 1a8 ,
     4   1p3e16.3 ,
     5   1a8 , 1a14 ,
     6   4a12 )
        write(19,*)
       end if
      end if
c ***********************************************************************
c ***********************************************************************
c ***********************************************************************
