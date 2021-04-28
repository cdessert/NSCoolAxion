      do i=1,imax,2
       ostoke(i)=stoke(i)
       stoke(i)=nstoke(i)
       obf_r(i)=bf_r(i)
       obf_t(i)=bf_t(i)
       bf_r(i)=nbf_r(i)
       bf_t(i)=nbf_t(i)
       obfield2(i)=bfield2(i)
       bfield2(i)=nbfield2(i)
       ocurrent(i)=current(i)
       current(i)=ncurrent(i)
      end do
