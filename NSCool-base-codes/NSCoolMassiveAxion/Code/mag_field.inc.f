c      real*16 stoke,ostoke,nstoke,onstoke
      dimension bf_r(0:isize),obf_r(0:isize),nbf_r(0:isize)
      dimension bf_t(0:isize),obf_t(0:isize),nbf_t(0:isize)
      dimension bfield2(0:isize),obfield2(0:isize),nbfield2(0:isize)
      dimension current(0:isize),ocurrent(0:isize),ncurrent(0:isize)
      dimension stoke(0:isize),ostoke(0:isize),nstoke(0:isize),
     1          onstoke(0:isize)
      dimension force_r(0:isize),force_t(0:isize)
      dimension sigma(0:isize)
      common /b_field/
     4                stoke,ostoke,nstoke,onstoke,
     1                bf_r,obf_r,nbf_r,bf_t,obf_t,nbf_t,
     2                bfield2,obfield2,nbfield2,
     3                current,ocurrent,ncurrent,
     5                bfield0,start_b_diffusion,qimp,ifield

