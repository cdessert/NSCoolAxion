      function emissivity(B,t,pf)
         Implicit None

         real*8 :: t,B,pf,emissivity
         integer indb,indt,indp
         integer,parameter :: h_dimb = 3
         integer,parameter :: h_dimp = 81
         integer,parameter :: h_dimt = 101
         real*8,dimension(1:h_dimb,1:h_dimp,1:h_dimt) :: h_bfield
         real*8,dimension(1:h_dimb,1:h_dimp,1:h_dimt) :: h_temp
         real*8,dimension(1:h_dimb,1:h_dimp,1:h_dimt) :: h_pfermi
         real*8,dimension(1:h_dimb,1:h_dimp,1:h_dimt) :: h_emissiv
         real*8 :: h_b0,h_b1,h_t0,h_t1,h_p0,h_p1
         common/h_tables/h_bfield,h_temp,h_pfermi,h_emissiv,
     1                h_b0,h_b1,h_t0,h_t1,h_p0,h_p1

         if ((B>1.5d13).and.(B<2.5d13)) then
            indb = 1
         else if ((B>2.5d13).and.(B<3d13)) then
            indb = 2
         else if ((B>6d13).and.(B<7d13)) then
            indb = 3
         end if
         indt = int((log10(t) -h_t0)/(h_t1-h_t0)*h_dimt + 1.5)
         indp = int((log10(pf)-h_p0)/(h_p1-h_p0)*h_dimp + 1.5)
         
         emissivity = h_emissiv(indb,indp,indt)
         return     
      end function



