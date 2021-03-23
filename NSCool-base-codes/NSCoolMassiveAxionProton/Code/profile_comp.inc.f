      dimension bar(0:isize),
     1          yneutr(0:isize),yprot(0:isize),
     2          yelect(0:isize),ymuon(0:isize),
     3          ylambda(0:isize),
     4          ysminus(0:isize),yszero(0:isize),ysplus(0:isize),
     5          yquarku(0:isize),yquarkd(0:isize),yquarks(0:isize),
     6          fhad(0:isize),
     7          theta_k(0:isize),theta_p(0:isize),
     8          a_cell(0:isize),a_ion(0:isize),z_ion(0:isize),
     9          v_ion(0:isize)
      common/profile_comp/bar,
     1                    yneutr,yprot,yelect,ymuon,
     2                    ylambda,ysminus,yszero,ysplus,
     3	 		  yquarku,yquarkd,yquarks,fhad,
     4                    theta_k,theta_p,
     4                    a_cell,a_ion,z_ion,v_ion,
     5                    istrange
 
c -----------------------------------------------------------------------
c theta_k and theta_p are the chiral angles for Kaon and Pion condensates
c
c z_ion  is the charge number of the nuclei
c a_ion  is the mass number of the nuclei
c a_cell is the number of nucleons per Wigner-Seitz cell 
c        (i.e., a_cell=a_ion + # dripped neutrons)
c -----------------------------------------------------------------------





