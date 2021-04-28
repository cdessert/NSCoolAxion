      dimension Tc_CFL(0:isize),Tc_3SC(0:isize),
     1 Tc_2SC(0:isize),Tc_u(0:isize),Tc_d(0:isize),Tc_s(0:isize)

      common/superfluid_quark/Tc_CFL,Tc_3SC,
     1                  Tc_2SC,Tc_u,Tc_d,Tc_s,
     2                  nb_3SC,nb_CFL,normal,
     3                  type_CFL,type_3SC,type_2SC,
     4                  type_u,type_d,type_s

c THREE DIFFERENT PAIRING PATTERNS:
c
c CFL PHASE: at nb >= nb_CFL
c Tc_CFL applies to u,d & s quarks of all colors
c No electrons are present
c
c 3SC PHASE: at nb_3SC <= nb < nb_CFL
c Tc_3SC applies to u,d & s quarks of all colors
c
c 2SC PHASE: at nb < nb_3SC
c Tc_2SC applies to u & d quarks of color 1 & 2
c Tc_u applies to u quarks of color 3
c Tc_d applies to d quarks of color 3
c Tc_s applies to s quarks of all colors
c
c If normal=1 then there is no pairing at all !
c
c


