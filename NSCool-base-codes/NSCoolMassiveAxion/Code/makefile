#--------------------------------------------------------------------------------
# For the Intel compiler ifort:
#COMPILER = ifort
#OPTIONS  = -O3 -m64
#OPTIONS  = -O3
#PLOT_LIB = -L/opt/pgplot/pgplot_ifort -lpgplot -L/usr/X11R6/lib -lX11
#PLOT_LIB = -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11
#--------------------------------------------------------------------------------
# For the GNU compiler gfortran:
#COMPILER = /sw/bin/gfortran
#OPTIONS  = -O3 -m64
#OPTIONS  = -O3
# The -Wall option prints out many "Warnings", good for debugging:
#OPTIONS  = -O3 -Wall
#PLOT_LIB = -L/sw/lib/pgplot -lpgplot -L/usr/X11R6/lib -lX11
#PLOT_LIB = -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -lcc_dynamic
#--------------------------------------------------------------------------------
# For the GNU compiler g95:
#COMPILER = /opt/local/bin/g95
COMPILER = gfortran -g
OPTIONS  = -fno-backslash -O3

#OPTIONS  = -fno-backslash -O3 -ftrace=full -Wall
#PLOT_LIB = -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
DEPENDS  = precool.o emissivity.o\
           conductivity.o conductivity_core.o conductivity_crust.o \
           opacity.o \
           neutrino.o neutrino_core.o neutrino_crust.o \
           heating.o spec_heat.o density.o boundary.o tc.o tc_Ioffe.o\
           eff_mass.o rotation.o accretion.o magnetic.o Tools.o
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# COOLING CODE: NSCool

NSCool.out: NSCool.o $(DEPENDS)
	$(COMPILER) $(OPTIONS) -o NSCool.out NSCool.o $(DEPENDS) $(PLOT_LIB)
	cp NSCool.out ..

NSCool.o: NSCool.f NSCool_READ.inc.f NSCool_OPEN.inc.f NSCool_PRINT.inc.f NSCool_INIT_TPROF.inc.f
	$(COMPILER) $(OPTIONS) -c NSCool.f
precool.o: precool.f
	$(COMPILER) $(OPTIONS) -c precool.f
emissivity.o: emissivity.f
	$(COMPILER) $(OPTIONS) -c emissivity.f
conductivity.o: conductivity.f
	$(COMPILER) $(OPTIONS) -c conductivity.f
conductivity_core.o: conductivity_core.f
	$(COMPILER) $(OPTIONS) -c conductivity_core.f
conductivity_crust.o: conductivity_crust.f
	$(COMPILER) $(OPTIONS) -c conductivity_crust.f
opacity.o: opacity.f
	$(COMPILER) $(OPTIONS) -c opacity.f
neutrino.o: neutrino.f
	$(COMPILER) $(OPTIONS) -c neutrino.f
neutrino_core.o: neutrino_core.f
	$(COMPILER) $(OPTIONS) -c neutrino_core.f
neutrino_crust.o: neutrino_crust.f
	$(COMPILER) $(OPTIONS) -c neutrino_crust.f
heating.o: heating.f
	$(COMPILER) $(OPTIONS) -c heating.f
spec_heat.o: spec_heat.f
	$(COMPILER) $(OPTIONS) -c spec_heat.f
density.o: density.f
	$(COMPILER) $(OPTIONS) -c density.f
boundary.o: boundary.f
	$(COMPILER) $(OPTIONS) -c boundary.f
tc.o: tc.f
	$(COMPILER) $(OPTIONS) -c tc.f
tc_Ioffe.o: tc_Ioffe.f
	$(COMPILER) $(OPTIONS) -c tc_Ioffe.f
eff_mass.o: eff_mass.f
	$(COMPILER) $(OPTIONS) -c eff_mass.f
rotation.o: rotation.f
	$(COMPILER) $(OPTIONS) -c rotation.f
accretion.o: accretion.f
	$(COMPILER) $(OPTIONS) -c accretion.f
magnetic.o: magnetic.f
	$(COMPILER) $(OPTIONS) -c magnetic.f
Tools.o: Tools.f
	$(COMPILER) $(OPTIONS) -c Tools.f
