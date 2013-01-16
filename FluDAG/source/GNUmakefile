          ########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 

# Modify MOAB to point to the local MOAB installation
# MOAB = /data/opt/dagmc/moab
MOAB = ${HOME}/dagmc_bld/MOAB
FLUDAG = ${HOME}/DAGMC/FluDAG
FLUDAG_SOURCE = ${FLUDAG}/source

# FLUPRO is an environment variable that should be set to the directory of the FLUKA installation
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp

# Local includes
FLUDAG_INCLUDES = -I${FLUDAG_SOURCE}/include

# MOAB_INCLUDES uses the ${MOAB} var to simplify compilation on different systems
MOAB_INCLUDES  = -I${MOAB}/include
MOAB_LIBDIR    =   ${MOAB}/lib

# MOAB_FLAGS determine what libraries we wish to include. Note: netcdf does not seem to be needed
MOAB_FLAGS = -lMOAB -ldagmc
MOAB_LIBS_LINK = -L${MOAB_LIBDIR} -lm ${MOAB_FLAGS} 


CPPFLAGS = ${MOAB_INCLUDES} ${FLUDAG_INCLUDES}
# CXXFLAGS = -g
# Modify this to the name of the directory from which myfluka is run
RUNDIR = ./userini_example

default: all

all: myfluka

# Create a cpp main that links in fluka libraries to get flukam_
mainDag: obj/mainDag.o obj/fludagW.o
	$(CXX) -o $@ $< $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lgfortran
        
mainDag.o : cpp/mainDag.cpp
	cd cpp
	g++ -c $< $(CPPFLAGS) 
	mv -f $@ obj
	cd ..

fludagW.o : cpp/fludagW.cpp
	cd cpp
	g++ -c $(CPPFLAGS) $<
	mv -f $@ obj
	cd ..
 
# Create a fortrain main linking object files compiled from fortran, C and C++
mainf: 
	cd obj; make main_f.o fluka_funcs.o; cd ..
	gfortran -o $@ obj/main_f.o  obj/fluka_funcs.o $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++

myfluka: obj/userini.o obj/fluka_funcs.o
	${FLUPRO}/flutil/lfluka -o myfluka -a MOAB -a dagmc -m fluka obj/userini.o -L ${FLUPRO} -L ${MOAB_LIBDIR} obj/fluka_funcs.o
	cp myfluka ${RUNDIR}

userini.o : for/userini.f
	cd obj; make userini.o

fluka_funcs.o : cpp/fluka_funcs.o
	cd obj; make fluka_funcs.o

obj: 
	make -f obj/GNUmakefile all

test.stl : test.h5m
	mbconvert $< $@

clean:
	rm obj/*.o
	rm ./myfluka

print-%:
	@echo $* = $($*)
