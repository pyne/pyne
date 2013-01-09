          ########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 

# Modify MOAB to point to the local MOAB installation
# MOAB = /data/opt/dagmc/moab
MOAB = /filespace/people/z/zachman/dagmc_bld/MOAB

# FLUPRO is an environment variable that should be set to the directory of the FLUKA installation
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp

# Local includes
FLUDAG_INCLUDES = -I./include

# MOAB_INCLUDES uses the ${MOAB} var to simplify compilation on different systems
MOAB_INCLUDES  = -I${MOAB}/include
MOAB_LIBDIR    =   ${MOAB}/lib

# MOAB_FLAGS determine what libraries we wish to include. Note: netcdf does not seem to be needed
MOAB_FLAGS = -lMOAB -ldagmc
MOAB_LIBS_LINK = -L${MOAB_LIBDIR} -lm ${MOAB_FLAGS} 


FSRCS = 

CCSRCS = 

OBJS = $(FSRCS:.f=.o) $(CCSRCS:.c=.o)
PRGS = $(FSRCS:.f=) $(CCSRCS:.c=)

.SUFFIXES:
.SUFFIXES: .f .c .o

.f.o:
	cd for
	gfortran -c $<
	mv -f $@ obj
	cd .. ;ls

.c.o: 
	g++ -c $<

CPPFLAGS = ${MOAB_INCLUDES} ${FLUDAG_INCLUDES}
CXXFLAGS = -g

default: all

all: mainDag mainf

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
mainf: obj/main_f.o  obj/fluka_funcs.o
	gfortran -o $@ obj/main_f.o  obj/fluka_funcs.o $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++

main_f.o: for/main_f.f
	cd for
	gfortran -c $<
	mv -f $@ obj
	cd .. ;ls

simple.o: for/simple.f 
	cd for
	gfortran -c $<
	mv -f simple.o obj
	cd ..; ls

test.o: c/test.c
	cd c
	gcc -c $<
	mv -f $@ obj
	cd ..; ls

fluka_funcs.o: cpp/fluka_funcs.cpp
	cd cpp
	g++ -c $(CPPFLAGS) $< 
	mv -f $@ obj
	cd ..; ls

dagmcinit.o: cpp/dagmcinit.cpp
	cd cpp
	g++ -c $(CPPFLAGS) $< 
	mv -f $@ obj
	cd ..; ls


clean:
	rm obj/*.o
	rm ./mainDag
	rm ./mainf

print-%:
	@echo $* = $($*)
