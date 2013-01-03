          ########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 
# Since we will be testing on various systems it makes sense to make this a variable for easy changes
MOAB = /data/opt/dagmc/moab
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp
# MOAB_INCLUDES now uses the ${MOAB} var to simplify compilation on different systems
MOAB_INCLUDES  = -I${MOAB}/include
MOAB_LIBDIR    =   ${MOAB}/lib
# MOAB_FLAGS determine what libraries we wish to include (In my case I have no netcdf)
MOAB_FLAGS = -lMOAB  -ldagmc
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

CPPFLAGS = ${MOAB_INCLUDES}
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
mainf: obj/main_f.o obj/simple.o obj/test.o obj/fludagW.o
	gfortran -o $@ obj/main_f.o obj/simple.o obj/test.o obj/fludagW.o $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++

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

clean:
	rm obj/*.o
	rm ./mainDag
	rm ./mainf

print-%:
	@echo $* = $($*)
