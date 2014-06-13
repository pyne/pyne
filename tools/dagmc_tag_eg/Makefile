include ${MOAB_DIR}/lib/moab.make

CXXFLAGS += ${MOAB_CXXFLAGS} -g 
LD_FLAGS = -g
CPPFLAGS += ${MOAB_INCLUDES} -g
CFLAGS   += ${MOAB_CFLAGS} -g


all: dagmc_add_density

clean:
	rm *.o dagmc_add_density

dagmc_add_density: dagmc_add_density.o 
	$(CXX) $(LD_FLAGS) -o dagmc_add_density $^ ${MOAB_LIBS_LINK}

