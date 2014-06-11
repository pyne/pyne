
libs += -lhdf5_hl
libs += -lhdf5
libs += -Li386-linux-gnu

source = main.cc pyne.cpp data.cpp nucname.cpp rxname.cpp

footest: $(source)
	g++ $(libs) -o $@ $^
