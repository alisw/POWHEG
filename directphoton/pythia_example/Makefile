ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC

PYTHIA8FLAGS=$(shell pythia8-config --cxxflags)
PYTHIA8LIBS=-L$(shell pythia8-config --libdir)  -L$(shell pythia8-config --libdir)/archive -lpythia8 -ldl -lstdc++

ROOT=$(shell root-config --cflags --libs)

example: example.cpp
	$(CXX) $(CXX_COMMON) -D _CPPPWHGHOOKS -o $@ $+ $(PYTHIA8FLAGS) $(PYTHIA8LIBS) $(ROOT)

example_QEDQCDHooks: example.cpp
	$(CXX) $(CXX_COMMON) -o $@  $+ $(PYTHIA8FLAGS) $(PYTHIA8LIBS) $(ROOT)

.PHONY : clean
clean:
	rm -f example
	rm -f example_QEDQCDHooks
