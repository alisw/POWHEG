#-*- Makefile -*-
## COMPILER: gfortran, ifort
COMPILER=gfortran
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf6 # for pPb
#PDF=lhapdf
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
#ANALYSIS=default
ANALYSIS=none
## For static linking uncomment the following
#STATIC= -static
# set to yes to enable compression
WITHZLIB=yes
#
ifeq ("$(COMPILER)","gfortran")	
FMISMATCH=
GFORTRAN_VERION = $(shell gfortran --version | head -1 | rev | cut -d " " -f1 | rev)
version_tuple := $(subst ., ,$(GFORTRAN_VERION:%=%))
GFORTRAN_VERION_MAJOR := $(word 1,$(version_tuple))
GFORTRAN_VERION_MINOR := $(word 2,$(version_tuple))
GFORTRAN_VERION_PATCH := $(word 3,$(version_tuple))
ifeq ($(shell test $(GFORTRAN_VERION_MAJOR) -ge 10; echo $$?), 0)
FMISMATCH=-fallow-argument-mismatch
endif
F77= gfortran -fno-automatic $(FMISMATCH)
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
##FFLAGS= -Wall -Wimplicit-interface
FFLAGS= -ffixed-line-length-132 -cpp
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -pg
ifdef DEBUG
OPT=-O0
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
endif
endif


ifeq ("$(COMPILER)","ifort")
F77 = ifort -save  -extend_source
CXX = icpc
LIBS = -limf
FFLAGS =  -check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
endif

SYSTEM=$(shell uname -s)
STDLIBCPP=-lstdc++
ifeq ("$(SYSTEM)", "Darwin")
STDLIBCPP=-lc++
endif

OBJ=obj-$(COMPILER)



PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:$(OBJ)/

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include

ifeq ("$(WITHZLIB)","yes")
LIBS+=-lz
endif

ifeq ("$(WITHZLIB)","dummy")
LIBZDUMMY=zlibdummy.o
endif


## LoopTools
LTAPP=$(shell which lt)
LT=$(subst bin/lt,,$(LTAPP))
LTFLAGS=-I$(LT)/include -L$(LT)/lib64 -looptools
LIBS+= -L$(LT)/lib64 -looptools

FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE2) $(LTFLAGS)


INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

## hendrik
#ifeq ("$(PDF)","lhapdf6")
LHAPDF_CONFIG=lhapdf-config
FJCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
PDFPACK=lhapdf6if.o lhapdf6ifcc.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF $(STDLIBCPP)
LIBS+=$(LIBSLHAPDF)
#endif


#ifeq ("$(PDF)","lhapdf")
#LHAPDF_CONFIG=lhapdf-config
#PDFPACK=lhapdfif.o
#LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
#ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH=/usr/lib/gcc/x86_64-redhat-linux/4.4
# LIBSTDCPP=/usr/lib/gcc/i386-redhat-linux/4.3.0/
# LIBSLHAPDF +=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
#endif
#LIBS+=$(LIBSLHAPDF)
#else
#PDFPACK=mlmpdfif.o hvqpdfpho.o
#endif

PWHGANAL=pwhg_analysis-dummy.o  pwhg_bookhist-multi.o 

ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STDLIBCPP)
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis.o pwhg_bookhist-multi.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

########## 
# This section is needed to compile and link routines relevant for madgraph:
# add libdhelas3.2a to the list of linked libraries:

# clear and reassign VPATH (needs also directories relevant to link madgraph):
vpath
VPATH=./:../:$(OBJ)/

%.o: %.f $(INCLUDE) | $(OBJ)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJ)
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc | $(OBJ)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o $(PWHGANAL) 

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	reshufflemoms.o setlocalscales.o validflav.o mint_upb.o  \
	pwhgreweight.o opencount.o ubprojections.o \
	sigcollsoft.o sigvirtual.o pwhg_io_interface.o rwl_weightlists.o\
	rwl_setup_param_weights.o lhefread.o printrwgtsuda.o $(PDFPACK) $(USER) $(FPEOBJ)

# target to generate LHEF output

$(shell ../svnversion/svnversion.sh>/dev/null)

pwhg_main: $(PWHG)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to analyze LHEF output
LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o	      \
	opencount.o powheginput.o pwhg_bookhist-multi.o lhefread.o  \
	newunit.o pwhg_analysis_driver.o $(FPEOBJ)            \
	$(PWHGANAL) 

lhef_analysis: $(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o powheginput.o \
	pwhg_bookhist-multi.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL)

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG)) $(LIBSFASTJET) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	pwhg_bookhist-multi.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBSFASTJET) $(STATIC) -o $@

# target to cleanup
.PHONY: clean
clean:
	rm -f $(OBJ)/*.o pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef 

$(OBJ):
	mkdir -p $(OBJ)
