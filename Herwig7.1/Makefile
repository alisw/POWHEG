#-*- Makefile -*-

######################################################
#This directory must be copied inside the process dir
#Remeber to set the analysis file
pwhg_analysis=pwhg_analysis.o 
#####################################################

COMPILER=gfortran
ff=$(COMPILER)
cc=gcc
cxx=g++

OBJ=../Herwig7.1/obj-$(COMPILER)
OBJDIR:=$(OBJ)

#Choose Analysis: dummy, process specific
## default analysis may require FASTJET package, that has to be installed separately (see below)

herwig-config=$(shell which herwig-config)
herwigdir=$(shell herwig-config --prefix)
libherwig=-L$(herwigdir)/lib/Herwig -l:HwShower.so

## set the HepMC and ThePEG installation directories
hepmcdir=$(herwigdir)
libhepmc=-L$(hepmcdir)/lib -lHepMC  

thepegdir=$(herwigdir)
libthepeg=-L$(thepegdir)/lib/ThePEG -lThePEG -l:LesHouches.so 

libextra=-lstdc++
librpath=-Wl,-rpath,$(hepmcdir)/lib:$(thepegdir)/lib/ThePEG

POWHEGDIR=$(PWD)/../../
POWHEGHERWIGDIR=$(PWD)


ifeq ("$(ff)","gfortran")	
ff=gfortran -fno-automatic -ffixed-line-length-none -fbounds-check 
endif

## For debugging uncomment the following
#DEBUG= -ggdb
ifdef DEBUG
OPTF=-O0  
OPTC=-O0  -Warray-bounds -Wno-div-by-zero -Waddress -Wall -W  -Wodr  -Wno-builtin-declaration-mismatch
else 
OPTF=-O2
OPTC=-O2
endif



INCLUDES=-I$(PWD) -I$(PWD)/.. -I$(hepmcdir)/include -I$(thepegdir)/include -I$(herwigdir)/include -I$(POWHEGHERWIGDIR) -I$(POWHEGDIR) -I$(POWHEGDIR)/include
FLAGSF=-fPIC $(OPTF) $(DEBUG) $(INCLUDES)
FLAGSC=-fPIC $(OPTC) $(DEBUG) $(INCLUDES) 
FF=$(ff) $(FLAGSF)
CC=$(cc)  $(FLAGSC)
CXX=$(cxx)  $(FLAGSC) -fpermissive
#-fpermissive
LIBS=$(libhepmc) $(libthepeg) $(libextra) $(librpath) $(libherwig)

VPATH=$(PWD):$(PWD)/..:$(POWHEGHERWIGDIR):$(POWHEGDIR)

FASTJET_CONFIG=$(shell which fastjet-config)
##WWarn: this fastjet MUST BE the one used to configure ThePeg!! 
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o multi_plot.o
PWHGANAL+=fastjetfortran.o pwhg_bookhist-multi-extra.o
PWHGANAL+=$(pwhg_analysis)




%.o: %.f | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.f90 | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJDIR) 
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc | $(OBJDIR)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS) -std=gnu++11

powhegAnalysis=powhegAnalysis.o herwig7_interface.o \
        boostrot.o bra_ket_subroutines.o\
	powheginput.o newunit.o utils.o \
        pwhg_io_interface.o $(PWHGANAL) $(FPEOBJ)

powhegReader=powhegLesHouchesFileReader.o opencount.o 

powhegHerwig=$(powhegAnalysis) $(powhegReader) 


powhegHerwig.so:$(powhegHerwig)
	$(FF) -shared $(patsubst %,$(OBJ)/%,$(powhegHerwig)) $(LIBSFASTJET) $(LIBS) -o $@ -std=gnu++11

clean:
	-rm -f $(OBJ)/*.o $(OBJ)/*.so *.so *.o

# target to generate object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

##########################################################################
