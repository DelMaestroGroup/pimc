#PIMC Makefile

CXX     ?= g++
LD      = $(CXX)
UNAME   = $(shell uname -s)

# Determine the comiler toolset, gcc or intel
ifeq ($(findstring g++,$(CXX)), g++)
TOOLSET = gcc
else ifeq ($(findstring icpc,$(CXX)), icpc)
TOOLSET = intel
else ifeq ($(findstring c++,$(CXX)), c++)
TOOLSET = clang
endif

#Number of dimensions to compile for
ndim = 3
DIM  = -D NDIM=$(ndim)

#Optimizations used: debug, basic and strict are valid
opts = basic

#If a user wants to override variables often, they can create
#there own section and declare a preset
preset   = none
OVERRIDE = $(preset)

# Get the svn version number into the executable
SVNDEV := -D'SVN_VERSION="$(shell svnversion -n .)"'

####################################################################
####################################################################
##OS Variables

###################################################################
#Linux
ifeq ($(UNAME), Linux)

codedir = $$HOME/local
CODEDIR = $(codedir)

#g++
ifeq ($(TOOLSET), gcc)
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
OPTS   = -Wall -fno-math-errno -O3

LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_system

#icpc
else ifeq ($(TOOLSET), intel)
DEBUG  = -D PIMC_DEBUG -debug -g
LDEBUG = -lblitz
OPTS   = -Wall -fast -fno-math-errno

LDFLAGS = -limf -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_system

endif #gcc, elseif intel
#Linux end#########################################################

###################################################################
#OS X
else ifeq ($(UNAME),Darwin)

codedir = /usr/local
CODEDIR = $(codedir)

#gcc
ifeq ($(TOOLSET), gcc)
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz

ifeq ($(opts), basic)
OPTS = -Wall -O3 -mtune=native  #-DNDEBUG
else ifeq ($(opts), strict)
OPTS = -Wall -O3 -W -Wshadow -fno-common -ansi -pedantic -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums
endif #basic, elseif strict

BOOSTVER ?= -gcc42-mt-1_49
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER)

#intel
else ifeq ($(TOOLSET), intel)
DEBUG  = -D PIMC_DEBUG -debug -g
LDEBUG = -lblitz -lboost
OPTS   = -Wall -fast -fno-math-errno

BOOSTVER ?= -il-mt-1_49
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER)

#clang
else ifeq ($(TOOLSET), clang)
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz

ifeq ($(opts), basic)
OPTS = -std=c++11 -stdlib=libc++ -Wall -O3 -mtune=native -Wno-deprecated-register
else ifeq ($(opts), strict)
OPTS = -Wall -O3 -W -Wshadow -fno-common -ansi -pedantic -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums
endif #basic, elseif strict

LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER)
#PROG = pimc_clang.e
endif #gcc, elseif intel elseif clang
#OS X end##########################################################

endif #linux, elseif osx
##OS Variables end##################################################

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include

####################################################################
####################################################################
##Overrides
ifeq ($(OVERRIDE), none)# skips remaining ifelse statements

###################################################################
#System Sharcnet
else ifeq ($(OVERRIDE), sharcnet)
OPTS = -Wall -TENV:simd_zmask=OFF -TENV:simd_imask=OFF -TENV:simd_omask=OFF -O3 -fno-math-errno
CODEDIR = /work/agdelma/local

ifeq ($(TOOLSET), intel)
LDFLAGS += -limf
OPTS += -vec-report0 -wd981
endif

CXXFLAGS  = $(OPTS) $(DIM) $(DEBUG) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_system
#Sharcnet end######################################################

###################################################################
#Westgrid
else ifeq ($(OVERRIDE), westgrid)
CXX = icpc
LD = icpc
OPTS = -O3 -ipo -axSSE4.1,axAVX  #-w -vec-report0 -opt-report0
CODEDIR = $$HOME/local
CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include # $(DEBUG)
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -limf
#Westgrid end######################################################

###################################################################
# System g++ macbook 
else ifeq ($(OVERRIDE), macbook)

CODEDIR = /usr/local
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
CXX = g++-mp-4.8
BOOSTVER = -gcc48-mt-1_55

ifeq ($(opts), basic)
OPTS = -Wall -O3 -std=c++11 -mtune=native  -Wno-unused-local-typedefs
else ifeq ($(opts), strict)
OPTS = -Wall -O3 -W -Wshadow -fno-common -ansi -pedantic -Wconversion -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums
endif #basic, elseif strict

CXXFLAGS  = $(OPTS) $(DIM) $(DEBUG) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER)
#macbook end
######################################################

endif # sharcnet, elseif westgrid, elseif macbook
##User Overrides end################################################

####################################################################
####################################################################
##Linking and Compiling Variables
RM     = /bin/rm -f
PROG  ?= pimc.e
SOURCE = pdrive.cpp pimc.cpp constants.cpp container.cpp path.cpp worm.cpp action.cpp potential.cpp move.cpp estimator.cpp lookuptable.cpp communicator.cpp setup.cpp wavefunction.cpp cmc.cpp
OBJS   = $(SOURCE:.cpp=.o)

COMPILE_PCH  = $(CXX) $(CXXFLAGS) $(SVNDEV) -pipe
COMPILE_WPCH = $(COMPILE_PCH) --include common.h
LINK         = $(LD) $(OBJS) -pipe $(LDFLAGS)

#Possible fix for older gccs, try without first
#COMPILE_WPCH += -fpch-preprocess

# -------------------------------------------------------------------------------
all: release

release: $(PROG)

debug: COMPILE_WPCH += $(DEBUG)
debug: LINK += $(LDEBUG)
debug: $(PROG)

# Link Objects
$(PROG): $(OBJS)
	$(LINK) -o $(PROG)

# Compile Objects
communicator.o: communicator.cpp common.h.gch constants.h
	$(COMPILE_WPCH) -c communicator.cpp

lookuptable.o: lookuptable.cpp common.h.gch constants.h container.h communicator.h path.h
	$(COMPILE_WPCH) -c lookuptable.cpp

constants.o: constants.cpp common.h.gch
	$(COMPILE_WPCH) -c constants.cpp

container.o: container.cpp common.h.gch constants.h
	$(COMPILE_WPCH) -c container.cpp

estimator.o: estimator.cpp common.h.gch path.h action.h potential.h communicator.h
	$(COMPILE_WPCH) -c estimator.cpp

potential.o: potential.cpp common.h.gch constants.h communicator.h path.h lookuptable.h
	$(COMPILE_WPCH) -c potential.cpp

pdrive.o: pdrive.cpp common.h.gch constants.h wavefunction.h container.h path.h potential.h action.h pimc.h lookuptable.h communicator.h setup.h cmc.h worm.h move.h
	$(COMPILE_WPCH) -c pdrive.cpp

wavefunction.o: wavefunction.cpp common.h.gch path.h 
	$(COMPILE_WPCH) -c wavefunction.cpp

action.o: action.cpp common.h.gch constants.h path.h potential.h lookuptable.h wavefunction.h
	$(COMPILE_WPCH) -c action.cpp

setup.o: setup.cpp common.h.gch constants.h communicator.h container.h potential.h action.h move.h estimator.h
	$(COMPILE_WPCH) -c setup.cpp

move.o: move.cpp common.h.gch path.h action.h lookuptable.h communicator.h
	$(COMPILE_WPCH) -c move.cpp

path.o: path.cpp common.h.gch constants.h container.h worm.h lookuptable.h communicator.h
	$(COMPILE_WPCH) -c path.cpp

pimc.o: pimc.cpp common.h.gch communicator.h estimator.h path.h lookuptable.h
	$(COMPILE_WPCH) -c pimc.cpp

worm.o: worm.cpp common.h.gch constants.h path.h
	$(COMPILE_WPCH) -c worm.cpp

cmc.o: cmc.cpp common.h.gch constants.h communicator.h potential.h container.h
	$(COMPILE_WPCH) -c cmc.cpp

# Precompile headers
#common.h.gch: common.h
#	$(COMPILE_PCH) -c common.cpp

# Precompile headers
common.h.gch: common.h
	$(COMPILE_PCH) -c common.h

# -------------------------------------------------------------------------------

clean:
	$(RM) $(PROG) $(OBJS) common.h.gch

print-%:
	@echo '$*=$($*)'
