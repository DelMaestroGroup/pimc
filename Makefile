#PIMC Makefile

$(info System compiler $$CXX is set to ${CXX})

CXX     ?= c++
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

$(info System $$TOOLSET is set to ${TOOLSET})


#Number of dimensions to compile for
ndim = 3
DIM  = -D NDIM=$(ndim)

#Optimizations used: debug, basic and strict are valid
opts ?= basic

#If a user wants to override variables often, they can create
#there own section and declare a preset
preset  ?= none

ifeq (,$(findstring none,$(preset)))
$(info You have specified $$preset ${preset})
endif

# Get the svn version number into the executable
SVNDEV := -D'SVN_VERSION="$(shell svnversion -n .)"'

####################################################################
####################################################################
## Preset Overrides
ifneq ($(preset), none) # skips remaining ifelse statements

###################################################################
#Westgrid
ifeq ($(preset), westgrid)
CXX = icpc
LD = icpc
OPTS = -O3 -ipo -axSSE4.1,axAVX  #-w -vec-report0 -opt-report0
CODEDIR = $$HOME/local
CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include # $(DEBUG)
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -limf
endif # westgrid
#Westgrid end######################################################

###################################################################
# System g++ 801 macbook 
ifeq ($(preset), 801)
CODEDIR = $$HOME/local
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
CXX = /opt/local/bin/g++-mp-6
BOOSTVER = 1_68 
BOOSTCOMP = gcc65-mt-x64

ifeq ($(opts), basic)
OPTS = -Wall -O3 -std=c++17 -mtune=native -Wno-unused-local-typedefs
else ifeq ($(opts), strict)
# OPTS = -std=c++17 -Wall -g -W -Wextra -Wshadow -fno-common -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fsanitize=address -fno-omit-frame-pointer -Wconversion -Wno-c++11-extensions -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-unused-parameter -Wno-ambiguous-member-template 
OPTS = -std=c++17 -Wall -Wextra -g -pedantic
endif #basic, elseif strict

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options-$(BOOSTCOMP)-$(BOOSTVER) -lboost_filesystem-$(BOOSTCOMP)-$(BOOSTVER) 
endif # 801 macbook end
######################################################

###################################################################
# System c++ 801p macbook 
ifeq ($(preset), 801p)

CODEDIR = $$HOME/local
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
BOOSTVER = 

ifeq ($(opts), basic)
OPTS = -std=c++17 -Wall -O3 -mtune=native -Wno-deprecated-declarations
else ifeq ($(opts), strict)
OPTS = -std=c++17 -Wall -Wextra -g -pedantic
endif #basic, elseif strict

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER)
endif # 801p macbook 
######################################################

###################################################################
# System c++ 801p macbook 
ifeq ($(preset), 801ps)

CODEDIR = $$HOME/local
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
BOOST = 1_73
BOOSTVER = -xgcc42-mt-x64-$(BOOST)

ifeq ($(opts), basic)
OPTS = -std=c++17 -Wall -O3 -mtune=native -Wno-deprecated-declarations
else ifeq ($(opts), strict)
OPTS = -std=c++17 -Wall -Wextra -g -pedantic
endif #basic, elseif strict

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include -I$(CODEDIR)/include/boost-$(BOOST)
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER) -lboost_serialization$(BOOSTVER)
endif # 801ps macbook 
######################################################

###################################################################
# System g++ local_target
ifeq ($(preset), local_target)

CODEDIR = $$HOME/local
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
BOOSTVER = 

ifeq ($(opts), basic)
OPTS = -std=c++17 -Wall -O3 -mtune=native -Wno-deprecated-declarations
else ifeq ($(opts), strict)
OPTS = -std=c++17 -Wall -Wextra -g -pedantic
endif #basic, elseif strict

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER) -lboost_serialization$(BOOSTVER)
endif # local_target 
######################################################

###################################################################
# System VACC
ifeq ($(preset), vacc)

CODEDIR = $$HOME/local

OPTS = -std=c++17 -Wall -O3 -mtune=native -Wno-deprecated-declarations 

BOOSTVER = -gcc73-mt-x64-1_68
DEBUG = -D PIMC_DEBUG -g
LDEBUG = -lblitz

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER)  -lboost_system$(BOOSTVER) 
endif # vacc
######################################################

##User Overrides end################################################

else

####################################################################
# Here we try to guess values based on some enironment variables
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
OPTS   = -Wall -fno-math-errno -O3 -std=c++17

LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_serialization

#icpc
else ifeq ($(TOOLSET), intel)
DEBUG  = -D PIMC_DEBUG -debug -g
LDEBUG = -lblitz
OPTS   = -Wall -fast -fno-math-errno

LDFLAGS = -limf -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem
endif #gcc, elseif intel
#Linux end#########################################################

###################################################################
#OS X
else ifeq ($(UNAME),Darwin)

codedir = $$HOME/local
CODEDIR = $(codedir)

#gcc
ifeq ($(TOOLSET), gcc)
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz

ifeq ($(opts), basic)
OPTS = -std=c++17 -Wall -O3 -mtune=native -Wno-deprecated-declarations #-Wshadow  #-DNDEBUG
else ifeq ($(opts), strict)
# OPTS = -std=c++11 -Wall -g -W -Wextra -Wshadow -fno-common -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fsanitize=address -fno-omit-frame-pointer -Wconversion -Wno-c++11-extensions -Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-unused-parameter -Wno-ambiguous-member-template 
OPTS = -std=c++17 -Wall -Wextra -g -pedantic
endif #basic, elseif strict

BOOSTVER ?= -gcc55-mt-x64-1_67 
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER)

#intel
else ifeq ($(TOOLSET), intel)
DEBUG  = -D PIMC_DEBUG -debug -g
LDEBUG = -lblitz -lboost
OPTS   = -Wall -fast -fno-math-errno

BOOSTVER ?= -il-mt-1_49
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER)

#clang
else ifeq ($(TOOLSET), clang)
DEBUG  = -D PIMC_DEBUG -g
LDEBUG = -lblitz
BOOSTVER ?= -xgcc42-mt-x64-1_68

ifeq ($(opts), basic)
OPTS = -std=c++17 -stdlib=libc++ -Wall -O3 -Wno-deprecated-register -Wshadow -Wno-parentheses
else ifeq ($(opts), strict)
# OPTS = -Wall -std=c++11 -stdlib=libc++ -O3 -W -Wshadow -fno-common -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fsanitize=address -fno-omit-frame-pointer -Wconversion -Wno-c++11-extensions -Wno-shorten-64-to-32 -ferror-limit=1000 -Wno-sign-conversion -Wno-unused-parameter -Wno-ambiguous-member-template
OPTS = -std=c++17 -stdlib=libc++ -W -Wall -Wextra -g -pedantic
else ifeq ($(opts),debug)
OPTS = -std=c++17 -stdlib=libc++ -Wall -Wno-deprecated-register
endif #basic, elseif strict
#-fsanitize=address -fno-omit-frame-pointer
#-Wconversion

LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER)
endif #gcc, elseif intel elseif clang
#OS X end##########################################################

endif #linux, elseif osx
##OS Variables end##################################################

CXXFLAGS  = $(OPTS) $(DIM) -I$(CODEDIR)/include

endif # preset else


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

debug: COMPILE_PCH += $(DEBUG)
debug: LINK += $(LDEBUG)
debug: $(PROG)

pigs: COMPILE_PCH += -DPIGS
pigs: PROG = pigs.e
pigs: $(PROG)

pigsdebug: COMPILE_PCH += -DPIGS $(DEBUG)
pigsdebug: LINK += $(LDEBUG)
pigsdebug: PROG = pigs.e
pigsdebug: $(PROG)

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

estimator.o: estimator.cpp common.h.gch path.h action.h potential.h communicator.h factory.h
	$(COMPILE_WPCH) -c estimator.cpp

potential.o: potential.cpp common.h.gch constants.h communicator.h path.h lookuptable.h
	$(COMPILE_WPCH) -c potential.cpp

pdrive.o: pdrive.cpp common.h.gch constants.h wavefunction.h container.h path.h potential.h action.h pimc.h lookuptable.h communicator.h setup.h cmc.h worm.h move.h
	$(COMPILE_WPCH) -c pdrive.cpp

wavefunction.o: wavefunction.cpp common.h.gch path.h 
	$(COMPILE_WPCH) -c wavefunction.cpp

action.o: action.cpp common.h.gch constants.h path.h potential.h lookuptable.h wavefunction.h
	$(COMPILE_WPCH) -c action.cpp

setup.o: setup.cpp common.h.gch constants.h communicator.h container.h potential.h action.h move.h estimator.h factory.h
	$(COMPILE_WPCH) -c setup.cpp

move.o: move.cpp common.h.gch path.h action.h lookuptable.h communicator.h factory.h
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
common.h.gch: common.h
	$(COMPILE_PCH) -c common.h

# -------------------------------------------------------------------------------

clean:
	$(RM) pimc.e pigs.e $(OBJS) common.h.gch

print-%:
	@echo '$*=$($*)'
