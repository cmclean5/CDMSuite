# Aaron Clauset
# Makefile  Oct2003
# feel free to hack this to pieces

#### local macros
# remove without fussing about it
RM = /bin/rm -f

# compiler name and flags
CCC = g++
CCFLAGS = -O3 -fopenmp -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -Wno-deprecated -g3 

### local program information
EXEC=run
SOURCES=communityDetection.C 

### intermediate objects
OBJECTS = $(SOURCES: .cc=.o)

### libraries
LIBS += 

### configuration options
CONFIG_OPTIONS=USEOMP

SET_CONFIG_OPTIONS=$(filter-out -%,$(CONFIG_OPTIONS))
UNSET_CONFIG_OPTIONS=$(patsubst -%,%,$(filter -%,$(CONFIG_OPTIONS)))
ALL_CONFIG_OPTIONS=$(SET_CONFIG_OPTIONS) $(UNSET_CONFIG_OPTIONS)

### Turn config options into make variables.
$(foreach cfg,$(SET_CONFIG_OPTIONS),$(eval $(cfg)=1))
$(foreach cfg,$(UNSET_CONFIG_OPTIONS),$(eval $(cfg)=))

### Make sure none of the options are set to anything except 1 or blank.
### Using "make OPTION=0" doesn't work, since "0" is set, you need "make OPTION=".
$(foreach cfg,$(ALL_CONFIG_OPTIONS), \
    $(if $(patsubst %1,%,$(value $(cfg))), \
        $(error Use "$(cfg)=1" OR "$(cfg)=" not "$(cfg)=$(value $(cfg))")))

### Turn them into tool flags (-D).
TOOL_DEFINES+=$(foreach cfg,$(ALL_CONFIG_OPTIONS),$(if $(value $(cfg)),-D$(cfg),-U$(cfg)))

### Configuration options are turned into make variables above.
ifdef USEOMP
$(info Building with option-useomp)
endif

### headers
HEADERS = nr3.h ran.h ludcmp.h eigen_sym.h Headers.h edge.h node.h Helper.h readInputFile.h

### targets, dependencies and actions
$(EXEC): $(OBJECTS) Makefile
	$(LINK.cc) $(TOOL_DEFINES) $(CCFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

### sort out dependencies
depend:
	makedepend $(INCLUDES) $(HEADERS) $(SOURCES) $(LIBS)

clean: 
	$(RM) $(EXEC) *.o *~