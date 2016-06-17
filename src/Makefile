# Copyright 2016 Mischa Schirmer
#
# This file is part of NEBULAR.
#
# NEBULAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# NEBULAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEBULAR.  If not, see <http://www.gnu.org/licenses/>.

# Which platfrom are we running on?
platform = 
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	platform = $(UNAME)
else
ifeq ($(UNAME), Darwin)
	platform = $(UNAME)
endif
endif

# Exit if platform is not Linux or Darwin
ifndef platform
  $(error Unsupported platform "$(platform)". Must be Linux or Darwin)
endif

# Include and library dirs
INCLUDE_DIRS := /usr/include/ /usr/local/include/ ../include
LIBRARY_DIRS := /usr/lib /usr/local/lib
LIBRARIES    := gsl gslcblas

# For Darwin
#ifeq ($(os),Darwin)
#  INCLUDE_DIRS += /usr/X11R6/include/X11
#  LIBRARY_DIRS += /usr/X11R6/lib
#endif

# add the include and library (dirs) to the makefile variables
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LDFLAGS  += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir))
LDFLAGS  += $(foreach library,$(LIBRARIES),-l$(library))

# compiler flags
CXXFLAGS += -std=c++0x -Wall -Wextra -pedantic -Wundef -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -O3 -g

# For Darwin
#ifeq ($(os),Darwin)
#  CPPFLAGS := -iquote -I. $(CPPFLAGS)
#  # and for MacOs 10.9 or newer, which default to libc++, we need to do this
#  NEWFLAG = $(shell sw_vers -productVersion | awk 'BEGIN{FS="."} {if ($$1$$2>=109) print "yes"; else print "no"}')
#  ifeq ($(NEWFLAG),yes)
#    CXXFLAGS := -stdlib=libstdc++ $(CXXFLAGS)
#    LDFLAGS  := -stdlib=libstdc++ $(LDFLAGS)
#  endif
#endif

vpath %.h ../include
HEADERS=functions.h emissionline.h freefree.h freebound.h twophoton.h spectrum.h HHeMix.h
SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
OBJALL=$(filter-out nebular.o,$(OBJECTS))
EXEC=nebular

all: $(EXEC)

$(EXEC): $(OBJECTS) $(HEADERS)
	$(CXX) -o $@ $@.o $(OBJALL) $(LDFLAGS) $(CXXFLAGS)

.PHONY : all clean

clean :
	rm -f nebular *.o
