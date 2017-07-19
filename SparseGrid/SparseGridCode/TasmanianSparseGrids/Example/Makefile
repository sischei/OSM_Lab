#
# Check the compile options here
# this is an example of how to link to
# the tasmanian library
#
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
COMPILE_OPTIONS = -fopenmp
else
COMPILE_OPTIONS = -Wno-unknown-pragmas
endif

CC = g++
OPTC = -O3 -ffast-math -mtune=native $(COMPILE_OPTIONS) -fPIC -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic
LIBTSG = ../libtasmaniansparsegrid.a
ITSG = -I../

all: example

example: example.cpp
	$(CC) $(OPTC) $(ITSG) example.cpp -o example $(LIBTSG)
	
clean:
	rm -fr example
