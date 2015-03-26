# Compiler definitions
# the compiler: gcc for C program
CC = gcc
# copiler flags:
#  -02 -g for optimization/debugging
CFLAGS = -Wall -Wextra -O2 -g
# library flags:
#  -lm: link <math.h> library 
LFLAGS = -lm

# Get config.h vars
N1D=$(shell awk <config.h '$$2 == "N1D" { print $$3; }')
InterpolationOrder=$(shell awk <config.h '$$2 == "InterpolationOrder" { print $$3; }')
Nt=$(shell awk <config.h '$$2 == "Nt" { print $$3; }')

# Make this executable
executable_name=dg_n$(N1D)_p$(InterpolationOrder)_t$(Nt)

# Default target (convention)
.PHONY: all
all: $(executable_name)

# Sources and headers
srcs = dg.c
headers = dg.h

# Objects from sources
#  patsubst: function that replace pattern in text
objs = $(patsubst %.c,%.o,$(srcs))

# Rule to compile sources
#  for rule manual, refer to gnu.org about GNU make
#  -c: compile source files without linking
#  automatic variables: $^ for names of all prerequisites, $@ for name of target
%.o: %.c
	$(CC) $(CFLAGS) -c $^

# Link executable
$(executable_name): $(objs) $(headers)
	$(CC) $(LFLAGS) $(objs) -o $@
	rm -f dg
	ln -s $@ dg

# Clean executables
#  RM = "rm -f"
clean:
	rm -f *.o
	rm -f dg_n*_p*_t*
	rm -f dg
	rm all.pbs
