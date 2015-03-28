# Compiler definitions
# the compiler: gcc for C program
CC = gcc
# copiler flags:
#  -02 -g: optimization/debugging
#  -Wall -Wextra suppresed for now
CFLAGS = -O2 -g
# library flags:
#  -lm: link <math.h> library 
LFLAGS = -lm

# Get config.h vars
#  use awk to print string
#  $$2 means "the value (leading $)" of "the second field ($2)"
N1D=$(shell awk <config.h '$$2 == "N1D" { print $$3; }')
InterpolationOrder=$(shell awk <config.h '$$2 == "InterpolationOrder" { print $$3; }')
Nt=$(shell awk <config.h '$$2 == "Nt" { print $$3; }')

# Make this executable
executable_name=dg_n$(N1D)_p$(InterpolationOrder)_t$(Nt)

# Default target (convention)
#  use of phony target; 
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
#  compiler executable
#  remove old dg link
#  link dg symbollically to executable
#  LFLAGS should be after objs
$(executable_name): $(objs) $(headers)
	$(CC) $(objs) $(LFLAGS) -o $@
	rm -f dg
	ln -s $@ dg

# Clean executables
#  RM = "rm -f"
clean:
	rm -f *.o
	rm -f dg_n*_p*_t*
	rm -f dg
	rm all.pbs
