# Default target (this line must go first)
all: dg

# Compiler definitions
# the compiler: gcc for C program
CC = gcc
# copiler flags:
#  -02 -g for optimization/debugging
CFLAGS = -O2 -g
# library flags:
#  -lm: link <math.h> library 
LFLAGS = -lm

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
dg: $(objs) $(headers)
	$(CC) $(objs) $(LFLAGS) -o $@

# Clean executables
#  RM = "rm -f"
clean:
	$(RM) dg *.o
