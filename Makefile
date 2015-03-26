# Default target (this line must go first)
all: dg

# Compiler definitions
CC=gcc
CFLAGS=-O2 -g
LFLAGS=-lm

# Sources and headers
srcs=dg.c
headers=dg.h

# Objects from sources
objs=$(patsubst %.c,%.o,$(srcs))

# Rule to compile sources
%.o: %.c
	$(CC) $(CFLAGS) -c $^

# Link executable
dg: $(objs) $(headers)
	$(CC) $(LFLAGS) $(objs) -o $@

clean:
	rm -f $@
