# DGadvection

Discontinuous Galerkin finite element solver for 2D uniform advection.

# Building

Macros needed to compile are written by configure to config.h. Once
these are created, the default 'make' target will build the executable
with a descriptive name and will add a link to it with the name 'dg'.

```
$ ./configure 8 1 100 # Configure for 8 elements, p1, using 100 timesteps
$ make                # Build executable (named dg_n8_p1_t100)
```
