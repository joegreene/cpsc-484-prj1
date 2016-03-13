# Makefile

objects := gmath_test
CFLAGS := -Wall -std=c++11
LDFLAGS := -Wall
CC := clang++

all: $(objects)

clean:
	-rm -f $(objects) *.o

.PHONY: all clean

gmath_test: gmath_test.cc gmath.hh
	$(CC) $(CFLAGS) $< -o $@
