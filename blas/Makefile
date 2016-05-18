CC=g++
CFLAGS=-O2 -ansi -pedantic -DFORTRAN_TRAILING_UNDERSCORE -I/System/Library/Frameworks/Kernel.framework/Versions/A/Headers/

## for Mac OS X
LIBS=-framework Accelerate 

## For Linux
#LIBS=-llapack -lblas -lg2c -lm

example: example.cpp
	$(CC) -o $@ $@.cpp $(CFLAGS) $(LIBS)

clean:
	rm example 2>&-
