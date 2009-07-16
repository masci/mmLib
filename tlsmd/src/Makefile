# Makefile for the tlsmdmodule
## Note: This is for the 64-bit compile. Christoph Champ, 2008-02-15
## Updated for 2009.1, 2009-07-10
TARGET    =tlsmdmodule.so
INCLUDES  =-I/usr/include/python2.6/
DEFINES   =-D_UNIX
DEFAULT_LIB_INSTALL_PATH = /home/tlsmd/tlsmd/src
CFLAGS 	  =-fPIC -g -O3 -Wall
LDFLAGS   = -lgfortran
CPP       =g++
CC 	  =gcc

LAPACK	= -llapack -lblas

MINPACK_ROOT  =/usr/local/lib
MINPACK       =$(MINPACK_ROOT)/libminpack.a

SOURCE 	= structure.cpp dgesdd.cpp tls_model.cpp tls_model_nl.cpp tls_model_engine.cpp tlsmdmodule.cpp
OBJ	= structure.o   dgesdd.o   tls_model.o   tls_model_nl.o   tls_model_engine.o   tlsmdmodule.o

all: $(TARGET) 

%.o: %.cpp
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CPP) -shared -fPIC $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) $(MINPACK) $(LAPACK) -lm

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) ../bin/tlsmdlib/.
