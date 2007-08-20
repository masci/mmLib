TARGET    =tlsmdmodule.so
INCLUDES  =-I$(HOME)/local/include/python2.4
DEFINES   =-D_UNIX
CFLAGS 	  =-fPIC -g -O3 -Wall
CPP       =g++
CC 	  =gcc

LAPACK_ROOT   =$(HOME)/build/LAPACK
LAPACK        =$(LAPACK_ROOT)/lapack_LINUX.a $(LAPACK_ROOT)/blas_LINUX.a
ATLAS         =-L$(HOME)/ATLAS/lib/Linux_P4SSE2/ -llapack -lf77blas -lcblas -latlas

MINPACK_ROOT  =$(HOME)/build/minpack
MINPACK       =$(MINPACK_ROOT)/libminpack.a

SOURCE 	= structure.cpp dgesdd.cpp tls_model.cpp tls_model_nl.cpp tls_model_engine.cpp tlsmdmodule.cpp
OBJ	= structure.o   dgesdd.o   tls_model.o   tls_model_nl.o   tls_model_engine.o   tlsmdmodule.o

all: $(TARGET) 

%.o: %.cpp
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CPP) -shared -fPIC $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) $(MINPACK) $(LAPACK) libg2c.a -lm

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) ../bin/tlsmdlib/.
