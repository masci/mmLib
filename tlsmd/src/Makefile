TARGET    =lineartls.so
INCLUDES  =-I/home/jpaint/local/include/python2.4
DEFINES   =-D_UNIX
CFLAGS 	  =-fPIC -g -O3 -Wall
CPP       =g++
CC 	  =gcc

LAPACK    = /home/jpaint/src/LAPACK/lapack_LINUX.a /home/jpaint/src/LAPACK/blas_LINUX.a
ATLAS     = -L/home/jpaint/build/ATLAS/lib/Linux_P4SSE2/ -llapack -lf77blas -lcblas -latlas
MINPACK   = /home/jpaint/build/minpack/libminpack.a

SOURCE 	= structure.cpp dgesdd.cpp tls_model.cpp tls_model_nl.cpp tls_model_engine.cpp tlsmdmodule.cpp
OBJ	= structure.o   dgesdd.o   tls_model.o   tls_model_nl.o   tls_model_engine.o   tlsmdmodule.o

all: $(TARGET) 

%.o: %.cpp
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CPP) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) $(MINPACK) $(LAPACK) -lg2c -lm

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) ../bin/tlsmdlib/.
