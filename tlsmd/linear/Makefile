TARGET    =lineartls.so
INCLUDES  =-I/home/jpaint/local/include/python2.4
DEFINES   =-D_UNIX
CFLAGS 	  =-fPIC -g -O3 -Wall
CC 	  =gcc

LAPACK    = /home/jpaint/src/LAPACK/lapack_LINUX.a /home/jpaint/src/LAPACK/blas_LINUX.a


SOURCE 	= lineartls.c
OBJ	= lineartls.o

all: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CC) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) $(LAPACK) -lg2c -lm

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) ../bin/tlsmdlib/.
