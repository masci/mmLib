TARGET    =tlsmdmodule.so
INCLUDES  =-I/home/jpaint/local/include/python2.4
DEFINES   =-D_UNIX
CFLAGS 	  =-g -march=pentium4 -O3 -Wall
CC 	  =gcc

SOURCE 	= tlsmdmodule.c
OBJ	= tlsmdmodule.o

all: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CC) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) /home/jpaint/build/LAPACK/lapack_LINUX.a /home/jpaint/build/LAPACK/blas_LINUX.a -lg2c -lm

clean:
	rm -f *.o $(TARGET)
