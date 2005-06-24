TARGET    =tlsmdmodule.so
INCLUDES  =-I/home/jpaint/local/include/python2.4
DEFINES   =-D_UNIX
CFLAGS 	  =-g -march=pentium4 -O2
CC 	  =gcc

SOURCE 	= tlsmdmodule.c
OBJ	= tlsmdmodule.o

all: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c -o $@ $<

$(TARGET): $(OBJ)
	$(CC) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS) minpack/minpack.a

clean:
	rm -f *.o $(TARGET)
