TARGET  =pdbmodule.so
INCLUDES=
CFLAGS 	=-g -O6 -march=pentium4 -I/opt/include/python2.3 -D_UNIX
LFLAGS 	=
CC 	=gcc

SOURCE 	= pdbmodule.c
OBJ	= pdbmodule.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS)

clean:
	rm -f *.o $(TARGET)
