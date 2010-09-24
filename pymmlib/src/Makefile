TARGET	= pdbmodule.so
INCLUDES= 
CFLAGS	= -g -I/usr/include/python2.5 -D_UNIX -fPIC
LFLAGS	= 
CC 	= gcc

SOURCE	= pdbmodule.c
OBJ	= pdbmodule.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -shared $(CLFAGS) -o $(TARGET) $(OBJ) $(LFLAGS)

clean:
	rm -f *.o $(TARGET)
