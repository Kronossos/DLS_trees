
CPPFLAGS=-g 
CC=g++
#LDFLAGS=-O2
#CPPFLAGS=-g


OBJ=dlsgen2.o dlslib.o 

all:   dlsgen2
dlslib.o: dlslib.cpp dlslib.h 
dlsgen2.o: dlsgen2.cpp dlslib.o dlslib.h 

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

dlsgen2: $(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ) 

dlsgen: $(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ) 


clean:
	rm -f dlsgen *.o
