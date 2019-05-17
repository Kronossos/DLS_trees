TARGET = tt 
OBJ =  tools.o rtree.o urtree.o
CFLAGS = -Wall -c 
CC = g++ 
LFLAGS =  -Wall 

all: tt urec treeman

rtree.o : rtree.h rtree.cpp
urtree.o : urtree.h urtree.cpp

%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ $<

urec : $(OBJ) urec.o
	$(CC) $(LFLAGS) -o $@ $(OBJ) $@.o

treeman : $(OBJ) treeman.o
	$(CC) $(LFLAGS) -o $@ $(OBJ) $@.o

tt : $(OBJ) tt.o 
	$(CC) $(LFLAGS)  -o  $@ $(OBJ) $@.o

clean :
	rm -f *.o $(TARGET) *.old *~ x *.log

tgz : 
	tar czvf urec.tgz *.cpp *.h Makefile README

bin:
	tar czvf ttbin.tgz urec tt README smpnniheur.sh


