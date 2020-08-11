CC = g++
CFLAG = -c
CCFLAG =-lm -o 
OBJECT1=check1d.o 
OBJECT2=Vector2d.o checkTriangle.o

all : $(OBJECT1) $(OBJECT2)
	$(CC) $(CCFLAG) check1d.out $(OBJECT1)
	$(CC) $(CCFLAG) checkTriangle.out $(OBJECT2)
	strip check1d.out
	strip checkTriangle.out

check1d.o : check1d.cpp
	$(CC) check1d.cpp $(CFLAG)

checkTriangle.o : checkTriangle.cpp
	$(CC) checkTriangle.cpp $(CFLAG)

Vector2d.o : ../Vector2d.cpp ../Vector2d.h
	$(CC) ../Vector2d.cpp $(CFLAG)

clean : 
	rm *.o *.out
