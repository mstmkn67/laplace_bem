CC = g++
CFLAG = -c -O3 
CCFLAG =  -lm -llapack
OBJECTS=checkLapack.o LapackFunctions.o

all : $(OBJECTS)
	$(CC) -o checkLapack  $(OBJECTS)  $(CCFLAG)

checkLapack.o : checkLapack.cpp 
	$(CC) checkLapack.cpp $(CFLAG)

LapackFunctions.o : LapackFunctions.cpp LapackFunctions.h
	$(CC) LapackFunctions.cpp $(CFLAG)

clean : 
	rm *.o checkLapack
