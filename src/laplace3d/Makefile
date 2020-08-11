CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o laplace3d  -L $(PF_FILES)/lib/linux_64 -lplatform -llapack -lblas

OBJECT= Timer.o gourmain.o cmdlineopts.o \
	bicgstab.o  jacobi_bicgstab.o non_bicgstab.o ilu_bicgstab.o \
	LapackFunctions.o \
	Mesh3d.o PotentialGreenFunc.o \
	LaplaceBEM3d.o Laplace3d.o Laplace3d_main.o 

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : ../common/Timer.cpp ../common/Timer.h
	$(CC) ../common/Timer.cpp $(CFLAG)

gourmain.o : ../common/udf/gourmain.cpp ../common/udf/gourmain.h
	$(CC) ../common/udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : ../common/udf/cmdlineopts.cpp ../common/udf/cmdlineopts.h
	$(CC) ../common/udf/cmdlineopts.cpp $(CFLAG)

bicgstab.o : ../common/bicgstab/bicgstab.h ../common/bicgstab/bicgstab.cpp
	$(CC) ../common/bicgstab/bicgstab.cpp $(CFLAG)

jacobi_bicgstab.o : ../common/bicgstab/jacobi_bicgstab.h ../common/bicgstab/jacobi_bicgstab.cpp
	$(CC) ../common/bicgstab/jacobi_bicgstab.cpp $(CFLAG)

non_bicgstab.o : ../common/bicgstab/non_bicgstab.h ../common/bicgstab/non_bicgstab.cpp
	$(CC) ../common/bicgstab/non_bicgstab.cpp $(CFLAG)

ilu_bicgstab.o : ../common/bicgstab/ilu_bicgstab.h ../common/bicgstab/ilu_bicgstab.cpp
	$(CC) ../common/bicgstab/ilu_bicgstab.cpp $(CFLAG)

LapackFunctions.o : ../common/lapack/LapackFunctions.h ../common/lapack/LapackFunctions.cpp
	$(CC) ../common/lapack/LapackFunctions.cpp $(CFLAG)

Mesh3d.o : ../common/mesh/Mesh3d.cpp ../common/mesh/Mesh3d.h
	$(CC) ../common/mesh/Mesh3d.cpp $(CFLAG)

PotentialGreenFunc.o : ../common/GreenFunc/PotentialGreenFunc.cpp ../common/GreenFunc/PotentialGreenFunc.h
	$(CC) ../common/GreenFunc/PotentialGreenFunc.cpp $(CFLAG)

LaplaceBEM3d.o : LaplaceBEM3d.h LaplaceBEM3d.cpp
	$(CC) LaplaceBEM3d.cpp $(CFLAG)

Laplace3d.o : Laplace3d.h Laplace3d.cpp
	$(CC) Laplace3d.cpp $(CFLAG)

Laplace3d_main.o : Laplace3d_main.cpp
	$(CC) Laplace3d_main.cpp $(CFLAG)

clean : 
	rm -f *.o laplace3d

install :
	strip laplace3d
	cp laplace3d ../../bin
