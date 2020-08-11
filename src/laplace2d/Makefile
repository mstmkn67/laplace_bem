CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o laplace2d  -L $(PF_FILES)/lib/linux_64 -lplatform -llapack -lblas

OBJECT= Timer.o gourmain.o cmdlineopts.o \
	bicgstab.o  jacobi_bicgstab.o non_bicgstab.o ilu_bicgstab.o \
	LapackFunctions.o \
	Mesh_c.o PotentialGreenFunc_c.o \
	LaplaceBEM2d.o Laplace2d.o Laplace2d_main.o 

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

Mesh_c.o : ../common/mesh/Mesh_c.cpp ../common/mesh/Mesh_c.h
	$(CC) ../common/mesh/Mesh_c.cpp $(CFLAG)

PotentialGreenFunc_c.o : ../common/GreenFunc/PotentialGreenFunc_c.cpp ../common/GreenFunc/PotentialGreenFunc_c.h
	$(CC) ../common/GreenFunc/PotentialGreenFunc_c.cpp $(CFLAG)

LaplaceBEM2d.o : LaplaceBEM2d.h LaplaceBEM2d.cpp
	$(CC) LaplaceBEM2d.cpp $(CFLAG)

Laplace2d.o : Laplace2d.h Laplace2d.cpp
	$(CC) Laplace2d.cpp $(CFLAG)

Laplace2d_main.o : Laplace2d_main.cpp
	$(CC) Laplace2d_main.cpp $(CFLAG)

clean : 
	rm -f *.o laplace2d

install :
	strip laplace2d
	cp laplace2d ../../bin
