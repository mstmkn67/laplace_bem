CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o laplace_fmm3d  -L $(PF_FILES)/lib/linux_64 -lplatform -llapack -lblas

OBJECT= Timer.o gourmain.o cmdlineopts.o \
	bicgstab.o fmm_non_bicgstab.o fmm_jacobi_bicgstab.o \
	fmm_block_bicgstab.o LapackFunctions.o \
	LaplaceFMMTree3d.o FMMTree3d.o Cell3d.o Mesh3d.o Moment3d.o \
	PotentialGreenFunc.o Vector3d.o \
	LaplaceBEM_FMM3d.o LaplaceFMM3d.o LaplaceFMM3d_main.o 

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : ../common/Timer.cpp ../common/Timer.h
	$(CC) ../common/Timer.cpp $(CFLAG)

gourmain.o : ../common/udf/gourmain.cpp ../common/udf/gourmain.h
	$(CC) ../common/udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : ../common/udf/cmdlineopts.cpp ../common/udf/cmdlineopts.h
	$(CC) ../common/udf/cmdlineopts.cpp $(CFLAG)

bicgstab.o : ../common/bicgstab/bicgstab.cpp ../common/bicgstab/bicgstab.h
	$(CC) ../common/bicgstab/bicgstab.cpp $(CFLAG)

fmm_non_bicgstab.o : fmm_non_bicgstab.cpp fmm_non_bicgstab.h
	$(CC) fmm_non_bicgstab.cpp $(CFLAG)

fmm_jacobi_bicgstab.o : fmm_jacobi_bicgstab.cpp fmm_jacobi_bicgstab.h
	$(CC) fmm_jacobi_bicgstab.cpp $(CFLAG)

fmm_block_bicgstab.o : fmm_block_bicgstab.cpp fmm_block_bicgstab.h
	$(CC) fmm_block_bicgstab.cpp $(CFLAG)

LapackFunctions.o : ../common/lapack/LapackFunctions.cpp ../common/lapack/LapackFunctions.h
	$(CC) ../common/lapack/LapackFunctions.cpp $(CFLAG)

FMMTree3d.o : ../common/mesh/FMMTree3d.cpp ../common/mesh/FMMTree3d.h
	$(CC) ../common/mesh/FMMTree3d.cpp $(CFLAG)

LaplaceFMMTree3d.o : ../common/mesh/LaplaceFMMTree3d.cpp ../common/mesh/LaplaceFMMTree3d.h
	$(CC) ../common/mesh/LaplaceFMMTree3d.cpp $(CFLAG)

Cell3d.o : ../common/mesh/Cell3d.cpp ../common/mesh/Cell3d.h
	$(CC) ../common/mesh/Cell3d.cpp $(CFLAG)

Mesh3d.o : ../common/mesh/Mesh3d.cpp ../common/mesh/Mesh3d.h
	$(CC) ../common/mesh/Mesh3d.cpp $(CFLAG)

Moment3d.o : ../common/mesh/Moment3d.cpp ../common/mesh/Moment3d.h
	$(CC) ../common/mesh/Moment3d.cpp $(CFLAG)

PotentialGreenFunc.o : ../common/GreenFunc/PotentialGreenFunc.cpp ../common/GreenFunc/PotentialGreenFunc.h
	$(CC) ../common/GreenFunc/PotentialGreenFunc.cpp $(CFLAG)

Vector3d.o : ../common/Vector3d.h ../common/Vector3d.cpp
	$(CC) ../common/Vector3d.cpp $(CFLAG)

LaplaceBEM_FMM3d.o : LaplaceBEM_FMM3d.h LaplaceBEM_FMM3d.cpp
	$(CC) LaplaceBEM_FMM3d.cpp $(CFLAG)

LaplaceFMM3d.o : LaplaceFMM3d.h LaplaceFMM3d.cpp
	$(CC) LaplaceFMM3d.cpp $(CFLAG)

LaplaceFMM3d_main.o : LaplaceFMM3d_main.cpp
	$(CC) LaplaceFMM3d_main.cpp $(CFLAG)

clean : 
	rm -f *.o laplace_fmm3d

install :
	strip laplace_fmm3d
	cp laplace_fmm3d ../../bin
