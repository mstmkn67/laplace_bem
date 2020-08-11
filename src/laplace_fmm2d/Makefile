CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o laplace_fmm2d  -L $(PF_FILES)/lib/linux_64 -lplatform -llapack -lblas

OBJECT= Timer.o gourmain.o cmdlineopts.o \
	bicgstab.o fmm_non_bicgstab.o fmm_jacobi_bicgstab.o \
	fmm_block_bicgstab.o LapackFunctions.o \
	FMMTree_c.o LaplaceFMMTree_c.o Cell_c.o Mesh_c.o Moment_c.o \
	PotentialGreenFunc_c.o \
	LaplaceBEM_FMM2d.o LaplaceFMM2d.o LaplaceFMM2d_main.o

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

FMMTree_c.o : ../common/mesh/FMMTree_c.cpp ../common/mesh/FMMTree_c.h
	$(CC) ../common/mesh/FMMTree_c.cpp $(CFLAG)

LaplaceFMMTree_c.o : ../common/mesh/LaplaceFMMTree_c.cpp ../common/mesh/LaplaceFMMTree_c.h
	$(CC) ../common/mesh/LaplaceFMMTree_c.cpp $(CFLAG)

Cell_c.o : ../common/mesh/Cell_c.cpp ../common/mesh/Cell_c.h
	$(CC) ../common/mesh/Cell_c.cpp $(CFLAG)

Mesh_c.o : ../common/mesh/Mesh_c.cpp ../common/mesh/Mesh_c.h
	$(CC) ../common/mesh/Mesh_c.cpp $(CFLAG)

Moment_c.o : ../common/mesh/Moment_c.cpp ../common/mesh/Moment_c.h
	$(CC) ../common/mesh/Moment_c.cpp $(CFLAG)

PotentialGreenFunc_c.o : ../common/GreenFunc/PotentialGreenFunc_c.cpp ../common/GreenFunc/PotentialGreenFunc_c.h
	$(CC) ../common/GreenFunc/PotentialGreenFunc_c.cpp $(CFLAG)

LaplaceBEM_FMM2d.o : LaplaceBEM_FMM2d.h LaplaceBEM_FMM2d.cpp
	$(CC) LaplaceBEM_FMM2d.cpp $(CFLAG)

LaplaceFMM2d.o : LaplaceFMM2d.h LaplaceFMM2d.cpp
	$(CC) LaplaceFMM2d.cpp $(CFLAG)

LaplaceFMM2d_main.o : LaplaceFMM2d_main.cpp
	$(CC) LaplaceFMM2d_main.cpp $(CFLAG)

clean : 
	rm -f *.o laplace_fmm2d

install :
	strip laplace_fmm2d
	cp laplace_fmm2d ../../bin
