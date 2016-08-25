all: aout

aout: mod_afssh.o AFSSH.o
	gfortran -o aout mod_afssh.o AFSSH.o -O2 ~/lib/libcommon.a ~/lapack-3.6.0/liblapack.a ~/BLAS/libblas_LINUX.a

%.o: %.f90
	gfortran -c $< -I/data/home/ambjain/lib

clean:
	rm *.o aout

