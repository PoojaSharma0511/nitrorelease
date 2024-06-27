MKLROOT=/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl

all: aout

aout: mod_md.o md.o
	ifort -o aout mod_md.o md.o -O3

%.o: %.f90
	ifort -c $< -O3

quick:
	gfortran -o aout mod_md.f90 md.f90

clean:
	rm *.o aout

