# fast gcmc make file


# Any optimisation flags?
OPTIM = -O3

# warning flags or any other compiler options?
#FLAGS = -Wall -Wno-conversion
FLAGS = 
# debugging
DEBUG =
# gfortran debugging
#DEBUG = -g -fbacktrace #-ffpe-trap=zero,overflow,underflow
# intel fortran debugging
#DEBUG = -debug extended -g -check all -traceback

#FC=gfortran
#FC=mpifort
FC=ifort

OBJNUM = mpmodule.f90 mpfuna.f90 mpfunbq.f90 mpfunc.f90 mpfund.f90 mpfune.f90 mpfunf.f90 mpfungq1.f90 second.f90

test: $(OBJNUM)
	$(FC) $(OPTIM) -c $^ $(FLAGS) $(DEBUG) 

clean:
	rm -f *.o *.mod 

#%.o: %.f90	${FC} ${FFLAGS} ${DEBUG} -c $< -o $@

