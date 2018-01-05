# fast gcmc make file


# Any optimisation flags?
# cray agressive optimization
#OPTIM = -O aggress
OPTIM =

# warning flags or any other compiler options?
#FLAGS = -Wall -Wno-conversion
FLAGS = 
# debugging
DEBUG =
# gfortran debugging
#DEBUG = -g -fbacktrace #-ffpe-trap=zero,overflow,underflow
# intel fortran debugging
#DEBUG = -debug extended -g -check all -traceback

SERIAL=ser
PARALLEL=par
BNLIBDIR=mpfun2015
.SUFFIXES: .f .f90 .o


OBJPAR = parse_module.o basic_comms.o utility_pack.o error.o \
		 vdw_module.o vdw_terms.o nlist_builders.o \
		 ewald_module.o ewald.o flex_module.o readinputs.o \
		 mc_moves.o wang_landau_module.o wang_landau.o gcmc.o

OBJSER = parse_module.o setup_module.o serial.o utility_pack.o \
		 error.o vdw_module.o vdw_terms.o nlist_builders.o \
		 ewald_module.o ewald.o flex_module.o readinputs.o \
		 mc_moves.o wang_landau_module.o wang_landau.o gcmc.o

# must choose a version to compile
all:
	@echo "[1m>> Please specify the target:[0m"
	@echo
	@echo "[1;32m * [0mintel"
	@echo "[1;32m * [0mcray"
	@echo "[1;32m * [0mgfortran"
	@echo "[1;32m * [0mparallel"
	@echo 
	@echo "NUMERICAL LIBRARY:"
	@echo 
	@echo "[1;32m * [0mmpfun"
	@echo

# parameters for different versions
gfortran:
	${MAKE} FC="gfortran" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc_gfort.x" ${SERIAL}
cray:
	${MAKE} FC="ftn -f fixed" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc-cray.x" ${SERIAL}
intel:
	${MAKE} FC="ifort" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc.x" ${SERIAL}
parallel:
	${MAKE} FC="mpifort" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc-par.x" ${PARALLEL}
mpfun: 
	for i in {1..2}; do\
	  $(MAKE) -C $(BNLIBDIR);\
	done 
	
# stuff that does the compilations
ser: ${OBJSER}
	${FC} ${FFLAGS} ${DEBUG} -o ${EXE} ${OBJSER}

par: ${OBJPAR}
	${FC} ${FFLAGS} ${DEBUG} -o ${EXE} ${OBJPAR}

clean:
	rm -f *.o *.mod

clean-all:
	$(MAKE) -C $(BNLIBDIR) clean 
	rm -f *.o *.mod 

%.o: ${BNLIBDIR}/%.f90
	${FC} ${FFLAGS} ${DEBUG} -c $< -o $@

%.o: %.f 
	${FC} ${FFLAGS} ${DEBUG} -c $< -o $@

