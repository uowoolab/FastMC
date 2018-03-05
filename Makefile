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
# The arbitrary precision number library courtesy of David Bailey.
# www.davidhbailey.com/dhbpapers/mpfun2015.pdf
# the directory is a local copy provided for the user's convenience.
# If one wishes to compile their own version, then change BNLIBDIR
# to the full path of their preferred version.
BNLIBDIR=mpfun2015
DOS=ddos
.SUFFIXES: .f .o

OBJPAR = parse_module.o basic_comms.o utility_pack.o error.o \
		 vdw_module.o vdw_terms.o nlist_builders.o \
		 ewald_module.o ewald.o flex_module.o readinputs.o \
		 mc_moves.o wang_landau_module.o wang_landau.o gcmc.o

OBJSER = parse_module.o setup_module.o serial.o utility_pack.o \
		 error.o vdw_module.o vdw_terms.o nlist_builders.o \
		 ewald_module.o ewald.o flex_module.o readinputs.o \
		 mc_moves.o wang_landau_module.o wang_landau.o gcmc.o

LIBOBJ = $(BNLIBDIR)/mpmodule.o \
	 $(BNLIBDIR)/mpfuna.o \
	 $(BNLIBDIR)/mpfunbq.o \
	 $(BNLIBDIR)/mpfunc.o \
	 $(BNLIBDIR)/mpfund.o \
	 $(BNLIBDIR)/mpfune.o \
	 $(BNLIBDIR)/mpfunf.o \
	 $(BNLIBDIR)/mpfungq1.o \
	 $(BNLIBDIR)/second.o

OBJDOS = parse_module.o setup_module.o serial.o utility_pack.o \
	 error.o vdw_module.o vdw_terms.o nlist_builders.o \
	 ewald_module.o ewald.o readinputs.o \
	 mc_moves.o wang_landau_module.o dos_to_isotherm.o\

# must choose a version to compile
all:
	@echo "[1m>> Please specify the target:[0m"
	@echo
	@echo "[1;32m * [0mintel"
	@echo "[1;32m * [0mcray"
	@echo "[1;32m * [0mgfortran"
	@echo "[1;32m * [0mparallel-gfort"
	@echo "[1;32m * [0mparallel-cray"
	@echo "[1;32m * [0mparallel-intel"
	@echo 
	@echo "NUMERICAL LIBRARY:"
	@echo 
	@echo "[1;32m * [0mmpfun"
	@echo
	@echo "TOOLS:"
	@echo 
	@echo "[1;32m * [0mdos_tool"
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
parallel-gfort:
	${MAKE} FC="mpifort" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc-par.x" ${PARALLEL}
parallel-intel:
	${MAKE} FC="mpiifort" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc-par.x" ${PARALLEL}
parallel-cray:
	${MAKE} FC="ftn -f fixed" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="gcmc-cray-par.x" ${PARALLEL}
mpfun: 
	for i in {1..2}; do\
	  $(MAKE) -C $(BNLIBDIR);\
	done 
dos_tool:
	${MAKE} FC="gfortran" \
	FFLAGS="${OPTIM} ${FLAGS}" EXE="dos_to_isotherm" $(DOS)

# stuff that does the compilations
ser: ${OBJSER}
	${FC} ${FFLAGS} ${DEBUG} -o ${EXE} $(LIBOBJ) $(OBJSER)

par: ${OBJPAR}
	${FC} ${FFLAGS} ${DEBUG} -o ${EXE} $(LIBOBJ) $(OBJPAR)

clean:
	rm -f *.o *.mod

clean-all:
	$(MAKE) -C $(BNLIBDIR) clean 
	rm -f *.o *.mod 

ddos: $(OBJDOS)
	$(FC) $(FFLAGS) $(DEBUG) -o $(EXE) $(LIBOBJ) $(OBJDOS) 

%.o: %.f 
	${FC} -I$(BNLIBDIR) ${FFLAGS} ${DEBUG} -c $< -o $@

