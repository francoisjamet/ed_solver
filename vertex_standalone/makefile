
#########################################################################
  SHELL       =  /bin/bash
  FCOMPL      =  mpif90
  FFLAG       =  -O3  -fopenmp -ffree-line-length-2000 -ffree-line-length-4000 -lstdc++

  MODFILES    = -I/home/jamet/src/ed/ed_cweber/ed_solver/MODFILES
  LIB         = -L/home/jamet/src/ed/ed_cweber/ed_solver/ -ldmftsolverED -lsplines
  LAPACK      = -L${HOME}/MYLIBRARY/LIBBLAS -L${HOME}/MYLIBRARY/LIBLAPACK -llapack -lblas
#########################################################################

prog: ${FILES}
	@echo " ........... compile : static / dynamic linking ........... "
	${FCOMPL} ${FFLAG} -o chiloc.out Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK}


#############################
#############################
#############################

	@echo " ...................... done .............................. "

#########################################################################
