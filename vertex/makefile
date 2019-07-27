
#########################################################################
  SHELL       =  /bin/bash
  FCOMPL      =  mpif90
  FFLAG       =  -O3  -fopenmp -ffree-line-length-2000 -ffree-line-length-4000 -lstdc++

  MODFILES    = -I../ed_solver/MODFILES
  LIB         = -L../ed_solver/ -ldmftsolverED -lsplines
  LAPACK      = -L../ed_solver/LIBBLAS -L../ed_solver/LIBLAPACK -llapack -lblas
#########################################################################

prog: ${FILES}
	@echo " ........... compile : static / dynamic linking ........... "
	${FCOMPL} ${FFLAG} -o chiloc.out Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK}


#############################
#############################
#############################

	@echo " ...................... done .............................. "

#########################################################################
