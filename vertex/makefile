include ../Makefile.in
ED_PATH  = ../ed_solver/
MODFILES    = -I${ED_PATH}/MODFILES
LIB         = -L${ED_PATH} -ldmftsolverED -lsplines



prog: ${FILES}
	@echo " ........... compile : static / dynamic linking ........... "
	${FCOMPL} ${FFLAG} -o dmft_chiloc Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK}


#############################
#############################
#############################

	@echo " ...................... done .............................. "
install :
	mkdir -p ../bin
	cp dmft_chiloc ../bin/.
	cp ../utils/vertex.py ../bin/readvertex
	cp ../utils/path.py ../bin/omega_path
	chmod +x ../bin/readvertex ../bin/omega_path

#########################################################################
