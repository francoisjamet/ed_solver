include ../Makefile.in
ED_PATH  = ../ed_solver/
MODFILES    = -I${ED_PATH}/MODFILES
LIB         = -L${ED_PATH} -ldmftsolverED -lsplines



prog: ${FILES}
	@echo " ........... compile : static / dynamic linking ........... "
	${FCOMPL} ${FFLAG} -o dmft_chiloc_1 Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK} -DM_NORB=1 -DM_NORB_C=1 -cpp
	${FCOMPL} ${FFLAG} -o dmft_chiloc_2 Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK} -DM_NORB=2 -DM_NORB_C=2 -cpp
	${FCOMPL} ${FFLAG} -o dmft_chiloc_3 Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK} -DM_NORB=3 -DM_NORB_C=3 -cpp
	${FCOMPL} ${FFLAG} -o dmft_chiloc_4 Chiloc_four_leg_rout.f90  Chiloc_four_leg_main.f90  ${MODFILES} ${LIB} ${LAPACK} -DM_NORB=4 -DM_NORB_C=4 -cpp


#############################
#############################
#############################

	@echo " ...................... done .............................. "
install :
	mkdir -p ../bin
	cp dmft_chiloc_* ../bin/.
	cp ../utils/vertex.py ../bin/readvertex
	cp ../utils/path.py ../bin/omega_path
	chmod +x ../bin/readvertex ../bin/omega_path

#########################################################################
