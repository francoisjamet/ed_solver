  SHELL  = /bin/bash
  FCOMPL = mpif90 -ffree-line-length-4000 -lstdc++    #gfortran-mp-6
  C++    = g++
  prep   = -cpp -D__complex -fopenmp -DOPENMP_MPI -D_arpack
  FFLAG  = -O3 -fopenmp
  FFLAG2 = -O3
  RANLIB = ranlib
  ARPACK  = -L./LIB_ARPACK -L./LIBLAPACK -L./LIBBLAS -lparpack -larpack


ROUTINES =  core_functions_fortran.f90 openmpmod.f90 GeneralVariables.f90 specFunc.f90 core_functions_def.f90 MPI_routines.f90 random.f90 Strings5.f90 Strings.f90 namelist_module.f90 core_functions_initsim_closesim.f90 splines.f90 mask_class.f90 masked_mat_class.f90 quantum_algebra.f90 quantum_fermion_ket_class.f90 quantum_readable_vec.f90 quantum_fermion_hilbert_class.f90 quantum_fermion_sector_class.f90 quantum_sector_class.f90 frequency_axis.f90 Global_Variables.f90 masked_matrix_real_comp.f90 rcvector_class.f90 rcmatrix_class.f90 conj_grad_civelli.f90  tools_minimization_wrapping.f90 liste3.f90 eigen_class.f90 eigen_sector_class.f90 impurity_class.f90 correl_class.f90 bath_class.f90 bath_class_vec.f90 bath_class_hybrid.f90 quantum_algebra_Anderson_Impurity_class.f90 quantum_algebra_Anderson_Impurity_class2.f90 HAIM_offdiag_arrays.f90 HAIMsz.f90 HAIMupdo.f90 density_matrix.f90 overlap.f90 apply_C.f90 apply_NS.f90 apply_P.f90 H_class.f90 matrices_tridiag.f90 Lanczos_gpu.f90 Lanczos_plain.f90 ED_ARPACK.f90 Block_Lanczos.f90 green_class.f90 green_class_compute_dynamic.f90 green_class_compute_symmetric_combination.f90 green_class_computeAA.f90 green_class_computeAB.f90 correlations.f90 Lanczos_Cullum_mult.f90 Lanczos_Cullum_lesub.f90 Lanczos_Cullum_hleval.f90 Lanczos_Cullum_leval.f90 Lanczos_Cullum_hlevec.f90 Lanczos_Cullum_levec.f90 Lanczos_Cullum_wrapper.f90  Chiloc_four_leg.f90 solver.f90 dmft_ed_solver.f90

ROUTINES2 = eigen_sector_class.f90

MODFILES= -I./MODFILES

comp:  ${ROUTINES}
	@echo "compile"
	$(FCOMPL) -c ${prof} $(FFLAG) ${prep} ${ROUTINES} ${MODFILES}
	@echo "done...."
	@echo "now compile C++ files..."
	$(C++) -c memory_address_functions.cc
	@mv *.mod MODFILES/

comp2:  ${ROUTINES2}
	@echo "compile2"
	$(FCOMPL) -c $(FFLAG2) ${prof} ${prep} ${ROUTINES2} ${MODFILES}
	@echo "set of files..."
	@mv *.mod MODFILES/

build: comp comp2 finish

prog:
	cp splines/libsplines.a .
	${FCOMPL} ${prep} -o dmft_solver.out  ${FFLAG} memory_address_functions.o ed_prog.f90 -L. -L${HOME}//LIBBLAS -L${HOME}//LIBLAPACK -ldmftsolverED -lsplines -I./MODFILES $(ARPACK)  -llapack -lblas

finish:
	cp splines/libsplines.a .
	ar ru libdmftsolverED.a `ls *.o | sort | uniq`
	${RANLIB} libdmftsolverED.a
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

clean:
	rm -f *.o
	rm -f *.a
	rm -f MODFILES/*.mod
cleanall : clean
	pushd LIBBLAS && make clean && 	popd
	pushd LIBLAPACK && 	make clean &&	popd
	pushd LIB_ARPACK && 	make clean &&	popd
