#!/bin/bash
pushd ed_solver
make cleanall
pushd splines
make build
popd
pushd LIBBLAS
make all -j $1
popd
pushd LIBLAPACK
make lib -j $1
popd
pushd LIB_ARPACK
make all -j $1
popd
make build
make prog
popd
pushd vertex
make prog
popd