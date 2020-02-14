#!/bin/sh 
#Sample shell script to configure, compile, and submit a supported simulation
#M.Huang@PNNL, 03/28/2017
#Last updated on 09/07/2018

#Note: need to load the following modules to compile on constance
# module load intel/15.0.1  
# module load mkl/15.0.1
# module load netcdf/4.3.2
# module load mvapich2/2.1

export BASE_DIR=/pic/projects/landuq/huan565/ngeet

#define the directories, files, cases
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export INPUTDATA_DIR=${BASE_DIR}/inputdata
#export SFC=ifort

#go to the directory
cd ${CESM_SRC_DIR}/cime/tools/mapping/gen_domain_files/src

# prepare the environ for compilation
export INC_NETCDF=/share/apps/netcdf/4.3.2/intel/15.0.1/include
export LIB_NETCDF=/share/apps/netcdf/4.3.2/intel/15.0.1/lib
export NETCDF_PATH=/share/apps/netcdf/4.3.2/intel/15.0.1
export COMPILER=intel
export NETCDF_HOME=${NETCDF_PATH}
export MKLROOT=${MKLROOT}

#compile gen_domain
../../../configure --machine constance --machines-dir ${CESM_SRC_DIR}/cime/config/cesm/machines --compiler intel --macros-format Makefile --mpilib intelmpi
gmake 
