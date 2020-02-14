#!/bin/sh 
#M.Huang@PNNL, 03/28/2017
#Last updated on 09/07/2018 

#Note: need to load the following modules to compile on constance
# module load intel/15.0.1  
# module load mkl/15.0.1
# module load netcdf/4.3.2
# module load mvapich2/2.1

#define the directories, files, cases
export BASE_DIR=/pic/projects/landuq/huan565/ngeet
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export USER_INPUTDATA_DIR=${BASE_DIR}/inputdata/user_inputdata
export GRIDNAME=1x1pt_km83
export OCNDOM=${GRIDNAME}
export ATMDOM=${GRIDNAME}
echo $OCNDOM
echo $ATMDOM
export CDATE=190627
export MAPFILE=${USER_INPUTDATA_DIR}/${GRIDNAME}/maps/map_${GRIDNAME}_noocean_to_${GRIDNAME}_nomask_aave_da_${CDATE}.nc

#go to the directory
cd ${CESM_SRC_DIR}/cime/tools/mapping/gen_domain_files

# prepare the environ for compilation
export INC_NETCDF=/share/apps/netcdf/4.3.2/intel/15.0.1/include
export LIB_NETCDF=/share/apps/netcdf/4.3.2/intel/15.0.1/lib
export NETCDF_PATH=/share/apps/netcdf/4.3.2/intel/15.0.1
export COMPILER=intel
export NETCDF_HOME=${NETCDF_PATH}
export MKLROOT=${MKLROOT}

#generate the domain file
./gen_domain -m ${MAPFILE} -o ${OCNDOM} -l ${ATMDOM}
mv ${CESM_SRC_DIR}/cime/tools/mapping/gen_domain_files/domain*.nc ${USER_INPUTDATA_DIR}/${GRIDNAME}/
