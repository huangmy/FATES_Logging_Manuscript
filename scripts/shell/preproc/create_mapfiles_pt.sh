#!/bin/sh -l
#Sample shell script to create mapping files for a user-defined grid
#M.Huang@PNNL, 08/11/2017
#Modified to use on constance on 09/09/2018

#Note: need to load the following modules to compile on constance
# module load intel/15.0.1  
# module load mkl/15.0.1
# module load netcdf/4.3.2
# module load mvapich2/2.1

#define the directories, files, cases
export BASE_DIR=/pic/projects/landuq/huan565/ngeet
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export INPUTDATA_DIR=${BASE_DIR}/inputdata
export GRIDNAME=1x1pt_km83
export OCNDOM=${GRIDNAME}
export ATMDOM=${GRIDNAME}
echo $OCNDOM
echo $ATMDOM
export CDATE=c190627
export GRIDFILE=${CESM_SRC_DIR}/tools/mkmapgrids/SCRIPgrid_${GRIDNAME}_noocean_${CDATE}.nc
export ESMFBIN_PATH=/pic/projects/climate/huan565/esmf/install/netcdf/bin/binO/Linux.pgi.64.mvapich2.default
export REGRID_PROC=1
#define the directories, files, cases
export CSMDATA=/pic/projects/climate/huan565/inputdata/user_inputdata
export REGRID_PROC=1
#export interactive=NO

# Change to the mkmapdata directory
cd ${CESM_SRC_DIR}/tools/mkmapdata

#create mapping files
./mkmapdata.sh -f ${GRIDFILE} -r ${GRIDNAME} -t regional -v
