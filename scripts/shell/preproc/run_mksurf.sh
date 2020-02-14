#!/bin/sh 
#Sample shell script to configure, compile, and submit a supported simulation
#M.Huang@PNNL, 03/28/2017
#Last modified on 06/27/2019 for constance
#

export BASE_DIR=/pic/projects/landuq/huan565/ngeet

#define the directories, files, cases
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export INPUTDATA_DIR=${BASE_DIR}/inputdata
export CESM_RES=1x1pt_km83
export CYYMMDD=c190628

# Change to the makesurfdata_map src directory
cp -f ${BASE_DIR}/scripts/preproc/user_namelists/mksurfdata_map.namelist_1x1pt_km83 ${CESM_SRC_DIR}/tools/mksurfdata_map
cd ${CESM_SRC_DIR}/tools/mksurfdata_map

# prepare the environ for compilation
export INC_NETCDF=/opt/cray/pe/netcdf/4.4.1.1.3/INTEL/16.0/include
export LIB_NETCDF=/opt/cray/pe/netcdf/4.4.1.1.3/INTEL/16.0/lib


# Running case :
./mksurfdata_map < mksurfdata_map.namelist_1x1pt_km83

#move the surfdata file to user_inputdata directory
mv  ${CESM_SRC_DIR}/tools/mksurfdata_map/surfdata_${CESM_RES}_78pfts_simyr2000_${CYYMMDD}.nc ${INPUTDATA_DIR}/user_inputdata/1x1pt_km83

#delete the namelist file
rm -f ${CESM_SRC_DIR}/tools/mksurfdata_map/mksurfdata_map.namelist_1x1pt_km83
