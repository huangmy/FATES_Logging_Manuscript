#!/bin/sh 
#Sample shell script to modify fates param files
#M.Huang@PNNL, 03/28/2017
#Last modified on 07/01/2019 for constance
#

export BASE_DIR=/pic/projects/landuq/huan565/ngeet

#define the directories, files, cases
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export INPUTDATA_DIR=${BASE_DIR}/inputdata
export CESM_RES=1x1pt_km83
export CYYMMDD=c190707
export EXPERIMENT=intact

# Change to the makesurfdata_map src directory
cd ${CESM_SRC_DIR}/src/fates/tools

# create a new fates param file with two tropical evergreen PFTs:
python tmp.py --pft-indices=1,1 --fin=${INPUTDATA_DIR}/user_inputdata/1x1pt_km83/fates_params_default_c190705.nc --fout=${INPUTDATA_DIR}/user_inputdata/1x1pt_km83/fates_params_default_2troppftclones_${CYYMMDD}.nc

#modify the pft-specific parameters
#python modify_fates_paramfile.py --fin=${INPUTDATA_DIR}/user_inputdata/1x1pt_km83/fates_params_api.7.3.0_2troppftclones_${CYYMMDD}.nc --fout=${INPUTDATA_DIR}/user_inputdata/1x1pt_km83/fates_params_api.7.3.0_2troppftclones_${CYYMMDD .${EXPERIMENT}.nc --pft=1 :
#create fates parameter files for the logging experiments
#mv  ${CESM_SRC_DIR}/tools/mksurfdata_map/surfdata_${CESM_RES}_78pfts_simyr2000_${CYYMMDD}.nc ${INPUTDATA_DIR}/user_inputdata/1x1pt_km83

