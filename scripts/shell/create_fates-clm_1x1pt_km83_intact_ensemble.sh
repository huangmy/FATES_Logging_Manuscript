#!/bin/sh
#############################################################################################################################
#  Purpose: spinning up fates-CLM at a single-point site (KM83) on corip1 
#  Maoyi Huang
#  Date: 29 September 2017
#  Modified to run on constance, PNNL Institutional Computing, 06/24/2019
#  Modified for ensemble simuluation following ckoven's bci script, 07/11/2019
#
##############################################################################################################################

export CCSMUSER=huan565
export BASE_DIR=/pic/projects/landuq/${CCSMUSER}/ngeet
export PROJECT_DIR=/pic/scratch/${CCSMUSER}
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export INPUTDATA_DIR=${BASE_DIR}/inputdata
export CESM_INPUTDATA_DIR=${INPUTDATA_DIR}/cesm_inputdata
export ARCHIVE_DIR=${PROJECT_DIR}/cesm_archive
export CIME_OUTPUT_ROOT=${PROJECT_DIR}/csmruns

export CESM_COMPSET=I2000Clm50FatesGs
export RES=CLM_USRDAT
export CESM_CASE_DIR=${BASE_DIR}/CESM_cases

export CLM_USRDAT_NAME=1x1pt_km83
export site_forcing=1x1pt_km67
export experiment=intact
export CESM_CASE_NAME_SPINUP=fates-clm5_${CLM_USRDAT_NAME}_${experiment}_spinup_${CESM_COMPSET}
export CESM_CASE_NAME=fates-clm5_${CLM_USRDAT_NAME}_${experiment}_${CESM_COMPSET}
export ARCHIVE_DIR=${PROJECT_DIR}/cesm_archive
export ninst=20    #number of ensemble members
export spnp_yr=2251 #year for the initial condition
#export fpar_fates=fates_params_default_2troppftclones_c190707.${experiment}.nc

# check the model versions
cd ${CESM_SRC_DIR}
GITHASH1=`git log -n 1 --format=%h`
cd src/fates
GITHASH2=`git log -n 1 --format=%h`

# delete the old case
rm -r ${CESM_CASE_DIR}/${CESM_CASE_NAME}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create soft links for CESM inputdata
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++ Link forcing data file, BCI
mkdir -p ${CESM_INPUTDATA_DIR}/atm/datm7/${CLM_USRDAT_NAME}/CLM1PT_data
#rm -f ${CESM_INPUTDATA_DIR}/atm/datm7/${CLM_USRDAT_NAME}/CLM1PT_data/*.nc
cd ${INPUTDATA_DIR}/user_inputdata/${site_forcing}/CLM1PT_data/  # 2001-2011
ls *.nc | while read filename
do
        echo ${filename}
        ln -s ${INPUTDATA_DIR}/user_inputdata/${site_forcing}/CLM1PT_data/${filename} ${CESM_INPUTDATA_DIR}/atm/datm7/${CLM_USRDAT_NAME}/CLM1PT_data/${filename}
done
cd ${BASE_DIR}/scripts/shell

mkdir -p ${CESM_INPUTDATA_DIR}/share/domains/domain.clm
rm -rf ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${CLM_USRDAT_NAME}_navy.nc
ln -s ${INPUTDATA_DIR}/user_inputdata/${CLM_USRDAT_NAME}/domain.lnd.${CLM_USRDAT_NAME}_navy.nc ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${CLM_USRDAT_NAME}_navy.nc
rm -rf ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${site_forcing}_navy.nc
#ln -s ${INPUTDATA_DIR}/user_inputdata/${site_forcing}/domain.lnd.${site_forcing}_navy.nc ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${site_forcing}_navy.nc

mkdir -p ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map
rm -rf ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/surfdata_${CLM_USRDAT_NAME}_simyr2000.nc
ln -s ${INPUTDATA_DIR}/user_inputdata/${CLM_USRDAT_NAME}/surfdata_${CLM_USRDAT_NAME}_simyr2000_c190701.nc ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/surfdata_${CLM_USRDAT_NAME}_simyr2000.nc

#modifiy the machine file
cp -f  ${BASE_DIR}/scripts/shell/user_machines/*constance*.xml ${CESM_SRC_DIR}/cime/config/cesm/machines

# Change to directory
cd ${CESM_SRC_DIR}/cime/scripts

# Creating case with command :
./create_newcase --case ${CESM_CASE_DIR}/${CESM_CASE_NAME} --res CLM_USRDAT --compset ${CESM_COMPSET} --mach constance --ninst=$ninst --multi-drive --run-unsupported

#+++ Configuring case :
cd ${CESM_CASE_DIR}/${CESM_CASE_NAME}

#+++ Modifying : env_batch.xml, if debugging
./xmlchange  --file env_batch.xml --id JOB_QUEUE --val "slurm" --force
./xmlchange  --file env_batch.xml --id JOB_WALLCLOCK_TIME --val "24:00:00"

#+++ Modifying : env_build.xml
./xmlchange  --file env_build.xml --id CIME_OUTPUT_ROOT --val ${CIME_OUTPUT_ROOT}

# Modifying : env_run.xml
./xmlchange --file env_run.xml --id CLM_USRDAT_NAME --val ${CLM_USRDAT_NAME}
./xmlchange --file env_run.xml --id CASESTR --val ${CESM_CASE_NAME}
./xmlchange --file env_run.xml --id ATM_DOMAIN_FILE --val domain.lnd.${CLM_USRDAT_NAME}_navy.nc
./xmlchange --file env_run.xml --id LND_DOMAIN_FILE --val domain.lnd.${CLM_USRDAT_NAME}_navy.nc
./xmlchange --file env_run.xml --id ATM_DOMAIN_PATH --val "\$DIN_LOC_ROOT/share/domains/domain.clm"
./xmlchange --file env_run.xml --id LND_DOMAIN_PATH --val "\$DIN_LOC_ROOT/share/domains/domain.clm"
./xmlchange --file env_run.xml --id DATM_MODE --val CLM1PT
./xmlchange --file env_run.xml --id ATM_NCPL --val 24
./xmlchange --file env_run.xml --id STOP_N --val 100
./xmlchange --file env_run.xml --id STOP_OPTION --val nyears
./xmlchange --file env_run.xml --id RUN_STARTDATE --val '2001-01-01'
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_START --val 2001
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val 2011
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_ALIGN --val 2001
./xmlchange --file env_run.xml --id REST_OPTION --val nyears
./xmlchange --file env_run.xml --id REST_N --val 1
./xmlchange --file env_run.xml --id DOUT_S_ROOT --val ${ARCHIVE_DIR}/${CESM_CASE_NAME}
./xmlchange --file env_run.xml --id DIN_LOC_ROOT --val ${CESM_INPUTDATA_DIR}
./xmlchange --file env_run.xml --id DIN_LOC_ROOT_CLMFORC --val "\$DIN_LOC_ROOT/atm/datm7"
./xmlchange --file env_run.xml --id DOUT_S_SAVE_INTERIM_RESTART_FILES --val TRUE
#./xmlchange --file env_run.xml --id CLM_FORCE_COLDSTART --val on
./xmlchange --file env_run.xml --id NINST_LND --val $ninst
./xmlchange --file env_run.xml --id ROOTPE_LND --val 1
#./xmlchange --file env_run.xml --id CLM_BLDNML_OPTS --val "-bgc ed -no-megan"

#setup and compile the case
./case.setup

export fsurdat=${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/surfdata_${CLM_USRDAT_NAME}_simyr2000.nc
for x  in `seq 1 1 $ninst`; do
  expstr=$(printf %04d $x)
  echo $expstr
  export fpar_fates=fates_params_default_2troppftclones_c190707.${experiment}.ensm${expstr}.nc
  export fates_paramfile=${INPUTDATA_DIR}/user_inputdata/${CLM_USRDAT_NAME}/${fpar_fates}
  export finidat=/pic/scratch/huan565/cesm_archive/${CESM_CASE_NAME_SPINUP}/rest/${spnp_yr}-01-01-00000/${CESM_CASE_NAME_SPINUP}.clm2_${expstr}.r.${spnp_yr}-01-01-00000.nc
  cat > user_nl_clm_${expstr} <<EOF
  co2_ppmv = 400.0
  finidat = '${finidat} '
  fates_paramfile = '${fates_paramfile}'  
  fsurdat = '${fsurdat}'
  hist_empty_htapes = .false.
  use_luna = .false.
  use_fates= .true.
  use_fates_logging = .false.
  hist_mfilt = 12, 10
  hist_nhtfrq = 0, -8760
  hist_fincl2 = 'ZSTAR_BY_AGE','RECRUITMENT','PFTbiomass','PATCH_AREA_BY_AGE','NPLANT_UNDERSTORY_SCPF','NPLANT_SCPF','NPLANT_SCAG','NPLANT_CANOPY_SCPF','MORTALITY_UNDERSTORY_SCPF','MORTALITY_CANOPY_SCPF','MORTALITY','M8_SCPF','M7_SCPF','M6_SCPF','M5_SCPF','M4_SCPF','M3_SCPF','M2_SCPF','M1_SCPF','LEAF_HEIGHT_DIST','LAI_UNDERSTORY_SCLS','LAI_CANOPY_SCLS','LAI_BY_AGE','FUEL_MOISTURE_NFSC','DDBH_UNDERSTORY_SCPF','DDBH_CANOPY_SCPF','CROWNAREA_CAN','CANOPY_HEIGHT_DIST','CANOPY_AREA_BY_AGE','BA_SCPF','GROWTHFLUX_SCPF','GROWTHFLUX_FUSION_SCPF'
EOF
done

# Modify datm streams
cat >> user_nl_datm << EOF
  taxmode = 'cycle','cycle'
EOF

#build the case
./case.build

#+++ Modify datm streams
#cp -f  ${CESM_CASE_DIR}/${CESM_CASE_NAME}/CaseDocs/datm.streams.txt.CLM1PT.CLM_USRDAT ${CESM_CASE_DIR}/${CESM_CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT
#chmod +rw ${CESM_CASE_DIR}/${CESM_CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT

./case.submit

