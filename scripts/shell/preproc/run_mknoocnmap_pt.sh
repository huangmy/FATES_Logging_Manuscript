#!/bin/sh
#Created by M.Huang on 09/09/2018

export GRIDNAME=1x1pt_km83
export CDATE=`date +%y%m%d`
export BASE_DIR=/pic/projects/landuq/huan565/ngeet
export CESM_SRC_DIR=${BASE_DIR}/clm5.0
export USER_INPUTDATA_DIR=${BASE_DIR}/inputdata/user_inputdata


# Create the SCRIP grid file for the single point location and
# create a unity mapping file for it
export lat=-3.018
export lon=305.00

cd ${CESM_SRC_DIR}/tools/mkmapdata
./mknoocnmap.pl -p ${lat},${lon} -n ${GRIDNAME}

# Move the created files to inputdata directory
mkdir ${USER_INPUTDATA_DIR}/${GRIDNAME}
mkdir ${USER_INPUTDATA_DIR}/${GRIDNAME}/maps
mkdir ${USER_INPUTDATA_DIR}/${GRIDNAME}/grids

# Check the data filename
ls ${CESM_SRC_DIR}/tools/mkmapdata/map_${GRIDNAME}_noocean_to_${GRIDNAME}_nomask_aave_da_${CDATE}.nc
ls ${CESM_SRC_DIR}/tools/mkmapgrids/SCRIPgrid_${GRIDNAME}_nomask_c${CDATE}.nc
mv ${CESM_SRC_DIR}/tools/mkmapdata/map_${GRIDNAME}_noocean_to_${GRIDNAME}_nomask_aave_da_${CDATE}.nc ${USER_INPUTDATA_DIR}/${GRIDNAME}/maps/
mv ${CESM_SRC_DIR}/tools/mkmapgrids/SCRIPgrid_${GRIDNAME}_nomask_c${CDATE}.nc ${USER_INPUTDATA_DIR}/${GRIDNAME}/grids/

