# FATES_Logging_Manuscript
Steps to check out codes, and driving/preprocessing/analysis scripts used in the FATES Logging Manuscript

## Repository structure
---scripts 

## Tutorial to configure the CLM(FATES)-logging experiments on CONSTANCE
We provide detailed notes on running the logging experiments on PNNL's CONSTANCE cluster.

### Download script  repository
    setenv BASE_DIR <dir-of-choice>
    cd $BASE_DIR
    git clone git@github.com:huangmy/FATES_Logging_Manuscript.git
  
### Check out CLM and FATES codes, please check http://www.cesm.ucar.edu/models/cesm2.0/land and https://fates-docs.readthedocs.io/en/latest/index.html for CLM5 and FATES documentations

#### CLM codes with fates api capatible with https://github.com/NGEET/fates/releases/tag/sci.1.27.2_api.7.3.0
    cd $BASE_DIR
    git clone -b fates_next_api https://github.com/ESCOMP/ctsm.git clm5.0
    setenv CLM_SRC_DIR $BASE_DIR/clm5.0
    cd $CLM_SRC_DIR

#### Update CLM externals to the fates branch used
    cp -f $BASE_DIR/scripts/shell/user_machines/Externals_CLM.cfg .

#### Check out FATES codes, specific branch used is https://github.com/huangmy/fates/tree/fix_export_frac, with minor update from M. Huang to fates release tag sci.1.27.2_api.7.3.0. The update has been integrated to newer fates releases
    ./manage_externals/checkout_externals

### Downlaod input and output datasets
    Click the link: https://drive.google.com/drive/folders/1ufnxkE2aw6bueLG7LttR-N8UrStLgCqZ?usp=sharing
    Unzip and untar inputdata.tar.gz and outputdata.tar.gz in $BASE_DIR

### Configure a user_defined single point CLM5 simulation
    cd $BASE_DIR/scripts/shell
#### Spin up the intact forest
    bash create_fates-clm_1x1pt_km83_intact_spinup.sh
#### Intact Simulation
    bash create_fates-clm_1x1pt_km83_intact_sim.sh 
#### Logging Experiments
    bash create_fates-clm_1x1pt_km83_logging.sh

### Extract and analyze your simulations
    cd $BASE_DIR/scripts/ncl/postproc
#### Extract ouputs
    ncl extract_outputs.ncl
#### Plotting carbon cycle variables
    ncl plot_Carbonflux_km83_bg.ncl
    ncl plot_Carbonpools_km83_bg.ncl
#### Plotting energy and water variables
    ncl plot_energywater_km83_bg.ncl
#### Plotting forest structure variables
    ncl check_PFTcomposition_sim_ts4sc_bg.ncl 

## Who do I talk to?
   Maoyi Huang (huangmy at gmail.com)

## Reference
Huang, M., Xu, Y., Longo, M., Keller, M., Knox, R., Koven, C., and Fisher, R.: Assessing impacts of selective logging on water, energy, and carbon budgets and ecosystem dynamics in Amazon forests using the Functionally Assembled Terrestrial Ecosystem Simulator, Biogeosciences Discuss., https://doi.org/10.5194/bg-2019-129, in review, 2019.


## Acknowledgment
This research was supported by The Next-Generation Ecosystem Experiments – Tropics project through the Terrestrial Ecosystem Science (TES) program within US Department of Energy’s Office of Biological and Environmental Research (BER).
