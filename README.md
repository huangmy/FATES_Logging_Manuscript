# FATES_Logging_Manuscript
steps to check out codes, and driving/preprocessing/analysis scripts used in the FATES Logging Manuscript

## Repository structure
---scripts 

## Tutorial to configure the CLM(FATES)-logging experiments on CONSTANCE
We provide detailed notes on running the logging experiments on PNNL's CONSTANCE cluster.

### Download script and data repository
    setenv BASE_DIR <dir-of-choice>
    cd $BASE_DIR
    git clone git@github.com:huangmy/FATES_Logging_Manuscript.git
  
### Download CLM code, please check http://www.cesm.ucar.edu/models/cesm2.0/land for CLM5 documentation
    cd $BASE_DIR
    git clone -b release-clm5.0 https://github.com/ESCOMP/ctsm.git clm5.0
    setenv CLM_SRC_DIR $BASE_DIR/clm5.0
    cd $CLM_SRC_DIR
    ./manage_externals/checkout_externals

### Download FATES code, please check for FATES documentation

### Downlaod input and output datasets
    Click the link: https://drive.google.com/drive/folders/1ufnxkE2aw6bueLG7LttR-N8UrStLgCqZ?usp=sharing
    Unzip and untar inputdata.tar.gz and outputdata.tar.gz in $BASE_DIR

### Configure a user_defined single point CLM5 simulation
    cd $BASE_DIR/scripts/shell
    bash create_fates-clm_1x1pt_km83_intact_spinup.sh
    bash create_fates-clm_1x1pt_km83_intact_sim.sh 
    bash create_fates-clm_1x1pt_km83_logging.sh

### Extract and analyze your simulations
    cd $BASE_DIR/scripts/ncl/postproc
    ncl extract_outputs.ncl
    ncl plot_Carbonflux_km83_bg.ncl
    ncl plot_Carbonpools_km83_bg.ncl
    ncl plot_energywater_km83_bg.ncl

## Who do I talk to?
    maoyi.huang at pnnl.gov

## Reference
Huang, M., Xu, Y., Longo, M., Keller, M., Knox, R., Koven, C., and Fisher, R.: Assessing impacts of selective logging on water, energy, and carbon budgets and ecosystem dynamics in Amazon forests using the Functionally Assembled Terrestrial Ecosystem Simulator, Biogeosciences Discuss., https://doi.org/10.5194/bg-2019-129, in review, 2019.


## Acknowledgment
This research was supported by The Next-Generation Ecosystem Experiments – Tropics project through the Terrestrial Ecosystem Science (TES) program within US Department of Energy’s Office of Biological and Environmental Research (BER).
