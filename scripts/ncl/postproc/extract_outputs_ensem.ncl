;******************************************************************
;  Compare simulated energy fluxes with observations 
;    Author: maoyi.huang@pnnl
;    Date :  11/30/2017
;******************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
 
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

;define variables
 site         = "km83"
 year_start   = 2980
 year_end     = 2990
 path_archive = "/pic/scratch/huan565/cesm_archive/"
 path_o = "./outputs_spinup/"
 compset      = "I2000Clm50FatesGs"
 nv      = 90
 exphdr  = "intact.ensm"
 expnum  = ispan(1, nv, 1) 
 npft         = 2;  2 PFTs
 ny           = year_end - year_start +1 ;# of nyear
 nm           = ny * 12
  
 ;****************************************************
 ; Modify variables for comparison
  ;carbon cycle variables
  YEAR = new(nm, integer)
  MON  = new(nm, integer)
  ER   = new((/nv, nm/), float); ecosystem respiration, autotrophic + heterotrophic (gC/m^2/s)
  HR   = new((/nv, nm/), float); heterotropic respiration (gC/m^2/s)
  AR   = new((/nv, nm/), float); autotrophic respiration (gC/m^2/s)
  AGB  = new((/nv, nm/), float); above ground biomass gC/m2
  AGB_pft  = new((/nv,nm,npft/), float); above ground biomass for PFTs gC/m2
  ED_biomass   = new((/nv, nm/), float); Total biomass (gC m-2)
  ED_LITTER_STOCK_COL  = new((/nv, nm/), float); total ED litter carbon at the column level (gC/m^2)
  CWD   = new((/nv, nm/), float); total CWD carbon at the column level(gC/m^2)
  TOTSOMC_1m = new((/nv, nm/), float); total soil organic matter carbon to 1 meter depth (gC/m^2)
  TOTSOMC = new((/nv, nm/), float); total soil organic matter carbon (gC/m^2)
  ED_bleaf = new((/nv, nm/), float); leaf biomas
  NEE  = new((/nv, nm/), float); net ecosystem exchange (gC/m^2/s)
  GPP  = new((/nv, nm/), float); gross primary production (gC/m^2/s)
  GPP_CANOPY = new((/nv, nm/), float); gross primary production of canopy plants (gC/m^2/s)
  GPP_UNDERSTORY = new((/nv, nm/), float); gross primary production of understory plants (gC/m^2/s)
  NPP  = new((/nv, nm/), float); net primary production (gC/m^2/s)
  ; energy variablesa
  LH   = new((/nv, nm/), float); latent heat flux, w/m^2
  SH   = new((/nv, nm/), float); sensible heat flux, w/m^2
  GH   = new((/nv, nm/), float); ground(soil) heat flux, w/m^2
  Rn   = new((/nv, nm/), float); net radiation, w/m^2 
  LAI  = new((/nv, nm/), float); leaf area index
  SOILWATER_10CM = new((/nv, nm/), float); volumetric soil moisture, first 10cm (kg/m2)
  BTRAN = new((/nv, nm/), float); transpiration beta factor (unitless)
  TWS = new((/nv, nm/), float); total water storage (mm)
  TSOI_10CM = new((/nv, nm/), float); soil temperature in top 10cm of soil (K)
  fcev = new((/nv, nm/), float); canopy evaporation, w/m^2
  fctr = new((/nv, nm/), float); canopy transpiration, w/m^2
  fgev = new((/nv, nm/), float); ground evaporation, w/m^2
  fsds = new((/nv, nm/), float); downward shortwave radiation, w/m^2
  flds = new((/nv, nm/), float); downward longwave radiation, w/m^2
  fire = new((/nv, nm/), float); emitted longwave radiation, w/m^2
  fsr = new((/nv, nm/), float); reflected shortwave radiaton, w/m^2
  fsa = new((/nv, nm/), float); absorbed shortwave radiation, w/m^2
  RAIN = new((/nv, nm/), float); atmospheric rain, mm/s
  QOVER = new((/nv, nm/), float); surface runoff, mm/s
  QINFL = new((/nv, nm/), float); infiltration, mm/s
  QDRIP = new((/nv, nm/), float); throughfall, mm/s
  QDRAI = new((/nv, nm/), float); sub-surface drainage, mm/s
  QCHARGE = new((/nv, nm/), float); aquifer recharge rate, mm/s
  QINTR  = new((/nv, nm/), float); interception, mm/s
  QVEGE = new((/nv, nm/), float);  canopy evaporation, mm/s
  QVEGT = new((/nv, nm/), float);  canopy transpiration, mm/s
  QSOIL = new((/nv, nm/), float);  Ground evaporation, mm/s
	
;retrieve simulations
;**********************************************************
;assign FATES values to the plotting variable, UTC time
;**********************************************************
 factor = 60*60*24; conversion factor from gC/m^2/s to gC/m^2/day or from mm/s ->mm/day
 filename     = new((/nv/), string)
 plot_time    = ispan(0,ny-1,1)
 
do i=0,nv-1
   casename = "fates-clm5_1x1pt"+"_" + site + "_" + exphdr + sprinti("%0.4i", expnum(i)) + "_spinup_" + compset
   do iy=0,ny-1
      do im=2,2
          YEAR(iy*12:(iy+1)*12-1) = iy + year_start
          MON(iy*12:(iy+1)*12-1)  = ispan(1,12,1)
         ;read file`
          filename(i)=path_archive + casename + "/lnd/hist/"+ casename + ".clm2.h0."+ sprinti("%0.4i", iy+year_start) + "-"+ sprinti("%0.2i", im) + "-01-00000.nc"
          print(filename(i))
          c = addfiles(filename(i),"r")

          ;read value from  FATES history
          fcev(i, iy*12:(iy+1)*12-1) = c[:]->FCEV(:,  0)
          fctr(i, iy*12:(iy+1)*12-1) = c[:]->FCTR(:,  0)
          fgev(i, iy*12:(iy+1)*12-1) = c[:]->FGEV(:,  0)
          flds(i, iy*12:(iy+1)*12-1) = c[:]->FLDS(:,  0)
          fire(i, iy*12:(iy+1)*12-1) = c[:]->FIRE(:,  0)
          fsr(i,  iy*12:(iy+1)*12-1)  = c[:]->FSR(:,  0)
          fsds(i, iy*12:(iy+1)*12-1) = c[:]->FSDS(:,  0)
          fsa(i,  iy*12:(iy+1)*12-1)  = c[:]->FSA(:,  0)
          LH(i,  iy*12:(iy+1)*12-1)   = c[:]->FCEV(:,  0) + c[:]->FCTR(:,  0) + c[:]->FGEV(:,  0)
          GH(i, iy*12:(iy+1)*12-1)   = c[:]->FGR(:,  0)
          SH(i, iy*12:(iy+1)*12-1)   = c[:]->FSH(:,  0)
          Rn(i, iy*12:(iy+1)*12-1)   = fsds(i, iy*12+im-1)-fsr(i, iy*12+im-1)+flds(i, iy*12+im-1)-fire(i, iy*12+im-1)
          GPP(i, iy*12:(iy+1)*12-1)   = c[:]->GPP(:,  0)*factor
          GPP_CANOPY(i, iy*12:(iy+1)*12-1)   = c[:]->GPP_CANOPY(:,  0)*factor
          GPP_UNDERSTORY(i, iy*12:(iy+1)*12-1) = c[:]->GPP_UNDERSTORY(:,  0)*factor
          NPP(i, iy*12:(iy+1)*12-1)   = c[:]->NPP(:,  0)*factor
          NEE(i, iy*12:(iy+1)*12-1)   = -c[:]->NEP(:,  0)*factor
          CWD(i, iy*12:(iy+1)*12-1)   = c[:]->CWD_STOCK_COL(:,  0)
          ED_biomass(i, iy*12:(iy+1)*12-1)   = c[:]->ED_biomass(:,  0)
          ED_LITTER_STOCK_COL(i, iy*12:(iy+1)*12-1)   = c[:]->ED_LITTER_STOCK_COL(:, 0)
          ED_bleaf(i, iy*12:(iy+1)*12-1)   = sum(c[:]->PFTleafbiomass(:,:, 0))
          AGB(i, iy*12:(iy+1)*12-1)   = c[:]->AGB(:,0)
          AGB_pft(i,iy*12:(iy+1)*12-1,:)   = 0.6*(c[:]->PFTbiomass(:,:,0) - c[:]->PFTleafbiomass(:,:,0)*2) + c[:]->PFTleafbiomass(:,:,0)
          HR(i, iy*12:(iy+1)*12-1)   = c[:]->HR(:,  0)*factor
          AR(i, iy*12:(iy+1)*12-1)   = c[:]->AR(:,  0)*factor
          ER(i, iy*12:(iy+1)*12-1)   = (c[:]->HR(:,  0) + c[:]->AR(:,  0))*factor
          SOILWATER_10CM(i, iy*12:(iy+1)*12-1)   = c[:]->SOILWATER_10CM(:,  0)
          BTRAN(i, iy*12:(iy+1)*12-1)   = c[:]->BTRAN(:,  0)
          TWS(i, iy*12:(iy+1)*12-1)   = c[:]->TWS(:,  0)
          TSOI_10CM(i, iy*12:(iy+1)*12-1)   = c[:]->TSOI_10CM(:,  0)
          TOTSOMC_1m(i, iy*12:(iy+1)*12-1)   = c[:]->TOTSOMC_1m(:,  0)
          TOTSOMC(i, iy*12:(iy+1)*12-1)   = c[:]->TOTSOMC(:,  0)
          LAI(i, iy*12:(iy+1)*12-1)   = c[:]->TLAI(:,  0)
          RAIN(i, iy*12:(iy+1)*12-1)   = c[:]->RAIN(:,  0)*factor
          QOVER(i, iy*12:(iy+1)*12-1)   = c[:]->QOVER(:,  0)*factor
          QINFL(i, iy*12:(iy+1)*12-1)   = c[:]->QINFL(:,  0)*factor
          QDRIP(i, iy*12:(iy+1)*12-1)   = c[:]->QDRIP(:,  0)*factor
          QDRAI(i, iy*12:(iy+1)*12-1)   = c[:]->QDRAI(:,  0)*factor
          QCHARGE(i, iy*12:(iy+1)*12-1)   = c[:]->QCHARGE(:,  0)*factor
          QINTR(i, iy*12:(iy+1)*12-1)   = c[:]->QINTR(:,  0)*factor
          QVEGE(i, iy*12:(iy+1)*12-1)   = c[:]->QVEGE(:,  0)*factor
          QSOIL(i, iy*12:(iy+1)*12-1)   = c[:]->QSOIL(:,  0)*factor
          QVEGT(i, iy*12:(iy+1)*12-1)   = c[:]->QVEGT(:,  0)*factor
          delete(c)
      end do; mon   
   end do; year
   delete(casename)
end do; nv
	
;*************************************************
; write variables into csv files
;*************************************************
lines = new(nm+1, string)
lines(0)=(/"  yr,mn,    intact01,    intact02,    intact03,    intact04,    intact05,    intact06,    intact07,    intact08,    intact09,    intact10,    intact11,    intact12,    intact13,    intact14,    intact15,    intact16,    intact17,    intact18,    intact19,    intact20,    RILlow,   RILhigh,     CLlow,    CLhigh"/); header

fnmo = path_o+"CWD_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",CWD(i,:))
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"AGB_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",AGB(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"AGB_early_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",AGB_pft(i,:,0)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"AGB_late_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",AGB_pft(i,:,1))
end do 
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"ED_biomass_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm)+ "," + sprintf("%10.2f",ED_biomass(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"ED_LITTER_STOCK_COL_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",ED_LITTER_STOCK_COL(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"ED_bleaf_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",ED_bleaf(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"TOTSOMC_1m_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",TOTSOMC_1m(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"TOTSOMC_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",TOTSOMC(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"NEE_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",NEE(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"GPP_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",GPP(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"GPP_CANOPY_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",GPP_CANOPY(i,:))
end do 
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"GPP_UNDERSTORY_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",GPP_UNDERSTORY(i,:)) 
end do 
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"NPP_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",NPP(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"ER_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",ER(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"HR_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",HR(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"AR_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",AR(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"LH_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",LH(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"SH_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",SH(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fcev_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fcev(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fgev_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fgev(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fctr_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fctr(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"flds_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",flds(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fsds_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fsds(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fsr_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fsr(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fsa_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fsa(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"fire_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",fire(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"Rn_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",Rn(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"GH_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",GH(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"LAI_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",LAI(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"SOILWATER_10CM_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",SOILWATER_10CM(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"BTRAN_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",BTRAN(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"TWS_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",TWS(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"TSOI_10CM_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",TSOI_10CM(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"RAIN_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",RAIN(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QOVER_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QOVER(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QINFL_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QINFL(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QDRIP_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QDRIP(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QDRAI_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QDRAI(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QCHARGE_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QCHARGE(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QINTR_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QINTR(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QVEGE_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QVEGE(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QVEGT_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QVEGT(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

fnmo = path_o+"QSOIL_"+site+".csv"
system("rm -rf " + fnmo)
lines(1:nm) = sprinti("%0.4i",YEAR) + "," + sprinti("%0.2i",MON)
do i=0, nv-1
   lines(1:nm) = lines(1:nm) + "," + sprintf("%10.2f",QSOIL(i,:)) 
end do
asciiwrite(fnmo,lines)
delete(fnmo)

end
