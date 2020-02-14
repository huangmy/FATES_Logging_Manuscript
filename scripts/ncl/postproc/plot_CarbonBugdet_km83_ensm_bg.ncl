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
 experiment   = "km83"
 year_start   = 2001
 year_end     = 2011
 path_o       = "./outputs_spinup/"
 fnm_obs      = "./inputs/km83_obs_Macros_filter_v3.csv"
 fnm_obs2     = "./inputs/km83_Miller_etal2011.csv"
 fnm_obs67    = "./inputs/km67_obs_Macos_filter_GPP_NEE_v3.csv"
 fnm_obs67_2  = "./inputs/km67_Miller_etal2011.csv"
 compset      = "I1PTCLM45"
 nsim         = 90
 nv           = nsim+4; sim+obs83+obs67+obs83_2+obs67_2
 iexp         = ispan(1, nsim,1)
 npft         = 2;  2 PFTs
 ny           = year_end - year_start +1 ;# of nyear
 nm           = ny * 12
  
 ;****************************************************
 ; Modify variables for comparison
  ;carbon cycle variables
  FillValue = -999.
  ER   = new((/nv, nm/), float, FillValue); ecosystem respiration, autotrophic + heterotrophic (gC/m^2/day)
  HR   = new((/nv, nm/), float, FillValue); heterotropic respiration (gC/m^2/day)
  AR   = new((/nv, nm/), float, FillValue); autotrophic respiration (gC/m^2/day)
  AGB  = new((/nv, nm/), float, FillValue); above ground biomass gC/m2
  AGB_early = new((/nv, nm/), float, FillValue); above ground biomass for early succession PFT (gC/m2)
  AGB_late  = new((/nv, nm/), float, FillValue); above ground biomass for late succession PFT (gC/m2)
  ED_biomass   = new((/nv, nm/), float, FillValue); Total biomass (gC m-2)
  CWD   = new((/nv, nm/), float, FillValue); total CWD carbon at the column level(gC/m^2)
  TOTSOMC_1m = new((/nv, nm/), float, FillValue); total soil organic matter carbon to 1 meter depth (gC/m^2)
  ED_bleaf = new((/nv, nm/), float, FillValue); factor limiting ground evap 
  NEE  = new((/nv, nm/), float, FillValue); net ecosystem exchange (gC/m^2/day)
  GPP  = new((/nv, nm/), float, FillValue); gross primary production (gC/m^2/day)
  GPP_CANOPY = new((/nv, nm/), float, FillValue); gross primary production of canopy trees (gC/m^2/day)
  GPP_UNDERSTORY = new((/nv, nm/), float, FillValue); gross primary production of understory trees (gC/m^2/day)
  NPP  = new((/nv, nm/), float, FillValue); net primary production (gC/m^2/day)
  ; energy variablesa
  LH   = new((/nv, nm/), float, FillValue); latent heat flux, w/m^2
  SH   = new((/nv, nm/), float, FillValue); sensible heat flux, w/m^2
  GH   = new((/nv, nm/), float, FillValue); ground(soil) heat flux, w/m^2
  Rn   = new((/nv, nm/), float, FillValue); net radiation, w/m^2 
  LAI  = new((/nv, nm/), float, FillValue); leaf area index
  SOILWATER_10CM = new((/nv, nm/), float, FillValue); volumetric soil moisture, first 10cm (kg/m2)
  TSOI_10CM = new((/nv, nm/), float, FillValue); soil temperature in top 10cm of soil (K)
  fcev = new((/nv, nm/), float, FillValue); canopy evaporation, w/m^2
  fctr = new((/nv, nm/), float, FillValue); canopy transpiration, w/m^2
  fgev = new((/nv, nm/), float, FillValue); ground evaporation, w/m^2
  fsds = new((/nv, nm/), float, FillValue); downward shortwave radiation, w/m^2
  flds = new((/nv, nm/), float, FillValue); downward longwave radiation, w/m^2
  fire = new((/nv, nm/), float, FillValue); emitted longwave radiation, w/m^2
  fsr = new((/nv, nm/), float, FillValue); reflected shortwave radiaton, w/m^2
  fsa = new((/nv, nm/), float, FillValue); absorbed shortwave radiation, w/m^2
  GPPann = new((/nv, ny/), float); annual GPP MgC/m^2/yr-1
  GPPann_can = new((/nv, ny/), float); annual GPP in the canopy layer MgC/m^2/yr-1
  GPPann_und = new((/nv, ny/), float); annual GPP in the understory MgC/m^2/yr-1
  ERann = new((/nv, ny/), float);  annual ER MgC/m^2/yr-1
  NEEann = new((/nv, ny/), float); annual NEE MgC/m^2/yr-1  
  NPPann = new((/nv, ny/), float); annual NPP MgC/m^2/yr-1
  HRann = new((/nv, ny/), float);  annual HR MgC/m^2/yr-1
  ARann = new((/nv, ny/), float);  annual AR MgC/m^2/yr-1
  YEAR = new(ny, integer)
  ;confidence interval
  CI_Rn  = new((/nv, nm/), float, FillValue); confidence interval of net radiation, w/m^2
  CI_SH  = new((/nv, nm/), float, FillValue); confidence interval of sensible heat flux, w/m^2
  CI_LH  = new((/nv, nm/), float, FillValue); confidence interval of latent heat flux, w/m^2
  CI_NEE = new((/nv, nm/), float, FillValue); confidence interval of NEE, gC/m^2/day 
  CI_GPP = new((/nv, nm/), float, FillValue); confidence interval of GPP, gC/m^2/day
  CI_ER  = new((/nv, nm/), float, FillValue); confidence interval of ER, gC/m^2/day

;retrieve obs
 factor = 1000.*60.*60.*24.; conversion factor from KgC/m^2/s to gC/m^2/day
 delim = ","
 ;need to get rid of obs from previous year starting from line 253 in the file fnm_obs
 nls=253
 lines_obs = asciiread(fnm_obs,-1,"string")
 printVarSummary(lines_obs)
 n_obs = 36
 GPP(nsim,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),6,delim))*factor
 NEE(nsim,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),5,delim))*factor
 ER(nsim,0:n_obs-1)  = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),7,delim))*factor

 lines_obs67 = asciiread(fnm_obs67,-1,"string")
 printVarSummary(lines_obs67)
 n_obs = 132
 GPP(nsim+1,0:n_obs-1) = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),4,delim))*factor
 NEE(nsim+1,0:n_obs-1) = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),3,delim))*factor
 ER(nsim+1,0:n_obs-1)  = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),5,delim))*factor

 ; starting from line 2 for Miller et al. obs for km83
 nls2=1 
 lines_obs2 = asciiread(fnm_obs2,-1,"string")
 n_obs2 = 47
 printVarSummary(lines_obs2)
 LAI(nsim+2,0:n_obs2-1)       = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),3,delim))
 TSOI_10CM(nsim+2,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),20,delim))
 SOILWATER_10CM(nsim+2,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),23,delim))
 Rn(nsim+2,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),10,delim))
 GH(nsim+2,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),26,delim)) 
 SH(nsim+2,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),28,delim))
 LH(nsim+2,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),29,delim))
 NEE(nsim+2,0:n_obs2-1)    = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),33,delim))
 ER(nsim+2,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),34,delim)) 
 GPP(nsim+2,0:n_obs2-1)    = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),37,delim))
 CI_Rn(nsim+2,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),40,delim))
 CI_SH(nsim+2,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),46,delim))
 CI_LH(nsim+2,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),47,delim))
 CI_NEE(nsim+2,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),50,delim))
 CI_ER(nsim+2,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),51,delim))
 CI_GPP(nsim+2,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),52,delim))

 ; starting from line 2 for Miller et al. obs for km67
 ; the 2011 CO2 flux at km67 is problematic, not showning here
 nls2=12
 lines_obs67_2 = asciiread(fnm_obs67_2,-1,"string")
 n_obs2 = 49
 printVarSummary(lines_obs67_2)
; factor_E2LH = 2.5*18 ; mmol/m2/s -> W/m2
; factor_Fc   = 12*1e-6 *24*3600  ; umol/m2/s -> gC/m2/day
 Rn(nsim+3,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),6,delim))
 SH(nsim+3,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),14,delim))
 LH(nsim+3,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),16,delim))
 NEE(nsim+3,nls2-1:nls2+n_obs2-2)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),18,delim))
 ER(nsim+3,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),20,delim))
 GPP(nsim+3,nls2-1:nls2+n_obs2-2)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),23,delim))
 CI_Rn(nsim+3,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),26,delim))
 CI_SH(nsim+3,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),31,delim))
 CI_LH(nsim+3,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),32,delim))
 CI_NEE(nsim+3,nls2-1:nls2+n_obs2-2) = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),33,delim))
 CI_ER(nsim+3,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),36,delim))
 CI_GPP(nsim+3,nls2-1:nls2+n_obs2-2) = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),38,delim))

;retrieve simulations
varlist = (/"ER", "HR", "AR", "AGB", "ED_biomass", "CWD", "TOTSOMC_1m","NEE","GPP","LAI","ED_bleaf","NPP","AGB_early","AGB_late", "GPP_CANOPY", "GPP_UNDERSTORY"/)
nvar= 16
ncol = nsim+2
printVarSummary(ER)
do i=0,nvar-1
   fnm = path_o+varlist(i)+"_"+experiment+".csv"
   print(fnm)
   ;lines = asciiread(fnm,-1,"string")
   ;printVarSummary(lines)
    data = readAsciiTable(fnm, ncol, "float", 1)
    printVarSummary(data)
   do j=0,nsim-1
      if (i.eq.0) then
         ER(j,:) = data(:,j+2)  ;gC/m2/day 
      end if
      if (i.eq.1) then
         HR(j,:) = data(:,j+2)  ;gC/m2/day
      end if
      if (i.eq.2) then
         AR(j,:) = data(:,j+2)  ;gC/m2/day
      end if
      if (i.eq.3) then
         AGB(j,:) = data(:,j+2)/1000.  ;kgC/m2
      end if
      if (i.eq.4) then
         ED_biomass(j,:) = data(:,j+2)/1000.  ;kgC/m2
      end if
      if (i.eq.5) then
         CWD(j,:) = data(:,j+2)/1000.  ;kgC/m2
      end if
      if (i.eq.6) then
         TOTSOMC_1m(j,:) = data(:,j+2)/1000.  ;kgC/m2
      end if
      if (i.eq.7) then
         NEE(j,:) = data(:,j+2)  ;gC/m2/day
      end if
      if (i.eq.8) then
         GPP(j,:) = data(:,j+2)  ;gC/m2/day
      end if
      if (i.eq.9) then
         LAI(j,:) = data(:,j+2)  ;gC/m2/day
      end if
      if (i.eq.10) then
         ED_bleaf(j,:) = data(:,j+2)  ;gC/m2
      end if
      if (i.eq.11) then
         NPP(j,:) = data(:,j+2)  ;gC/m2
      end if
      if (i.eq.12) then
         AGB_early(j,:) = data(:,j+2)/1000  ;kgC/m2
      end if
      if (i.eq.13) then
         AGB_late(j,:) = data(:,j+2)/1000  ;kgC/m2
      end if
      if (i.eq.14) then
         GPP_CANOPY(j,:) = data(:,j+2)/1000  ;kgC/m2
      end if
      if (i.eq.15) then
         GPP_UNDERSTORY(j,:) = data(:,j+2)/1000  ;kgC/m2
      end if
   end do; nv 
end do; nvar

;calculate mean, std, median and ranges of the simulated variables
GPPmean      = dim_avg_n(GPP(0:nsim-1,:),0)
NEEmean      = dim_avg_n(NEE(0:nsim-1,:),0)
ERmean       = dim_avg_n(ER(0:nsim-1,:),0)
NPPmean      = dim_avg_n(NPP(0:nsim-1,:),0)
HRmean       = dim_avg_n(HR(0:nsim-1,:),0)
ARmean       = dim_avg_n(AR(0:nsim-1,:),0)
AGBmean      = dim_avg_n(AGB(0:nsim-1,:),0)
CWDmean      = dim_avg_n(CWD(0:nsim-1,:),0)
AGB_earlymean = dim_avg_n(AGB_early(0:nsim-1,:),0)
AGB_latemean  = dim_avg_n(AGB_late(0:nsim-1,:),0)

GPPstddev    = dim_stddev_n(GPP(0:nsim-1,:),0)
NEEstddev    = dim_stddev_n(NEE(0:nsim-1,:),0)
ERstddev     = dim_stddev_n(ER(0:nsim-1,:),0)
NPPstddev    = dim_stddev_n(NPP(0:nsim-1,:),0)
HRstddev     = dim_stddev_n(HR(0:nsim-1,:),0)
ARstddev     = dim_stddev_n(AR(0:nsim-1,:),0)
AGBstddev    = dim_stddev_n(AGB(0:nsim-1,:),0)
CWDstddev    = dim_stddev_n(CWD(0:nsim-1,:),0)
AGB_earlystddev = dim_stddev_n(AGB_early(0:nsim-1,:),0)
AGB_latestddev  = dim_stddev_n(AGB_late(0:nsim-1,:),0)


GPPmedian   = dim_median_n(GPP(0:nsim-1,:),0)
NEEmedian   = dim_median_n(NEE(0:nsim-1,:),0)
ERmedian    = dim_median_n(ER(0:nsim-1,:),0)
NPPmedian   = dim_median_n(NPP(0:nsim-1,:),0)
HRmedian    = dim_median_n(HR(0:nsim-1,:),0)
ARmedian    = dim_median_n(AR(0:nsim-1,:),0)
AGBmedian   = dim_median_n(AGB(0:nsim-1,:),0)
CWDmedian   = dim_median_n(CWD(0:nsim-1,:),0)
AGB_earlymedian = dim_median_n(AGB_early(0:nsim-1,:),0)
AGB_latemedian  = dim_median_n(AGB_late(0:nsim-1,:),0)

GPPmax      = dim_max_n(GPP(0:nsim-1,:),0)
NEEmax      = dim_max_n(NEE(0:nsim-1,:),0)
ERmax       = dim_max_n(ER(0:nsim-1,:),0)
NPPmax      = dim_max_n(NPP(0:nsim-1,:),0)
HRmax       = dim_max_n(HR(0:nsim-1,:),0)
ARmax       = dim_max_n(AR(0:nsim-1,:),0)
AGBmax      = dim_max_n(AGB(0:nsim-1,:),0)
CWDmax      = dim_max_n(CWD(0:nsim-1,:),0)
AGB_earlymax = dim_max_n(AGB_early(0:nsim-1,:),0)
AGB_latemax  = dim_max_n(AGB_late(0:nsim-1,:),0)

GPPmin      = dim_min_n(GPP(0:nsim-1,:),0)
NEEmin      = dim_min_n(NEE(0:nsim-1,:),0)
ERmin       = dim_min_n(ER(0:nsim-1,:),0)
NPPmin      = dim_min_n(NPP(0:nsim-1,:),0)
HRmin       = dim_min_n(HR(0:nsim-1,:),0)
ARmin       = dim_min_n(AR(0:nsim-1,:),0)
AGBmin      = dim_min_n(AGB(0:nsim-1,:),0)
CWDmin      = dim_min_n(CWD(0:nsim-1,:),0)
AGB_earlymin = dim_min_n(AGB_early(0:nsim-1,:),0)
AGB_latemin  = dim_min_n(AGB_late(0:nsim-1,:),0)

; calculate annual statistics
factor2 = 365./100.; gC m-2 d-1 -> MgC ha-1 yr-1 
do iy = 0, ny-1
   GPPann(:,iy) = dim_avg(GPP(:, iy*12:(iy+1)*12-1))*factor2
   GPPann_can(:,iy) = dim_avg(GPP_CANOPY(:, iy*12:(iy+1)*12-1))*factor2
   GPPann_und(:,iy) = dim_avg(GPP_UNDERSTORY(:, iy*12:(iy+1)*12-1))*factor2
   NEEann(:,iy) = dim_avg(NEE(:, iy*12:(iy+1)*12-1))*factor2
   ERann(:,iy)  = dim_avg(ER(:, iy*12:(iy+1)*12-1))*factor2
   NPPann(:,iy)  = dim_avg(NPP(:, iy*12:(iy+1)*12-1))*factor2
   HRann(:,iy)  = dim_avg(HR(:, iy*12:(iy+1)*12-1))*factor2
   ARann(:,iy)  = dim_avg(AR(:, iy*12:(iy+1)*12-1))*factor2
end do
GPPpre6m = dim_avg(GPP(:, 3:8))*factor2
NEEpre6m = dim_avg(NEE(:, 3:8))*factor2
ERpre6m = dim_avg(ER(:, 3:8))*factor2

;write to ascii files
lines = new(ny+1, string)
lines(0)=(/"  yr,    intact,    RILlow,   RILhigh,     CLlow,    CLhigh"/); header
YEAR = ispan(year_start,year_end,1)
fmtf   = ny + "f10.2,"
opt = True
opt@title  = "simulations"
opt@tspace = 22

; output the annual statistics
fnmo = path_o+"GPPann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(GPPann(0:nsim-1,:), fmtf, opt)

fnmo = path_o+"GPPann_can"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(GPPann_can, fmtf, opt)

fnmo = path_o+"GPPann_und"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(GPPann_und, fmtf, opt)

fnmo = path_o+"NEEann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(NEEann, fmtf, opt)

fnmo = path_o+"ERann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(ERann, fmtf, opt)

fnmo = path_o+"NPPann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(NPPann, fmtf, opt)

fnmo = path_o+"HRann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(HRann, fmtf, opt)

fnmo = path_o+"ARann"+experiment+".csv"
opt@fout = fnmo
system("rm -rf " + fnmo)
write_matrix(ARann, fmtf, opt)

;calculate statistics for annual variables
GPPannmean   = dim_avg_n(GPPann(0:nsim-1,:),0)
NEEannmean   = dim_avg_n(NEEann(0:nsim-1,:),0)
ERannmean    = dim_avg_n(ERann(0:nsim-1,:),0)
NPPannmean   = dim_avg_n(NPPann(0:nsim-1,:),0)
HRannmean    = dim_avg_n(HRann(0:nsim-1,:),0)
ARannmean    = dim_avg_n(ARann(0:nsim-1,:),0)

GPPannmedian   = dim_median_n(GPPann(0:nsim-1,:),0)
NEEannmedian   = dim_median_n(NEEann(0:nsim-1,:),0)
ERannmedian    = dim_median_n(ERann(0:nsim-1,:),0)
NPPannmedian   = dim_median_n(NPPann(0:nsim-1,:),0)
HRannmedian    = dim_median_n(HRann(0:nsim-1,:),0)
ARannmedian    = dim_median_n(ARann(0:nsim-1,:),0)

GPPannmax      = dim_max_n(GPPann(0:nsim-1,:),0)
NEEannmax      = dim_max_n(NEEann(0:nsim-1,:),0)
ERannmax       = dim_max_n(ERann(0:nsim-1,:),0)
NPPannmax      = dim_max_n(NPPann(0:nsim-1,:),0)
HRannmax       = dim_max_n(HRann(0:nsim-1,:),0)
ARannmax       = dim_max_n(ARann(0:nsim-1,:),0)

GPPannmin      = dim_min_n(GPPann(0:nsim-1,:),0)
NEEannmin      = dim_min_n(NEEann(0:nsim-1,:),0)
ERannmin       = dim_min_n(ERann(0:nsim-1,:),0)
NPPannmin      = dim_min_n(NPPann(0:nsim-1,:),0)
HRannmin       = dim_min_n(HRann(0:nsim-1,:),0)
ARannmin       = dim_min_n(ARann(0:nsim-1,:),0)


;************************************************
; create plots
;************************************************
; define variables for the plots
  wks = gsn_open_wks("pdf",path_o+"CarbonFluxes"+experiment)
  gsn_define_colormap(wks,"gui_default")        ; select color map
  plot = new(6, graphic)
  wks1 = gsn_open_wks("pdf",path_o+"CarbonPools"+experiment)
  gsn_define_colormap(wks1,"gui_default")        ; select color map
  plot1 = new(6, graphic)

  ;plotting properties 1
  xint=12   
  xlb=0
  xub=10*xint-1
  res = True
  res@vpWidthF = 0.65                           ; change aspect ratio of plot
  res@vpHeightF = 0.2                           ;
  res@tiXAxisString = "Year"                    ; x-axis title
  res@trXMinF  = xlb
  res@trXMaxF  = xub
  res@tmXBMinorPerMajor = 0
  res@tmXBMode = "Explicit"
  res@tmXBTickEndF = xub
  res@tmXBTickSpacingF = 12
  res@tmXBTickStartF = 0
  res@tmXBLabels = (/2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011/)
  res@tmXBValues = (/xlb, xlb+xint, xlb+2*xint,xlb+3*xint,xlb+4*xint,xlb+5*xint,xlb+6*xint,xlb+7*xint,xlb+8*xint,xlb+9*xint,xlb+10*xint/)
  res@xyExplicitLegendLabels=(/"Mean","Median","Obs67"/)
  res@xyDashPatterns     = (/0, 0, 0/)
  res@xyLineThicknesses  = (/2.0, 2.0, 2.0/)
  res@xyLineColors       = (/"green","blue","red"/)

  ;plotting properties 2
  xint1 = 1
  xlb1  = 0
  xub1  = ny-1
  res1 = True
  res1@vpWidthF = 0.65                           ; change aspect ratio of plot
  res1@vpHeightF = 0.2                           ;
  res1@tiXAxisString = "Year"                    ; x-axis title
  res1@trXMinF  = xlb1
  res1@trXMaxF  = xub1
  res1@tmXBMinorPerMajor = 0
  res1@tmXBMode = "Explicit"
  res1@tmXBTickEndF = xub1
  res1@tmXBTickSpacingF = 1
  res1@tmXBTickStartF = 0
  res1@tmXBLabels = (/2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011/)
  res1@tmXBValues = (/xlb1, xlb1+xint1, xlb1+2*xint1,xlb1+3*xint1,xlb1+4*xint1,xlb1+5*xint1,xlb1+6*xint1,xlb1+7*xint1,xlb1+8*xint1,xlb1+9*xint1,xlb1+10*xint1/)
  res1@xyExplicitLegendLabels=(/"Mean","Median","estimate"/)
  res1@xyDashPatterns     = (/0, 0, 0/)
  res1@xyLineThicknesses  = (/2.0, 2.0, 2.0/)
  res1@xyLineColors       = (/"green","blue","red"/)
   
; define variables for the plots
  ;plotting ranges
  plot_time    = ispan(0,nm-1,1)
  ;index for observations at km67
  iobs67_miller = nsim+3

; Gross Primary Production
  ;add CI from Miller et al
  res@trYMinF  = 0
  res@trYMaxF  = 12
  xp  = new((/2*nm/), float)
  yp  = new((/2*nm/), float)
  xs  = new((/2*nm/), float)
  ys  = new((/2*nm/), float)
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_GPP(iobs67_miller,k)
     yp(k)                 = GPP(iobs67_miller,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*nm-1-k)          = plot_time(k)
     yp(2*nm-1-k)          = GPP(iobs67_miller,k) - dx
  end do
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                 = plot_time(k)
     xs(2*nm-1-k)          = plot_time(k)
     dsu                   = GPPmax(k)-GPPmedian(k)
     dsl                   = GPPmin(k)-GPPmedian(k)
     ys(k)                 = GPPmedian(k) + dsu
     ys(2*nm-1-k)          = GPPmedian(k) + dsl
  end do
  plot(0) = gsn_csm_xy(wks, plot_time(xlb:xub), (/GPPmean(xlb:xub), GPPmedian(xlb:xub), GPP(iobs67_miller,xlb:xub)/), res  )   ; create plot
  setvalues plot(0)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(a) Gross Primary Production (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.1                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
  end setvalues
  gsres1                   = True                        ; poly res
  gsres1@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres1@gsFillColor       = "pink"                      ; color chosen
  dummy1 = gsn_add_polygon (wks,plot(0),xp,yp,gsres1)
  gsres2                   = True                        ; poly res
  gsres2@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres2@gsFillColor       = "skyblue"                      ; color chosen
  dummy2 = gsn_add_polygon (wks,plot(0),xs,ys,gsres2)

; Ecosystem respiration
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_ER(iobs67_miller,k)
     yp(k)                 = ER(iobs67_miller,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*nm-1-k) = plot_time(k)
     yp(2*nm-1-k) = ER(iobs67_miller,k) - dx
  end do
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                 = plot_time(k)
     xs(2*nm-1-k)          = plot_time(k)
     dsu                   = ERmax(k)-ERmedian(k)
     dsl                   = ERmin(k)-ERmedian(k)
     ys(k)                 = ERmedian(k) + dsu
     ys(2*nm-1-k)          = ERmedian(k) + dsl
  end do
  plot(1) = gsn_csm_xy(wks, plot_time(xlb:xub), (/ERmean(xlb:xub), ERmedian(xlb:xub), ER(iobs67_miller,xlb:xub)/), res  )   ; create plot
  setvalues plot(1)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(b) Ecosystem respiration (intact)"
  end setvalues
  gsres3                   = True                        ; poly res
  gsres3@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres3@gsFillColor       = "pink"                      ; color chosen
  dummy3 = gsn_add_polygon (wks,plot(1),xp,yp,gsres3)
  gsres4                   = True                        ; poly res
  gsres4@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres4@gsFillColor       = "skyblue"                      ; color chosen
  dummy4 = gsn_add_polygon (wks,plot(1),xs,ys,gsres4)

; Net ecoystem exchange
  res@trYMinF  = -3
  res@trYMaxF  = 3
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_NEE(iobs67_miller,k)
     yp(k)                 = NEE(iobs67_miller,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*nm-1-k)          = plot_time(k)
     yp(2*nm-1-k)          = NEE(iobs67_miller,k) - dx
  end do
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                 = plot_time(k)
     xs(2*nm-1-k)          = plot_time(k)
     dsu                   = NEEmax(k)-NEEmedian(k)
     dsl                   = NEEmin(k)-NEEmedian(k)
     ys(k)                 = NEEmedian(k) + dsu
     ys(2*nm-1-k)          = NEEmedian(k) + dsl
  end do
  plot(2) = gsn_csm_xy(wks, plot_time(xlb:xub), (/NEEmean(xlb:xub), NEEmedian(xlb:xub), NEE(iobs67_miller,xlb:xub)/) , res  )   ; create plot
  setvalues plot(2)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(c) Net Ecosystem Exchange (intact)"
  end setvalues
  gsres5                   = True                        ; poly res
  gsres5@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres5@gsFillColor       = "pink"                      ; color chosen
  dummy5 = gsn_add_polygon (wks,plot(2),xp,yp,gsres5)
  gsres6                   = True                        ; poly res
  gsres6@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres6@gsFillColor       = "skyblue"                      ; color chosen
  dummy6 = gsn_add_polygon (wks,plot(2),xs,ys,gsres6)


; Net ecosystem production
  ;plotting ranges
  plot_time1    = ispan(0,ny-1,1)
  xs1  = new((/2*ny/), float)
  ys1  = new((/2*ny/), float)

  res1@trYMinF  = 1
  res1@trYMaxF  = 30
  NPPannobs = 9.5+ ispan(0,ny-1,1)*0.0000001
  ;add range for the simulations
  do k=xlb1,xub1
     xs1(k)                  = plot_time1(k)
     xs1(2*ny-1-k)           = plot_time1(k)
     dsu                     = NPPannmax(k)-NPPannmedian(k)
     dsl                     = NPPannmin(k)-NPPannmedian(k)
     ys1(k)                  = NPPannmedian(k) + dsu
     ys1(2*ny-1-k)           = NPPannmedian(k) + dsl
  end do
  plot(3) = gsn_csm_xy(wks, plot_time1(xlb1:xub1), (/NPPannmean(xlb1:xub1),NPPannmedian(xlb1:xub1), NPPannobs(xlb1:xub1)/), res1  )   ; create plot
  setvalues plot(3) 
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(d) Net Primary Production (intact)"
  end setvalues
  gsres7                   = True                        ; poly res
  gsres7@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres7@gsFillColor       = "skyblue"                      ; color chosen
  dummy7 = gsn_add_polygon (wks,plot(3),xs1,ys1,gsres7)

; heteotropic respiration
  res1@trYMinF  = 0
  res1@trYMaxF  = 25
  HRannobs = 8.9+ ispan(0,ny-1,1)*0.0000001
  ;add range for the simulations
  do k=xlb1,xub1
     xs1(k)                   = plot_time1(k)
     xs1(2*ny-1-k)            = plot_time1(k)
     dsu                      = HRannmax(k)-HRannmedian(k)
     dsl                      = HRannmin(k)-HRannmedian(k)
     ys1(k)                   = HRannmedian(k) + dsu
     ys1(2*ny-1-k)            = HRannmedian(k) + dsl
  end do
  plot(4) = gsn_csm_xy(wks, plot_time1(xlb1:xub1), (/HRannmean(xlb1:xub1),HRannmedian(xlb1:xub1), HRannobs(xlb1:xub1)/), res1  )   ; create plot
  setvalues plot(4)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(e) heteotropic respiration (intact)"
  end setvalues
  gsres8                   = True                        ; poly res
  gsres8@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres8@gsFillColor       = "skyblue"                      ; color chosen
  dummy8 = gsn_add_polygon (wks,plot(4),xs1,ys1,gsres8)

; Autotropic respiration
  res1@trYMinF  = 0
  res1@trYMaxF  = 30
  ARannobs = 23.1+ ispan(0,ny-1,1)*0.0000001
  ;add range for the simulations
  do k=xlb1,xub1
     xs1(k)                   = plot_time1(k)
     xs1(2*ny-1-k)            = plot_time1(k)
     dsu                      = ARannmax(k)-ARannmedian(k)
     dsl                      = ARannmin(k)-ARannmedian(k)
     ys1(k)                   = ARannmedian(k) + dsu
     ys1(2*ny-1-k)            = ARannmedian(k) + dsl
  end do
  plot(5) = gsn_csm_xy(wks, plot_time1(xlb1:xub1), (/ARannmean(xlb1:xub1),ARannmedian(xlb1:xub1), ARannobs(xlb1:xub1)/), res1  )   ; create plot
  setvalues plot(5)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(f) autotropic respiration (intact)"
  end setvalues
  gsres9                   = True                        ; poly res
  gsres9@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres9@gsFillColor       = "skyblue"                      ; color chosen
  dummy9 = gsn_add_polygon (wks,plot(5),xs1,ys1,gsres9)

; Above ground biomass
  res@trYMinF  = 0
  res@trYMaxF  = 30
  AGBobs = 16.5 + ispan(0,nm-1,1)*0.00000001
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = AGBmax(k)-AGBmedian(k)
     dsl                     = AGBmin(k)-AGBmedian(k)
     ys(k)                   = AGBmedian(k) + dsu
     ys(2*nm-1-k)            = AGBmedian(k) + dsl
  end do
  plot1(0) = gsn_csm_xy(wks1, plot_time(xlb:xub), (/AGBmean(xlb:xub),AGBmedian(xlb:xub),AGBobs(xlb:xub)/), res  )   ; create plot
  setvalues plot1(0)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(a) Above Ground Biomass (intact)"
  end setvalues
  gsres10                  = True                        ; poly res
  gsres10@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres10@gsFillColor       = "skyblue"                      ; color chosen
  dummy10 = gsn_add_polygon (wks1,plot1(0),xs,ys,gsres10)

; Coarse woody debris
  res@trYMinF  = 0
  res@trYMaxF  = 20
  CWDobs = 5.84 + ispan(0,nm-1,1)*0.00000001
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = CWDmax(k)-CWDmedian(k)
     dsl                     = CWDmin(k)-CWDmedian(k)
     ys(k)                   = CWDmedian(k) + dsu
     ys(2*nm-1-k)            = CWDmedian(k) + dsl
  end do
  plot1(1) = gsn_csm_xy(wks1, plot_time(xlb:xub), (/CWDmean(xlb:xub),CWDmedian(xlb:xub),CWDobs(xlb:xub)/), res  )   ; create plot
  setvalues plot1(1)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(b) Coarse Woody Debris (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.1                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
  end setvalues
  gsres11                  = True                        ; poly res
  gsres11@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres11@gsFillColor       = "skyblue"                      ; color chosen
  dummy11 = gsn_add_polygon (wks1,plot1(1),xs,ys,gsres11)

; Above ground biomass, early
  AGB_earlyobs = 7.6 + ispan(0,nm-1,1)*0.00000001
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = AGB_earlymax(k)-AGB_earlymedian(k)
     dsl                     = AGB_earlymin(k)-AGB_earlymedian(k)
     ys(k)                   = AGB_earlymedian(k) + dsu
     ys(2*nm-1-k)            = AGB_earlymedian(k) + dsl
  end do
  plot1(2) = gsn_csm_xy(wks1, plot_time(xlb:xub), (/AGB_earlymean(xlb:xub),AGB_earlymedian(xlb:xub),AGB_earlyobs(xlb:xub)/), res  )   ; create plot
  setvalues plot1(2)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(c) Above Ground Biomass, Early (intact)"
  end setvalues
  gsres12                  = True                        ; poly res
  gsres12@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres12@gsFillColor       = "skyblue"                      ; color chosen
  dummy12 = gsn_add_polygon (wks1,plot1(2),xs,ys,gsres12)

; Above ground biomass, late
  AGB_earlyobs = 8.9 + ispan(0,nm-1,1)*0.00000001
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = AGB_earlymax(k)-AGB_earlymedian(k)
     dsl                     = AGB_earlymin(k)-AGB_earlymedian(k)
     ys(k)                   = AGB_earlymedian(k) + dsu
     ys(2*nm-1-k)            = AGB_earlymedian(k) + dsl
  end do
  plot1(3) = gsn_csm_xy(wks1, plot_time(xlb:xub), (/AGB_earlymean(xlb:xub),AGB_earlymedian(xlb:xub),AGB_earlyobs(xlb:xub)/), res  )   ; create plot
  setvalues plot1(3)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(d) Above Ground Biomass, Late (intact)"
  end setvalues
  gsres13                   = True                        ; poly res
  gsres13@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres13@gsFillColor       = "skyblue"                      ; color chosen
  dummy13 = gsn_add_polygon (wks1,plot1(3),xs,ys,gsres13)



;*************************************************
; create panel
;*************************************************
  resP = True   ; modify the panel plot
  resP@gsnMaximize = True           ; use full page
  gsn_panel(wks, plot, (/4, 2/), resP)
  gsn_panel(wks1, plot1, (/3, 2/), resP)

end
