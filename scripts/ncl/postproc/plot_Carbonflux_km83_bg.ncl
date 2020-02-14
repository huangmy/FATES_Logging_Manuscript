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
 year_start   = 2001
 year_end     = 2099
 path_o       = "./outputs/"
 fnm_obs      = "./inputs/km83_obs_Macros_filter_v3.csv"
 fnm_obs2     = "./inputs/km83_Miller_etal2011.csv"
 fnm_obs67    = "./inputs/km67_obs_Macos_filter_GPP_NEE_v3.csv"
 fnm_obs67_2  = "./inputs/km67_Miller_etal2011.csv"
 compset      = "I1PTCLM45"
 experiments  = (/"intact","RILlow","RILhigh","CLlow","CLhigh"/)
 nv           = 5+4; sim+obs83+obs67+obs83_2+obs67_2
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
 GPP(5,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),6,delim))*factor
 NEE(5,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),5,delim))*factor
 ER(5,0:n_obs-1)  = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),7,delim))*factor

 lines_obs67 = asciiread(fnm_obs67,-1,"string")
 printVarSummary(lines_obs67)
 n_obs = 132
 GPP(6,0:n_obs-1) = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),4,delim))*factor
 NEE(6,0:n_obs-1) = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),3,delim))*factor
 ER(6,0:n_obs-1)  = tofloat(str_get_field(lines_obs67(nls:nls+n_obs-1),5,delim))*factor

 ; starting from line 2 for Miller et al. obs for km83
 nls2=1 
 lines_obs2 = asciiread(fnm_obs2,-1,"string")
 n_obs2 = 47
 printVarSummary(lines_obs2)
 LAI(7,0:n_obs2-1)       = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),3,delim))
 TSOI_10CM(7,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),20,delim))
 SOILWATER_10CM(7,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),23,delim))
 Rn(7,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),10,delim))
 GH(7,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),26,delim)) 
 SH(7,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),28,delim))
 LH(7,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),29,delim))
 NEE(7,0:n_obs2-1)    = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),33,delim))
 ER(7,0:n_obs2-1)     = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),34,delim)) 
 GPP(7,0:n_obs2-1)    = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),37,delim))
 CI_Rn(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),40,delim))
 CI_SH(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),46,delim))
 CI_LH(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),47,delim))
 CI_NEE(7,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),50,delim))
 CI_ER(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),51,delim))
 CI_GPP(7,0:n_obs2-1) = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),52,delim))

 ; starting from line 2 for Miller et al. obs for km67
 ; the 2011 CO2 flux at km67 is problematic, not showning here
 nls2=12
 lines_obs67_2 = asciiread(fnm_obs67_2,-1,"string")
 n_obs2 = 49
 printVarSummary(lines_obs67_2)
; factor_E2LH = 2.5*18 ; mmol/m2/s -> W/m2
; factor_Fc   = 12*1e-6 *24*3600  ; umol/m2/s -> gC/m2/day
 Rn(8,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),6,delim))
 SH(8,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),14,delim))
 LH(8,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),16,delim))
 NEE(8,nls2-1:nls2+n_obs2-2)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),18,delim))
 ER(8,nls2-1:nls2+n_obs2-2)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),20,delim))
 GPP(8,nls2-1:nls2+n_obs2-2)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),23,delim))
 CI_Rn(8,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),26,delim))
 CI_SH(8,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),31,delim))
 CI_LH(8,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),32,delim))
 CI_NEE(8,nls2-1:nls2+n_obs2-2) = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),33,delim))
 CI_ER(8,nls2-1:nls2+n_obs2-2)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),36,delim))
 CI_GPP(8,nls2-1:nls2+n_obs2-2) = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),38,delim))

;retrieve simulations
varlist = (/"ER", "HR", "AR", "AGB", "ED_biomass", "CWD", "TOTSOMC_1m","NEE","GPP","LAI","ED_bleaf","NPP","AGB_early","AGB_late", "GPP_CANOPY", "GPP_UNDERSTORY"/)
nvar= 16
ncol = nv-2
do i=0,nvar-1
   fnm = path_o+varlist(i)+"_"+site+".csv"
   print(fnm)
   ;lines = asciiread(fnm,-1,"string")
   ;printVarSummary(lines)
    data = readAsciiTable(fnm, ncol, "float", 1)
   do j=0,nv-5
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
         AGB_early(j,:) = data(:,j+2)  ;gC/m2
      end if
      if (i.eq.13) then
         AGB_late(j,:) = data(:,j+2)  ;gC/m2
      end if
      if (i.eq.14) then
         GPP_CANOPY(j,:) = data(:,j+2)  ;gC/m2
      end if
      if (i.eq.15) then
         GPP_UNDERSTORY(j,:) = data(:,j+2)  ;gC/m2
      end if
   end do; nv 
end do; nvar

;calculate statistics
GPPmean   = dim_avg(GPP(:,0:n_obs-1))
NEEmean   = dim_avg(NEE(:,0:n_obs-1))
ERmean   = dim_avg(ER(:,0:n_obs-1))
AGBmean  = dim_avg(AGB(:,0:n_obs-1))
CWDmean  = dim_avg(CWD(:,0:n_obs-1))
GPPstddev = dim_stddev(GPP(:,0:n_obs-1))
NEEstddev = dim_stddev(NEE(:,0:n_obs-1))
ERstddev = dim_stddev(ER(:,0:n_obs-1))
AGBstddev = dim_stddev(AGB(:,0:n_obs-1))
CWDstddev = dim_stddev(CWD(:,0:n_obs-1))
GPPrmsd   = dim_rmsd(GPP(1,0:n_obs-1), GPP(5,0:n_obs-1))
NEErmsd   = dim_rmsd(NEE(1,0:n_obs-1), NEE(5,0:n_obs-1))
ERrmsd   = dim_rmsd(ER(1,0:n_obs-1), ER(5,0:n_obs-1))
;GPPrl     = avg(abs((GPP(1,0:n_obs-1)-GPP(5,0:n_obs-1))/GPP(5,0:n_obs-1)))
;NEErl     = avg(abs((NEE(1,0:n_obs-1)-NEE(5,0:n_obs-1))/NEE(5,0:n_obs-1)))
;ERrl     = avg(abs((ER(1,0:n_obs-1)-ER(5,0:n_obs-1))/ER(5,0:n_obs-1)))
;print((/GPPmean(1), GPPstddev(1), GPPmean(5), GPPstddev(5), GPPrmsd, GPPrl/))
;print((/NEEmean(1), NEEstddev(1), NEEmean(5), NEEstddev(5), NEErmsd, NEErl/))
;print((/ERmean(1), ERstddev(1), ERmean(5), ERstddev(5), ERrmsd, ERrl/))
;print((/AGBmean(1), AGBstddev(1), AGBmean(5), AGBstddev(5)/))
;print((/CWDmean(1), CWDstddev(1), CWDmean(5), CWDstddev(5)/))

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
lines = new(ny+1, string)
lines(0)=(/"  yr,    intact,    RILlow,   RILhigh,     CLlow,    CLhigh"/); header
YEAR = ispan(year_start,year_end,1)
fnmo = path_o+"GPPann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",GPPann(0,:)) + "," + sprintf("%10.2f",GPPann(1,:))+ "," + sprintf("%10.2f",GPPann(2,:)) + "," + sprintf("%10.2f",GPPann(3,:)) + "," + sprintf("%10.2f",GPPann(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"GPPann_can"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",GPPann_can(0,:)) + "," + sprintf("%10.2f",GPPann_can(1,:))+ "," + sprintf("%10.2f",GPPann_can(2,:)) + "," + sprintf("%10.2f",GPPann_can(3,:)) + "," + sprintf("%10.2f",GPPann_can(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"GPPann_und"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",GPPann_und(0,:)) + "," + sprintf("%10.2f",GPPann_und(1,:))+ "," + sprintf("%10.2f",GPPann_und(2,:)) + "," + sprintf("%10.2f",GPPann_und(3,:)) + "," + sprintf("%10.2f",GPPann_und(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"NEEann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",NEEann(0,:)) + "," + sprintf("%10.2f",NEEann(1,:))+ "," + sprintf("%10.2f",NEEann(2,:)) + "," + sprintf("%10.2f",NEEann(3,:)) + "," + sprintf("%10.2f",NEEann(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"ERann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",ERann(0,:)) + "," + sprintf("%10.2f",ERann(1,:))+ "," + sprintf("%10.2f",ERann(2,:)) + "," + sprintf("%10.2f",ERann(3,:)) + "," + sprintf("%10.2f",ERann(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"NPPann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",NPPann(0,:)) + "," + sprintf("%10.2f",NPPann(1,:))+ "," + sprintf("%10.2f",NPPann(2,:)) + "," + sprintf("%10.2f",NPPann(3,:)) + "," + sprintf("%10.2f",NPPann(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"HRann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",HRann(0,:)) + "," + sprintf("%10.2f",HRann(1,:))+ "," + sprintf("%10.2f",HRann(2,:)) + "," + sprintf("%10.2f",HRann(3,:)) + "," + sprintf("%10.2f",HRann(4,:))
asciiwrite(fnmo,lines)
fnmo = path_o+"ARann"+site+".csv"
system("rm -rf " + fnmo)
lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",ARann(0,:)) + "," + sprintf("%10.2f",ARann(1,:))+ "," + sprintf("%10.2f",ARann(2,:)) + "," + sprintf("%10.2f",ARann(3,:)) + "," + sprintf("%10.2f",ARann(4,:))
asciiwrite(fnmo,lines)

;************************************************
; create plots
;************************************************
  ;wks1 = gsn_open_wks("pdf",path_o+"CarbonPools"+site)
  ;gsn_define_colormap(wks1,"gui_default")        ; select color map
  wks2 = gsn_open_wks("pdf",path_o+"CarbonFluxes"+site)
  gsn_define_colormap(wks2,"gui_default")        ; select color map
; define variables for the plots
  ;plot1 = new(6, graphic) 
  plot2 = new(12, graphic)

  xlb=0
  xub=10*12-1
  xint=12
  plot_time    = ispan(0,nm-1,1)
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
  res@xyExplicitLegendLabels=(/"Intact","Obs67"/)
  res@xyMarkLineModes =  (/"Lines","Lines"/)   ; Three cases of sensitivity analysis
  res@xyDashPattern     = (/0, 0/)
  res@xyLineThicknesses = (/2.0, 1.0/)
  ;res@gsnXRefLine = 8
  ;res@gsnXRefLineColor = (/"black"/)
  res@xyLineColors      = (/"blue","red"/)
  res@xyMarker = 3
  res@xyMarkerColor =(/"red"/)
  res@xyMarkerSizeF = 0.005

  res1 = res
  delete(res1@xyExplicitLegendLabels)
  delete(res1@xyMarkLineModes)
  delete(res1@xyLineThicknesses)
  delete(res1@xyLineColors)
  delete(res1@xyDashPattern)
  res1@xyExplicitLegendLabels=(/"RIL~B~low","RIL~B~high","CL~B~low","CL~B~high","Obs83"/)
  res1@xyMarkLineModes   =  (/"Lines","Lines","Lines","Lines","Lines"/)
  res1@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0, 1.0/)
  res1@xyLineColors      = (/"royalblue3","skyblue","green","orange","red"/)
  res1@xyDashPatterns     = (/0, 0, 0, 0, 0/)
  res1@gsnXRefLine = 8.
  res1@gsnXRefLineColor = (/"black"/)
  res1@gsnXRefLineDashPattern = 1

  res2 = res
  xlb1 = 0
  xub1 = 10-1
  xint1 = 1
  delete(res2@tmXBValues)
  delete(res2@trXMinF)
  delete(res2@trXMaxF)
  delete(res2@tmXBTickEndF)
  delete(res2@tmXBTickSpacingF)
  res2@trXMinF = xlb1
  res2@trXMaxF = xub1
  res2@tmXBTickEndF = xub1
  res2@tmXBTickSpacingF = 1
  res2@tmXBValues = (/xlb, xlb+xint1, xlb+2*xint1,xlb+3*xint1,xlb+4*xint1,xlb+5*xint1,xlb+6*xint1,xlb+7*xint1,xlb+8*xint1,xlb+9*xint1,xlb+10*xint1/)

  res3 = res1
  delete(res3@tmXBValues)
  delete(res3@trXMinF)
  delete(res3@trXMaxF)
  delete(res3@tmXBTickEndF)
  delete(res3@tmXBTickSpacingF)
  res3@trXMinF = xlb1
  res3@trXMaxF = xub1
  res3@tmXBTickEndF = xub1
  res3@tmXBTickSpacingF = 1
  res3@tmXBValues = (/xlb, xlb+xint1, xlb+2*xint1,xlb+3*xint1,xlb+4*xint1,xlb+5*xint1,xlb+6*xint1,xlb+7*xint1,xlb+8*xint1,xlb+9*xint1,xlb+10*xint1/)
  res3@gsnXRefLine = 0.75

; plotting simulations
  ipleft=(/0, 8/)
  ipright=(/1,2,3,4,7/)


;************************************************
;  define a polygon centered the width of 2 sigma
;************************************************
  ;add CI from Miller et al
  res@trYMinF  = 0
  res@trYMaxF  = 12
  xp  = new((/2*(xub-xlb-1)/), float)
  yp  = new((/2*(xub-xlb-1)/), float)
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_GPP(8,k)
     yp(k)                 = GPP(8,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = GPP(8,k) - dx
  end do
  plot2(0) = gsn_csm_xy(wks2, plot_time(xlb:xub), GPP(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(0)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(a) Gross Primary Production (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.08                ; height of legend.
    "pmLegendParallelPosF":   0.75                ; move units right
    "pmLegendOrthogonalPosF": -1.3              ; more neg = down
  end setvalues
  gsres1                   = True                        ; poly res
  gsres1@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres1@gsFillColor       = "pink"                      ; color chosen
  dummy1 = gsn_add_polygon (wks2,plot2(0),xp,yp,gsres1)

  res1@trYMinF  = 0
  res1@trYMaxF  = 12
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_GPP(7,k)
     yp(k)                 = GPP(7,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = GPP(7,k) - dx
  end do  
  plot2(1) = gsn_csm_xy(wks2, plot_time(xlb:xub), GPP(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(1)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(b) Gross Primary Production (logged)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.08                 ; height of legend.
    "pmLegendParallelPosF":   0.75                ; move units right
    "pmLegendOrthogonalPosF": -1.3              ; more neg = down
  end setvalues
  gsres2                   = True                        ; poly res
  gsres2@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres2@gsFillColor       = "pink"                      ; color chosen
  dummy2 = gsn_add_polygon (wks2,plot2(1),xp,yp,gsres2)
  ;draw(plot2(1))                                            ; draw
  ;frame(wks2)

  do k=xlb,xub
     dx                    = CI_ER(8,k)
     yp(k)                 = ER(8,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = ER(8,k) - dx
  end do
  plot2(2) = gsn_csm_xy(wks2, plot_time(xlb:xub), ER(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(2)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(c) Ecosystem respiration (intact)"
  end setvalues
  gsres3                   = True                        ; poly res
  gsres3@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres3@gsFillColor       = "pink"                      ; color chosen
  dummy3 = gsn_add_polygon (wks2,plot2(2),xp,yp,gsres3)

  do k=xlb,xub
     dx                    = CI_ER(7,k)
     yp(k)                 = ER(7,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = ER(7,k) - dx
  end do
  plot2(3) = gsn_csm_xy(wks2, plot_time(xlb:xub), ER(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(3)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(d) Ecosystem respiration (logged)"
  end setvalues
  gsres4                   = True                        ; poly res
  gsres4@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres4@gsFillColor       = "pink"                      ; color chosen
  dummy4 = gsn_add_polygon (wks2,plot2(3),xp,yp,gsres4)
  ;draw(plot2(3))                                            ; draw
  ;frame(wks2)  

  res@trYMinF  = -3
  res@trYMaxF  = 3
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_NEE(8,k)
     yp(k)                 = NEE(8,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = NEE(8,k) - dx
  end do
  plot2(4) = gsn_csm_xy(wks2, plot_time(xlb:xub), NEE(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(4)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(e) Net Ecosystem Exchange (intact)"
  end setvalues
  gsres5                   = True                        ; poly res
  gsres5@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres5@gsFillColor       = "pink"                      ; color chosen
  dummy5 = gsn_add_polygon (wks2,plot2(4),xp,yp,gsres5)

  do k=xlb,xub
     dx                    = CI_NEE(7,k)
     yp(k)                 = NEE(7,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = NEE(7,k) - dx
  end do
  res1@trYMinF  = -3
  res1@trYMaxF  = 3
  plot2(5) = gsn_csm_xy(wks2, plot_time(xlb:xub), NEE(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(5)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(f) Net Ecosystem Exchange (logged)"
  end setvalues
  gsres6                   = True                        ; poly res
  gsres6@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres6@gsFillColor       = "pink"                      ; color chosen
  dummy6 = gsn_add_polygon (wks2,plot2(5),xp,yp,gsres6)
  ;draw(plot2(5))
  ;frame(wks2)                                            ; advance frame

; Net ecosystem production
  ;plotting ranges
  plot_time1    = ispan(0,ny-1,1)
  res2@trYMinF  = 1
  res2@trYMaxF  = 30
  res2@gsnYRefLine = 9.5
  res2@gsnYRefLineColor = (/"red"/)
  res2@gsnYRefLineDashPattern = (/1/)
  plot2(6) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), NPPann(ipleft,xlb1:xub1), res2  )   ; create plot
  setvalues plot2(6)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(g) Net Primary Production (intact)"
  end setvalues
  
  res3@trYMinF  = 0
  res3@trYMaxF  = 30
  res3@gsnYRefLine = 9.8
  res3@gsnYRefLineColor = (/"red"/)
  res3@gsnYRefLineDashPattern = (/1/)
  plot2(7) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), NPPann(ipright,xlb1:xub1), res3  )   ; create plot
  setvalues plot2(7)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(h) Net Primary Production (logged)"
  end setvalues

  res2@trYMinF  = 0
  res2@trYMaxF  = 25
  res2@gsnYRefLine = 8.9
  plot2(8) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), HRann(ipleft,xlb1:xub1), res2  )   ; create plot
  setvalues plot2(8)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(i) Heterotrophic respiration (intact)"
  end setvalues

  res3@trYMinF  = 0
  res3@trYMaxF  = 25
  res3@gsnYRefLine = 10.4
  plot2(9) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), HRann(ipright,xlb1:xub1), res3  )   ; create plot
  setvalues plot2(9)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(j) Heterotrophic respiration (logged)"
  end setvalues

  res2@trYMinF  = 0
  res2@trYMaxF  = 30
  res2@gsnYRefLine = 23.1
  plot2(10) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), ARann(ipleft,xlb1:xub1), res2  )   ; create plot
  setvalues plot2(10)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(k) Autotrophic respiration (intact)"
  end setvalues

  res3@trYMinF  = 0
  res3@trYMaxF  = 30
  res3@gsnYRefLine = 20.1
  plot2(11) = gsn_csm_xy(wks2, plot_time1(xlb1:xub1), ARann(ipright,xlb1:xub1), res3  )   ; create plot
  setvalues plot2(11)
    "tiYAxisString" : "MgC ha~S~-1~N~yr~S~-1"   ; y-axis title
    "tiMainString": "(i) Autotrophic respiration (logged)"
  end setvalues

;*************************************************
; create panel
;*************************************************
  resP = True   ; modify the panel plot
  resP@gsnMaximize = True           ; use full page
  ;gsn_panel(wks1, plot1, (/3, 2/), resP)
  gsn_panel(wks2, plot2, (/6, 2/), resP)

end
