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
 fnm_obs      = "./ed2_logging_tapajos/km83_obs_Macros_filter_v3.csv"
 fnm_obs67      = "./ed2_logging_tapajos/km67_obs_Macos_filter_GPP_NEE_v3.csv"
 compset      = "I2000Clm50FatesGs" 
 experiments  = (/"intact","RILlow","RILhigh","CLlow","CLhigh"/)
 nv           = 5+2; sim+obs83+obs67
 npft         = 2;  2 PFTs
 ny           = year_end - year_start +1 ;# of nyear
 nm           = ny * 12
  
 ;****************************************************
 ; Modify variables for comparison
  ;carbon cycle variables
  ER   = new((/nv, nm/), float); ecosystem respiration, autotrophic + heterotrophic (gC/m^2/day)
  HR   = new((/nv, nm/), float); heterotropic respiration (gC/m^2/day)
  AR   = new((/nv, nm/), float); autotrophic respiration (gC/m^2/day)
  AGB  = new((/nv, nm/), float); above ground biomass gC/m2
  AGB_early = new((/nv, nm/), float); above ground biomass for early succession PFT (gC/m2)
  AGB_late  = new((/nv, nm/), float); above ground biomass for late succession PFT (gC/m2)
  ED_biomass   = new((/nv, nm/), float); Total biomass (gC m-2)
  CWD   = new((/nv, nm/), float); total CWD carbon at the column level(gC/m^2)
  TOTSOMC_1m = new((/nv, nm/), float); total soil organic matter carbon to 1 meter depth (gC/m^2)
  ED_bleaf = new((/nv, nm/), float); factor limiting ground evap 
  NEE  = new((/nv, nm/), float); net ecosystem exchange (gC/m^2/day)
  GPP  = new((/nv, nm/), float); gross primary production (gC/m^2/day)
  GPP_CANOPY = new((/nv, nm/), float); gross primary production of canopy trees (gC/m^2/day)
  GPP_UNDERSTORY = new((/nv, nm/), float); gross primary production of understory trees (gC/m^2/day)
  NPP  = new((/nv, nm/), float); net primary production (gC/m^2/day)
  ; energy variablesa
  LH   = new((/nv, nm/), float); latent heat flux, w/m^2
  SH   = new((/nv, nm/), float); sensible heat flux, w/m^2
  GH   = new((/nv, nm/), float); ground(soil) heat flux, w/m^2
  Rn   = new((/nv, nm/), float); net radiation, w/m^2 
  LAI  = new((/nv, nm/), float); leaf area index
  SOILWATER_10CM = new((/nv, nm/), float); volumetric soil moisture, first 10cm (kg/m2)
  TSOI_10CM = new((/nv, nm/), float); soil temperature in top 10cm of soil (K)
  fcev = new((/nv, nm/), float); canopy evaporation, w/m^2
  fctr = new((/nv, nm/), float); canopy transpiration, w/m^2
  fgev = new((/nv, nm/), float); ground evaporation, w/m^2
  fsds = new((/nv, nm/), float); downward shortwave radiation, w/m^2
  flds = new((/nv, nm/), float); downward longwave radiation, w/m^2
  fire = new((/nv, nm/), float); emitted longwave radiation, w/m^2
  fsr = new((/nv, nm/), float); reflected shortwave radiaton, w/m^2
  fsa = new((/nv, nm/), float); absorbed shortwave radiation, w/m^2
  GPPann = new((/nv, ny/), float); annual GPP MgC/m^2/yr-1
  GPPann_can = new((/nv, ny/), float); annual GPP in the canopy layer MgC/m^2/yr-1
  GPPann_und = new((/nv, ny/), float); annual GPP in the understory MgC/m^2/yr-1
  ERann = new((/nv, ny/), float);  annual ER MgC/m^2/yr-1
  NEEann = new((/nv, ny/), float); annual NEE MgC/m^2/yr-1  
  NPPann = new((/nv, ny/), float); annual NPP MgC/m^2/yr-1
  HRann = new((/nv, ny/), float);  annual HR MgC/m^2/yr-1
  ARann = new((/nv, ny/), float);  annual AR MgC/m^2/yr-1
  YEAR = new(ny, integer)

;retrieve obs
 factor = 1000.*60.*60.*24.; conversion factor from KgC/m^2/s to gC/m^2/day
 delim = ","
 ;need to get rid of obs from previous year starting from line 253
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
	
;retrieve simulations
varlist = (/"ER", "HR", "AR", "AGB", "ED_biomass", "CWD", "TOTSOMC_1m","NEE","GPP","LAI","ED_bleaf","NPP","AGB_early","AGB_late", "GPP_CANOPY", "GPP_UNDERSTORY"/)
nvar= 16
ncol = nv
do i=0,nvar-1
   fnm = path_o+varlist(i)+"_"+site+".csv"
   print(fnm)
   ;lines = asciiread(fnm,-1,"string")
   ;printVarSummary(lines)
    data = readAsciiTable(fnm, ncol, "float", 1)
   do j=0,nv-3
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
         AGB_early(j,:) = data(:,j+2)/1000.  ;KgC/m2
      end if
      if (i.eq.13) then
         AGB_late(j,:) = data(:,j+2)/1000.  ;KgC/m2
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
  wks2 = gsn_open_wks("pdf",path_o+"CarbonPools"+site)
  gsn_define_colormap(wks2,"gui_default")        ; select color map
; define variables for the plots
  ;plot1 = new(6, graphic) 
  plot2 = new(8, graphic)

  xlb=0
  xub=nm-1
  xint=10*12
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
  res@tmXBValues = (/xlb, xlb+1*xint, xlb+2*xint,xlb+3*xint,xlb+4*xint,xlb+5*xint,xlb+6*xint,xlb+7*xint,xlb+8*xint,xlb+9*xint,xlb+nm/)
  res@tmXBLabels = (/2001,2011,2021,2031,2041,2051,2061,2071,2081,2091,2101/)
  ;res@xyExplicitLegendLabels=(/"Intact","RIL~B~low","RIL~B~high","CL~B~low","CL~B~high","Obs"/)
  res@xyExplicitLegendLabels=(/"Intact","Obs67"/)
  res@xyMarkLineModes =  (/"Lines","Markers"/)   ; Three cases of sensitivity analysis
  res@xyDashPatterns     = (/0/)
  res@xyLineThicknesses = (/2.0/)
  ;res@gsnXRefLine = 8
  ;res@gsnXRefLineColor = (/"black"/)
  res@xyLineColors      = (/"blue"/)
  res@xyMarker = 3
  res@xyMarkerColor =(/"red"/)
  res@xyMarkerSizeF = 0.005
  res@xyDashPattern     = 0
  res@tiMainFontHeightF = 0.023
  ;res@gsnYRefLineDashPattern = 1
  ;plotting ranges
  plot_time    = ispan(0,nm-1,1)

;plot2
  res1 = res
  delete(res1@xyExplicitLegendLabels)
  delete(res1@xyMarkLineModes)
  delete(res1@xyLineThicknesses)
  delete(res1@xyLineColors)
  delete(res1@xyDashPatterns)
  res1@xyExplicitLegendLabels=(/"RIL~B~low","RIL~B~high","CL~B~low","CL~B~high","Obs83"/)
  res1@xyMarkLineModes =  (/"Lines","Lines","Lines","Lines","Markers"/)
  res1@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0,1.0/)
  res1@xyLineColors = (/"royalblue3","skyblue","green","orange","red"/)
  res1@xyDashPatterns     = (/0, 0, 0, 0, 0/)
  res1@xyMarker = 4
  res1@gsnXRefLine = 8
  res1@gsnXRefLineColor = (/"black"/)
  res1@gsnXRefLineDashPattern     = 1

  ipleft=(/0/)
  ipright=(/1,2,3,4/)
  res@trYMinF  = 0
  res@trYMaxF  = 20
  res@gsnYRefLine = 16.5
  res@gsnYRefLineColor = (/"red"/)
  print(AGB(ipleft, 1:12))
  plot2(0) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(0)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(a) Above Ground Biomass (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Bottom"                ; Change location o
    "pmLegendWidthF"     :    0.15                 ; Change width and
    "pmLegendHeightF"    :    0.08                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
    "pmLegendOrthogonalPosF": -1.0                ; more neg = down
  end setvalues

  res1@trYMinF  = 0
  res1@trYMaxF  = 20
  res1@gsnYRefLine = 14.7
  res1@gsnYRefLineColor = (/"red"/)
  plot2(1) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(1)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(b) Above Ground Biomass (logged)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Bottom"                ; Change location o
    "pmLegendWidthF"     :    0.15                 ; Change width and
    "pmLegendHeightF"    :    0.08                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
    "pmLegendOrthogonalPosF": -1.0                ; more neg = down
  end setvalues
  res@trYMinF  = 0
  res@trYMaxF  = 12
  res@gsnYRefLine = 5.84
  res@gsnYRefLineColor = (/"red"/)
  plot2(2) = gsn_csm_xy(wks2, plot_time(xlb:xub), CWD(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(2)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(c) Coarse Woody Debris (intact)"
  end setvalues
  res1@trYMinF  = 0
  res1@trYMaxF  = 12
  res1@gsnYRefLine = 7.44
  res1@gsnYRefLineColor = (/"red"/)
  plot2(3) = gsn_csm_xy(wks2, plot_time(xlb:xub), CWD(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(3)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(d) Coarse Woody Debris (logged)"
  end setvalues
  ;delete(res@gsnYRefLine)
  ;delete(res1@gsnYRefLine) 
  res@trYMinF  = 0
  res@trYMaxF  = 10
  res@gsnYRefLine = 7.6
  plot2(4) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB_early(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(4)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(e) Above Ground Biomass, Early (intact)"
  end setvalues
  res1@trYMinF  = 0
  res1@trYMaxF  = 10 
  res1@gsnYRefLine = 6.8
  plot2(5) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB_early(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(5)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(f) Above Ground Biomass, Early (logged)"
  end setvalues
  res@trYMinF  = 0
  res@trYMaxF  = 18
  res@gsnYRefLine = 8.9
  plot2(6) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB_late(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(6)
    "tiYAxisString" : "KgC m~S~-2"   ; y-axis title
    "tiMainString": "(g) Above Ground Biomass, Late (intact)"
  end setvalues
  res1@trYMinF  = 0
  res1@trYMaxF  = 18
  res1@gsnYRefLine = 7.9
  plot2(7) = gsn_csm_xy(wks2, plot_time(xlb:xub), AGB_late(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(7)
    "tiYAxisString" : "gC m~S~-2~N~d~S~-1"   ; y-axis title
    "tiMainString": "(h) Above Ground Biomass, Late (logged)"
  end setvalues

; create panel
;*************************************************
  resP = True   ; modify the panel plot
  ;resP@txString = "State variables" ; add common title
  ;resP@gsnMaximize = True           ; use full page

  gsn_panel(wks2, plot2, (/4, 2/), resP)

end
