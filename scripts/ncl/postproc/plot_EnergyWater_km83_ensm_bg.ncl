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
 fnm_obs      = "./inputs/km83_obs_Macros_filter_v3"
 fnm_obs2     = "./inputs/km83_Miller_etal2011.csv"
 fnm_obs67    = "./inputs/km67_obs_Macos"
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
  FillValue = -999.
  ; energy variablesa
  LH   = new((/nv, nm/), float, FillValue); latent heat flux, w/m^2
  SH   = new((/nv, nm/), float, FillValue); sensible heat flux, w/m^2
  GH   = new((/nv, nm/), float, FillValue); ground(soil) heat flux, w/m^2
  Rn   = new((/nv, nm/), float, FillValue); net radiation, w/m^2 
  LAI  = new((/nv, nm/), float, FillValue); leaf area index
  SOILWATER_10CM = new((/nv, nm/), float, FillValue); volumetric soil moisture, first 10cm (kg/m2)
  BTRAN = new((/nv, nm/), float, FillValue); transpiration beta factor
  TWS = new((/nv, nm/), float, FillValue); total water storage, mm
  TSOI_10CM = new((/nv, nm/), float, FillValue); soil temperature in top 10cm of soil (K)
  fcev = new((/nv, nm/), float, FillValue); canopy evaporation, w/m^2
  fctr = new((/nv, nm/), float, FillValue); canopy transpiration, w/m^2
  fgev = new((/nv, nm/), float, FillValue); ground evaporation, w/m^2
  fsds = new((/nv, nm/), float, FillValue); downward shortwave radiation, w/m^2
  flds = new((/nv, nm/), float, FillValue); downward longwave radiation, w/m^2
  fire = new((/nv, nm/), float, FillValue); emitted longwave radiation, w/m^2
  fsr = new((/nv, nm/), float, FillValue); reflected shortwave radiaton, w/m^2
  fsa = new((/nv, nm/), float, FillValue); absorbed shortwave radiation, w/m^2
  RAIN  = new((/nv, nm/), float, FillValue); atmospheric rain,  mm/day
  QOVER = new((/nv, nm/), float, FillValue); surface runoff, mm/day
  QINFL = new((/nv, nm/), float, FillValue); infiltration, mm/day
  QDRIP = new((/nv, nm/), float, FillValue); throughfall, mm/day
  QDRAI = new((/nv, nm/), float, FillValue); sub-surface drainage, mm/day
  QRUNOFF = new((/nv, nm/), float, FillValue); total runoff, mm/day
  QCHARGE = new((/nv, nm/), float, FillValue); aquifer recharge rate, mm/day
  QINTR  = new((/nv, nm/), float, FillValue); interception, mm/day
  QVEGE = new((/nv, nm/), float, FillValue);  canopy evaporation, mm/day
  QVEGT = new((/nv, nm/), float, FillValue);  canopy transpiration, mm/day
  QSOIL = new((/nv, nm/), float, FillValue);  Ground evaporation, mm/day
  ET = new((/nv, nm/), float, FillValue);  Evapotranspiration, mm/day

  YEAR = new(ny, integer)
  ;confidence interval
  CI_Rn  = new((/nv, nm/), float, FillValue); confidence interval of net radiation, w/m^2
  CI_SH  = new((/nv, nm/), float, FillValue); confidence interval of sensible heat flux, w/m^2
  CI_LH  = new((/nv, nm/), float, FillValue); confidence interval of latent heat flux, w/m^2

;retrieve obs
 factor = 1; mm/day -> mm/day
 delim = ","
 ;need to get rid of obs from previous year starting from line 253
 nls=253
 lines_obs = asciiread(fnm_obs+".csv",-1,"string")
 printVarSummary(lines_obs)
 n_obs = 36
 SH(5,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),4,delim))
 LH(5,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),3,delim))
 Rn(5,0:n_obs-1)  = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),8,delim))
 ET(5,0:n_obs-1)  = LH(5,0:n_obs-1)/28.94; converting from LH
 lines_obs = asciiread(fnm_obs+"_water.csv",-1,"string")
 SOILWATER_10CM(5,0:n_obs-1)  = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),3,delim))

 ; km83
 n_obs = 132
 lines_obs = asciiread(fnm_obs67+"_filter_LH_v3.csv",-1,"string")
 LH(6,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),3,delim))
 lines_obs = asciiread(fnm_obs67+"_filter_SH_v3.csv",-1,"string")
 SH(6,0:n_obs-1) = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),3,delim))
 lines_obs = asciiread(fnm_obs67+"_filter_Rn_v3.csv",-1,"string")
 Rn(6,0:n_obs-1)  = tofloat(str_get_field(lines_obs(nls:nls+n_obs-1),3,delim))
 ET(6,0:n_obs-1)  = LH(1,0:n_obs-1)/28.94; converting from LH

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
 CI_Rn(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),40,delim))
 CI_SH(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),46,delim))
 CI_LH(7,0:n_obs2-1)  = tofloat(str_get_field(lines_obs2(nls2:nls2+n_obs2-1),47,delim))

 ; starting from line 2 for Miller et al. obs for km67
 ; the 2011 CO2 flux at km67 is problematic, not showning here
 nls2=  1
 lines_obs67_2 = asciiread(fnm_obs67_2,-1,"string")
 n_obs2 = 60
 printVarSummary(lines_obs67_2)
 factor_E2LH = 2.5*18 ; mmol/m2/s -> W/m2
 LAI(nsim+3,0:n_obs2-1)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),3,delim))
 Rn(nsim+3,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),6,delim))
 SH(nsim+3,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),14,delim))
 LH(nsim+3,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),16,delim))*factor_E2LH
 LH(nsim+3,0:11) = -999.
 LH(nsim+3, 47) = -999.
 CI_Rn(nsim+3,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),26,delim))
 CI_SH(nsim+3,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),31,delim))
 CI_LH(nsim+3,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),32,delim))*factor_E2LH
 CI_LH(nsim+3,0:11) = -999.
 CI_LH(nsim+3, 47) = -999.

;retrieve simulations
varlist = (/"Rn", "SH", "LH", "GH", "SOILWATER_10CM", "TSOI_10CM", "LAI","fcev","fctr","fgev","QOVER","QINFL","QDRIP","QDRAI","QCHARGE","QINTR","QVEGE","QVEGT","QSOIL","RAIN","TWS","BTRAN", "fsr","fire"/)
nvar= 24
ncol = nsim+2
do i=0,nvar-1
   fnm = path_o+varlist(i)+"_"+experiment+".csv"
   data = readAsciiTable(fnm, ncol, "float", 1)
   do j=0,nsim-1
      if (i.eq.0) then
         Rn(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.1) then
         SH(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.2) then
         LH(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.3) then
         GH(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.4) then
         SOILWATER_10CM(j,:) = data(:,j+2)/100.
         print(SOILWATER_10CM(j,:))
      end if
      if (i.eq.5) then
         TSOI_10CM(j,:) = data(:,j+2)-273.15  ;degree C
      end if
      if (i.eq.6) then
         LAI(j,:) = data(:,j+2)  ;m^2/m^2
      end if
      if (i.eq.7) then
         fcev(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.8) then
         fctr(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.9) then
         fgev(j,:) = data(:,j+2)  ;w/m^2
      end if
      if (i.eq.10) then
         QOVER(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.11) then
         QINFL(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.12) then
         QDRIP(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.13) then
         QDRAI(j,:) = data(:,j+2)*factor  ;mm/day
         QRUNOFF(j,:) = QOVER(j,:)+QDRAI(j,:)
      end if
      if (i.eq.14) then
         QCHARGE(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.15) then
         QINTR(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.16) then
         QVEGE(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.17) then
         QVEGT(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.18) then
         QSOIL(j,:) = data(:,j+2)*factor  ;mm/day
         ET(j,:) = QVEGE(j,:) + QVEGT(j,:) + QSOIL(j,:)
      end if
      if (i.eq.19) then
         RAIN(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.20) then
         TWS(j,:) = data(:,j+2)*factor  ;mm/day
      end if
      if (i.eq.21) then
         BTRAN(j,:) = data(:,j+2)       ;unitless
      end if
      if (i.eq.22) then
         fsr(j,:) = data(:,j+2)       ;unitless
      end if
      if (i.eq.23) then
         fire(j,:) = data(:,j+2)       ;unitless
      end if
   end do; nv
   delete(data)
end do; nvar

;calculate mean, std, median and ranges of the simulated variables
Rnmean                   = dim_avg_n(Rn(0:nsim-1,:),0)
SHmean                   = dim_avg_n(SH(0:nsim-1,:),0)
LHmean                   = dim_avg_n(LH(0:nsim-1,:),0)
LAImean                  = dim_avg_n(LAI(0:nsim-1,:),0)
SOILWATER_10CMmean       = dim_avg_n(SOILWATER_10CM(0:nsim-1,:),0)
BTRANmean                = dim_avg_n(BTRAN(0:nsim-1,:),0)
TSOI_10CMmean            = dim_avg_n(TSOI_10CM(0:nsim-1,:),0)

Rnstddev                 = dim_stddev_n(Rn(0:nsim-1,:),0)
SHstddev                 = dim_stddev_n(SH(0:nsim-1,:),0)
LHstddev                 = dim_stddev_n(LH(0:nsim-1,:),0)
LAIstddev                = dim_stddev_n(LAI(0:nsim-1,:),0)
SOILWATER_10CMstddev     = dim_stddev_n(SOILWATER_10CM(0:nsim-1,:),0)
BTRANstddev              = dim_stddev_n(BTRAN(0:nsim-1,:),0)
TSOI_10CMstddev          = dim_stddev_n(TSOI_10CM(0:nsim-1,:),0)

Rnmedian                 = dim_median_n(Rn(0:nsim-1,:),0)
SHmedian                 = dim_median_n(SH(0:nsim-1,:),0)
LHmedian                 = dim_median_n(LH(0:nsim-1,:),0)
LAImedian                = dim_median_n(LAI(0:nsim-1,:),0)
SOILWATER_10CMmedian     = dim_median_n(SOILWATER_10CM(0:nsim-1,:),0)
BTRANmedian              = dim_median_n(BTRAN(0:nsim-1,:),0)
TSOI_10CMmedian          = dim_median_n(TSOI_10CM(0:nsim-1,:),0)

Rnmax                    = dim_max_n(Rn(0:nsim-1,:),0)
SHmax                    = dim_max_n(SH(0:nsim-1,:),0)
LHmax                    = dim_max_n(LH(0:nsim-1,:),0)
LAImax                   = dim_max_n(LAI(0:nsim-1,:),0)
SOILWATER_10CMmax        = dim_max_n(SOILWATER_10CM(0:nsim-1,:),0)
BTRANmax                 = dim_max_n(BTRAN(0:nsim-1,:),0)
TSOI_10CMmax             = dim_max_n(TSOI_10CM(0:nsim-1,:),0)

Rnmin                    = dim_min_n(Rn(0:nsim-1,:),0)
SHmin                    = dim_min_n(SH(0:nsim-1,:),0)
LHmin                    = dim_min_n(LH(0:nsim-1,:),0)
LAImin                   = dim_min_n(LAI(0:nsim-1,:),0)
SOILWATER_10CMmin        = dim_min_n(SOILWATER_10CM(0:nsim-1,:),0)
BTRANmin                 = dim_min_n(BTRAN(0:nsim-1,:),0)
TSOI_10CMmin             = dim_min_n(TSOI_10CM(0:nsim-1,:),0)

;************************************************
; create plots
;************************************************
; define variables for the plots
  wks = gsn_open_wks("pdf",path_o+"WaterEnergy"+experiment)
  gsn_define_colormap(wks,"gui_default")        ; select color map
  plot = new(6, graphic)

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

; define variables for the plots
  ;plotting ranges
  plot_time    = ispan(0,nm-1,1)
  ;index for observations at km67
  iobs67_miller = nsim+3

  ;add CI from Miller et al
  res@trYMinF  = 0.
  res@trYMaxF  = 180.
  xp  = new((/2*nm/), float)
  yp  = new((/2*nm/), float)
  xs  = new((/2*nm/), float)
  ys  = new((/2*nm/), float)
; Sensible heat
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_SH(iobs67_miller,k)
     yp(k)                 = SH(iobs67_miller,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*nm-1-k)          = plot_time(k)
     yp(2*nm-1-k)          = SH(iobs67_miller,k) - dx
  end do
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                 = plot_time(k)
     xs(2*nm-1-k)          = plot_time(k)
     dsu                   = SHmax(k) - SHmedian(k)
     dsl                   = SHmin(k) - SHmedian(k)
     ys(k)                 = SHmedian(k) + dsu
     ys(2*nm-1-k)          = SHmedian(k) + dsl
  end do
  plot(0) = gsn_csm_xy(wks, plot_time(xlb:xub), (/SHmean(xlb:xub), SHmedian(xlb:xub), SH(iobs67_miller,xlb:xub)/), res  )   ; create plot
  setvalues plot(0)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(a) Sensible heat flux (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.1                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
  end setvalues
  gsres3                   = True                        ; poly res
  gsres3@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres3@gsFillColor       = "pink"                      ; color chosen
  dummy3 = gsn_add_polygon (wks,plot(0),xp,yp,gsres3)
  gsres4                   = True                        ; poly res
  gsres4@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres4@gsFillColor       = "skyblue"                      ; color chosen
  dummy4 = gsn_add_polygon (wks,plot(0),xs,ys,gsres4)

; Latent heat flux
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_LH(iobs67_miller,k)
     yp(k)                 = LH(iobs67_miller,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*nm-1-k)          = plot_time(k)
     yp(2*nm-1-k)          = LH(iobs67_miller,k) - dx
  end do
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                 = plot_time(k)
     xs(2*nm-1-k)          = plot_time(k)
     dsu                   = LHmax(k)-LHmedian(k)
     dsl                   = LHmin(k)-LHmedian(k)
     ys(k)                 = LHmedian(k) + dsu
     ys(2*nm-1-k)          = LHmedian(k) + dsl
  end do
  plot(1) = gsn_csm_xy(wks, plot_time(xlb:xub), (/LHmean(xlb:xub), LHmedian(xlb:xub), LH(iobs67_miller,xlb:xub)/) , res  )   ; create plot
  setvalues plot(1)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(b) Latent heat flux (intact)"
  end setvalues
  gsres5                   = True                        ; poly res
  gsres5@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres5@gsFillColor       = "pink"                      ; color chosen
  dummy5 = gsn_add_polygon (wks,plot(1),xp,yp,gsres5)
  gsres6                   = True                        ; poly res
  gsres6@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres6@gsFillColor       = "skyblue"                      ; color chosen
  dummy6 = gsn_add_polygon (wks,plot(1),xs,ys,gsres6)


; Leaf Area Index
  res@trYMinF  = 0.
  res@trYMaxF  = 10.
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                  = plot_time(k)
     xs(2*nm-1-k)           = plot_time(k)
     dsu                    = LAImax(k)-LAImedian(k)
     dsl                    = LAImin(k)-LAImedian(k)
     ys(k)                  = LAImedian(k) + dsu
     ys(2*nm-1-k)           = LAImedian(k) + dsl
  end do
  plot(2) = gsn_csm_xy(wks, plot_time(xlb:xub), (/LAImean(xlb:xub), LAImedian(xlb:xub), LAI(iobs67_miller,xlb:xub)/), res  )   ; create plot
  setvalues plot(2) 
    "tiYAxisString" : "m~S~2~N~m~S~-2"   ; y-axis title
    "tiMainString": "(c) Leaf area index (intact)"
  end setvalues
  gsres7                   = True                        ; poly res
  gsres7@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres7@gsFillColor       = "skyblue"                      ; color chosen
  dummy7 = gsn_add_polygon (wks,plot(2),xs,ys,gsres7)

; Soil Moisture
  res@trYMinF  = 0.2
  res@trYMaxF  = 0.5
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = SOILWATER_10CMmax(k)-SOILWATER_10CMmedian(k)
     dsl                     = SOILWATER_10CMmin(k)-SOILWATER_10CMmedian(k)
     ys(k)                   = SOILWATER_10CMmedian(k) + dsu
     ys(2*nm-1-k)            = SOILWATER_10CMmedian(k) + dsl
  end do
  plot(3) = gsn_csm_xy(wks, plot_time(xlb:xub), (/SOILWATER_10CMmean(xlb:xub),SOILWATER_10CMmedian(xlb:xub)/), res  )   ; create plot
  setvalues plot(3)
    "tiYAxisString" : "m~S~3~N~m~S~-3"   ; y-axis title
    "tiMainString": "(d) Soil water in top 10 cm (intact)"
  end setvalues
  gsres8                   = True                        ; poly res
  gsres8@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres8@gsFillColor       = "skyblue"                      ; color chosen
  dummy8 = gsn_add_polygon (wks,plot(3),xs,ys,gsres8)

; Autotropic respiration
  res@trYMinF  = 0.
  res@trYMaxF  = 1.
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = BTRANmax(k)-BTRANmedian(k)
     dsl                     = BTRANmin(k)-BTRANmedian(k)
     ys(k)                   = BTRANmedian(k) + dsu
     ys(2*nm-1-k)            = BTRANmedian(k) + dsl
  end do
  plot(4) = gsn_csm_xy(wks, plot_time(xlb:xub), (/BTRANmean(xlb:xub),BTRANmedian(xlb:xub)/), res  )   ; create plot
  setvalues plot(4)
    "tiYAxisString" : "unitless"   ; y-axis title
    "tiMainString": "(e) Transpiration beta factor (intact)"
  end setvalues
  gsres9                   = True                        ; poly res
  gsres9@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres9@gsFillColor       = "skyblue"                      ; color chosen
  dummy9 = gsn_add_polygon (wks,plot(4),xs,ys,gsres9)

; Soil temperature
  res@trYMinF  = 20 
  res@trYMaxF  = 35
  ;add range for the simulations
  do k=xlb,xub
     xs(k)                   = plot_time(k)
     xs(2*nm-1-k)            = plot_time(k)
     dsu                     = TSOI_10CMmax(k)-TSOI_10CMmedian(k)
     dsl                     = TSOI_10CMmin(k)-TSOI_10CMmedian(k)
     ys(k)                   = TSOI_10CMmedian(k) + dsu
     ys(2*nm-1-k)            = TSOI_10CMmedian(k) + dsl
  end do
  plot(5) = gsn_csm_xy(wks, plot_time(xlb:xub), (/TSOI_10CMmean(xlb:xub),TSOI_10CMmedian(xlb:xub),TSOI_10CM(iobs67_miller,xlb:xub)/), res  )   ; create plot
  setvalues plot(5)
    "tiYAxisString" : "degree C"   ; y-axis title
    "tiMainString": "(f) 10-cm soil temperature (intact)"
  end setvalues
  gsres10                  = True                        ; poly res
  gsres10@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres10@gsFillColor       = "skyblue"                      ; color chosen
  dummy10 = gsn_add_polygon (wks,plot(5),xs,ys,gsres10)

;*************************************************
; create panel
;*************************************************
  resP = True   ; modify the panel plot
  resP@gsnMaximize = True           ; use full page
  gsn_panel(wks, plot, (/3, 2/), resP)

end
