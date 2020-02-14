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
 fnm_obs      = "./inputs/km83_obs_Macros_filter_v3"
 fnm_obs2     = "./inputs/km83_Miller_etal2011.csv"
 fnm_obs67    = "./inputs/km67_obs_Macos"
 fnm_obs67_2  = "./inputs/km67_Miller_etal2011.csv"
 compset      = "I1PTCLM45"
 experiments  = (/"intact","RILlow","RILhigh","CLlow","CLhigh"/)
 nv           = 5+4; sim+obs83+obs67+obs83_2+obs67_2
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
  RAIN  = new((/nv, nm/), float, FillValue); atmospheric rain,  w/m^2
  SOILWATER_10CM = new((/nv, nm/), float, FillValue); volumetric soil moisture, first 10cm (kg/m2)
  BTRAN = new((/nv, nm/), float, FillValue); transpiration beta factor
  TWS = new((/nv, nm/), float, FillValue); total water storage, mm
  TSOI_10CM = new((/nv, nm/), float, FillValue); soil temperature in top 10cm of soil (K)
  LAI  = new((/nv, nm/), float, FillValue); leaf area index
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
  ;confidence interval
  CI_Rn  = new((/nv, nm/), float, FillValue); confidence interval of net radiation, w/m^2
  CI_SH  = new((/nv, nm/), float, FillValue); confidence interval of sensible heat flux, w/m^2
  CI_LH  = new((/nv, nm/), float, FillValue); confidence interval of latent heat flux, w/m^2
  CI_NEE = new((/nv, nm/), float, FillValue); confidence interval of NEE, gC/m^2/day
  CI_GPP = new((/nv, nm/), float, FillValue); confidence interval of GPP, gC/m^2/day
  CI_ER  = new((/nv, nm/), float, FillValue); confidence interval of ER, gC/m^2/day

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
	
;starting from line 2 for Miller et al. obs
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
 nls2=1
 lines_obs67_2 = asciiread(fnm_obs67_2,-1,"string")
 n_obs2 = 60
 printVarSummary(lines_obs67_2)
 factor_E2LH = 2.5*18 ; mmol/m2/s -> W/m2
; factor_Fc   = 12*1e-6 *24*3600  ; umol/m2/s -> gC/m2/day
 LAI(8,0:n_obs2-1)    = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),3,delim))
 Rn(8,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),6,delim))
 SH(8,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),14,delim))
 LH(8,0:n_obs2-1)     = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),16,delim))*factor_E2LH
 LH(8,0:11) = -999.
 LH(8, 47) = -999.
 CI_Rn(8,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),26,delim))
 CI_SH(8,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),31,delim))
 CI_LH(8,0:n_obs2-1)  = tofloat(str_get_field(lines_obs67_2(nls2:nls2+n_obs2-1),32,delim))*factor_E2LH
 CI_LH(8,0:11) = -999.
 CI_LH(8, 47) = -999.

;retrieve simulations
varlist = (/"Rn", "SH", "LH", "GH", "SOILWATER_10CM", "TSOI_10CM", "LAI","fcev","fctr","fgev","QOVER","QINFL","QDRIP","QDRAI","QCHARGE","QINTR","QVEGE","QVEGT","QSOIL","RAIN","TWS","BTRAN", "fsr","fire"/)
nvar= 24
ncol = nv-2
do i=0,nvar-1
   fnm = path_o+varlist(i)+"_"+site+".csv"
   print(fnm)
   data = readAsciiTable(fnm, ncol, "float", 1)
   do j=0,nv-5
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

Rnmean=dim_avg(Rn(:,0:n_obs2-1))
SHmean=dim_avg(SH(:,0:n_obs2-1))
LHmean=dim_avg(LH(:,0:n_obs2-1))
LHstddev = dim_stddev(LH(:,0:n_obs2-1))
SHstddev = dim_stddev(SH(:,0:n_obs2-1))
Rnstddev = dim_stddev(Rn(:,0:n_obs2-1))
CI_Rnmean=dim_avg(CI_Rn(:,0:n_obs2-1))
CI_SHmean=dim_avg(CI_SH(:,0:n_obs2-1))
CI_LHmean=dim_avg(CI_LH(:,0:n_obs2-1))
print((/Rnmean(0), Rnstddev(0),Rnmean(1),Rnstddev(1), Rnmean(2), Rnstddev(2),Rnmean(3),Rnstddev(3),Rnmean(4),Rnstddev(4),Rnmean(8), Rnstddev(8), Rnmean(7), Rnstddev(7), CI_Rnmean(8), CI_Rnmean(7)/))
print((/LHmean(0), LHstddev(0),LHmean(1),LHstddev(1), LHmean(2), LHstddev(2),LHmean(3),LHstddev(3),LHmean(4),LHstddev(4),LHmean(8), LHstddev(8), LHmean(7), LHstddev(7), CI_LHmean(8), CI_LHmean(7)/))
print((/SHmean(0), SHstddev(0),SHmean(1),SHstddev(1), SHmean(2), SHstddev(2),SHmean(3),SHstddev(3),SHmean(4),SHstddev(4),SHmean(8), SHstddev(8), SHmean(7), SHstddev(7), CI_SHmean(8), CI_SHmean(7)/))

;************************************************
; create plots
;************************************************
  wks2 = gsn_open_wks("pdf",path_o+"WaterEnergy"+site)
  gsn_define_colormap(wks2,"gui_default")        ; select color map
; define variables for the plots
  plot2 = new(12, graphic)

  xlb=0
  xub=10*12-1
  xint=12
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
  res@xyExplicitLegendLabels=(/"Intact","Obs67","Obs67"/)
  res@xyMarkLineModes =  (/"Lines","Lines","Markers"/)   ; Three cases of sensitivity analysis
  res@xyDashPattern     = (/0,0/)
  res@xyLineThicknesses = (/2.0,1.0/)
  ;res@gsnXRefLine = 8
  ;res@gsnXRefLineColor = (/"black"/)
  ;res@gsnXRefLineDashPattern = 1
  res@xyLineColors      = (/"blue","red"/)
  res@xyMarker = 3
  res@xyMarkerColor =(/"red"/)
  res@xyMarkerSizeF = 0.002
; define variables for the plots
  ;plotting ranges
  plot_time    = ispan(0,nm-1,1)

;plot2
  res1 = res
  delete(res1@xyExplicitLegendLabels)
  delete(res1@xyMarkLineModes)
  delete(res1@xyLineThicknesses)
  delete(res1@xyLineColors)
  res1@xyExplicitLegendLabels=(/"RIL~B~low","RIL~B~high","CL~B~low","CL~B~high","Obs83"/)
  res1@xyMarkLineModes =  (/"Lines","Lines","Lines","Lines","Lines","Markers"/)
  res1@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0, 1.0/)
  res1@xyLineColors = (/"royalblue3","skyblue","green","orange","red"/)
  res1@xyDashPatterns     = (/0, 0, 0, 0, 0/)
  res1@xyMarker = 3
  res1@xyMarkerColor =(/"red"/)
  res@xyMarkerSizeF = 0.005
  res1@gsnXRefLine = 8
  res1@gsnXRefLineColor = (/"black"/)
  res1@gsnXRefLineDashPattern = 1


  ipleft=(/0,8/)
  ipright=(/1,2,3,4,7/)
  res@trYMinF  = 0.
  res@trYMaxF  = 180.
  xp  = new((/2*(xub-xlb-1)/), float)
  yp  = new((/2*(xub-xlb-1)/), float)
  res@trYMinF  = 0.
  res@trYMaxF  = 180.
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_SH(8,k)
     yp(k)                 = SH(8,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = SH(8,k) - dx
  end do
  plot2(0) = gsn_csm_xy(wks2, plot_time(xlb:xub), SH(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(0)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(a) Sensible heat flux (intact)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.09                 ; height of legend.
    "pmLegendOrthogonalPosF": -0.8               ; more neg = down
  end setvalues
  gsres3                   = True                        ; poly res
  gsres3@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres3@gsFillColor       = "pink"                      ; color chosen
  dummy3 = gsn_add_polygon (wks2,plot2(0),xp,yp,gsres3)

  res1@trYMinF  = 0.
  res1@trYMaxF  = 180.
  do k=xlb,xub
     dx                    = CI_SH(7,k)
     yp(k)                 = SH(7,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = SH(7,k) - dx
  end do
  plot2(1) = gsn_csm_xy(wks2, plot_time(xlb:xub), SH(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(1)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(b) Sensible heat flux (logged)"
    "pmLegendDisplayMode":   "Always"             ; turn on legend
    "pmLegendSide"       :   "Top"                ; Change location o
    "pmLegendWidthF"     :    0.2                 ; Change width and
    "pmLegendHeightF"    :    0.09                 ; height of legend.
    "pmLegendParallelPosF":   0.80                ; move units right
    "pmLegendOrthogonalPosF": -0.8               ; more neg = down
  end setvalues
  gsres4                   = True                        ; poly res
  gsres4@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres4@gsFillColor       = "pink"                      ; color chosen
  dummy4 = gsn_add_polygon (wks2,plot2(1),xp,yp,gsres4)

  res@trYMinF  = 0.
  res@trYMaxF  = 180.
  ;add CI from Miller et al
  do k=xlb,xub
     dx                    = CI_LH(8,k)
     yp(k)                 = LH(8,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = LH(8,k) - dx
  end do
  plot2(2) = gsn_csm_xy(wks2, plot_time(xlb:xub), LH(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(2)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(c) Latent heat flux (intact)"
  end setvalues
  gsres5                   = True                        ; poly res
  gsres5@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres5@gsFillColor       = "pink"                      ; color chosen
  dummy5 = gsn_add_polygon (wks2,plot2(2),xp,yp,gsres5)

  res1@trYMinF  = 0.
  res1@trYMaxF  = 180.
  do k=xlb,xub  
     dx                    = CI_LH(7,k)
     yp(k)                 = LH(7,k) + dx
     xp(k)                 = plot_time(k)
     xp(2*(xub-xlb-1)-1-k) = plot_time(k)
     yp(2*(xub-xlb-1)-1-k) = LH(7,k) - dx
  end do
  plot2(3) = gsn_csm_xy(wks2, plot_time(xlb:xub), LH(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(3)
    "tiYAxisString" : "W m~S~-2"   ; y-axis title
    "tiMainString": "(d)  Latent heat flux (logged)"
  end setvalues
  gsres6                   = True                        ; poly res
  gsres6@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres6@gsFillColor       = "pink"                      ; color chosen
  dummy6 = gsn_add_polygon (wks2,plot2(3),xp,yp,gsres6)

  res@trYMinF  = 0.
  res@trYMaxF  = 6.
  delete(ipleft)
  ipleft=(/0,8,8/)
  ;res@gsnYRefLine = 4.8
  ;res@gsnYRefLineColor = (/"red"/)
  ;res@gsnYRefLineDashPattern = 1
  plot2(4) = gsn_csm_xy(wks2, plot_time(xlb:xub), LAI(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(4)
    "tiYAxisString" : "m~S~2~N~m~S~-2"   ; y-axis title
    "tiMainString": "(e) Leaf area index (intact)"
  end setvalues

  res1@trYMinF  = 1.
  res1@trYMaxF  = 6.
  ;res1@gsnYRefLine = 5.5
  ;res1@gsnYRefLineColor = (/"red"/)
  ;res1@gsnYRefLineDashPattern = 1
  plot2(5) = gsn_csm_xy(wks2, plot_time(xlb:xub), LAI(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(5)
    "tiYAxisString" : "m~S~2~N~m~S~-2"   ; y-axis title
    "tiMainString": "(f) Leaf area index (logged)"
  end setvalues
  ;delete(res@gsnYRefLine)
  ;delete(res1@gsnYRefLine)

  res@trYMinF  = 0.2
  res@trYMaxF  = 0.5
  plot2(6) = gsn_csm_xy(wks2, plot_time(xlb:xub), SOILWATER_10CM(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(6)
    "tiYAxisString" : "m~S~3~N~m~S~-3"   ; y-axis title
    "tiMainString": "(g) Soil water in top 10 cm (intact)"
  end setvalues

  res1@trYMinF  = 0.2
  res1@trYMaxF  = 0.5
  plot2(7) = gsn_csm_xy(wks2, plot_time(xlb:xub), SOILWATER_10CM(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(7)
    "tiYAxisString" : "m~S~3~N~m~S~-3"   ; y-axis title
    "tiMainString": "(h) Soil water in top 10 cm (logged)"
  end setvalues

  res@trYMinF  = 0.
  res@trYMaxF  = 1. 
  plot2(8) = gsn_csm_xy(wks2, plot_time(xlb:xub), BTRAN(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(8)
    "tiYAxisString" : "unitless"   ; y-axis title
    "tiMainString": "(i) Transpiration beta factor (intact)"
  end setvalues

  res1@trYMinF  = 0.
  res1@trYMaxF  = 1.
  plot2(9) = gsn_csm_xy(wks2, plot_time(xlb:xub), BTRAN(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(9)
    "tiYAxisString" : "unitless"   ; y-axis title
    "tiMainString": "(j) Transpiration beta factor (logged)"
  end setvalues

  res@trYMinF  = 20.
  res@trYMaxF  = 35.
  plot2(10) = gsn_csm_xy(wks2, plot_time(xlb:xub), TSOI_10CM(ipleft,xlb:xub), res  )   ; create plot
  setvalues plot2(10)
    "tiYAxisString" : "degree C"   ; y-axis title
    "tiMainString": "(k) 10-cm soil temperature (intact)"
  end setvalues

  res1@trYMinF  = 20.
  res1@trYMaxF  = 35.
  plot2(11) = gsn_csm_xy(wks2, plot_time(xlb:xub), TSOI_10CM(ipright,xlb:xub), res1  )   ; create plot
  setvalues plot2(11)
    "tiYAxisString" : "degree C"   ; y-axis title
    "tiMainString": "(l) 10-cm soil temperature (logged)"
  end setvalues

; create panel
;*************************************************
  resP = True   ; modify the panel plot
  ;resP@txString = "State variables" ; add common title
  resP@gsnMaximize = True          ; use full page

  gsn_panel(wks2, plot2, (/6, 2/), resP)

end
