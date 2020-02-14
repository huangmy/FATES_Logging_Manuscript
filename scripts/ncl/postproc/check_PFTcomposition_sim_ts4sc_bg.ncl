;******************************************************************
;  Compare simulated energy fluxes with observations 
;    Author: maoyi.huang@pnnl
;    Date :  10/11/2017
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
 path_archive = "/pic/scratch/huan565/cesm_archive/"
 path_o       = "./outputs/"
 compset      = "I2000Clm50FatesGs"
 ;experiments   = (/"intact","RILlow","RILhigh","CLlow","CLhigh"/)
 experiments   = (/"intact.ensm0063_sim", "intact.ensm0063.RILlow", "intact.ensm0063.RILhigh", "intact.ensm0063.CLlow", "intact.ensm0063.CLhigh"/)
 nv           = 5
 npft         = 2;  2 PFTs
 nint         = 1
 ny           = (year_end - year_start +1)/nint ;# of nyear
 nscpfp       = 4; size classes for plottting
 nscpf        = 13; size classes in fates
 ; size class lower bound 
 ;    0.000 5.000 10.00 15.00 20.00 30.00 40.00 50.00 60.00 70.00  80.00  90.00 100.0
 
 ;****************************************************
 ; Modify variables for comparison
  basal_area_early  = new((/nscpfp, nv, ny/), float); basal area  m2/ha
  basal_area_late   = new((/nscpfp, nv, ny/), float); basal area  m2/ha
  
  stem_density_early = new((/nscpfp, nv, ny/), float); stem density N/ha
  stem_density_late = new((/nscpfp, nv, ny/), float); stem density N/ha
  stem_density_canopy = new((/nscpfp, nv, ny/), float); stem density N/ha
  stem_density_understory = new((/nscpfp, nv, ny/), float); stem density N/ha
  ;stem_density_late_canopy = new((/nscpfp, nv, ny/), float); stem density N/ha
  ;stem_density_late_understory = new((/nscpfp, nv, ny/), float); stem density N/ha

  ddbh_early   = new((/nscpfp, nv, ny/), float); diameter growth increment by pft/size, cm/yr/ha
  ddbh_late    = new((/nscpfp, nv, ny/), float); diameter growth increment by pft/size, cm/yr/ha

  mortality_canopy_early   = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by pft/size, N/ha/yr
  mortality_canopy_late   = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by pft/size, N/ha/yr

  mortality_understory_early   = new((/nscpfp, nv, ny/), float); total mortality of understory plants by pft/size, N/ha/yr
  mortality_understory_late   = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by pft/size, N/ha/yr

  mortality_early = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by size
  mortality_late   = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by size 	
  mortality_canopy = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by size
  mortality_understory = new((/nscpfp, nv, ny/), float); total mortality of canopy plants by size
;retrieve simulations
;**********************************************************
;assign FATES values to the plotting variable, UTC time
;**********************************************************

 filename     = new((/nv/), string)

;size class lower and upper boundaries 
sz_l = (/0, 2, 5, 7/); 
sz_u = (/1, 4, 6, 12/); 
im = 8 ; only plot for the beginning of the year

do iscpfp =0, nscpfp-1
do i=0,nv-1
   casename = "fates-clm5_1x1pt"+"_" + site + "_" + experiments(i) + "_" + compset
   do iy=0,ny-1
      ;read file
       filename(i)=path_archive + casename + "/lnd/hist/"+ casename + ".clm2.h0."+ sprinti("%0.4i", iy+year_start) + "-"+ sprinti("%0.2i", 2) + "-01-00000.nc"
       print(filename(i))
       c = addfiles(filename(i),"r")
      ;read value from  FATES history
       basal_area_early(iscpfp,i,iy)   = sum(c[:]->BA_SCPF(im,sz_l(iscpfp):sz_u(iscpfp), 0 ) )
       basal_area_late(iscpfp,i,iy)    = sum(c[:]->BA_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf), 0 ))
       stem_density_early(iscpfp,i,iy)  = sum(c[:]->NPLANT_SCPF(im,sz_l(iscpfp):sz_u(iscpfp) , 0 ))
       stem_density_late(iscpfp,i,iy)   = sum(c[:]->NPLANT_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf),0)) 
       stem_density_canopy(iscpfp,i,iy)   = sum(c[:]->NPLANT_CANOPY_SCLS(im,sz_l(iscpfp):sz_u(iscpfp) , 0 ) )
       stem_density_understory(iscpfp,i,iy)   = sum(c[:]->NPLANT_UNDERSTORY_SCLS(im,sz_l(iscpfp):sz_u(iscpfp) , 0 ) )
       ddbh_early(iscpfp,i,iy)   = sum(c[:]->DDBH_SCPF(im,sz_l(iscpfp):sz_u(iscpfp),0 ) )/stem_density_early(iscpfp,i,iy)
       ddbh_late(iscpfp,i,iy)    = sum(c[:]->DDBH_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf), 0 ) )/stem_density_late(iscpfp,i,iy)
       if (stem_density_early(iscpfp,i,iy).gt.0.5) then
          mortality_early(iscpfp,i,iy)   = 100*(sum(c[:]->MORTALITY_CANOPY_SCPF(im,sz_l(iscpfp):sz_u(iscpfp) ,0 ) )+sum(c[:]->MORTALITY_UNDERSTORY_SCPF(im,sz_l(iscpfp):sz_u(iscpfp) , 0 ) )) / stem_density_early(iscpfp,i,iy)
       end if
       if (stem_density_late(iscpfp,i,iy).gt.0.5) then
          mortality_late(iscpfp,i,iy)   = 100*(sum(c[:]->MORTALITY_CANOPY_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf), 0 ) )+sum(c[:]->MORTALITY_UNDERSTORY_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf),0 ) ))/ stem_density_late(iscpfp,i,iy)
       end if
       if (stem_density_canopy(iscpfp,i,iy).gt.0.5) then
          mortality_canopy(iscpfp,i,iy)    = 100*(sum(c[:]->MORTALITY_CANOPY_SCPF(im,sz_l(iscpfp):sz_u(iscpfp) ,0 ) )+sum(c[:]->MORTALITY_CANOPY_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf),0 ) ) )/ stem_density_canopy(iscpfp,i,iy)
       end if
       if (stem_density_understory(iscpfp,i,iy).gt.0.5) then
          mortality_understory(iscpfp,i,iy)   = 100*(sum(c[:]->MORTALITY_UNDERSTORY_SCPF(im,sz_l(iscpfp):sz_u(iscpfp) , 0 ) + sum(c[:]->MORTALITY_UNDERSTORY_SCPF(im,(sz_l(iscpfp)+nscpf):(sz_u(iscpfp)+nscpf), 0 ) )))/stem_density_understory(iscpfp,i,iy)
       end if
       if (mortality_early(iscpfp,i,iy).gt.1000) then
          mortality_early(iscpfp,i,iy) = mortality_early@_FillValue
       end if
       if (mortality_late(iscpfp,i,iy).gt.1000) then
          mortality_late(iscpfp,i,iy) = mortality_canopy_late@_FillValue
       end if
       if (mortality_understory(iscpfp,i,iy).gt.500) then
          mortality_understory(iscpfp,i,iy) = mortality_understory@_FillValue
       end if
       if (mortality_canopy(iscpfp,i,iy).gt.500) then
          mortality_canopy(iscpfp,i,iy) = mortality_canopy@_FillValue
       end if
       delete(c)
   end do; year
   delete(casename)
end do; nv
end do; pft


;*************************************************
; write variables into csv files
;*************************************************
lines = new(ny+1, string)
lines(0)=(/"  yr,    0-10,     10-30,     30-50,     >=50"/); header
YEAR = ispan(year_start,year_end,1)

do i=0,nv-1
   fnmo = path_o+"BA_early"+experiments(i)+site+".csv"
   system("rm -rf " + fnmo)
   lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",basal_area_early(0,i,:)) + "," + sprintf("%10.2f",basal_area_early(1,i,:))+ "," + sprintf("%10.2f", basal_area_early(2,i,:)) + "," + sprintf("%10.2f",basal_area_early(3,i,:)) 
   asciiwrite(fnmo,lines)
   fnmo = path_o+"BA_late"+experiments(i)+site+".csv"
   system("rm -rf " + fnmo)
   lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",basal_area_late(0,i,:)) + "," + sprintf("%10.2f",basal_area_late(1,i,:))+ "," + sprintf("%10.2f", basal_area_late(2,i,:)) + "," + sprintf("%10.2f",basal_area_late(3,i,:))
   asciiwrite(fnmo,lines)
   fnmo = path_o+"stem_density_early"+experiments(i)+site+".csv"
   system("rm -rf " + fnmo)
   lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",stem_density_early(0,i,:)) + "," + sprintf("%10.2f",stem_density_early(1,i,:))+ "," + sprintf("%10.2f", stem_density_early(2,i,:)) + "," + sprintf("%10.2f",stem_density_early(3,i,:))
   asciiwrite(fnmo,lines)
   fnmo = path_o+"stem_density_late"+experiments(i)+site+".csv"
   system("rm -rf " + fnmo)
   lines(1:ny) = sprinti("%0.4i",YEAR) + "," + sprintf("%10.2f",stem_density_late(0,i,:)) + "," + sprintf("%10.2f",stem_density_late(1,i,:))+ "," + sprintf("%10.2f", stem_density_late(2,i,:)) + "," + sprintf("%10.2f",stem_density_late(3,i,:))
   asciiwrite(fnmo,lines)
end do
	
;************************************************
; create plots
;************************************************
  plot_time    = ispan(0,ny*nint-1,nint)
  nmp = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)
  nmpLeft   = (/"(a)","(c)","(e)","(g)","(i)"/)
  nmpRight  = (/"(b)","(d)","(f)","(h)","(j)"/)
  scpf    = (/"0-10","10-30", "30-50",">=50"/)
  iexp = ispan(0,nv-1,1)

; define variables for the plots
  xlb=0
  xub=ny*nint
  xint=10
  res = True
  res@vpWidthF = 0.65                           ; change aspect ratio of plot
  res@vpHeightF = 0.2                           ;
  res@tiXAxisString = "Year"                    ; x-axis title
  res@trXMinF  = xlb
  res@trXMaxF  = xub
  res@tmXBMinorPerMajor = 0
  res@tmXBMode = "Explicit"
  res@tmXBTickEndF = xub
  res@tmXBTickSpacingF = xint
  res@tmXBTickStartF = 0
  res@tmXBValues = (/xlb, xlb+xint, xlb+2*xint,xlb+3*xint,xlb+4*xint,xlb+5*xint,xlb+6*xint,xlb+7*xint,xlb+8*xint,xlb+9*xint,xlb+ny*nint/)
  res@tmXBLabels = (/2001,2011,2021,2031,2041,2051,2061,2071,2081,2091,2101/)
  res@xyLineColors      = (/"blue", "royalblue3","skyblue","green","orange"/)  ;1st one will be ignored
  ;res@xyLineColors      = (/"green","orange"/)
  res@xyMarkers = 16
  res@xyMarkerColor =(/"red"/)
  res@xyMarkerSizeF = 0.005
  res@xyDashPattern     = 1
  res@xyExplicitLegendLabels=(/"Intact","RIL~B~low","RIL~B~high","CL~B~low","CL~B~high"/)
  ;res@xyExplicitLegendLabels=(/"RIL~B~low","RIL~B~high","CL~B~low","CL~B~high"/)
  res@xyMarkLineModes =  (/"Lines","Lines","Lines","Lines", "Lines"/)   ; Three cases of sensitivity analysis
  ;res@xyExplicitLegendLabels=(/"Early","Late"/)
  ;res@xyMarkLineModes =  (/"Lines","Lines"/)
  res@xyDashPattern     = 1
  res@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0, 2.0/)
  ;res@xyLineThicknesses = (/2.0, 2.0/)
  res@gsnXRefLine = 0.75
  res@gsnXRefLineColor = (/"black"/)
  res@gsnXRefLineDashPattern = 1
  res@trYMinF  = 0.
  res@trYMaxF  = 50.

  print("plot1")
; plot1
   plot1 = new(10, graphic)   ;
   wks1 = gsn_open_wks("pdf",path_o+"FATES_BA_sim_4sc"+site)
   gsn_define_colormap(wks1,"gui_default")        ; select color map
   ;Basal Area from inventory
   BAearly_obs = (/-999.,2.5,2.6,6.7/)
   BAlogearly_obs = (/-999.,2.2,2.4,5.8/)
   BAlate_obs  = (/-999.,1.7,1,6,5.9/)
   BAloglate_obs = (/-999.,1.7,1.6,5.1/)

   do iscpfp=0,nscpfp-1
      if (iscpfp.ge.1) then
         res@gsnYRefLine = (/BAearly_obs(iscpfp),BAlogearly_obs(iscpfp)/)
         res@gsnYRefLineColor = (/"red","cyan"/)
         res@gsnYRefLineDashPatterns = (/0,1/)
      end if
      if (iscpfp.ge.1) then
         res@trYMinF  = 0
         res@trYMaxF  = 10
      else
         res@trYMinF  = 0
         res@trYMaxF  = 0.2
      end if
      plot1(iscpfp*2) = gsn_csm_xy(wks1, plot_time, basal_area_early(iscpfp,iexp,:), res  )   ; create plot
      setvalues plot1(iscpfp*2)
        "tiYAxisString" : "m~S~2~N~ha~S~-1"   ; y-axis title
        ;"tiYAxisString" : "%" ; y-axis title
        "tiMainString": nmpLeft(iscpfp)+" Basal area,"+scpf(iscpfp)+"cm, early"
      end setvalues
      delete(res@gsnYRefLine)
      if (iscpfp.eq.0) then
         res@trYMinF  = 0
         res@trYMaxF  = 3
      end if
      if (iscpfp.eq.1.or.iscpfp.eq.2) then
         res@trYMinF  = 0
         res@trYMaxF  = 10
      end if
      if (iscpfp.eq.3) then
         res@trYMinF  = 0
         res@trYMaxF  = 45
      end if
      if (iscpfp.ge.1) then
         res@gsnYRefLine =  (/BAlate_obs(iscpfp),BAloglate_obs(iscpfp)/)
      end if
      plot1(iscpfp*2+1) = gsn_csm_xy(wks1, plot_time, basal_area_late(iscpfp,iexp,:), res  )   ;
      setvalues plot1(iscpfp*2+1)
        "tiYAxisString" : "m~S~2~N~ha~S~-1"   ; y-axis title
        ;"tiYAxisString" : "%" ; y-axis title
        "tiMainString": nmpRight(iscpfp)+" Basal area,"+scpf(iscpfp)+"cm, late"
      end setvalues
      delete(res@gsnYRefLine)
   end do; nscpfp
      delete(res@trYMinF)
      delete(res@trYMaxF)
      delete(res@gsnYRefLine)
      delete(res@gsnYRefLineColor)
      setvalues plot1(0)
           "pmLegendDisplayMode":   "Always"             ; turn on legend
           "pmLegendSide"       :   "Top"                ; Change location o
           "pmLegendWidthF"     :    0.16                 ; Change width and
           "pmLegendHeightF"    :    0.08                ; height of legend.
           "pmLegendParallelPosF":   0.82                ; move units right
           "pmLegendOrthogonalPosF": -0.80               ; more neg = down
      end setvalues

  print("plot2")
   plot2 = new(10, graphic)   ;
   wks2 = gsn_open_wks("pdf",path_o+"FATES_SD_sim_4sc"+site)
   SD_obs         = (/-999,399,30,30/)
   pctSDearly_obs  = (/-999,0.58,0.59,0.54/)
   SDlog_obs      =(/-999,396,29,18/)
   pctSDlogearly_obs  = (/-999,0.58,0.59,0.32/)
   gsn_define_colormap(wks2,"gui_default")        ; select color map
   do iscpfp=0,nscpfp-1
      if (iscpfp.ge.1) then
         res@gsnYRefLine = (/SD_obs(iscpfp),SDlog_obs(iscpfp)/)
         res@gsnYRefLineColor = (/"red","cyan"/)
         res@gsnYRefLineDashPatterns = (/0,1/)
      end if
      if (iscpfp.eq.1) then
         res@trYMinF  = 0
         res@trYMaxF  = 450
      end if
      if (iscpfp.eq.2) then
         res@trYMinF  = 0
         res@trYMaxF  = 100
      end if
      if (iscpfp.eq.3) then
         res@trYMinF  = 0
         res@trYMaxF  = 100
      end if
      plot2(iscpfp*2) = gsn_csm_xy(wks2, plot_time, stem_density_early(iscpfp,iexp,:)+stem_density_late(iscpfp,iexp,:),res )
      setvalues plot2(iscpfp*2)
        "tiYAxisString" : "N~N~ha~S~-1"   ; y-axis title
        ;"tiYAxisString" : "%" ; y-axis title
        "tiMainString": nmpLeft(iscpfp)+" Stem density,"+scpf(iscpfp)+"cm"
      end setvalues
      if (iscpfp.ge.1) then
          delete(res@gsnYRefLine)
      end if
      if (iscpfp.ge.1) then
          res@gsnYRefLine = (/100 * pctSDearly_obs(iscpfp),100 * pctSDlogearly_obs(iscpfp)/)
      end if
      if (iscpfp.ge.2) then
         res@trYMinF  = 0
         res@trYMaxF  = 70
      else
         res@trYMinF  = 0
         res@trYMaxF  = 70
      end if
      plot2(iscpfp*2+1) = gsn_csm_xy(wks2, plot_time, 100*stem_density_early(iscpfp,iexp,:)/(stem_density_early(iscpfp,iexp,:)+stem_density_late(iscpfp,iexp,:)),res )
     setvalues plot2(iscpfp*2+1)
       ;"tiYAxisString" : "N~N~ha~S~-1"   ; y-axis title
        "tiYAxisString" : "%" ; y-axis title
        "tiMainString": nmpRight(iscpfp)+" Fraction of early successional,"+scpf(iscpfp)+"cm"
      end setvalues
      if (iscpfp.ge.1) then
         delete(res@gsnYRefLine)
      end if
   end do; nscpfp
      setvalues plot2(1)
           "pmLegendDisplayMode":   "Always"             ; turn on legend
           "pmLegendSide"       :   "Top"                ; Change location o
           "pmLegendWidthF"     :    0.2                 ; Change width and
           "pmLegendHeightF"    :    0.1                 ; height of legend.
           "pmLegendParallelPosF":   0.82                ; move units right
           "pmLegendOrthogonalPosF": -0.98               ; more neg = down
      end setvalues

   delete(res@trYMinF)
   delete(res@trYMaxF)
   delete(res@gsnYRefLine)
   delete(res@gsnYRefLineColor)

  print("plot3")
   plot3 = new(10, graphic)   ;
   wks3 = gsn_open_wks("pdf",path_o+"FATES_MTearlylate_sim_4sc"+site)
   gsn_define_colormap(wks3,"gui_default")        ; select color map
   do iscpfp=0,nscpfp-1
      if (iscpfp.eq.0) then
          res@trYMinF  = 0.
          res@trYMaxF  = 20.
      else
          res@trYMinF  = 0.
          res@trYMaxF  = 10.
      end if
      x_runave = runave(mortality_early(iscpfp,iexp,:),5,0)
      plot3(iscpfp*2) = gsn_csm_xy(wks3, plot_time, x_runave, res  )   ; create plot
      setvalues plot3(iscpfp*2)
        ;"tiYAxisString" : "N ha~S~-1~N~yr~S~-1"   ; y-axis title
        "tiYAxisString" : "% yr~S~-1"   ; y-axis title
        "tiMainString": nmpLeft(iscpfp)+" Mortality,"+scpf(iscpfp)+"cm, early"
      end setvalues
      ;delete(res@trYMinF)
      ;delete(res@trYMaxF)
      x_runave = runave(mortality_late(iscpfp,iexp,:),5,0)
      plot3(iscpfp*2+1) = gsn_csm_xy(wks3, plot_time, x_runave, res  )   ; create plot
      setvalues plot3(iscpfp*2+1)
        ;"tiYAxisString" : "N ha~S~-1~N~yr~S~-1"   ; y-axis title
        "tiYAxisString" : "% yr~S~-1"   ; y-axis title
        "tiMainString": nmpRight(iscpfp)+" Mortality,"+scpf(iscpfp)+"cm, late"
      end setvalues
   end do; nscpfp
      setvalues plot3(3)
           "pmLegendDisplayMode":   "Always"             ; turn on legend
           "pmLegendSide"       :   "Top"                ; Change location o
           "pmLegendWidthF"     :    0.2                 ; Change width and
           "pmLegendHeightF"    :    0.1                 ; height of legend.
           "pmLegendParallelPosF":   0.82                ; move units right
           "pmLegendOrthogonalPosF": -0.90               ; more neg = down
      end setvalues

  print("plot4")
   plot4 = new(10, graphic)   ;
   wks4 = gsn_open_wks("pdf",path_o+"FATES_MTcanunder_sim_4sc"+site)
   gsn_define_colormap(wks4,"gui_default")        ; select color map
   do iscpfp=0,nscpfp-1
      if (iscpfp.le.1) then
          res@trYMinF  = 0.
          res@trYMaxF  = 15.
      end if
      if (iscpfp.ge.2) then
          res@trYMinF  = 0.
          res@trYMaxF  = 15.
      end if
      x_runave = runave(mortality_canopy(iscpfp,iexp,:),5,0)
      plot4(iscpfp*2) = gsn_csm_xy(wks4, plot_time, x_runave, res  )   ; create plot
      setvalues plot4(iscpfp*2)
        ;"tiYAxisString" : "N ha~S~-1~N~yr~S~-1"   ; y-axis title
        "tiYAxisString" : "% yr~S~-1"   ; y-axis title
        "tiMainString": nmpLeft(iscpfp)+" Mortality,"+scpf(iscpfp)+"cm, canopy"
      end setvalues
      x_runave = runave(mortality_understory(iscpfp,iexp,:),5,0)
      plot4(iscpfp*2+1) = gsn_csm_xy(wks4, plot_time, x_runave, res  )   ; create plot
      setvalues plot4(iscpfp*2+1)
        ;"tiYAxisString" : "N ha~S~-1~N~yr~S~-1"   ; y-axis title
        "tiYAxisString" : "% yr~S~-1"   ; y-axis title
        "tiMainString": nmpRight(iscpfp)+" Mortality,"+scpf(iscpfp)+"cm, understory"
      end setvalues
      ;delete(res@trYMinF)
      ;delete(res@trYMaxF)
   end do; nscpfp
      setvalues plot4(0)
           "pmLegendDisplayMode":   "Always"             ; turn on legend
           "pmLegendSide"       :   "Top"                ; Change location o
           "pmLegendWidthF"     :    0.2                 ; Change width and
           "pmLegendHeightF"    :    0.1                 ; height of legend.
           "pmLegendParallelPosF":   0.82                ; move units right
           "pmLegendOrthogonalPosF": -0.90               ; more neg = down
      end setvalues

  print("plot5")
   plot5 = new(10, graphic)   ;
   wks5 = gsn_open_wks("pdf",path_o+"FATES_DDBH_sim_4sc"+site)
   gsn_define_colormap(wks5,"gui_default")        ; select color map
   do iscpfp=0,nscpfp-1
      if (iscpfp.eq.0) then
          res@trYMinF  = 0.0
          res@trYMaxF  = 0.5
      else
          res@trYMinF  = 0
          res@trYMaxF  = 3.0
      end if
      plot5(iscpfp*2) = gsn_csm_xy(wks5, plot_time, ddbh_early(iscpfp,iexp,:), res  )   ; create plot
      setvalues plot5(iscpfp*2)
        "tiYAxisString" : "cm yr~S~-1~N~ha~S~-1"   ; y-axis title
        "tiMainString": nmpLeft(iscpfp)+" Diameter increment,"+scpf(iscpfp)+"cm, early"
      end setvalues
      if (iscpfp.eq.0) then
          res@trYMinF  = 0.0
          res@trYMaxF  = 0.5
      else
          res@trYMinF  = 0
          res@trYMaxF  = 3.0
      end if
      plot5(iscpfp*2+1) = gsn_csm_xy(wks5, plot_time, ddbh_late(iscpfp,iexp,:), res  )   ; create plot
      setvalues plot5(iscpfp*2+1)
        "tiYAxisString" : "cm yr~S~-1~N~ha~S~-1"   ; y-axis title
        "tiMainString": nmpRight(iscpfp)+" Diameter increment,"+scpf(iscpfp)+"cm"
      end setvalues
   end do; nscpfp
      setvalues plot5(0)
           "pmLegendDisplayMode":   "Always"             ; turn on legend
           "pmLegendSide"       :   "Top"                ; Change location o
           "pmLegendWidthF"     :    0.16                 ; Change width and
           "pmLegendHeightF"    :    0.08                 ; height of legend.
           "pmLegendParallelPosF":   0.82                ; move units right
           "pmLegendOrthogonalPosF": -0.80               ; more neg = down
      end setvalues


   ;*************************************************
   ; create panel
   ;*************************************************
   resP = True   ; modify the panel plot
   ;resP@gsnMaximize = True           ; use full page
   gsn_panel(wks1, plot1, (/5, 2/), resP)
   gsn_panel(wks2, plot2, (/5, 2/), resP)
   gsn_panel(wks3, plot3, (/5, 2/), resP)
   gsn_panel(wks2, plot4, (/5, 2/), resP)
   gsn_panel(wks3, plot5, (/5, 2/), resP)

end

