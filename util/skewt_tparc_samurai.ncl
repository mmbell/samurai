 ;**************************************************
 ; skewt_tparc_samurai.ncl
 ; Plot dropsondes vs. Samurai-output
 ; Directory of samurai-data needs to be specified below.
 ; To define absolute x,y max/min bounds adapt lines 177 pp.
 ; in $NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl
 ;**************************************************
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 ;**************************************************

   begin

 ; --- Read Data; set to _FillValue------------;
     if (.not. isvar("dir")) then      ; is diri on command line?
       dir = "/home/quinting/Diplomarbeit/Daten/Dropsondendaten_P3/"
       print("CALL ncl dir=<PATH> fname=<file> skewt_tparc_ecmwf.ncl")       
       exit
     end if
     if (.not. isvar("fname")) then      ; is diri on command line?
       fname = "D20080919_012914_PQC.eol.Wwind"
       print("CALL ncl dir=<PATH> fname=<file> skewt_tparc_ecmwf.ncl") 
      exit
     end if
 plot_dir=dir+"/plots/"

 inc=10 ;average model data in a 2*inc+1 square around closest GP
        ;depends on grid-resolution

 system("mkdir "+plot_dir)
 diri=dir
 fili=fname
 nlvl = numAsciiRow(diri+fili)
 ncol = 19
 UNDEF=-999.0
; filename, number of columns, type of column, number of lines to skip 
; at beginning of file
 TestData = readAsciiTable(diri+fili,ncol,"float",14)

 print("WORKING on file: "+fili)
 
 ;=====================
 ;for radiosondes
 ;=====================
 ;p  = TestData(0:,0)
 ;tc = TestData(0:,2)
 ;tcd = TestData(0:,3)
 ;z = TestData(0:,1)
 ;wspd = TestData(0:,7)
 ;wdir = TestData(0:,6)
 ;lon = 140.13      ;<---- location of radiosonde
 ;lat = 36.05      ;<---- location of radiosonde

 ;=====================
 ;for dropsondes
 ;=====================
 p  = TestData(0:,4)
 prev = p(nlvl-15:0)
 tc = TestData(0:,5)
 tcrev = tc(nlvl-15:0)
 tcd = TestData(0:,6)
 tcdrev = tcd(nlvl-15:0)
 z = TestData(0:,13)
 zrev = z(nlvl-15:0)
 wspd = TestData(::,10)
 wspdrev = wspd(nlvl-15:0)/0.514444 ;3.6/1.852 (windspeed in knots)
 wdir = TestData(::,11)
 wdirrev = wdir(nlvl-15:0)
 u = TestData(::,8)
 urev = u(nlvl-15:0)
 v = TestData(::,9)
 vrev = v(nlvl-15:0)
 lon = TestData(::,14)
 lonrev = lon(nlvl-15:0)
 lat = TestData(::,15)
 latrev = lat(nlvl-15:0)

 sonde_date=str_get_cols(fname,1,8)
 sonde_time=str_get_cols(fname,10,15)
 


;========================================================================
; open and average the SAMURAI data
;========================================================================
  samurai_dir = "/media/Drive-n-go/MAC_Julian/samurai_400x400_0500_0800_ecmwf_ds_ft_sat_uverr5_rerr3_res5_n/"

     
; SAMURAI FILE name
  grb_name = samurai_dir + "samurai_XYZ_analysis.nc"
  grb = addfile(grb_name,"r")


; see what variables are available in this file
  names=getfilevarnames(grb)
  ;print(names)
; transfer data and coordinates	
  t_samurai = grb->T
  qv_samurai= grb->QV
  fi_samurai= grb->altitude
  u_samurai= grb->U*3.6/1.852
  v_samurai= grb->V*3.6/1.852
  rh_samurai= grb->RH
  lat2d   = grb->latitude
  lon2d   = grb->longitude
  levels  = grb->P

  dims_samurai = dimsizes(t_samurai)

  t_samurai@lat2d = lat2d
  t_samurai@lon2d = lon2d
  qv_samurai@lat2d = lat2d
  qv_samurai@lon2d = lon2d
  fi_samurai@lat2d = lat2d
  fi_samurai@lon2d = lon2d
  v_samurai@lat2d = lat2d
  v_samurai@lon2d = lon2d
  u_samurai@lat2d = lat2d
  u_samurai@lon2d = lon2d
  rh_samurai@lat2d = lat2d
  rh_samurai@lon2d = lon2d

;calculate dew point Temp.
; dewpoint temperature in K
  td_samurai=dewtemp_trh(t_samurai,rh_samurai)
  copy_VarCoords(t_samurai,td_samurai)

;find closest grid point to Drop location
  lon_av_drop=(min(lon)+max(lon))/2;
  lat_av_drop=(min(lat)+max(lat))/2;
  print("DROP  location: "+lon_av_drop+", "+lat_av_drop)
  dx=abs(lon2d(1)-lon2d(0))
  dy=abs(lat2d(1)-lat2d(0))
  i_GP=round((lon_av_drop-lon2d(0))/dx,3)
  j_GP=abs(round((lat_av_drop-lat2d(0))/dy,3))


  if(j_GP-inc.lt.0.or.j_GP+inc.gt.(dims_samurai(2)).or.i_GP-inc.lt.0.or.i_GP+inc.gt.max(dims_samurai(3))) then
    print("ERROR: dropsonde data out off analysis domain")
    exit
  end if

  print("MODEL location: "+lon2d(i_GP)+", "+lat2d(j_GP))

; calculate stdev in a square around closet model grid point
  tmp  = t_samurai (:,1:dims_samurai(1)-1,j_GP-inc:j_GP+inc,i_GP-inc:i_GP+inc)
  tmp2 = td_samurai(:,1:dims_samurai(1)-1,j_GP-inc:j_GP+inc,i_GP-inc:i_GP+inc)
  ump  = u_samurai(:,1:dims_samurai(1)-1,j_GP-inc:j_GP+inc,i_GP-inc:i_GP+inc)
  vmp  = v_samurai(:,1:dims_samurai(1)-1,j_GP-inc:j_GP+inc,i_GP-inc:i_GP+inc)
  dims=dimsizes(tmp)
  
  t_samurai_avg  = new(dims(1),float)
  td_samurai_avg = new(dims(1),float)
  u_samurai_avg  = new(dims(1),float)
  v_samurai_avg  = new(dims(1),float)
  t_samurai_std  = new(dims(1),float)
  td_samurai_std = new(dims(1),float)
  p_samurai      = new(dims(1),float)
  
  do k=0,dims(1)-2
    t_samurai_avg(k)  = avg(tmp(0,k,:,:))
    td_samurai_avg(k) = avg(tmp2(0,k,:,:)) 
    t_samurai_std(k)  = stddev(tmp(0,k,:,:))
    td_samurai_std(k) = stddev(tmp2(0,k,:,:))
    u_samurai_avg(k)  = ump(0,k,inc,inc)
    v_samurai_avg(k)  = vmp(0,k,inc,inc)
    p_samurai(k)      = 1013.25*(1-0.0065*fi_samurai(k+1)*1000/t_samurai_avg(k))^5.255 ;calculate pressure in box around model location
  end do

;========================================================================
 ;*************************
 ; create plot
 ;*************************
  wks_type = "ps" 
  wks_type@wkOrientation = "portrait"
  plot_name=fili+"_ncl"
  wks = gsn_open_wks(wks_type,plot_dir+plot_name+"_samurai")

   gsn_define_colormap(wks,"hlu_default")
 ; --- Create background skew-T; plot sounding ---------------------
 ; for further information visit: http://www.ncl.ucar.edu/Applications/skewt.shtml
    skewtOpts				= True
;   skewtOpts@DrawColAreaFill		= True    	; default is False
    skewtOpts@tiMainString		= sonde_date+sonde_time+" vs samurai-analysis"
    skewtOpts@tiMainFontHeightF  	= 0.02
    skewtOpts@DrawFahrenheit 		= False
    skewtOpts@DrawColLine 		= True
;    skewtOpts@DrawIsotherm      	= True                                       
;    skewtOpts@DrawIsobar        	= True      
;    skewtOpts@DrawMixRatio      	= True         
;    skewtOpts@DrawDryAdiabat    	= True  
;    skewtOpts@DrawMoistAdiabat  	= True  	; aka: saturation or pseudo adibat     
;    skewtOpts@DrawWind          	= False
    skewtOpts@DrawStandardAtm   	= False      
;    skewtOpts@DrawColLine       	= True      
;    skewtOpts@DrawColAreaFill   	= False       
;    skewtOpts@DrawFahrenheit    	= False 	; True  ; Fahrenheit "x" axis      
    skewtOpts@DrawHeightScale   	= False
    skewtOpts@DrawHeightScaleFt 	= False		; True  ;default is feet [otherwise km] 
    
   
    skewt_bkgd = skewT_BackGround_Original(wks, skewtOpts)
    print("Plot background")
    draw(skewt_bkgd)

; ===================
; the sounding

    dataOpts 				= True
    dataOpts@PrintZ 			= True
    dataOpts@Wthin 			= 20 		;reduces number of windbarbs
    dataOpts@xpWind    			= 33.      	; new location for winds [default 45] in celsius
    dataOpts@WspdWdir 			= True		; use u,v instead of wspd, wdir
    dataOpts@HspdHdir			= True
    dataOpts@DrawFahrenheit   		= False    	; default is True
    dataOpts@ThermoInfo       		= False    	; do not print thermodynamic info (CAPE, CIN)
    dataOpts@offsetThermoInfo 		= 0.011 	;offset for Thermoinfo default 0
    dataOpts@colCape			="orange"
    dataOpts@linePatternCape 		= 1  
    dataOpts@colTemperature 		= "gray35"
    print("Plot sounding")
    plot = skewT_PlotData(wks,skewt_bkgd,prev,tcrev,tcdrev,zrev,wspdrev, wdirrev,dataOpts)
    draw(plot)

; ====================
; the model profile

    dataOpts_samurai 			= True
    dataOpts_samurai@PrintZ 		= True
    dataOpts@WspdWdir			= False;
    dataOpts@PlotWindH 			= False
    dataOpts_samurai@Wthin 		= 0 		;reduces number of windbarbs
    dataOpts_samurai@DrawFahrenheit   	= False    	; default is True
    dataOpts_samurai@WspdWdir		= False 	; use u,v instead of wspd, wdir
    dataOpts_samurai@HspdHdir		= False
    dataOpts_samurai@ThermoInfo       	= False    	; print thermodynamic info
    dataOpts_samurai@offsetThermoInfo 	= -0.007 	;offset for Thermoinfo default 0
    dataOpts_samurai@linePatternTemperature = 2 	; change to dashed (1=solid)
    dataOpts_samurai@linePatternDewPt	=11 		; change to dashed (1=solid)
    dataOpts_samurai@linePatternCape	=15 		; change to dashed (1=solid)
    dataOpts_samurai@xpWind    		= 31.      	; new location for winds [default 45] in celsius
    dataOpts_samurai@colThermoInfo 	= "Red"
    dataOpts_samurai@colWindP 		= "Red"
    dataOpts_samurai@colDewPt 		= "Blue"
    dataOpts_samurai@lineThicknessTemperature=2.0
    dataOpts_samurai@lineThicknessDewPt	= 2.0
    print("Plot model profile")
    plot_samurai = skewT_PlotData(wks,skewt_bkgd,p_samurai(0:dims(1)-1),t_samurai_avg(0:dims(1)-1)-273.15,td_samurai_avg(0:dims(1)-1)-273.15,fi_samurai(0:dims(1)-1),u_samurai_avg(0:dims(1)-1),v_samurai_avg(0:dims(1)-1),dataOpts_samurai)
    draw(plot_samurai)

; Plot stdevs 
    print("Plot model profile stddevs")
    dataOpts_samurai@lineThicknessTemperature=1.0
    dataOpts_samurai@lineThicknessDewPt=1.0
    dataOpts_samurai@PlotWindP=False
    dataOpts_samurai@ThermoInfo= False

    plot_samurai_s1 = skewT_PlotData(wks,skewt_bkgd,p_samurai(0:dims(1)-1),t_samurai_avg(0:dims(1)-1)-273.15+t_samurai_std(0:dims(1)-1),td_samurai_avg(0:dims(1)-1)-273.15+td_samurai_std(0:dims(1)-1),fi_samurai(0:dims(1)-1),u_samurai_avg(0:dims(1)-1),v_samurai_avg(0:dims(1)-1), dataOpts_samurai)
    draw(plot_samurai_s1)

    plot_samurai_s2 = skewT_PlotData(wks,skewt_bkgd,p_samurai(0:dims(1)-1),t_samurai_avg(0:dims(1)-1)-273.15-t_samurai_std(0:dims(1)-1),td_samurai_avg(0:dims(1)-1)-273.15-td_samurai_std(0:dims(1)-1),fi_samurai(0:dims(1)-1),u_samurai_avg(0:dims(1)-1),v_samurai_avg(0:dims(1)-1), dataOpts_samurai)
    draw(plot_samurai_s2)
    frame(wks)
    system("mogrify -density 150 -trim -format png "+plot_dir+plot_name+"_samurai."+wks_type)


end
