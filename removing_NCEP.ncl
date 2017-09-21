;**************************************************************
;IBTRACS: Storm Cyclogenisis
;  Function: remove TC signals from NCEP/DOE reanalysis
;            using KBR method (Kurihara et al., 1993,1995,MWR)
;  Created by TaoFENG@NJU, 12 Sep 2017
;**************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../../tcRemovingKBR/tcRemoving.ncl"


;==========================================================
begin

        iplot=1
        lev=850       ; level for relocating
        yrStart=1979
        yrEnd=2016
        mostart=6
        monthEnd=10
        rad = 12 ; KBR:12 degree; Gaus:500 km

        diri     = "./"
        diro     = "./"
        fili     = "/Users/fengtao/Working/Data_proc/ibt/Allstorms.ibtracs_all.v03r10.nc"

        ;==========================================================
        print ("Start year: "+ yrStart)
        print ("End year: "+ yrEnd  )
        f        = addfile (fili, "r")
        names = getfilevarnames(f)
        ;print(names)

        ; read Tc data
        basin     = f->genesis_basin
        tracktype = f->track_type
        stormYear = f->season

        ;basin2d=ndtooned(conform_dims(dimsizes(totaltime),basin,0))
        ;track2d=ndtooned(conform_dims(dimsizes(totaltime),tracktype,0))
        ;stmyr2d=ndtooned(conform_dims(dimsizes(totaltime),stormYear,0))

        tcsn      = f->storm_sn
        numobs    = f->numObs
        time2d    = f->source_time
        lat2d = short2flt(f->lat_for_mapping)
        lon2d = short2flt(f->lon_for_mapping)
        nature2d= f->nature_for_mapping

        basin2d=conform_dims(dimsizes(time2d),basin,0)
        track2d=conform_dims(dimsizes(time2d),tracktype,0)
        stmyr2d=conform_dims(dimsizes(time2d),stormYear,0)

        lat1d    = ndtooned(   lat2d)
        lon1d    = ndtooned(   lon2d)
        time1d   = ndtooned(  time2d)
        nature1d = ndtooned(nature2d)
        basin1d  = ndtooned( basin2d)
        track1d  = ndtooned( track2d)
        stmyr1d  = ndtooned( stmyr2d)
        lon1d    = lonReorder(lon1d)  ; reorder Long

        time1d@units = "days since 1858-11-17 00:00:00"
        time1d   = cd_convert( time1d, "hours since 1800-01-01 00:00") ; convert date

        ;==========================================================
        datadir = "/Volumes/MyDrive/Datasets/NCEP2_daily/"
        outdir = "/Volumes/MyDrive/Datasets/NCEP2_daily_tcremoved_KBR"

        do yr=yrStart,yrEnd
                timeUnits = "hours since 1800-01-01 00:00:00"
                minTime = ut_inv_calendar( yr, 01, 01, 00, 0, 0, timeUnits, 0 )
                maxTime = ut_inv_calendar( yr, 12, 31, 18, 0, 0, timeUnits, 0 )

                tcind = ind((basin1d.eq.2) \           ; WNP
                .and.(track1d.eq.0) \           ; Not split track
                .and.(time1d.ge.minTime) \   ; In current year
                .and.(time1d.le.maxTime) \   
                .and.nature1d.ne.3 \              ; No disturbances
                .and.nature1d.ne.2)               ; No extratropical
                numTcs=dimsizes(tcind)
                print("============================================================")
                print("Tc record Numbers:" + numTcs)
                print("maxTime = " + maxTime)
                print("minTime = " + minTime)
                print(max(time1d(tcind)))
                print(min(time1d(tcind)))
                timeyr=time1d(tcind)
                latyr = lat1d(tcind)
                lonyr = lon1d(tcind)

                ; read nc data
                filename=outdir+"/uwnd."+yr+".nc"
                system("cp -f "+datadir+"/uwnd."+yr+".nc"+ " "+ filename)

                infile    = addfile(filename,"r")
                flevel    = infile->level
                flat      = infile->lat(::-1)
                flon      = infile->lon
                ftime     = infile->time
                ftimebnds = infile->time_bnds
                uwnd     = short2flt(infile->uwnd)
                uwnd     = uwnd(:,:,::-1,:)

                filename1=outdir+"/vwnd."+yr+".nc"
                system("cp -f "+datadir+"/vwnd."+yr+".nc"+ " "+ filename1)

                infile1   = addfile(filename1,"r")
                flevel1   = infile1->level
                flat1     = infile1->lat(::-1)
                flon1     = infile1->lon
                ftime1    = infile1->time
                ftimebnds1= infile1->time_bnds
                vwnd    = short2flt(infile1->vwnd)
                vwnd    = vwnd (:,:,::-1,:)
                tempu=infile->uwnd
                tempv=infile1->vwnd

                ; var initialization
                u_dis=uwnd
                v_dis=vwnd
                u_bas=u_dis
                v_bas=v_dis
                u_env=u_dis
                v_env=v_dis
                u_tc =u_dis
                v_tc =v_dis

                ftimehr   = cd_convert( ftime, "hours since 1800-01-01 00:00") ; convert date
                print("File:"+filename)
                print("Total Tcs: " + numTcs)
                uwndin = u_env
                vwndin = v_env

                do i=0,numTcs-1
                        time_char = cd_calendar(timeyr(i), 0)
                        time_hour = tointeger(time_char(:,3))
                        if(time_hour.eq.0)  ; for ncep only
                                time_utc  = cd_calendar(timeyr(i), 3)
                                time_jul  = tointeger(timeyr(i))
                                centerLatRel=latyr(i)
                                centerLonRel=lonyr(i)  
                                plat=new(25,float)
                                plon=new(25,float)

                                print("Time:"+ time_utc + " " + time_jul +" LatLon: "+latyr(i)+ " "+lonyr(i)) 
                                ;varout({time_days(j)},:,:,:)=removeTcGaus(varin({time_days(j)},:,:,:),flat,flon,lat(i,j),lon(i,j),rad)
                                removeTcKBR(uwndin({time_jul},:,:,:),vwndin({time_jul},:,:,:),latyr(i),lonyr(i),rad,lev,u_bas({time_jul},:,:,:),v_bas({time_jul},:,:,:),u_dis({time_jul},:,:,:),v_dis({time_jul},:,:,:),u_env({time_jul},:,:,:),v_env({time_jul},:,:,:),u_tc({time_jul},:,:,:),v_tc({time_jul},:,:,:),centerLatRel,centerLonRel,plat,plon)
                                uwndin({time_jul},:,:,:) = u_env({time_jul},:,:,:)
                                vwndin({time_jul},:,:,:) = v_env({time_jul},:,:,:)

                                if(iplot.eq.1) then
                                        plotTcRemoving(uwnd({time_jul},{lev},:,:),vwnd({time_jul},{lev},:,:),u_bas({time_jul},{lev},:,:),v_bas({time_jul},{lev},:,:),u_dis({time_jul},{lev},:,:),v_dis({time_jul},{lev},:,:),u_env({time_jul},{lev},:,:),v_env({time_jul},{lev},:,:),u_tc({time_jul},{lev},:,:),v_tc({time_jul},{lev},:,:),time_utc,centerLatRel,centerLonRel,latyr(i),lonyr(i),plat,plon)
                                end if

                                dumpInfo(time_utc,latyr(i),lonyr(i),plat,plon)
                                
                                delete(time_utc)
                                delete(time_jul)
                        end if  ; for ncep only
                end do

                ; output to nc file
                copy_VarCoords(uwnd,u_env)
                copy_VarCoords(vwnd,v_env)
                u_env = u_env(:,:,::-1,:)
                v_env = v_env(:,:,::-1,:)
                flat = flat(::-1)
                
                opt=True
                opt@scale_factor = tempu@scale_factor
                opt@add_offset   = tempu@add_offset
                varoutshort=pack_values(u_env,"short",opt)
                copy_VarMeta(tempu,varoutshort)
                printVarSummary(varoutshort)
                outfile = outdir+"/uwnd."+yr+".nc"
                system("/bin/rm -f " + outfile)   ; remove any pre-existing file
                ncdf     = addfile(outfile ,"c")  ; open output netCDF file

                filedimdef(ncdf,"time",-1,True) 
                ncdf->level = flevel
                ncdf->lat   = flat
                ncdf->lon   = flon
                ncdf->time  = ftime
                ncdf->time_bnds = ftimebnds
                ncdf->uwnd  = varoutshort

                varoutshort2=pack_values(v_env,"short",opt)
                copy_VarMeta(tempv,varoutshort2)
                outfile1 = outdir+"/vwnd."+yr+".nc"
                system("/bin/rm -f " + outfile1)   ; remove any pre-existing file
                ncdf1     = addfile(outfile1 ,"c")  ; open output netCDF file

                filedimdef(ncdf1,"time",-1,True)
                ncdf1->level = flevel
                ncdf1->lat   = flat
                ncdf1->lon   = flon
                ncdf1->time  = ftime
                ncdf1->time_bnds = ftimebnds
                ncdf1->vwnd   = varoutshort2



               ; delete([/tcind,tcsn,numobs,lat,lon,time,nature,filename,infile,flevel,flat,flon,ftime,ftimebnds,numTcs,varoutshort,varoutshort2,ftime1,uwnd,vwnd,latyr,lonyr,u_bas,v_bas,u_dis,v_dis,u_env,v_env,u_tc,v_tc,plat,plon/])
               delete([/uwndin,vwndin,tempu,tempv,varoutshort,varoutshort2,tcind,timeyr,latyr,lonyr,ftime,ftimebnds,uwnd,ftime1,ftimebnds1,vwnd,u_dis,v_dis,u_bas,v_bas,u_env,v_env,u_tc,v_tc,ftimehr/])
        end do
end
