;**********************************************************
;IBTRACS: Storm Cyclogenisis 
;**********************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**********************************************************

undef("calGaus")
function calGaus(dis[*][*],rad)
local weight,a
begin
        c=0.5
        a=rad*sqrt(2*log(2))
        ;weight=1-exp(-(dis*dis)/(2*a*a))
        weight=1-exp(-(dis*dis)/(2*a*a))

        return(weight)
end

undef("removeTC")
function removeTC(var[*][*][*],lat,lon,centerlat,centerlon,rad)
local dlat,dlon,dis
begin
        dlat=lat-centerlat
        dlon=lon-centerlon
        dlat2d=conform_dims((/dimsizes(lat),dimsizes(lon)/),dlat,0)
        dlon2d=conform_dims((/dimsizes(lat),dimsizes(lon)/),dlon,1)
        dis=sqrt(dlat2d*dlat2d+dlon2d*dlon2d) * 110
        weight=calGaus(dis,rad)
        weight3d=conform_dims(dimsizes(var),weight,(/1,2/))
        var=var*weight3d
        return(var)
end


;==========================================================
begin

        yrStart=1980
        yrEnd=2016
        mostart=6
        monthEnd=10
        rad = 500
        var = "vwnd"

        diri     = "./"
        diro     = "./"
        fili     = "Allstorms.ibtracs_all.v03r09.nc"

        ;==========================================================
        print ("Start year: "+ yrStart)
        print ("End year: "+ yrEnd  )
        f        = addfile (diri+fili, "r")
        names = getfilevarnames(f)
        ;print(names)

        basin     = f->genesis_basin
        tracktype = f->track_type
        stormYear = f->season
        totaltime = f->source_time

        ;==========================================================

        ;==========================================================
        datadir = "/Volumes/MyDrive/Datasets/NCEP2_daily/"
        outdir = "/Volumes/MyDrive/Datasets/NCEP2_daily_tcremoved"
        do yr=yrStart,yrEnd
                ; read TC data
                tcind = ind((basin.eq.2).and.(tracktype.eq.0).and.(stormYear.eq.yr))
                tcsn = f->storm_sn(tcind,:)
                numobs=f->numObs(tcind)
                lat = short2flt(f->lat_for_mapping(tcind,:))
                lon = short2flt(f->lon_for_mapping(tcind,:))
                time= f->source_time(tcind,:)
                nature= f->nature_for_mapping(tcind,:)

                ; read nc data
                yr1=yr-1
                yr2=yr+1
                filename=outdir+"/"+var+"."+yr+".nc"
                filename2=outdir+"/"+var+"."+yr1+".nc"
                filename3=datadir+"/"+var+"."+yr2+".nc"
                ; read two years data for those TCs across the years
                if(yr.eq.yrStart) then
                        system("cp "+datadir+"/"+var+"."+yr1+".nc"+ " "+ filename2)
                        system("cp "+datadir+"/"+var+"."+yr+".nc"+ " "+ filename)
                end if
                infile=addfiles((/filename2,filename,filename3/),"r")
                flevel=infile[1]->level
                flat=infile[1]->lat
                flon=infile[1]->lon
                ftime=infile[1]->time
                ftime1=infile[0]->time
                ftime3=infile[2]->time
                print("lie")
                ftime4=infile[:]->time
                print("lie")
                ftimebnds=infile[1]->time_bnds
                ftimebnds2=infile[2]->time_bnds
                varin=short2flt(infile[:]->$var$)
                print("lie1;")

                print(varin&time)
                ftime5= cd_convert( ftime4, "hours since 1800-01-01 00:00") ; convert date
                ;varin&time(0:dimsizes(ftime2)-1)=ftime2
                varin&time=ftime5
                printVarSummary(varin&time)
                numtcs  = dimsizes(lat(:,0))
                print("File:"+filename)
                print("Total TCs: " + numtcs)

                varout=varin
                do i=0,numtcs-1
                        time_days = cd_convert( time(i,0:numobs(i)-1), "hours since 1800-01-01 00:00")
                        ;print("Time_day:"+time_days)
                        do j=0,dimsizes(time_days)-1
                                varout({time_days(j)},:,:,:)=removeTC(varin({time_days(j)},:,:,:),flat,flon,lat(i,j),lon(i,j),rad)
                        end do
                        delete(time_days)
                end do

                ; output to nc file
                varoutshort=floattoshort(varout({ftime},:,:,:))
                copy_VarMeta(varin({ftime},:,:,:),varoutshort)
                varoutshort&time=ftime
                outfile = outdir+"/"+var+"."+yr+".nc"
                system("/bin/rm -f " + outfile)   ; remove any pre-existing file
                ncdf     = addfile(outfile ,"c")  ; open output netCDF file

                filedimdef(ncdf,"time",-1,True) 
                ncdf->level = flevel
                ncdf->lat   = flat
                ncdf->lon   = flon
                ncdf->time  = ftime
                ncdf->time_bnds = ftimebnds
                ncdf->$var$  = varoutshort

                varoutshort2=floattoshort(varout({ftime3},:,:,:))
                copy_VarMeta(varin({ftime3},:,:,:),varoutshort2)
                varoutshort2&time=ftime3
                outfile = outdir+"/"+var+"."+yr2+".nc"
                system("/bin/rm -f " + outfile)   ; remove any pre-existing file
                ncdf     = addfile(outfile ,"c")  ; open output netCDF file
                
                filedimdef(ncdf,"time",-1,True)
                ncdf->level = flevel
                ncdf->lat   = flat
                ncdf->lon   = flon
                ncdf->time  = ftime3
                ncdf->time_bnds = ftimebnds2
                ncdf->$var$  = varoutshort2

                delete([/tcind,tcsn,numobs,lat,lon,time,nature,filename,infile, \
                flevel,flat,flon,ftime,ftimebnds,varin,varout,numtcs,varoutshort,\
                       ftime3,varoutshort2,ftimebnds2,ftime1,ftime4,ftime5/])
        end do
        ystartn1=yrStart-1
        system("rm -f "+outdir+"/"+var+"."+ystartn1+".nc")
end
