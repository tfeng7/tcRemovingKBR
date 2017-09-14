;**************************************************************
;IBTRACS: Storm Cyclogenisis
;  Function: remove TC signals from NCEP/DOE reanalysis
;            using KBR method (Kurihara et al., 1993,1995,MWR)
;  Created by TaoFENG@NJU, 12 Sep 2017
;**************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************************************

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

undef("removeTcGaus")
function removeTcGaus(var[*][*][*],lat,lon,centerlat,centerlon,rad)
local dlat,dlon,dis
begin
        dlat=lat-centerlanstallt
        dlon=lon-centerlon
        dlat2d=conform_dims((/dimsizes(lat),dimsizes(lon)/),dlat,0)
        dlon2d=conform_dims((/dimsizes(lat),dimsizes(lon)/),dlon,1)
        dis=sqrt(dlat2d*dlat2d+dlon2d*dlon2d) * 110
        weight=calGaus(dis,rad)
        weight3d=conform_dims(dimsizes(var),weight,(/1,2/))
        var=var*weight3d
        return(var)
end

undef("lonReorder")
function lonReorder(lon[*])
local i
begin
        reLon=where(lon.lt.0,(360+lon),lon)
        return(reLon)
end

undef("smoothKBR")
function smoothKBR(var[*][*][*],centerlat,centerlon,rad)
local k,m,numIteration,nlat,nlon
begin
        numIteration = 11
        pi = 3.1415926
        m = (/2,3,4,2,5,6,7,2,8,9,2/)
        lat=var&lat
        lon=var&lon
        ;varin=var(:,{centerlat-rad:centerlat+rad},{centerlon-rad:centerlon+rad})
        varin=var
        nlon=dimsizes(varin(0,0,:))
        nlat=dimsizes(varin(0,:,0))
        do j=0,numIteration-1
                k=1./2/(1-cos(2*pi/m(j)))
                varin(:,0:nlat-1,1:nlon-2)=varin(:,0:nlat-1,1:nlon-2) + k*(varin(:,0:nlat-1,0:nlon-3)+varin(:,0:nlat-1,2:nlon-1)-2*varin(:,0:nlat-1,1:nlon-2))
        end do
        do j=0,numIteration-1
                k=1/2/(1-cos(2*pi/m(j)))
                varin(:,1:nlat-2,0:nlon-1)=varin(:,1:nlat-2,0:nlon-1) + k*(varin(:,0:nlat-3,0:nlon-1)+varin(:,2:nlat-1,0:nlon-1)-2*varin(:,1:nlat-2,0:nlon-1))
        end do
        ;var(:,{centerlat-rad:centerlat+rad},{centerlon-rad:centerlon+rad}) = varin

        return(varin)
end

undef("creatPolorCoord")
function creatPolorCoord(centerlat,centerlon,nRadial,nAzimth,maxRadial,maxAzimth)
local  DEGTORAD
begin
        ;
        ; incrRadial:degree 0~6     0.2 deg
        ; incrAzimth:degree 0~360  15.0 deg
        ; nRadial:30
        ; nAzimth:24
        ;
        DEGTORAD = 0.017453292519943
        x        = fspan(0,maxRadial  ,nRadial) 
        y        = fspan(0,maxAzimth,nAzimth)  ; 0=360
        dimpolor = new((/nRadial,nAzimth/),"float")
        x2d      = conform(dimpolor,x,0)
        y2d      = conform(dimpolor,y,1)
        xarr     = centerlon + x2d  * cos(DEGTORAD * y2d     )
        yarr     = centerlat + x2d  * sin(DEGTORAD * y2d     )

        array    = (/yarr,xarr/)
        return(array)
end

undef("tcFirstGuess")
function tcFirstGuess(wspd[*][*],centerlat,centerlon)
local w,d,rad
begin
        rad=5.5 ; 11 * 11 degree box
        d1=200
        d2=400
        lat=wspd&lat
        lon=wspd&lon
        lat&lat=wspd&lat
        lon&lon=wspd&lon
        slat=lat({centerlat-rad:centerlat+rad})
        slon=lon({centerlon-rad:centerlon+rad})
        svar=wspd({centerlat-rad:centerlat+rad},{centerlon-rad:centerlon+rad})
        distance=new((/dimsizes(slat),dimsizes(slon)/),"float")
        slat2d=conform(svar,slat,0)
        slon2d=conform(svar,slon,1)
        dlat=slat2d-centerlat
        dlon=slon2d-centerlon
        distance=sqrt(dlat*dlat+dlon*dlon) * 110.0 ;!!! distances per degree should vary with latitude
        w2d = where(distance.le.d1,1,exp(-(distance-d1)*(distance-d1)/(d2*d2)))
        clat=sum(w2d*svar*slat2d)/sum(w2d*svar)
        clon=sum(w2d*svar*slon2d)/sum(w2d*svar)

        center=(/clat,clon/)
        return(center)
end

undef("calTanWspd")
function calTanWspd(u[*][*],v[*][*],centerlat,centerlon)
local i
begin
        pi=3.14159
        lat=u&lat
        lon=u&lon
        dlat=lat-centerlat
        dlon=lon-centerlon
        dlat2d=conform(u,dlat,0)
        dlon2d=conform(u,dlon,1)
        theta=atan2(dlat2d,dlon2d)
        wspd=sqrt(u*u+v*v)
        fcos=cos(theta)
        fsin=sin(theta)
        utan=v*fcos-u*fsin
        copy_VarMeta(u,utan)
        return(utan)
end

undef("GridtoPolar")
function GridtoPolar(varGrid[*][*],centerlat,centerlon,nRadial,nAzimth,maxRadial,maxAzimth)
local clat,clon
begin
        llarray = creatPolorCoord(centerlat,centerlon,nRadial,nAzimth,maxRadial,maxAzimth)
        ; interpolate tanWspd on Polor points
        ilat2d=llarray(0,:,:)
        ilon2d=llarray(1,:,:)
        varPolar2d =  rgrid2rcm_Wrap(varGrid&lat,varGrid&lon,varGrid, ilat2d, ilon2d, 1)
        varPolar2d@lat2d=ilat2d
        varPolar2d@lon2d=ilon2d

        ;Deprecated
        ;tanWspdPol =  linint2_points_Wrap(lon,lat,tanWspd, True, ilon1d, ilat1d, 1)
        ;tanWspdPol2d=onedtond(tanWspdPol,(/nRadial,nAzimth/)) 
        ;tanWspdPol2d@lat1d=ilat1d
        ;tanWspdPol2d@lon1d=ilon1d

        delete(llarray)
        delete(ilat2d)
        delete(ilon2d)
        return(varPolar2d)
end

undef("findTcBoundary")
function findTcBoundary(tanWspdPol2d[*][*],centerLatRel,centerLonRel,maxRadial,nRadial,rminDeg)
local clat,clon
begin
        ilat2d=tanWspdPol2d@lat2d
        ilon2d=tanWspdPol2d@lon2d
        tanWspdAveRadial  = dim_avg_n(tanWspdPol2d,1)
        radii=sqrt((ilat2d-centerLatRel)*(ilat2d-centerLatRel)+(ilon2d-centerLonRel)*(ilon2d-centerLonRel)) * 110. ;distance to center

        ; cal partial deviation of tanwind
        dvtan=center_finite_diff_n(tanWspdPol2d,radii,False,0,0) 
        dvdr=-dvtan/(maxRadial/nRadial*110.)

        print("Searching for rmax")
        rmin=floor(tofloat(rminDeg)*nRadial/maxRadial)
        rmax=maxind(tanWspdAveRadial(:)) ; find rmax from 1 degree with 0.1 degree increment
        rmax=where(rmax.gt.55,55,rmax)
        rs=max((/rmin,rmax/))
        ra=floattointeger(floor(1.5*rs))

        rf=new(24,integer)
        do i=0,dimsizes(tanWspdPol2d(0,:))-1
                rf(i)=maxind(tanWspdPol2d(:,i))
        end do

        print("Searching for grid points where match wind speed criterions")
        rfdim=new((/25/),"integer")
        do i=0,dimsizes(tanWspdPol2d(0,:))-1
                temp=ind((tanWspdPol2d(ra::,i).lt.6.and.dvdr(ra::,i).lt.4e-6).or.tanWspdPol2d(ra::,i).lt.3)
                temp=where(ismissing(temp),nRadial-1,temp)
                if(tanWspdPol2d(temp(0),i).lt.3.or.dimsizes(temp).eq.1) then
                        rfdim(i)=temp(0)+ra
                else
                        rfdim(i)=temp(1)+ra
                end if
                delete(temp)
        end do
        print("Check if the boundary grid point exeeded 1200km")
        rfdim1=where(rfdim.gt.110,110,rfdim)
        rfdim1=where(ismissing(rfdim1),110,rfdim1)

        print("Check wind consistency for extended boundary")
        rfdim2=floattointeger(floor(rfdim1*1.25))
        rfdim = rfdim2
        ; test if tanWspd became negative from rfdim to rfdim*1.25
        do i=0,dimsizes(tanWspdPol2d(0,:))-1
                temp=ind(tanWspdPol2d(rfdim1(i):rfdim2(i),i).lt.0)
                temp=where(ismissing(temp),rfdim2(i)-rfdim1(i),temp)
                rfdim(i)= temp(0) + rfdim1(i)
                delete(temp)
        end do

        plon=new((/25/),"float")
        plat=new((/25/),"float")
        do i=0,dimsizes(tanWspdPol2d(0,:))-1
                plon(i)=ilon2d(rfdim(i),i)
                plat(i)=ilat2d(rfdim(i),i)
        end do
        plon(24)=plon(0)
        plat(24)=plat(0)
        rfdim(24)=rfdim(0)

        return((/plat,plon,rfdim/))
end

undef("tcRelocating")
function tcRelocating(u[*][*],v[*][*],centerlat,centerlon,nRadial,nAzimth,maxRadial,maxAzimth,iFineGrid)
local clat,clon
begin
        wspd = sqrt(u*u+v*v)
        copy_VarMeta(u,wspd)
        lat=u&lat
        lon=u&lon
        lat&lat=wspd&lat
        lon&lon=wspd&lon

        center  = tcFirstGuess(wspd,centerlat,centerlon) ; use wspd to guess Tc center
        clat = center(0)
        clon = center(1)
        print("First guess center: "+clat+" "+ clon)
        ; calculate averaged tangential winds for each grid point inside the domain
        if(iFineGrid.eq.1) then
                print("Fine grid used.")
                gueLats = fspan(clat-3.5,clat+3.5,71) ; for locate TC center in finer grid
                gueLons = fspan(clon-3.5,clon+3.5,71)
        else
                gueLats = lat({clat-3.5:clat+3.5})
                gueLons = lon({clon-3.5:clon+3.5})
        end if

        tanWspdAve2d=new((/dimsizes(gueLats),dimsizes(gueLons)/),float)
        do i=0,dimsizes(gueLats)-1
                do j=0,dimsizes(gueLons)-1
                        print("Calculating averaged tangential winds on grid point:" + gueLats(i) + " " + gueLons(j))
                        tanWspd = calTanWspd(u,v, gueLats(i), gueLons(j))
                        tanWspdPol2d = GridtoPolar(tanWspd,gueLats(i), gueLons(j),nRadial,nAzimth,maxRadial,maxAzimth)
                        tanWspdAve2d(i,j)=avg(tanWspdPol2d)
                end do
        end do

        ; find the grid points with maximum averaged tangential winds
        a1D      = ndtooned(tanWspdAve2d)
        dsizes_a = dimsizes(tanWspdAve2d)
        indices  = ind_resolve(maxind(a1D),dsizes_a)
        clat=gueLats(indices(0,0))
        clon=gueLons(indices(0,1))
        return((/clat,clon/))

        delete(clat)
        delete(clon)
        delete(gueLats)
        delete(gueLons)
        delete(tanWspd)
        delete(llarray)
        delete(ilat1d)
        delete(ilon1d)
        delete(tanWspdPol)
        delete(tanWspdPol2d)
        delete(tanWspdAve)
        delete(tanWspdAve2d)
        delete(a1D)
        delete(dsizes_a)
        delete(indices)
end


undef("plotTcRemoving")
procedure plotTcRemoving(u[*][*],v[*][*],u_bas[*][*],v_bas[*][*],u_dis[*][*],v_dis[*][*],u_env[*][*],v_env[*][*],u_tc[*][*],v_tc[*][*],time,centerlat,centerlon,centerlatorg,centerlonorg,plat,plon)
local i 
begin
        wspd_org = sqrt(u*u+v*v)
        wspd_bas = sqrt(u_bas*u_bas+v_bas*v_bas)
        wspd_dis = sqrt(u_dis*u_dis+v_dis*v_dis)
        wspd_env = sqrt(u_env*u_env+v_env*v_env)
        wspd_tc  = sqrt(u_tc*u_tc+v_tc*v_tc)
        copy_VarMeta(u    ,wspd_org)
        copy_VarMeta(u_bas,wspd_bas)
        copy_VarMeta(u_dis,wspd_dis)
        copy_VarMeta(u_env,wspd_env)
        copy_VarMeta(u_tc ,wspd_tc )
        lat=u&lat
        lon=u&lon

        figname=tostring(time)+"_"+centerlat+"_"+centerlon
        wks   = gsn_open_wks ("eps",figname)              ; send graphics to PNG file
        ;gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
        gsn_define_colormap(wks,"WhBlGrYeRe")


        res             = True                              ; plot mods desired
        res@gsnFrame    = False                             ; don't advance yet
        res@gsnDraw     = False                             ; don't advance yet
        ;res@gsnMaximize = True

        res@mpMinLatF = min(centerlat) - 20
        res@mpMaxLatF = max(centerlat) + 20
        res@mpMinLonF = min(centerlon) - 20
        res@mpMaxLonF = max(centerlon) + 20
        ;res@tiMainString = "Locations of stations"
        res@mpOutlineBoundarySets = "Geophysical"
        res@gsnRightString = time

        cnres = True
        cnres@gsnFrame    = False                             ; don't advance yet
        cnres@gsnDraw     = False                             ; don't advance yet
        cnres@cnFillOn              = True     ; turn on color fill
        cnres@cnLinesOn             = False    ; turn of contour lines
        ;cnres@cnLevelSelectionMode  = "ExplicitLevels"
        ;cnres@cnLevels              = (/-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10/)
        cnres@cnLevelSelectionMode  = "ManualLevels"
        cnres@cnInfoLabelOn  = False
        cnres@cnMinLevelValF = 4
        cnres@cnMaxLevelValF = 20
        cnres@cnLevelSpacingF=1
        cnres@cnSmoothingOn = True
        cnres@gsnLeftString = ""
        cnres@gsnRightString = ""
        ;cnres@sfXArray    = ilon1d
        ;cnres@sfYArray    = ilat1d
        cnres@lbLabelBarOn = False

        vcres = True
        vcres@gsnFrame    = False                             ; don't advance yet
        vcres@gsnDraw     = False                             ; don't advance yet
        ;vcres@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
        vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
        vcres@vcRefLengthF            = 0.045            ; define length of vec ref
        ;vcres@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
        vcres@vcMinDistanceF          = 0.01
        vcres@vcRefAnnoPerimOn        = False
        vcres@gsnLeftString           = ""
        vcres@gsnRightString           = ""
        vcres@vcRefAnnoOrthogonalPosF = -0.1


        txres               = True
        txres@txFontHeightF = 0.050
        txres@txFontColor   = "orange"
        txres@txJust        = "CenterCenter"
        txres@txFontThicknessF = 1.0

        txres2=txres
        txres2@txFontColor        = "gray"
        txres2@txFontHeightF = 0.020

        xyRes = True
        xyRes@gsnDraw  = False
        xyRes@gsnFrame = False
        xyRes@xyMarkLineMode     = "MarkLines"
        xyRes@xyMonoDashPattern  = True
        xyRes@xyLineColor        = "orange"
        xyRes@xyLineThicknessF   = 4
        xyRes@xyMarker           = 1
        xyRes@xyMarkerSizeF      = 0.01
        xyRes@xyMarkerColor      = "black"
        xyRes@gsnLeftString           = ""
        xyRes@gsnRightString           = ""

        res@gsnLeftString = "Original"
        map1 = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
        cn1 = gsn_csm_contour(wks,wspd_org, cnres)
        vc1 = gsn_csm_vector(wks,u,v,vcres)
        cyc1= gsn_add_text( wks, map1, "~F35~m", centerlon   , centerlat   , txres )
        ll1 = gsn_add_text( wks, map1, "~F35~m", centerlonorg, centerlatorg, txres2)
        out1= gsn_csm_xy( wks, plon, plat, xyRes )

        overlay(map1,cn1)
        overlay(map1,vc1)
        overlay(map1,out1)

        res@gsnLeftString = "Envorionment"
        map2 = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
        cn2 = gsn_csm_contour(wks,wspd_env, cnres)
        vc2 = gsn_csm_vector(wks,u_env,v_env,vcres)
        cyc2= gsn_add_text( wks, map2, "~F35~m", centerlon   , centerlat   , txres )
        ll2 = gsn_add_text( wks, map2, "~F35~m", centerlonorg, centerlatorg, txres2)
        out2= gsn_csm_xy( wks, plon, plat, xyRes )

        overlay(map2,cn2)
        overlay(map2,vc2)
        overlay(map2,out2)

        res@gsnLeftString = "Tropical cyclone"
        map3 = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
        cn3 = gsn_csm_contour(wks,wspd_tc, cnres)
        vc3 = gsn_csm_vector(wks,u_tc,v_tc,vcres)
        cyc3= gsn_add_text( wks, map3, "~F35~m", centerlon   , centerlat   , txres )
        ll3 = gsn_add_text( wks, map3, "~F35~m", centerlonorg, centerlatorg, txres2)
        out3= gsn_csm_xy( wks, plon, plat, xyRes )

        overlay(map3,cn3)
        overlay(map3,vc3)
        overlay(map3,out3)

        res@gsnLeftString = "Basic"
        map4 = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
        cn4 = gsn_csm_contour(wks,wspd_bas, cnres)
        vc4 = gsn_csm_vector(wks,u_bas,v_bas,vcres)
        cyc4= gsn_add_text( wks, map4, "~F35~m", centerlon   , centerlat   , txres )
        ll4 = gsn_add_text( wks, map4, "~F35~m", centerlonorg, centerlatorg, txres2)
        out4= gsn_csm_xy( wks, plon, plat, xyRes )

        overlay(map4,cn4)
        overlay(map4,vc4)
        overlay(map4,out4)
        
        res@gsnLeftString = "Disturbance"
        map5 = gsn_csm_map(wks,res)                ; Create the map plot, don't draw
        cn5 = gsn_csm_contour(wks,wspd_dis, cnres)
        vc5 = gsn_csm_vector(wks,u_dis,v_dis,vcres)
        cyc5= gsn_add_text( wks, map5, "~F35~m", centerlon   , centerlat   , txres )
        ll5 = gsn_add_text( wks, map5, "~F35~m", centerlonorg, centerlatorg, txres2)
        out5= gsn_csm_xy( wks, plon, plat, xyRes )

        overlay(map5,cn5)
        overlay(map5,vc5)
        overlay(map5,out5)

        resP                     = True                ; modify the panel plot
        ;resP@txString            = "A plot with a common label bar"
        resP@gsnPanelLabelBar    = True                ; add common colorbar
        resP@lbLabelFontHeightF  = 0.007               ; make labels smaller
        plot=(/map1,map2,map3,map4,map5/)
        gsn_panel(wks,plot,(/2,3/),resP)               ; now draw as one plot

end

undef("removeTcKBR")
procedure removeTcKBR(u[*][*][*],v[*][*][*],centerlat,centerlon,rad,lev,u_bas[*][*][*],v_bas[*][*][*],u_dis[*][*][*],v_dis[*][*][*],u_env[*][*][*],v_env[*][*][*],u_tc[*][*][*],v_tc[*][*][*],centerLatRel,centerLonRel,plat,plon)
local i 
begin
        ; some initialization
        iRelocate=1   ; 0=relocating off, use latlon in best track dataset
        iFineGrid=0   ; 0=Find grid off, TC center is located at the original grid point
        nRadial = 30
        nAzimth = 24
        maxRadial = 6.0
        maxAzimth = 345.0
        rminDeg=3 ; degree, ignore the strong-wind region near TC center, 2 degree recommended, param added by ft
        missingvalue=u@_FillValue_original

        ; wind fields initialization
        u_org=u
        v_org=v
        u_bas = smoothKBR(u_org,centerlat,centerlon,rad)
        v_bas = smoothKBR(v_org,centerlat,centerlon,rad)
        u_dis = u-u_bas
        v_dis = v-v_bas
        copy_VarMeta(u,u_bas)
        copy_VarMeta(u,u_dis)
        copy_VarMeta(v,v_bas)
        copy_VarMeta(v,v_dis)

        u850  = u({lev},:,:)
        v850  = v({lev},:,:)


        lat=u&lat
        lon=u&lon
        lat&lat=u&lat
        lon&lon=u&lon

        if(iRelocate.eq.1) then
                print("Begin relocating using tangential wind")
                ; wind fields at 850 hPa for relocating
                u850 = u_org({lev},:,:) - u_bas({lev},:,:)
                v850 = v_org({lev},:,:) - v_bas({lev},:,:)
                relCenters=tcRelocating(u850,v850,centerlat,centerlon,nRadial,nAzimth,maxRadial,maxAzimth,iFineGrid)

                centerLatRel=relCenters(0)
                centerLonRel=relCenters(1)
                print("Relocated Tc center: "+ centerLatRel+ " " + centerLonRel)
        else
                print("Begin relocating using tangential wind")
                centerLatRel=centerlat
                centerLonRel=centerlon
        end if

        print("Calculate wind fields on grid point:"+centerLatRel+ " " + centerLonRel)
        nRadial = 180
        maxRadial = 18.0  ; 0.1 degree increment from 0 to 15 degree (approximately 1700 km)

        ; cal tang wind
        tanWspd = calTanWspd(u850,v850, centerLatRel, centerLonRel)

        ; interpolate tanWspd to Polor points
        print("Find the outermost boundary")
        tanWspdPol2d = GridtoPolar(tanWspd,centerLatRel,centerLonRel,nRadial,nAzimth,maxRadial,maxAzimth)
        ilat2d=tanWspdPol2d@lat2d
        ilon2d=tanWspdPol2d@lon2d

        ; find the outermost boundary of a TC
        llboundary=findTcBoundary(tanWspdPol2d,centerLatRel,centerLonRel,maxRadial,nRadial,rminDeg)
        plat=llboundary(0,:)
        plon=llboundary(1,:)
        pdim=tointeger(llboundary(2,:))

        ; construct none-hurricane component
        print("Set variables inside the boundary to missing value")
        llgrid=new((/73,144/),"float")

        u_none=u_dis
        v_none=v_dis
        do i=0,dimsizes(lat)-1
                do j=0,dimsizes(lon)-1
                        instate = gc_inout(lat(i), lon(j), plat, plon)
                        if(instate.eq.True) then
                                u_none(:,i,j)=missingvalue
                                v_none(:,i,j)=missingvalue
                        end if
                end do
        end do

        print("Set variables inside the boundary to 0")
        u_none(:,{centerLatRel},{centerLatRel}) = 0
        v_none(:,{centerLatRel},{centerLatRel}) = 0

        ; commented because tang wind will be not used
        ;tanWspdPolEnv=tanWspdPol2d
        ;do i=0,dimsizes(tanWspdPol2d(0,:))-1
        ;        tanWspdPolEnv(0:pdim(i),i)=tanWspdPolEnv@_FillValue
        ;end do
        ;tanWspdPolEnv(0,0)=0

        ; grid fill
        print("Construct non-hurricane component by using poission grid fill")
        guess     = 1                ; use zonal means
        is_cyclic = True             ; cyclic [global]
        nscan     = 1500             ; usually much less than this
        eps       = 1.e-2            ; variable dependent
        relc      = 0.6              ; relaxation coefficient
        opt       = 0                ; not used

        poisson_grid_fill( u_none, is_cyclic, guess, nscan, eps, relc, opt)
        poisson_grid_fill( v_none, is_cyclic, guess, nscan, eps, relc, opt)

        print("Construct hurricane and environment component")
        u_tc=u_dis-u_none
        v_tc=v_dis-v_none
        copy_VarMeta(u,u_tc)
        copy_VarMeta(v,v_tc)

        u_env=u_none+u_bas
        v_env=v_none+v_bas
        copy_VarMeta(u,u_env)
        copy_VarMeta(v,v_env)

end

;==========================================================
; Here is a sample to use tcRemoving on NCEP/DOE dataset
; remember to load this tcRemoving.ncl first
;==========================================================
; begin

;         iplot=1
;         lev=850       ; level for relocating
;         yrStart=1979
;         yrEnd=2015
;         mostart=6
;         monthEnd=10
;         rad = 12 ; KBR:12 degree; Gaus:500 km
;         var = "vwnd"

;         diri     = "./"
;         diro     = "./"
;         fili     = "Allstorms.ibtracs_all.v03r09.nc"

;         ;==========================================================
;         print ("Start year: "+ yrStart)
;         print ("End year: "+ yrEnd  )
;         f        = addfile (diri+fili, "r")
;         names = getfilevarnames(f)
;         ;print(names)

;         ; read Tc data
;         basin     = f->genesis_basin
;         tracktype = f->track_type
;         stormYear = f->season

;         ;basin2d=ndtooned(conform_dims(dimsizes(totaltime),basin,0))
;         ;track2d=ndtooned(conform_dims(dimsizes(totaltime),tracktype,0))
;         ;stmyr2d=ndtooned(conform_dims(dimsizes(totaltime),stormYear,0))

;         tcsn      = f->storm_sn
;         numobs    = f->numObs
;         time2d    = f->source_time
;         lat2d = short2flt(f->lat_for_mapping)
;         lon2d = short2flt(f->lon_for_mapping)
;         nature2d= f->nature_for_mapping

;         basin2d=conform_dims(dimsizes(time2d),basin,0)
;         track2d=conform_dims(dimsizes(time2d),tracktype,0)
;         stmyr2d=conform_dims(dimsizes(time2d),stormYear,0)

;         lat1d    = ndtooned(   lat2d)
;         lon1d    = ndtooned(   lon2d)
;         time1d   = ndtooned(  time2d)
;         nature1d = ndtooned(nature2d)
;         basin1d  = ndtooned( basin2d)
;         track1d  = ndtooned( track2d)
;         stmyr1d  = ndtooned( stmyr2d)
;         lon1d    = lonReorder(lon1d)  ; reorder Long

;         time1d@units = "days since 1858-11-17 00:00:00"
;         time1d   = cd_convert( time1d, "hours since 1800-01-01 00:00") ; convert date

;         ;==========================================================
;         datadir = "/Volumes/MyDrive/Datasets/NCEP2_daily/"
;         outdir = "/Volumes/MyDrive/Datasets/NCEP2_daily_tcremoved_kuri"

;         do yr=yrStart,yrEnd
;                 timeUnits = "hours since 1800-01-01 00:00:00"
;                 minTime = ut_inv_calendar( yr, 01, 01, 00, 0, 0, timeUnits, 0 )
;                 maxTime = ut_inv_calendar( yr, 12, 31, 18, 0, 0, timeUnits, 0 )

;                 tcind = ind((basin1d.eq.2) \           ; WNP
;                 .and.(track1d.eq.0) \           ; Not split track
;                 .and.(time1d.ge.minTime) \   ; In current year
;                 .and.(time1d.le.maxTime) \   
;                 .and.nature1d.ne.3 \              ; No disturbances
;                 .and.nature1d.ne.2)               ; No extratropical
;                 numTcs=dimsizes(tcind)
;                 print("Tc record Numbers:" + numTcs)
;                 print("maxTime = " + maxTime)
;                 print("minTime = " + minTime)
;                 print(max(time1d(tcind)))
;                 print(min(time1d(tcind)))
;                 timeyr=time1d(tcind)
;                 latyr = lat1d(tcind)
;                 lonyr = lon1d(tcind)

;                 ; read nc data
;                 filename=outdir+"/uwnd."+yr+".nc"
;                 system("cp -f "+datadir+"/uwnd."+yr+".nc"+ " "+ filename)

;                 infile    = addfile(filename,"r")
;                 flevel    = infile->level
;                 flat      = infile->lat(::-1)
;                 flon      = infile->lon
;                 ftime     = infile->time
;                 ftimebnds = infile->time_bnds
;                 uwnd     = short2flt(infile->uwnd)
;                 uwnd     = uwnd(:,:,::-1,:)

;                 filename1=outdir+"/vwnd."+yr+".nc"
;                 system("cp -f "+datadir+"/vwnd."+yr+".nc"+ " "+ filename1)

;                 infile1   = addfile(filename1,"r")
;                 flevel1   = infile1->level
;                 flat1     = infile1->lat(::-1)
;                 flon1     = infile1->lon
;                 ftime1    = infile1->time
;                 ftimebnds1= infile1->time_bnds
;                 vwnd    = short2flt(infile1->vwnd)
;                 vwnd    = vwnd (:,:,::-1,:)

;                 ; var initialization
;                 u_dis=uwnd
;                 v_dis=vwnd
;                 u_bas=u_dis
;                 v_bas=v_dis
;                 u_env=u_dis
;                 v_env=v_dis
;                 u_tc =u_dis
;                 v_tc =v_dis

;                 ftimehr   = cd_convert( ftime, "hours since 1800-01-01 00:00") ; convert date
;                 print("File:"+filename)
;                 print("Total Tcs: " + numTcs)

;                 do i=0,numTcs-1
;                         time_char = cd_calendar(timeyr(i), 0)
;                         time_hour = tointeger(time_char(:,3))
;                         if(time_hour.eq.0)  ; for ncep only
;                                 time_utc  = cd_calendar(timeyr(i), 3)
;                                 time_jul  = tointeger(timeyr(i))
;                                 centerLatRel=latyr(i)
;                                 centerLonRel=lonyr(i)  
;                                 plat=new(25,float)
;                                 plon=new(25,float)

;                                 print("Time:"+ time_utc + " " + time_jul +" LatLon: "+latyr(i)+ " "+lonyr(i)) 
;                                 ;varout({time_days(j)},:,:,:)=removeTcGaus(varin({time_days(j)},:,:,:),flat,flon,lat(i,j),lon(i,j),rad)
;                                 removeTcKBR(uwnd({time_jul},:,:,:),vwnd({time_jul},:,:,:),latyr(i),lonyr(i),rad,lev,u_bas({time_jul},:,:,:),v_bas({time_jul},:,:,:),u_dis({time_jul},:,:,:),v_dis({time_jul},:,:,:),u_env({time_jul},:,:,:),v_env({time_jul},:,:,:),u_tc({time_jul},:,:,:),v_tc({time_jul},:,:,:),centerLatRel,centerLonRel,plat,plon)

;                                 if(iplot.eq.1) then
;                                         plotTcRemoving(uwnd({time_jul},{lev},:,:),vwnd({time_jul},{lev},:,:),u_bas({time_jul},{lev},:,:),v_bas({time_jul},{lev},:,:),u_dis({time_jul},{lev},:,:),v_dis({time_jul},{lev},:,:),u_env({time_jul},{lev},:,:),v_env({time_jul},{lev},:,:),u_tc({time_jul},{lev},:,:),v_tc({time_jul},{lev},:,:),time_utc,centerLatRel,centerLonRel,latyr(i),lonyr(i),plat,plon)
;                                 end if
;                                 uwnd=u_env
;                                 vwnd=v_env
;                                 delete(time_utc)
;                                 delete(time_jul)
;                         end if  ; for ncep only
;                 end do

;                 ; output to nc file
;                 varoutshort=floattoshort(u_env(:,:,:,:))
;                 outfile = outdir+"/uwnd."+yr+".nc"
;                 system("/bin/rm -f " + outfile)   ; remove any pre-existing file
;                 ncdf     = addfile(outfile ,"c")  ; open output netCDF file

;                 filedimdef(ncdf,"time",-1,True) 
;                 ncdf->level = flevel
;                 ncdf->lat   = flat
;                 ncdf->lon   = flon
;                 ncdf->time  = ftime
;                 ncdf->time_bnds = ftimebnds
;                 ncdf->uwnd  = varoutshort

;                 varoutshort2=floattoshort(v_env(:,:,:,:))
;                 outfile1 = outdir+"/uwnd."+yr+".nc"
;                 system("/bin/rm -f " + outfile)   ; remove any pre-existing file
;                 ncdf1     = addfile(outfile1 ,"c")  ; open output netCDF file

;                 filedimdef(ncdf1,"time",-1,True)
;                 ncdf1->level = flevel
;                 ncdf1->lat   = flat
;                 ncdf1->lon   = flon
;                 ncdf1->time  = ftime
;                 ncdf1->time_bnds = ftimebnds
;                 ncdf1->vwnd   = varoutshort

;                ; delete([/tcind,tcsn,numobs,lat,lon,time,nature,filename,infile,flevel,flat,flon,ftime,ftimebnds,numTcs,varoutshort,varoutshort2,ftime1,uwnd,vwnd,latyr,lonyr,u_bas,v_bas,u_dis,v_dis,u_env,v_env,u_tc,v_tc,plat,plon/])
;                delete([/varoutshort,varoutshort2,tcind,timeyr,latyr,lonyr,ftime,ftimebnds,uwnd,ftime1,ftimebnds1,vwnd,u_dis,v_dis,u_bas,v_bas,u_env,v_env,u_tc,v_tc,ftimehr/])
;         end do
;         ystartn1=yrStart-1
;         system("rm -f "+outdir+"/"+var+"."+ystartn1+".nc")
; end
