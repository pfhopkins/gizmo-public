;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Load a snapshot into IDL                             ;;
;; -takes as inputes snapshot directory and number      ;;
;;  (routine originally written by TJ Cox, modified     ;;
;;   many times since then by Phil Hopkins)             ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;  use as follows:                                     ;;
;;                                                      ;;
;;      ok= fload_snapshot_bh(frun,0,[....])            ;;
;;                                                      ;;
;;      frun= full path to data directory.              ;;
;;            it assumes that all snapshots are         ;;
;;            called snapshot_###                       ;;
;;                                                      ;;
;;      ##= snapshot number that you want to load       ;;
;;          (don't need preceeding 0's)                 ;;
;;                                                      ;;
;;      ....= options, see below, but likely need       ;;
;;            to ask Phil about these.                  ;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;                                                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro check_field, var, varlabel

	print, "checking: "+varlabel, "   ",size(var)
	idx=where(finite(var) eq 0, n_i)
	if idx(0) ne -1 then begin
	    print, varlabel+' has '+strcompress(string(n_i),/remove_all)+' bad number(s)'
	    ;if n_i lt 10 then begin
	    ;	print, idx
	    ;	var(idx)= 1.0e-9
	    ;endif
	endif

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function readsnap, frun, num, $
			h0=h0, $
			ics=ics, $ 						;; the file is a gadget IC file
			arepo=arepo, $					;; the file is AREPO-format
			nopot_in_snap=nopot_in_snap, $	;; no potential (downgraded, since default is nopot)
			havepot_in_snap=havepot_in_snap, $ ;; force read of potential
			header_only=header_only, $		;; only read header info
			show_header=show_header, $		;; show header info
			direct_name_used=direct_name_used, $
			find_center=find_center, $
			center_on_bh=center_on_bh, $
			do_four=do_four, $				;; 
			skip_extra_gas_info=skip_extra_gas_info, $
			skip_bh=skip_bh,$
			SET_BASENAME=SET_BASENAME,$
			cosmological=cosmological

    COMMON GalaxyHeader, N,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal, $
		flag_cooling,numfiles,cosmocrap, $
		flag_multiphase, $
		flag_stellarage, $
		flag_snaphaspot,$
		flag_metals, la, $
		flag_stargens, flag_energydetail, flag_parentid, flag_starorig
    COMMON HaloData, xhalo,yhalo,zhalo,vxhalo,vyhalo,vzhalo, mhalo,dmid
    COMMON DiskData, xdisk,ydisk,zdisk,vxdisk,vydisk,vzdisk, mdisk,did
    COMMON OtherData, id, mass
    COMMON GasData, xgas,ygas,zgas,vxgas,vygas,vzgas,u,rho,volume,hsml,mgas,gid
    COMMON OldStarData, diskage, diskmetals, bulgeage, bulgemetals
    COMMON SfrData, mfs, sfr, stellage, gmetals, smetals
    COMMON MultiphaseData, mclouds
    COMMON CoolingData, nume, numh
    COMMON FeedbackData, tpu
    COMMON BulgeData, xbulge,ybulge,zbulge,vxbulge,vybulge,vzbulge,mbulge,bib
    COMMON NewStarData, xstars,ystars,zstars,vxstars,vystars,vzstars,mstars,nsid
    COMMON FileInfo, rundir, runnum, fname, exts
    COMMON PotData, pgas, phalo, pdisk, pbulge, pstars, pbh
    COMMON EnergyDetail, totradiated, totshocked, totfeedback
    COMMON ParentID, parentid
    COMMON StarOrig, origmass, orighsml
    COMMON Center, com, alternate_com, b_com
    COMMON BlackHoleData, xbh, ybh, zbh, vxbh, vybh, vzbh, mbh, bhmass, bhaccrate, bhid
;    COMMON Galaxy1, startid1, numpart1, gas_startid1, gas_numpart1
;    COMMON Galaxy2, startid2, numpart2, gas_startid2, gas_numpart2
;    COMMON SnapInfo, frun

	;; 
	;;
	;; All in GADGET units
	;;   length = h^-1 kpc = 1.43 kpc
	;;	 time = h^-1 Gyr = 1.43 Gyr
	;;   mass = 10^10 h^-1 msun = 1.43x10^10 msun
	;;   velocity = km/s
	;;   
	;; (we always use h=0.7)
	;;

    if not keyword_set(num) then num=0
    if not keyword_set(frun) then begin
	print, "  "
	print, "fload_snapshot, frun, num"
	return, -1
    endif

    exts='0000'
    exts=exts+strcompress(string(num),/remove_all)

    if num ge 1000 then do_four= 1

    ; does the snapshot have four, or three,
    ; numbers after it
    ; ---------------------------------------
    do_four_k=0
	 if keyword_set(do_four) OR num ge 1000 then do_four_k=1

    if do_four_k eq 0 then begin
	exts=strmid(exts,strlen(exts)-3,3)
    endif else begin
	exts=strmid(exts,strlen(exts)-4,4)
		print, "******************************"
		print, "******************************"
		print, "   fload_snapshot_bh "
		print, "          is  "
		print, "  momentarily changed"
		print, "******************************"
		print, "******************************"
    endelse

    if keyword_set(ics) then begin
	fname= frun
	if keyword_set(direct_name_used) then fname = direct_name_used
	result= frun
	if strlen(fname) gt 0 and strmid(fname,strlen(fname)-4,4) eq "hdf5" then goto, file_is_hdf5
	goto, foundfile
    endif

    basename='snap'
    basename='snapshot'
    if keyword_set(SET_BASENAME) then basename=SET_BASENAME

	MULTI_PART_FILE_KEY = 0

    ; HDF5 snapshot file?
    cmd= "/bin/ls "+frun+"/"+basename+"*_"+exts+".hdf5"
    spawn, cmd, result   &    print, "result= ", result
    if strlen(result) gt 0 and strmid(result,strlen(result)-4,4) eq "hdf5" then goto, file_is_hdf5
	;; also check for a multi-part hdf5 filename :: 
    cmd= "/bin/ls "+frun+"/"+basename+"*_"+exts+'.0'+".hdf5"
    spawn, cmd, result   &    print, "result= ", result
    if strlen(result) gt 0 and strmid(result,strlen(result)-4,4) eq "hdf5" then goto, file_is_hdf5


    ; regular snapshot file
    cmd= "/bin/ls "+frun+"/"+basename+"*_"+exts  
        if keyword_set(direct_name_used) then cmd= "/bin/ls "+frun+'/'+direct_name_used
    print, "spawning cmd= ", cmd
    spawn, cmd, result   &    print, "result= ", result
    if strlen(result) eq 0 then return, -1

    fname=strcompress(result,/remove_all)
    rundir=strcompress(frun,/remove_all)
    runnum=num



;=========================================================================================
;=========================================================================================


foundfile:

    if keyword_set(h0) then hubble= 0.7  
    ;if keyword_set(h0) then hubble=fload_cosmology('h')
    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedbacktp=0L
    npartTotal=lonarr(6)	
    flag_cooling=0L
    numfiles=0L
    cosmocrap=dblarr(4)
    flag_stellarage=0L
    flag_metals=0L
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4
    la=intarr(bytesleft/2)


;    if file_size(fname) le 0 then return, -1


    openr,1,fname,/f77_unformatted,ERROR=err,/SWAP_IF_BIG_ENDIAN
    if (err NE 0) then begin
	print, "  "
	print, !ERR_STRING
	print, "  "
	close, 1
	return, -1
    endif else begin
	print, "opening: ",fname
    endelse

    catch, error_status
    if error_status ne 0 then begin
	print, " "
	print, "Error: ", error_status
	print, "Error message: ", !ERROR_STATE.MSG
	close, 1
	catch, /cancel
	return, -1
    endif

    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedbacktp,npartTotal,flag_cooling,numfiles,cosmocrap,flag_stellarage,flag_metals,la

	if flag_metals gt 1 then begin
	  print, " Reading file in multiple metal species mode: N_metals = ",flag_metals
	endif

    if keyword_set(show_header) then begin
	print, "** HEADER **"
	print, "npart= ", npart
	print, "massarr= ", massarr
	print, "time= ", time
	print, "flag_sfr", flag_sfr
	print, "flag_feedbacktp", flag_feedbacktp
	print, "npartTotal", npartTotal
	print, "flag_cooling", flag_cooling
	print, "numfiles", numfiles
	print, "cosmocrap", cosmocrap
	print, "flag_stellarage", flag_stellarage
	print, "flag_metals", flag_metals
	print, "strlen(la)", strlen(la)
	print, "** ******* **"
    endif

    if time lt 0.0 or time gt 1.0e10 then begin
	print, " "
	print, " PROBLEM with time read from header, returning .... "
	print, " "
	return, -1
    endif

    ; volker isn't using a flag, so do it manually
    flag_snaphaspot= 0L
    if keyword_set(nopot_in_snap) then flag_snaphaspot= 0L 
    if keyword_set(havepot_in_snap) then flag_snaphaspot= 1L
    ;flag_snaphaspot= 0L
    flag_stargens= 1L
    flag_energydetail= 0L
    flag_parentid= 0L
    flag_starorig= 0L
    

    if keyword_set(header_only) then begin
	close, 1
	return, 1
    endif
   
    N=long(total(npart,/double))    ; could also do   total(npart,/integer)
    ;N=long(npart(0)+npart(1)+npart(2)+npart(3)+npart(4)+npart(5))
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
	Nwithmass= long(total(npart(ind),/double))
        mass=fltarr(Nwithmass)
        print, "Nwithmasses= ",Nwithmass
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    check_field, pos, 'pos'
    readu,1,vel
    check_field, vel, 'vel'
    readu,1,id
    check_field, id, 'id'
    if Nwithmass gt 0 then begin
      readu,1,mass
      check_field, mass, 'mass'
      ;print, "Nwithmasses= ",Nwithmass
    endif

    Ngas=long(npart(0))
    Nhalo=long(npart(1))
    Ndisk=long(npart(2))
    Nbulge=long(npart(3))
    Nstars=long(npart(4))
    Nbh=long(npart(5))

    N_baryons= Ngas + Ndisk + Nbulge + Nstars
    N_tot= total(npart)
    print, "Ntot= ", N_tot
    print, npart

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u
	check_field, u, 'u'
    endif

    if Ngas gt 0 then begin
        rho=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,rho
        check_field, rho, 'rho'
    endif

    if (Ngas gt 0) and keyword_set(arepo) then begin
        volume=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,volume
        check_field, volume, 'volume'
    endif

if keyword_set(skip_extra_gas_info) then goto, moveon

    if flag_cooling gt 0 and Ngas gt 0 then begin
	nume=fltarr(Ngas)
	numh=fltarr(Ngas)
	if not keyword_set(ics) then readu,1,nume
	check_field, nume, 'nume'
	if not keyword_set(ics) then readu,1,numh
	check_field, numh, 'numh'
    endif

    if Ngas gt 0 then begin
	hsml=fltarr(Ngas)
	if not keyword_set(ics) then readu,1,hsml
	check_field, hsml, 'hsml'
    endif

    if flag_sfr gt 0 and Ngas gt 0 then begin
        sfr=fltarr(Ngas)
        if not keyword_set(ics) then readu,1,sfr
	check_field, sfr, 'sfr'
	print, "sfr= ", total(sfr)
    endif

	if flag_sfr gt 0 then begin
	    if flag_stellarage gt 0 and Nstars gt 0 then begin
		    stellage= fltarr(Nstars)
		    if not keyword_set(ics) then readu,1,stellage
			check_field, stellage, 'stellage'
	    endif
	endif

    if eof(1) eq 1 then goto, moveon
    goto, moveon

    if flag_metals gt 0 then begin
	if (Ngas+Nstars) gt 0 then begin
                mets= fltarr(flag_metals,Ngas+Nstars)
		if not keyword_set(ics) then readu,1,mets
		if flag_metals gt 1 then begin
 	               dim_metal_tmp = size(mets,/dimensions)
 	               if (dim_metal_tmp(0) ne Ngas+Nstars) then mets=transpose(mets)
		endif
		check_field, mets, 'mets'
		gmetals= fltarr(Ngas,flag_metals)
		if Nstars gt 0 then smetals= fltarr(Nstars,flag_metals)
		if flag_metals gt 1 then begin
		  gmetals(*,*)= mets(0:Ngas-1,*)
		  if Nstars gt 0 then smetals(*,*)= mets(Ngas:Ngas+Nstars-1,*)
		endif else begin
                  gmetals(*)= mets(0:Ngas-1)
                  if Nstars gt 0 then smetals(*)= mets(Ngas:Ngas+Nstars-1)
		endelse
	endif
    endif

    if keyword_set(nopot_in_snap) OR not keyword_set(havepot_in_snap) then goto, moveon

    if flag_snaphaspot gt 0 then begin
        pot= fltarr(N)
        if not keyword_set(ics) then readu,1,pot
	check_field, pot, 'pot'
    endif

    if eof(1) eq 1 then goto, moveon

    if Nbh gt 0 and not keyword_set(skip_bh) then begin
        bhmass= fltarr(Nbh)
        readu, 1, bhmass
        check_field, bhmass, 'bhmass'
        bhaccrate= fltarr(Nbh)
        readu, 1, bhaccrate
        check_field, bhaccrate, 'bhaccrate'

	print, "bhmass= ", bhmass
	print, "bhaccrate= ", bhaccrate
    endif


moveon:
    close,1


    if keyword_set(h0) then begin
	;hubble= fload_cosmology('h0')
	hubble= 0.7
	print, "snapshot values corrected for h= ", hubble
	
    if keyword_set(cosmological) then begin
      ascale=time
      hubble=cosmocrap[3]
      hinv=1.d0/hubble
      print,' cosmological sim: a=',ascale,' h=',hubble
    endif
	if not keyword_set(cosmological) then time= time/hubble
        pos= pos/hubble
        mass= mass/hubble
		massarr= massarr/hubble
        rho=rho*hubble*hubble
        if keyword_set(arepo) then volume=volume/hubble/hubble/hubble
        ;mfs=mfs/hubble
        ;if flag_multiphase gt 0 then mclouds=mclouds/hubble
		hsml=hsml/hubble
	if not keyword_set(cosmological) then begin
	if flag_stellarage gt 0 and Nstars gt 0 then stellage= stellage/hubble
	endif
	;if flag_metals gt 0 then gmetals=gmetals/hubble
	;if flag_metals gt 0 and Nstars gt 0 then smetals=smetals/hubble
	;if flag_starorig gt 0 and Nstars gt 0 then origmass=origmass/hubble
	;if flag_starorig gt 0 and Nstars gt 0 then orighsml=orighsml/hubble

	;
	; potential and velocity should be fine

	if keyword_set(cosmological) then begin
	  ascale=time & hsml=hsml*ascale & pos=pos*ascale & rho=rho/(ascale*ascale*ascale) 
	endif
    endif

    ; ----------------------------------
    ;   now start parsing file
    ; ----------------------------------
    if Ngas gt 0 then begin
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)

        vxgas=fltarr(Ngas) &  vygas=fltarr(Ngas)  & vzgas=fltarr(Ngas)
        vxgas(*)=vel(0,0:Ngas-1)
        vygas(*)=vel(1,0:Ngas-1)
        vzgas(*)=vel(2,0:Ngas-1)

        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse

	if flag_snaphaspot gt 0 then begin
	    pgas=fltarr(Ngas)
	    pgas(*)=pot(0:Ngas-1)
	endif
    endif

    if Nhalo gt 0 then begin
        xhalo=fltarr(Nhalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo)
        xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)

        vxhalo=fltarr(Nhalo) &  vyhalo=fltarr(Nhalo) & vzhalo=fltarr(Nhalo)
        vxhalo(*)=vel(0,0+Ngas:Nhalo+Ngas-1)
        vyhalo(*)=vel(1,0+Ngas:Nhalo+Ngas-1)
        vzhalo(*)=vel(2,0+Ngas:Nhalo+Ngas-1)

        if massarr(1) eq 0 then begin
	    skip=0L
            for t=0,0 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mhalo(*)=mass(0+skip:Nhalo-1+skip)	
        endif else begin
            mhalo(*)= massarr(1)
	endelse

        if flag_snaphaspot gt 0 then begin
            phalo=fltarr(Nhalo)
            phalo(*)=pot(0+Ngas:Nhalo+Ngas-1)
        endif
    endif


    ; Disk
    ; ------
    mdisk= 0
    if Ndisk gt 0 then begin
        xdisk=fltarr(Ndisk) &  ydisk=fltarr(Ndisk) &  zdisk=fltarr(Ndisk) & mdisk=fltarr(Ndisk)
        xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)

        vxdisk=fltarr(Ndisk) &  vydisk=fltarr(Ndisk) &  vzdisk=fltarr(Ndisk)
        vxdisk(*)=vel(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vydisk(*)=vel(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vzdisk(*)=vel(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)

        if massarr(2) eq 0 then begin
	    skip=0L
            for t=0,1 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mdisk(*)=mass(0+skip:Ndisk-1+skip)	
        endif else begin
            mdisk(*)= massarr(2)
	endelse

        if flag_snaphaspot gt 0 then begin
            pdisk=fltarr(Ndisk)
            pdisk(*)=pot(Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        endif

	diskage= -5.0*randomu(seed,Ndisk)
	diskage= -1.0*randomu(seed,Ndisk)    ; changed this for Marijn
	if(flag_metals gt 0) then begin
		diskmetals= 10^(1.5*randomu(seed,Ndisk,flag_metals) - 2.0)*0.02
	endif

    endif



    ; Bulge
    ; -------
    mbulge= 0
    if Nbulge gt 0 then begin
        xbulge=fltarr(Nbulge) &  ybulge=fltarr(Nbulge) &  zbulge=fltarr(Nbulge) & mbulge=fltarr(Nbulge)
        xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)

        vxbulge=fltarr(Nbulge) &  vybulge=fltarr(Nbulge) &  vzbulge=fltarr(Nbulge)
        vxbulge(*)=vel(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vybulge(*)=vel(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vzbulge(*)=vel(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)

        if massarr(3) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mbulge(*)=mass(0+skip:Nbulge-1+skip)	
        endif else begin
            mbulge(*)= massarr(3)
	endelse

        if flag_snaphaspot gt 0 then begin
            pbulge=fltarr(Nbulge)
            pbulge(*)=pot(Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        endif

	bulgeage= -7.0 + 0.0*mbulge
	if(flag_metals gt 0) then begin
		bulgemetals= 0.001 + fltarr(Nbulge,flag_metals)
	endif
    endif


    if Nstars gt 0 then begin
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)

        vxstars=fltarr(Nstars) &  vystars=fltarr(Nstars) &  vzstars=fltarr(Nstars)
        vxstars(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        vystars(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        vzstars(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)

        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,3 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            mstars(*)= massarr(4)
	endelse

        if flag_snaphaspot gt 0 then begin
            pstars=fltarr(Nstars)
            pstars(*)=pot(Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        endif
    endif


    if Nbh gt 0 then begin
        xbh=fltarr(Nbh) &  ybh=fltarr(Nbh)  & zbh=fltarr(Nbh)  & mbh=fltarr(Nbh)
        xbh(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        ybh(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        zbh(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)

        vxbh=fltarr(Nbh) &  vybh=fltarr(Nbh) &  vzbh=fltarr(Nbh)
        vxbh(*)=vel(0,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        vybh(*)=vel(1,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        vzbh(*)=vel(2,Nhalo+Ndisk+Ngas+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)

        if massarr(5) eq 0 then begin
            skip=0L
            for t=0,4 do begin
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                        skip=skip + npart(t)
                endif
            endfor
            mbh(*)=mass(0+skip:Nbh-1+skip)
        endif else begin
            mbh(*)= massarr(5)
        endelse

        if flag_snaphaspot gt 0 then begin
            pbh=fltarr(Nbh)
            pbh(*)=pot(Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+Nbh-1)
        endif
    endif


    goto, calccenter




;=========================================================================================
;=========================================================================================

    file_is_hdf5:

    fname=strcompress(result,/remove_all)
	if keyword_set(direct_name_used) then fname = direct_name_used
    rundir=strcompress(frun,/remove_all)
    runnum=num

    if h5f_is_hdf5(fname) eq 0 then begin
	print, " "
	print, "Problem, not (or bad) HDF5 file: ", fname
	print, " "
	return, -1
    endif

    print, " "
    print, "Reading HDF5 file: ", fname
    print, " "

    ; open the file
    ; -------------------------
    ;print, h5_parse(fname)   ; prints a useful overview
    file_id= h5f_open(fname)


    ; open the header group
    ; -------------------------
    hdr_group_id= h5g_open(file_id,"Header")
    if hdr_group_id lt 0 then begin
	print, "Problem, bad header"
	return, -1
    endif

    hubble=1.d0
    if keyword_set(h0) then begin
      hubble=0.7d0
      print, "snapshot values corrected for h= ", hubble
    endif
    hinv=1.d0/hubble

    npart= h5a_read(h5a_open_name(hdr_group_id,"NumPart_ThisFile"))
    massarr= h5a_read(h5a_open_name(hdr_group_id,"MassTable"))
    time= h5a_read(h5a_open_name(hdr_group_id,"Time"))
    redshift= h5a_read(h5a_open_name(hdr_group_id,"Redshift"))
    flag_sfr= h5a_read(h5a_open_name(hdr_group_id,"Flag_Sfr"))
    flag_feedbacktp= h5a_read(h5a_open_name(hdr_group_id,"Flag_Feedback"))
    npartTotal= h5a_read(h5a_open_name(hdr_group_id,"NumPart_Total"))
    flag_cooling= h5a_read(h5a_open_name(hdr_group_id,"Flag_Cooling"))
    numfiles= h5a_read(h5a_open_name(hdr_group_id,"NumFilesPerSnapshot"))
    cosmocrap= fltarr(4)
    cosmocrap[0]= h5a_read(h5a_open_name(hdr_group_id,"BoxSize"))
    cosmocrap[1]= h5a_read(h5a_open_name(hdr_group_id,"Omega0"))
    cosmocrap[2]= h5a_read(h5a_open_name(hdr_group_id,"OmegaLambda"))
    cosmocrap[3]= h5a_read(h5a_open_name(hdr_group_id,"HubbleParam"))
    flag_stellarage= h5a_read(h5a_open_name(hdr_group_id,"Flag_StellarAge"))
    flag_metals= h5a_read(h5a_open_name(hdr_group_id,"Flag_Metals"))

    h5g_close, hdr_group_id
    
    
	ascale=1.0
    if keyword_set(cosmological) then begin
      ascale=time
      hubble=cosmocrap[3]
      hinv=1.d0/hubble
      print,' cosmological sim: a=',ascale,' h=',hubble
    endif
    massarr = hinv * massarr
	if not keyword_set(cosmological) then time=time*hinv


	if flag_metals gt 1 then begin
	  print, " Reading file in multiple metal species mode: N_metals = ",flag_metals
	endif

   ; process things a bit
   if time lt 0.0 or time gt 1.0e4 then begin
        print, " "
        print, " PROBLEM with time read from header, returning .... "
        print, " "
        return, -1
    endif

    ; volker isn't using a flag, so do these manually
    flag_stargens= 1L
    flag_energydetail= 0L
    flag_parentid= 0L
    flag_starorig= 0L


    if keyword_set(header_only) then begin
        close, 1
        return, 1
    endif

    N=long(total(npartTotal,/double))    ; could also do   total(npart,/integer)

    ind=where((npartTotal gt 0) and (massarr eq 0))
    if ind(0) ne -1 then begin
        Nwithmass= long(total(npart(ind),/double))
        mass=fltarr(Nwithmass)
        print, "Nwithmasses= ",Nwithmass
    endif else begin
        Nwithmass= 0
    endelse

    Ngas=long(npartTotal(0))
    Nhalo=long(npartTotal(1))
    Ndisk=long(npartTotal(2))
    Nbulge=long(npartTotal(3))
    Nstars=long(npartTotal(4))
    Nbh=long(npartTotal(5))

    N_baryons= Ngas + Ndisk + Nbulge + Nstars
    N_tot= total(npart)
    print, "Ntot= ", N_tot
    print, npart

	if MULTI_PART_FILE_KEY eq 0 then begin 
	npart_done = 0*npartTotal
	if npartTotal(0) gt 0 then begin
	  xgas=fltarr(npartTotal(0)) & ygas=xgas & zgas=xgas & vxgas=xgas & vygas=xgas & vzgas=xgas
	  mgas=xgas & u=xgas & rho=xgas & nume=xgas & numh=xgas & sfr=xgas & hsml=xgas & gid=lonarr(npartTotal(0))
	  if flag_metals gt 0 then gmetals=fltarr(npartTotal(0),flag_metals)
	endif
	if npartTotal(1) gt 0 then begin
	  xhalo=fltarr(npartTotal(1)) & yhalo=xhalo & zhalo=xhalo & vxhalo=xhalo & vyhalo=xhalo & vzhalo=xhalo
	  mhalo=xhalo & dmid=lonarr(npartTotal(1))
	endif
	if npartTotal(2) gt 0 then begin
	  xdisk=fltarr(npartTotal(2)) & ydisk=xdisk & zdisk=xdisk & vxdisk=xdisk & vydisk=xdisk & vzdisk=xdisk
	  mdisk=xdisk & did=lonarr(npartTotal(2))
	endif
	if npartTotal(3) gt 0 then begin
	  xbulge=fltarr(npartTotal(3)) & ybulge=xbulge & zbulge=xbulge & vxbulge=xbulge & vybulge=xbulge & vzbulge=xbulge
	  mbulge=xbulge & bid=lonarr(npartTotal(3))
	endif
	if npartTotal(4) gt 0 then begin
	  xstars=fltarr(npartTotal(4)) & ystars=xstars & zstars=xstars & vxstars=xstars & vystars=xstars & vzstars=xstars
	  mstars=xstars & nsid=lonarr(npartTotal(4))
	  stellage=xstars & if flag_metals gt 0 then smetals=fltarr(npartTotal(4),flag_metals)
	endif
	if npartTotal(5) gt 0 then begin
	  xbh=fltarr(npartTotal(5)) & ybh=xbh & zbh=xbh & vxbh=xbh & vybh=xbh & vzbh=xbh
	  mbh=xbh & bhid=lonarr(npartTotal(5))
	  bhmdot=xbh & bhmass=xbh
	endif
	endif


    ; now get the particle data
    ; ------------------------------

    ;  Gas
    ; -----
    if npart(0) gt 0 then begin
	group_id= h5g_open(file_id,"PartType0")
	xyz= hinv * ascale * h5d_read(h5d_open(group_id,"Coordinates"))
	xgas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[0,*])
	ygas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[1,*])
	zgas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxgas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[0,*])
	vygas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[1,*])
	vzgas[npart_done(0):npart_done(0)+npart(0)-1]= transpose(xyz[2,*])
	if massarr[0] le 0 then mgas[npart_done(0):npart_done(0)+npart(0)-1]= hinv*h5d_read(h5d_open(group_id,"Masses")) else mgas= 0.0*xgas + massarr[0]
	u[npart_done(0):npart_done(0)+npart(0)-1]= h5d_read(h5d_open(group_id,"InternalEnergy"))
	rho[npart_done(0):npart_done(0)+npart(0)-1]= hubble*hubble * h5d_read(h5d_open(group_id,"Density")) / (ascale*ascale*ascale)
	if flag_metals gt 0 then gmetals_t= h5d_read(h5d_open(group_id,"Metallicity"))
	  if flag_metals gt 1 then begin
		dim_metal_tmp = size(gmetals_t,/dimensions)
		if (dim_metal_tmp(0) ne npart(0)) then $
			gmetals[npart_done(0):npart_done(0)+npart(0)-1,*]=transpose(gmetals_t) $
		else $
			gmetals[npart_done(0):npart_done(0)+npart(0)-1,*]=gmetals_t
	  endif ;else begin
	  		;gmetals[npart_done(0):npart_done(0)+npart(0)-1,*]=gmetals_t
	  ;endelse
	if flag_cooling gt 0 then nume[npart_done(0):npart_done(0)+npart(0)-1]= h5d_read(h5d_open(group_id,"ElectronAbundance"))
	if flag_cooling gt 0 then numh[npart_done(0):npart_done(0)+npart(0)-1]= h5d_read(h5d_open(group_id,"NeutralHydrogenAbundance"))
	if flag_sfr gt 0 then sfr[npart_done(0):npart_done(0)+npart(0)-1]= h5d_read(h5d_open(group_id,"StarFormationRate"))
	hsml[npart_done(0):npart_done(0)+npart(0)-1]= hinv * h5d_read(h5d_open(group_id,"SmoothingLength")) * ascale
	;pgas= h5d_read(h5d_open(group_id,"Potential"))
	gid[npart_done(0):npart_done(0)+npart(0)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  DM
    ; -----
    if npart(1) gt 0 then begin
	group_id= h5g_open(file_id,"PartType1")
	xyz= hinv * h5d_read(h5d_open(group_id,"Coordinates")) * ascale
	xhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[0,*])
	yhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[1,*])
	zhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[0,*])
	vyhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[1,*])
	vzhalo[npart_done(1):npart_done(1)+npart(1)-1]= transpose(xyz[2,*])
	if massarr[1] le 0 then mhalo[npart_done(1):npart_done(1)+npart(1)-1]= hinv* h5d_read(h5d_open(group_id,"Masses")) else mhalo= 0.0*xhalo + massarr[1]
	;phalo= h5d_read(h5d_open(group_id,"Potential"))
	dmid[npart_done(1):npart_done(1)+npart(1)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  Disk
    ; -------
    if npart(2) gt 0 then begin
	group_id= h5g_open(file_id,"PartType2")
	xyz= hinv* h5d_read(h5d_open(group_id,"Coordinates")) * ascale
	xdisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[0,*])
	ydisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[1,*])
	zdisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxdisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[0,*])
	vydisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[1,*])
	vzdisk[npart_done(2):npart_done(2)+npart(2)-1]= transpose(xyz[2,*])
	if massarr[2] le 0 then mdisk[npart_done(2):npart_done(2)+npart(2)-1]= hinv* h5d_read(h5d_open(group_id,"Masses")) else mdisk= 0.0*xdisk + massarr[2]
	;pdisk= h5d_read(h5d_open(group_id,"Potential"))
	did[npart_done(2):npart_done(2)+npart(2)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))

	diskage= -5.0*randomu(seed,Ndisk)
	diskage= -1.0*randomu(seed,Ndisk)    ; changed this for Marijn
	if (flag_metals gt 0) then $
	diskmetals= 10^(1.5*randomu(seed,Ndisk,flag_metals) - 2.0)*0.02

	h5g_close, group_id
    endif


    ;  Bulge
    ; -------
    if npart(3) gt 0 then begin
	group_id= h5g_open(file_id,"PartType3")
	xyz= hinv* h5d_read(h5d_open(group_id,"Coordinates")) * ascale
	xbulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[0,*])
	ybulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[1,*])
	zbulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxbulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[0,*])
	vybulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[1,*])
	vzbulge[npart_done(3):npart_done(3)+npart(3)-1]= transpose(xyz[2,*])
	if massarr[3] le 0 then mbulge[npart_done(3):npart_done(3)+npart(3)-1]= hinv* h5d_read(h5d_open(group_id,"Masses")) else mbulge= 0.0*xbulge + massarr[3]
	;pbulge= h5d_read(h5d_open(group_id,"Potential"))
	bid[npart_done(3):npart_done(3)+npart(3)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))

	bulgeage= -7.0 + 0.0*mbulge
	if flag_metals gt 0 then $
		bulgemetals= 0.001 + fltarr(Nbulge,flag_metals)

	h5g_close, group_id
    endif


    ;  New Stars
    ; ------------
    if npart(4) gt 0 then begin
	group_id= h5g_open(file_id,"PartType4")
	xyz= hinv* h5d_read(h5d_open(group_id,"Coordinates")) * ascale
	xstars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[0,*])
	ystars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[1,*])
	zstars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[2,*])
	xyz= h5d_read(h5d_open(group_id,"Velocities"))
	vxstars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[0,*])
	vystars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[1,*])
	vzstars[npart_done(4):npart_done(4)+npart(4)-1]= transpose(xyz[2,*])
	if massarr[4] le 0 then mstars[npart_done(4):npart_done(4)+npart(4)-1]=hinv* h5d_read(h5d_open(group_id,"Masses")) else mstars= 0.0*xstars + massarr[4]
	if flag_metals gt 0 then begin
		smetals_t= h5d_read(h5d_open(group_id,"Metallicity"))
	  if flag_metals gt 1 then begin
		dim_metal_tmp = size(smetals_t,/dimensions)
		if (dim_metal_tmp(0) ne npart(4)) then $
			smetals[npart_done(4):npart_done(4)+npart(4)-1,*]=transpose(smetals_t) $
		else $
			smetals[npart_done(4):npart_done(4)+npart(4)-1,*]=smetals_t
	  endif else begin
	  		smetals[npart_done(4):npart_done(4)+npart(4)-1]=smetals_t
	  endelse
	endif
	if flag_sfr gt 0 and flag_stellarage gt 0 then begin
		stellage[npart_done(4):npart_done(4)+npart(4)-1]= h5d_read(h5d_open(group_id,"StellarFormationTime"))
		if not keyword_set(cosmological) then stellage=stellage*hinv
	endif
	;pstars= h5d_read(h5d_open(group_id,"Potential"))
	nsid[npart_done(4):npart_done(4)+npart(4)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))
	h5g_close, group_id
    endif


    ;  Black Holes
    ; --------------
    if npart(5) gt 0 then begin
        group_id= h5g_open(file_id,"PartType5")
        xyz= hinv* h5d_read(h5d_open(group_id,"Coordinates")) * ascale
        xbh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[0,*])
        ybh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[1,*])
        zbh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[2,*])
        xyz= h5d_read(h5d_open(group_id,"Velocities"))
        vxbh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[0,*])
        vybh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[1,*])
        vzbh[npart_done(5):npart_done(5)+npart(5)-1]= transpose(xyz[2,*])
        if massarr[5] le 0 then mbh[npart_done(5):npart_done(5)+npart(5)-1]= hinv* h5d_read(h5d_open(group_id,"Masses")) else mbh= 0.0*xbh + massarr[5]
        ;pbh= h5d_read(h5d_open(group_id,"Potential"))
        bhid[npart_done(5):npart_done(5)+npart(5)-1]= h5d_read(h5d_open(group_id,"ParticleIDs"))
		if not keyword_set(skip_bh) then begin
        	bhmdot[npart_done(5):npart_done(5)+npart(5)-1]= h5d_read(h5d_open(group_id,"BH_Mdot"))
        	bhmass[npart_done(5):npart_done(5)+npart(5)-1]= h5d_read(h5d_open(group_id,"BH_Mass"))
        endif
        h5g_close, group_id
    endif

	npart_done = npart_done + npart

    ; done, get out of here
    ; -------------------------
    h5f_close, file_id


	if (MULTI_PART_FILE_KEY lt numfiles-1) then begin
	 print,'processed file ',MULTI_PART_FILE_KEY,'of ',numfiles
	 MULTI_PART_FILE_KEY = MULTI_PART_FILE_KEY + 1

	  exts_pt = string(MULTI_PART_FILE_KEY,'(I1)')
	   if MULTI_PART_FILE_KEY ge 10 then exts_pt=string(MULTI_PART_FILE_KEY,'(I2)')
      cmd= "/bin/ls "+frun+"/"+basename+"*_"+exts+'.'+exts_pt+".hdf5"
      spawn, cmd, result   &    print, "result= ", result
      if strlen(result) gt 0 and strmid(result,strlen(result)-4,4) eq "hdf5" then goto, file_is_hdf5
	endif



    ;
    ; a bit of data processing
    ;
    ; combine ids
    id= [-1]
    if n_elements(gid) gt 0 then id= [id, gid]
    if n_elements(dmid) gt 0 then id= [id, dmid]
    if n_elements(did) gt 0 then id= [id, did]
    if n_elements(bid) gt 0 then id= [id, bid]
    if n_elements(nsid) gt 0 then id= [id, nsid]
    if n_elements(bhid) gt 0 then id= [id, bhid]
    id= id[1:*]


    goto, calccenter

;=========================================================================================
;=========================================================================================


calccenter:

	com=[0.,0.,0.]
    if keyword_set(center_on_bh) then if Nbh gt 0 then com=[xbh(0),ybh(0),zbh(0)]
    if keyword_set(find_center) then com=fload_center(1) ;; fload_center_felix(1)
    alternate_com=com

    return, 0

end


