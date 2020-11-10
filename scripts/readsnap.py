import numpy as np
import h5py as h5py
import os.path
## This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO ##



def readsnap(sdir,snum,ptype,
    snapshot_name='snapshot',
    extension='.hdf5',
    h0=0,cosmological=0,skip_bh=0,four_char=0,
    header_only=0,loud=0):
    '''
    This is a sub-routine designed to copy a GIZMO snapshot portion - specifically
    all the data corresponding to particles of a given type - into active memory in 
    a parent structure for use in python. The routine automatically handles multi-part 
    snapshot files for you (concatenating), and works with python2.x and python3.x, 
    and both GIZMO hdf5 and un-formatted binary outputs.

    Syntax:
      P = readsnap(sdir,snum,ptype,....)
      
      Here "P" is a structure which contains the data. The snapshot file[s] are opened,
      the data fully copied out, and the file closed. This attempts to copy 
      all the common data types, all together, into "P". Three things to note: 
      (1) the fields in P (visible by typing P.keys()) are not given the same names 
      as those in the raw snapshot, but 'shorthand' names, for convenience. you should
      look at the keys and be sure you know which associate with which files.
      (2) because of the full-copy approach into new-named fields, this will not 
      handle arbitrary new data types (it is impossible to handle these in full 
      generality with un-formatted binary, you need to code the file-order and byte 
      numbers for each new data structure). if you add new fields (with e.g. 
      additional physics modules) beyond what this code looks for, you need to 
      add code here, or use the more general 'load_from_snapshot.py' routine.
      (3) also because of the full copy strategy, this routine is much more expensive 
      in time and memory compared to 'load_from_snapshot.py'. use that if you want 
      a light-weight, more flexible reading option, and have HDF5 outputs.

      For example, after calling, the 'Coordinates' field from the snapshot is 
      accessible from the new structure P by calling P['p']. 'Velocities' as P['v'],
      'Masses' as P['m']. Fields specific to gas include 'Density' as P['rho'], 
      'InternalEnergy' as P['u'], and more.

      More details and examples are given in the GIZMO user guide.

    Arguments:               
      sdir: parent directory (string) of the snapshot file or immediate snapshot sub-directory 
            if it is a multi-part file.
            
      snum: number (int) of the snapshot. e.g. snapshot_001.hdf5 is '1'
            Note for multi-part files, this is just the number of the 'set', i.e. 
            if you have snapshot_001.N.hdf5, set this to '1', not 'N' or '1.N'

      ptype: element type (int) = 0[gas],1,2,3,4,5[meaning depends on simulation, see
             user guide for details]. if your chosen 'value' is in the file header, 
             this will be ignored
      
    Optional:
      header_only: the structure "P" will return the file header, instead of the 
        particle data. you can see the data in the header then by simply typing 
        P.keys() -- this contains data like the time of the snapshot, as 
        P['Time']. With this specific routine the header information is only 
        saved if you choose this option. Default 0/False, turn on by setting to 
        1 or True.
      
      cosmological: default 0/False: turn on (set to 1/True) to convert cosmological 
        co-moving units to physical units. will specifically convert Coordinates, 
        Masses, Velocities, Densities, Smoothing Lengths, and Times/Ages. If this 
        is on, you do not need to set 'h0' (this will force it to be set -also-), 
        but it does no harm to set it as well.

      h0: default 0/False: turn on (set to 1/True) for the code to use the value 
        of the hubble constant (h = H0/100 km/s/Mpc) saved in the snapshot to convert 
        from code units. Recall, units of time, length, and mass in the code are in 
        h^-1. So this on means your units are physical, with no "h" in them. 
        Will specifically convert Coordinates, Masses, Densities, Smoothing Lengths, 
        and Times/Ages. 

      skip_bh: default 0/False: turn on (set to 1/True) to skip black hole-specific 
        fields for particles of type 5 (use if your snapshot contains elements of type 5, 
        but the black hole physics modules were not actually used; otherwise you will 
        get an error).
      
      four_char: default numbering is that snapshots with numbers below 1000 have 
        three-digit numbers. if they were numbered with four digits (e.g. snapshot_0001), 
        set this to 1 or True (default 0/False)
        
      snapshot_name: default 'snapshot': the code will automatically try a number of 
        common snapshot and snapshot-directory prefixes. but it can't guess all of them, 
        especially if you use an unusual naming convention, e.g. naming your snapshots 
        'xyzBearsBeetsBattleStarGalactica_001.hdf5'. In that case set this to the 
        snapshot name prefix (e.g. 'xyzBearsBeetsBattleStarGalactica')
        
      extension: default 'hdf5': again like 'snapshot' set if you use a non-standard 
        extension (it checks multiply options like 'h5' and 'hdf5' and 'bin'). but 
        remember the file must actually be hdf5 format!

      loud: print additional checks as it reads, useful for debugging, 
        set to 1 or True if desired (default 0/False)
    


    '''


    if (ptype<0): return {'k':-1};
    if (ptype>5): return {'k':-1};

    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    if(fname=='NULL'): return {'k':-1}
    if(loud==1): print('loading file : '+fname)

    ## open file and parse its header information
    nL = 0 # initial particle point to start at 
    if(fname_ext=='.hdf5'):
        file = h5py.File(fname,'r') # Open hdf5 snapshot file
        header_topdict = file["Header"] # Load header dictionary (to parse below)
        header_toparse = header_topdict.attrs
    else:
        file = open(fname) # Open binary snapshot file
        header_toparse = load_gadget_format_binary_header(file)

    npart = header_toparse["NumPart_ThisFile"]
    massarr = header_toparse["MassTable"]
    time = header_toparse["Time"]
    redshift = header_toparse["Redshift"]
    flag_sfr = header_toparse["Flag_Sfr"]
    flag_feedbacktp = header_toparse["Flag_Feedback"]
    npartTotal = header_toparse["NumPart_Total"]
    flag_cooling = header_toparse["Flag_Cooling"]
    numfiles = header_toparse["NumFilesPerSnapshot"]
    boxsize = header_toparse["BoxSize"]
    hubble = header_toparse["HubbleParam"]
    flag_stellarage = header_toparse["Flag_StellarAge"]
    flag_metals = header_toparse["Flag_Metals"]
    print("npart_file: ",npart)
    print("npart_total:",npartTotal)

    hinv=1.
    if (h0==1):
        hinv=1./hubble
    ascale=1.
    if (cosmological==1):
        ascale=time
        hinv=1./hubble
    if (cosmological==0): 
        time*=hinv
    
    boxsize*=hinv*ascale
    if (npartTotal[ptype]<=0): file.close(); return {'k':-1};
    if (header_only==1): file.close(); return {'k':0,'time':time,
        'boxsize':boxsize,'hubble':hubble,'npart':npart,'npartTotal':npartTotal};

    # initialize variables to be read
    pos=np.zeros([npartTotal[ptype],3],dtype=np.float64)
    vel=np.copy(pos)
    ids=np.zeros([npartTotal[ptype]],dtype=long)
    mass=np.zeros([npartTotal[ptype]],dtype=np.float64)
    if (ptype==0):
        ugas=np.copy(mass)
        rho=np.copy(mass)
        hsml=np.copy(mass) 
        #if (flag_cooling>0): 
        nume=np.copy(mass)
        numh=np.copy(mass)
        #if (flag_sfr>0): 
        sfr=np.copy(mass)
        metal=np.copy(mass)
    if (ptype == 0 or ptype == 4) and (flag_metals > 0):
        metal=np.zeros([npartTotal[ptype],flag_metals],dtype=np.float64)
    if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0):
        stellage=np.copy(mass)
    if (ptype == 5) and (skip_bh == 0):
        bhmass=np.copy(mass)
        bhmdot=np.copy(mass)

    # loop over the snapshot parts to get the different data pieces
    for i_file in range(numfiles):
        if (numfiles>1):
            file.close()
            fname = fname_base+'.'+str(i_file)+fname_ext
            if(fname_ext=='.hdf5'):
                file = h5py.File(fname,'r') # Open hdf5 snapshot file
            else:
                file = open(fname) # Open binary snapshot file
                header_toparse = load_gadget_format_binary_header(file)
                
        if (fname_ext=='.hdf5'):
            input_struct = file
            npart = file["Header"].attrs["NumPart_ThisFile"]
            bname = "PartType"+str(ptype)+"/"
        else:
            npart = header_toparse['NumPart_ThisFile']
            input_struct = load_gadget_format_binary_particledat(file, header_toparse, ptype, skip_bh=skip_bh)
            bname = ''
            
        
        # now do the actual reading
        if(npart[ptype]>0):
            nR=nL + npart[ptype]
            pos[nL:nR,:]=input_struct[bname+"Coordinates"]
            vel[nL:nR,:]=input_struct[bname+"Velocities"]
            ids[nL:nR]=input_struct[bname+"ParticleIDs"]
            mass[nL:nR]=massarr[ptype]
            if (massarr[ptype] <= 0.):
                mass[nL:nR]=input_struct[bname+"Masses"]
            if (ptype==0):
                ugas[nL:nR]=input_struct[bname+"InternalEnergy"]
                rho[nL:nR]=input_struct[bname+"Density"]
                hsml[nL:nR]=input_struct[bname+"SmoothingLength"]
                if (flag_cooling > 0): 
                    nume[nL:nR]=input_struct[bname+"ElectronAbundance"]
                    numh[nL:nR]=input_struct[bname+"NeutralHydrogenAbundance"]
                if (flag_sfr > 0):
                    sfr[nL:nR]=input_struct[bname+"StarFormationRate"]
            if (ptype == 0 or ptype == 4) and (flag_metals > 0):
                metal_t=input_struct[bname+"Metallicity"]
                if (flag_metals > 1):
                    if (metal_t.shape[0] != npart[ptype]): 
                        metal_t=np.transpose(metal_t)
                else:
                    metal_t=np.reshape(np.array(metal_t),(np.array(metal_t).size,1))
                metal[nL:nR,:]=metal_t
            if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0):
                stellage[nL:nR]=input_struct[bname+"StellarFormationTime"]
            if (ptype == 5) and (skip_bh == 0):
                bhmass[nL:nR]=input_struct[bname+"BH_Mass"]
                bhmdot[nL:nR]=input_struct[bname+"BH_Mdot"]
            nL = nR # sets it for the next iteration	

	## correct to same ID as original gas particle for new stars, if bit-flip applied
    if ((np.min(ids)<0) | (np.max(ids)>1.e9)):
        bad = (ids < 0) | (ids > 1.e9)
        ids[bad] += (long(1) << 31)

    # do the cosmological conversions on final vectors as needed
    pos *= hinv*ascale # snapshot units are comoving
    mass *= hinv
    vel *= np.sqrt(ascale) # remember gizmo's (and gadget's) weird velocity units!
    if (ptype == 0):
        rho *= (hinv/((ascale*hinv)**3))
        hsml *= hinv*ascale
    if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0) and (cosmological == 0):
        stellage *= hinv
    if (ptype == 5) and (skip_bh == 0):
        bhmass *= hinv

    file.close();
    if (ptype == 0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'u':ugas,'rho':rho,'h':hsml,'ne':nume,'nh':numh,'sfr':sfr,'z':metal};
    if (ptype == 4):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'z':metal,'age':stellage}
    if (ptype == 5) and (skip_bh == 0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'mbh':bhmass,'mdot':bhmdot}
    return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids}



def check_if_filename_exists(sdir,snum,snapshot_name='snapshot',extension='.hdf5',four_char=0):
    for extension_touse in [extension,'.bin','']:
        fname=sdir+'/'+snapshot_name+'_'
        ext='00'+str(snum);
        if (snum>=10): ext='0'+str(snum)
        if (snum>=100): ext=str(snum)
        if (four_char==1): ext='0'+ext
        if (snum>=1000): ext=str(snum)
        fname+=ext
        fname_base=fname

        s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1];
        if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2];

        ## try several common notations for the directory/filename structure
        fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is it a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap(snapdir)' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+snapdir_specific+'_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory? (we assume this means multi-part files)
            fname_base=sdir+'/snapdir_'+ext+'/'+snapshot_name+'_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory AND named 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snapdir_'+ext+'/'+'snap_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## wow, still couldn't find it... ok, i'm going to give up!
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        fname_found = fname;
        fname_base_found = fname_base;
        fname_ext = extension_touse
        break; # filename does exist! 
    return fname_found, fname_base_found, fname_ext;



def load_gadget_format_binary_header(f):
    ### Read header.
    import array
    # Skip 4-byte integer at beginning of header block.
    f.read(4)
    # Number of particles of each type. 6*unsigned integer.
    Npart = array.array('I')
    Npart.fromfile(f, 6)
    # Mass of each particle type. If set to 0 for a type which is present, 
    # individual particle masses from the 'mass' block are used instead.
    # 6*double.
    Massarr = array.array('d')
    Massarr.fromfile(f, 6)
    # Expansion factor (or time, if non-cosmological sims) of output. 1*double. 
    a = array.array('d')
    a.fromfile(f, 1)
    a = a[0]
    # Redshift of output. Should satisfy z=1/a-1. 1*double.
    z = array.array('d')
    z.fromfile(f, 1)
    z = float(z[0])
    # Flag for star formation. 1*int.
    FlagSfr = array.array('i')
    FlagSfr.fromfile(f, 1)
    # Flag for feedback. 1*int.
    FlagFeedback = array.array('i')
    FlagFeedback.fromfile(f, 1)
    # Total number of particles of each type in the simulation. 6*int.
    Nall = array.array('i')
    Nall.fromfile(f, 6)
    # Flag for cooling. 1*int.
    FlagCooling = array.array('i')
    FlagCooling.fromfile(f, 1)
    # Number of files in each snapshot. 1*int.
    NumFiles = array.array('i')
    NumFiles.fromfile(f, 1)
    # Box size (comoving kpc/h). 1*double.
    BoxSize = array.array('d')
    BoxSize.fromfile(f, 1)
    # Matter density at z=0 in units of the critical density. 1*double.
    Omega_Matter = array.array('d')
    Omega_Matter.fromfile(f, 1)
    # Vacuum energy density at z=0 in units of the critical density. 1*double.
    Omega_Lambda = array.array('d')
    Omega_Lambda.fromfile(f, 1)
    # Hubble parameter h in units of 100 km s^-1 Mpc^-1. 1*double.
    h = array.array('d')
    h.fromfile(f, 1)
    h = float(h[0])
    # Creation times of stars. 1*int.
    FlagAge = array.array('i')
    FlagAge.fromfile(f, 1)
    # Flag for metallicity values. 1*int.
    FlagMetals = array.array('i')
    FlagMetals.fromfile(f, 1)

    # For simulations that use more than 2^32 particles, most significant word 
    # of 64-bit total particle numbers. Otherwise 0. 6*int.
    NallHW = array.array('i')
    NallHW.fromfile(f, 6)

    # Flag that initial conditions contain entropy instead of thermal energy
    # in the u block. 1*int.
    flag_entr_ics = array.array('i')
    flag_entr_ics.fromfile(f, 1)

    # Unused header space. Skip to particle positions.
    f.seek(4+256+4+4)

    return {'NumPart_ThisFile':Npart, 'MassTable':Massarr, 'Time':a, 'Redshift':z, \
    'Flag_Sfr':FlagSfr[0], 'Flag_Feedback':FlagFeedback[0], 'NumPart_Total':Nall, \
    'Flag_Cooling':FlagCooling[0], 'NumFilesPerSnapshot':NumFiles[0], 'BoxSize':BoxSize[0], \
    'Omega_Matter':OmegaMatter[0], 'Omega_Lambda':OmegaLambda[0], 'HubbleParam':h, \
    'Flag_StellarAge':FlagAge[0], 'Flag_Metals':FlagMetals[0], 'Nall_HW':NallHW, \
    'Flag_EntrICs':flag_entr_ics[0]}


def load_gadget_format_binary_particledat(f, header, ptype, skip_bh=0):
    ## load old format=1 style gadget-format binary snapshot files (unformatted fortran binary)
    import array
    gas_u=0.; gas_rho=0.; gas_ne=0.; gas_nhi=0.; gas_hsml=0.; gas_SFR=0.; star_age=0.; 
    zmet=0.; bh_mass=0.; bh_mdot=0.; mm=0.;
    Npart = header['NumPart_ThisFile']
    Massarr = header['MassTable']
    NpartTot = np.sum(Npart)
    NpartCum = np.cumsum(Npart)
    n0 = NpartCum[ptype] - Npart[ptype]
    n1 = NpartCum[ptype]
    
    ### particles positions. 3*Npart*float.
    pos = array.array('f')
    pos.fromfile(f, 3*NpartTot)
    pos = np.reshape(pos, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### particles velocities. 3*Npart*float.
    vel = array.array('f')
    vel.fromfile(f, 3*NpartTot)
    vel = np.reshape(vel, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### Particle IDs. # (Npart[0]+...+Npart[5])*int
    id = array.array('i')
    id.fromfile(f, NpartTot)
    id = np.array(id)
    f.read(4+4) # Read block size fields.
        
    ### Variable particle masses. 
    Npart_MassCode = np.copy(np.array(Npart))
    Npart=np.array(Npart)
    Npart_MassCode[(Npart <= 0) | (np.array(Massarr,dtype='d') > 0.0)] = long(0)
    NwithMass = np.sum(Npart_MassCode)
    mass = array.array('f')
    mass.fromfile(f, NwithMass)
    f.read(4+4) # Read block size fields.
    if (Massarr[ptype]==0.0):
        Npart_MassCode_Tot = np.cumsum(Npart_MassCode)
        mm = mass[Npart_MassCode_Tot[ptype]-Npart_MassCode[ptype]:Npart_MassCode_Tot[ptype]]

    if ((ptype == 0) | (ptype == 4) | (ptype == 5)):
        if (Npart[0]>0):
            ### Internal energy of gas particles ((km/s)^2).
            gas_u = array.array('f')
            gas_u.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.
            ### Density for the gas paraticles (units?).
            gas_rho = array.array('f')
            gas_rho.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.

            if (header['Flag_Cooling'] > 0):
                ### Electron number density for gas particles (fraction of n_H; can be >1).
                gas_ne = array.array('f')
                gas_ne.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.
                ### Neutral hydrogen number density for gas particles (fraction of n_H).
                gas_nhi = array.array('f')
                gas_nhi.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

            ### Smoothing length (kpc/h). ###
            gas_hsml = array.array('f')
            gas_hsml.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.

            if (header['Flag_Sfr'] > 0):
                ### Star formation rate (Msun/yr). ###
                gas_SFR = array.array('f')
                gas_SFR.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

        if (Npart[4]>0):
            if (header['Flag_Sfr'] > 0):
                if (header['Flag_StellarAge'] > 0):
                    ### Star formation time (in code units) or scale factor ###
                    star_age = array.array('f')
                    star_age.fromfile(f, Npart[4])
                    f.read(4+4) # Read block size fields.
        
        if (Npart[0]+Npart[4]>0):
            if (header['Flag_Metals'] > 0):
                ## Metallicity block (species tracked = Flag_Metals)
                if (Npart[0]>0):
                    gas_z = array.array('f')
                    gas_z.fromfile(f, header['Flag_Metals']*Npart[0])
                if (Npart[4]>0):
                    star_z = array.array('f')
                    star_z.fromfile(f, header['Flag_Metals']*Npart[4])
                f.read(4+4) # Read block size fields.
                if (ptype==0): zmet=np.reshape(gas_z,(-1,header['Flag_Metals']))
                if (ptype==4): zmet=np.reshape(star_z,(-1,header['Flag_Metals']))
        
        if (Npart[5]>0):
            if (skip_bh > 0):
                ## BH mass (same as code units, but this is the separately-tracked BH mass from particle mass)
                bh_mass = array.array('f')
                bh_mass.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
                ## BH accretion rate in snapshot
                bh_mdot = array.array('f')
                bh_mdot.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
    
    return {'Coordinates':pos[n0:n1,:], 'Velocities':vel[n0:n1,:], 'ParticleIDs':id[n0:n1], \
        'Masses':mm, 'Metallicity':zmet, 'StellarFormationTime':star_age, 'BH_Mass':bh_mass, \
        'BH_Mdot':bh_mdot, 'InternalEnergy':gas_u, 'Density':gas_rho, 'SmoothingLength':gas_hsml, \
        'ElectronAbundance':gas_ne, 'NeutralHydrogenAbundance':gas_nhi, 'StarFormationRate':gas_SFR}
