import h5py
import numpy
import os
## This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO ##


def load_from_snapshot(value,ptype,sdir,snum,particle_mask=numpy.zeros(0),axis_mask=numpy.zeros(0),
        units_to_physical=True,four_char=False,snapshot_name='snapshot',snapdir_name='snapdir',extension='.hdf5'):
    '''

    The routine 'load_from_snapshot' is designed to load quantities directly from GIZMO 
    snapshots in a robust manner, independent of the detailed information actually saved 
    in the snapshot. It is able to do this because of how HDF5 works, so it  --only-- 
    works for HDF5-format snapshots [For binary-format, you need to know -exactly- the 
    datasets saved and their order in the file, which means you cannot do this, and should
    use the 'readsnap.py' routine instead.]
    
    The routine automatically handles multi-part snapshot files for you (concatenating). 
    This should work with both python2.x and python3.x

    Syntax:
      loaded_value = load_from_snapshot(value,ptype,sdir,snum,....)

      For example, to load the coordinates of gas (type=0) elements in the file
      snapshot_001.hdf5 (snum=1) located in the active directory ('.'), just call
      xyz_coordinates = load_from_snapshot('Coordinates',0,'.',1)

      More details and examples are given in the GIZMO user guide.

    Arguments:
      value: the value to extract from the HDF5 file. this is a string with the same name 
             as in the HDF5 file. if you arent sure what those values might be, setting 
             value to 'keys' will return a list of all the HDF5 keys for the chosen 
             particle type, or 'header_keys' will return all the keys in the header.
             (example: 'Time' returns the simulation time in code units (single scalar). 
                       'Coordinates' will return the [x,y,z] coordinates in an [N,3] 
                       matrix for the N resolution elements of the chosen type)
  
      ptype: element type (int) = 0[gas],1,2,3,4,5[meaning depends on simulation, see
             user guide for details]. if your chosen 'value' is in the file header, 
             this will be ignored
             
      sdir: parent directory (string) of the snapshot file or immediate snapshot 
            sub-directory if it is a multi-part file
            
      snum: number (int) of the snapshot. e.g. snapshot_001.hdf5 is '1'
            Note for multi-part files, this is just the number of the 'set', i.e. 
            if you have snapshot_001.N.hdf5, set this to '1', not 'N' or '1.N'
      
    Optional:
      particle_mask: if set to a mask (boolean array), of length N where N is the number
        of elements of the desired ptype, will return only those elements
      
      axis_mask: if set to a mask (boolean array), return only the chosen -axis-. this
        is useful for some quantities like metallicity fields, with [N,X] dimensions 
        where X is large (lets you choose to read just one of the "X")
        
      units_to_physical: default 'True': code will auto-magically try to detect if the 
        simulation is cosmological by comparing time and redshift information in the 
        snapshot, and if so, convert units to physical. if you want default snapshot units
        set this to 'False'
        
      four_char: default numbering is that snapshots with numbers below 1000 have 
        three-digit numbers. if they were numbered with four digits (e.g. snapshot_0001), 
        set this to 'True' (default False)
        
      snapshot_name: default 'snapshot': the code will automatically try a number of 
        common snapshot and snapshot-directory prefixes. but it can't guess all of them, 
        especially if you use an unusual naming convention, e.g. naming your snapshots 
        'xyzBearsBeetsBattleStarGalactica_001.hdf5'. In that case set this to the 
        snapshot name prefix (e.g. 'xyzBearsBeetsBattleStarGalactica')
      
      snapdir_name: default 'snapdir': like 'snapshot_name', set this if you use a 
        non-standard prefix for snapshot subdirectories (directories holding multi-part
        snapshots pieces)
        
      extension: default 'hdf5': again like 'snapshot' set if you use a non-standard 
        extension (it checks multiply options like 'h5' and 'hdf5' and 'bin'). but 
        remember the file must actually be hdf5 format!

    '''



    # attempt to verify if a file with this name and directory path actually exists
    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,snapdir_name=snapdir_name,extension=extension,four_char=four_char)
    # if no valid file found, give up
    if(fname=='NULL'): 
        print('Could not find a valid file with this path/name/extension - please check these settings')
        return 0
    # check if file has the correct extension
    if(fname_ext!=extension): 
        print('File has the wrong extension, you specified ',extension,' but found ',fname_ext,' - please specify this if it is what you actually want')
        return 0
    # try to open the file
    try: 
        file = h5py.File(fname,'r') # Open hdf5 snapshot file
    except:
        print('Unexpected error: could not read hdf5 file ',fname,' . Please check the format, name, and path information is correct')
        return 0
        
    # try to parse the header
    try:
        header_toparse = file["Header"].attrs # Load header dictionary (to parse below)
    except:
        print('Was able to open the file but not the header, please check this is a valid GIZMO hdf5 file')
        file.close()
        return 0
    # check if desired value is contained in header -- if so just return it and exit
    if(value=='header_keys')|(value=='Header_Keys')|(value=='HEADER_KEYS')|(value=='headerkeys')|(value=='HeaderKeys')|(value=='HEADERKEYS')|((value=='keys' and not (ptype == 0 or ptype == 1 or ptype == 2 or ptype == 3 or ptype == 4 or ptype == 5))):
        q = header_toparse.keys()
        print('Returning list of keys from header, includes: ',q)
        file.close()
        return q
    if(value in header_toparse):
        q = header_toparse[value] # value contained in header, no need to go further
        file.close()
        return q

    # ok desired quantity is not in the header, so we need to go into the particle data

    # check that a valid particle type is specified
    if not (ptype == 0 or ptype == 1 or ptype == 2 or ptype == 3 or ptype == 4 or ptype == 5):
        print('Particle type needs to be an integer = 0,1,2,3,4,5. Returning 0')
        file.close()
        return 0
    # check that the header contains the expected data needed to parse the file
    if not ('NumFilesPerSnapshot' in header_toparse and 'NumPart_Total' in header_toparse
        and 'Time' in header_toparse and 'Redshift' in header_toparse 
        and 'HubbleParam' in header_toparse and 'NumPart_ThisFile' in header_toparse):
        print('Header appears to be missing critical information. Please check that this is a valid GIZMO hdf5 file')
        file.close()
        return 0
    # parse data needed for checking sub-files 
    numfiles = header_toparse["NumFilesPerSnapshot"]
    npartTotal = header_toparse["NumPart_Total"]
    if(npartTotal[ptype]<1): 
        print('No particles of designated type exist in this snapshot, returning 0')
        file.close()
        return 0
    # parse data needed for converting units [if necessary]
    if(units_to_physical):
        time = header_toparse["Time"]
        z = header_toparse["Redshift"]
        hubble = header_toparse["HubbleParam"]
        cosmological = False
        ascale = 1.0;
        # attempt to guess if this is a cosmological simulation from the agreement or lack thereof between time and redshift. note at t=1,z=0, even if non-cosmological, this won't do any harm
        if(numpy.abs(time*(1.+z)-1.) < 1.e-6): 
            cosmological=True; ascale=time;
    # close the initial header we are parsing
    file.close()
    
    # now loop over all snapshot segments to identify and extract the relevant particle data
    check_counter = 0
    for i_file in range(numfiles):
        # augment snapshot sub-set number
        if (numfiles>1): fname = fname_base+'.'+str(i_file)+fname_ext  
        # check for existence of file
        if(os.stat(fname).st_size>0):
            # exists, now try to read it
            try: 
                file = h5py.File(fname,'r') # Open hdf5 snapshot file
            except:
                print('Unexpected error: could not read hdf5 file ',fname,' . Please check the format, name, and path information is correct, and that this file is not corrupted')
                return 0
            # read in, now attempt to parse. first check for needed information on particle number
            npart = file["Header"].attrs["NumPart_ThisFile"]
            if(npart[ptype] > 1):
                # return particle key data, if requested
                if((value=='keys')|(value=='Keys')|(value=='KEYS')): 
                    q = file['PartType'+str(ptype)].keys()
                    print('Returning list of valid keys for this particle type: ',q)
                    file.close()
                    return q
                # check if requested data actually exists as a valid keyword in the file
                if not (value in file['PartType'+str(ptype)].keys()):
                    print('The value ',value,' given does not appear to exist in the file ',fname," . Please check that you have specified a valid keyword. You can run this routine with the value 'keys' to return a list of valid value keys. Returning 0")
                    file.close()
                    return 0
                # now actually read the data
                axis_mask = numpy.array(axis_mask)
                if(axis_mask.size > 0):
                    q_t = numpy.array(file['PartType'+str(ptype)+'/'+value+'/']).take(axis_mask,axis=1)
                else:
                    q_t = numpy.array(file['PartType'+str(ptype)+'/'+value+'/'])                    
                # check data has non-zero size
                if(q_t.size > 0): 
                    # if this is the first time we are actually reading it, parse it and determine the shape of the vector, to build the data container
                    if(check_counter == 0): 
                        qshape=numpy.array(q_t.shape); qshape[0]=0; q=numpy.zeros(qshape); check_counter+=1;
                    # add the data to our appropriately-shaped container, now
                    try:
                        q = numpy.concatenate([q,q_t],axis=0)
                    except:
                        print('Could not concatenate data for ',value,' in file ',fname,' . The format appears to be inconsistent across your snapshots or with the usual GIZMO conventions. Please check this is a valid GIZMO snapshot file.')
                        file.close()
                        return 0
            file.close()
        else:
            print('Expected file ',fname,' appears to be missing. Check if your snapshot has the complete data set here')
            
    # convert units if requested by the user. note this only does a few obvious units: there are many possible values here which cannot be anticipated!
    if(units_to_physical):
        hinv=1./hubble; rconv=ascale*hinv;
        if((value=='Coordinates')|(value=='SmoothingLength')): q*=rconv; # comoving length
        if(value=='Velocities'): q *= numpy.sqrt(ascale); # special comoving velocity units
        if((value=='Density')|(value=='Pressure')): q *= hinv/(rconv*rconv*rconv); # density = mass/comoving length^3
        if((value=='StellarFormationTime')&(cosmological==False)): q*=hinv; # time has h^-1 in non-cosmological runs
        if((value=='Masses')|('BH_Mass' in value)|(value=='CosmicRayEnergy')|(value=='PhotonEnergy')): q*=hinv; # mass x [no-h] units

    # return final value, if we have not already
    particle_mask=numpy.array(particle_mask)
    if(particle_mask.size > 0): q=q.take(particle_mask,axis=0)
    return q





def check_if_filename_exists(sdir,snum,snapshot_name='snapshot',snapdir_name='snapdir',extension='.hdf5',four_char=False):
    '''
    This subroutine attempts to check if a snapshot or snapshot directory with 
    valid GIZMO outputs exists. It will check several common conventions for 
    file and directory names, and extensions. 
    
    Input:
      sdir: parent directory of the snapshot file or immediate snapshot sub-directory 
            if it is a multi-part file. string.
            
      snum: number (int) of the snapshot. e.g. snapshot_001.hdf5 is '1'
    
    Optional:
      snapshot_name: default 'snapshot': the code will automatically try a number of 
        common snapshot and snapshot-directory prefixes. but it can't guess all of them, 
        especially if you use an unusual naming convention, e.g. naming your snapshots 
        'xyzBearsBeetsBattleStarGalactica_001.hdf5'. In that case set this to the 
        snapshot name prefix (e.g. 'xyzBearsBeetsBattleStarGalactica')
      
      snapdir_name: default 'snapdir': like 'snapshot_name', set this if you use a 
        non-standard prefix for snapshot subdirectories (directories holding multi-part
        snapshots pieces)
        
      extension: default 'hdf5': again like 'snapshot' set if you use a non-standard 
        extension (it checks multiply options like 'h5' and 'hdf5' and 'bin'). but 
        remember the file must actually be hdf5 format!

      four_char: default numbering is that snapshots with numbers below 1000 have 
        three-digit numbers. if they were numbered with four digits (e.g. snapshot_0001), 
        set this to 'True' (default False)        
    '''    
    
    # loop over possible extension names to try and check for valid files
    for extension_touse in [extension,'.h5','.bin','']:
        fname=sdir+'/'+snapshot_name+'_'
        
        # begin by identifying the snapshot extension with the file number
        ext='00'+str(snum);
        if (snum>=10): ext='0'+str(snum)
        if (snum>=100): ext=str(snum)
        if (four_char==True): ext='0'+ext
        if (snum>=1000): ext=str(snum)
        fname+=ext
        fname_base=fname

        # isolate the specific path up to the snapshot name, because we will try to append several different choices below
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
            fname_base=sdir+'/'+snapdir_name+'_'+ext+'/'+snapshot_name+'_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory AND named 'snap' instead of 'snapshot'?
            fname_base=sdir+'/'+snapdir_name+'_'+ext+'/'+'snap_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## wow, still couldn't find it... ok, i'm going to give up!
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        if(os.stat(fname).st_size <= 0):
            ## file exists but is null size, do not use
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        fname_found = fname;
        fname_base_found = fname_base;
        fname_ext = extension_touse
        break; # filename does exist! 
    return fname_found, fname_base_found, fname_ext;
