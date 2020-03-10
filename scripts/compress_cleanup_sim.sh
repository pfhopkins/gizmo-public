#!/bin/sh
##
## This is a script designed to 'clean up' GIZMO simulation directories after the 
##   simulation is complete, for purposes of archiving and analyzing the simulations.
##   It assumes certain directory and filename conventions, so you'll need to modify
##   it to match your conventions. But I (PFH) am including it in the GIZMO repo 
##   because the notes below contain some useful information about what you can delete
##   or compress after simulations are completed.
##
## Begin by checking for a dummy empty file called 'clean_me' which is required for 
##   a directory to be cleaned. this is to prevent accidentally working on active dirs
done=0
if [ -f "clean_me" ] 
then
done=1
fi
## Now clean up the directory. Assumes the source code is a folder named 'code' 
##   in your run-time directory, the snapshots are written to a sub-folder named 
##   'output', and the stdout/stderr files (named gizmo.out/gizmo.err) are written
##   into that same run-time directory. This is all good code practice so that all
##   the information from your run is in one place, completely recoverable
if [ $done -eq 1 ]
then
  ## delete stdout/stderr, log files. these are very large and only useful while running
  rm gizmo.*
  rm *.err
  rm *.log
  rm *.out
  ## ewald tables are only used in-run, they can be removed
  rm ewald*
  ## by default, I never copy the full spcool_tables directory into a run-folder, but 
  ##   use symlinks to a common directory, to save space. but if you did copy them, 
  ##   they take a lot of space for identical files. 
  rm -r spcool_tables
  ## check if the source code directory is already tarballed and zipped. if not, 
  ##   go into it and 'make clean' then tar and zip the source code directory. 
  ##   this is really important for archival purposes (e.g. backing up simulations). 
  ##   you always want to save the source code, of course. but because of the 
  ##   git/version-control system, even if the source code it small, it can easily
  ##   contain thousands of tiny hidden files from version control. this plays 
  ##   havok with disk systems if you have a lot of them and will get you into trouble.
  if [ ! -f ./code.tgz ]; then
  if [ -d "code" ]; then
    cd code
    make clean
    cd ..
    tar -czvf code.tgz code
    rm -r code
  fi
  fi
fi
## Remove any file labeled 'core' in the run or output directory. These are core-dumps
##   from certain types of errors on certain machines, only used for debugging.
rm core*
rm output/core*
## Touch the files to update their access times appropriately
#touch -a *
#touch -a */*
#touch -a */*/*
## Check that the folder 'output' exists, which this script assumes is where all 
##   the simulation outputs are being written (should be set in the parameter file)
if [ -d "output" ]
then
    cd output
if [ $done -eq 1 ]
then
    ## remove the 'timing' files in output. especially if you are running with extended
    ##   outputs, these files can be enormous (they are ascii, so poorly-compressed), 
    ##   and written every timestep in some cases. and they contain only information 
    ##   about optimization and timing and code checks, used in-run but not to archive
    rm timebin.txt
    rm timings.txt
    rm balance.txt
    rm cpu.txt
    rm energy.txt
    rm info.txt
    ## remove the restartfiles. since this is a complete flash of system memory 
    ##   (with a backup), it can easily be huge (bigger than all the snapshots), 
    ##   but is useless once runs are completed (and cannot be transferred between
    ##   systems or even slightly different compiled versions of the code.
    rm -r restartfiles
    ## some of the output files are useful, but can get incredibly large 
    ##   (depending on the details of the simulation), because they are written 
    ##   as ascii with a lot of repeating numbers. this includes all of the 'baryonic' 
    ##   diagnostic files like sfr.txt, HIIheating.txt, blackhole information, etc. 
    ##   (especially the blackhole_details files). these can be zipped to compress
    ##   them by typical factors of ~100, so we do that here if they exist.
    if [ -f "sfr.txt" ]
    then
      tar -czvf sfr.tgz sfr.txt
      rm sfr.txt
    fi
    if [ -f "HIIheating.txt" ]
    then
      tar -czvf HIIheating.tgz HIIheating.txt
      rm HIIheating.txt
    fi
    if [ -f "blackholes.txt" ]
    then
      tar -czvf blackholes.tgz blackholes.txt
      rm blackholes.txt
    fi
    if [ -d "blackhole_details" ]
    then
      tar -czvf blackhole_details.tgz blackhole_details
      rm -r blackhole_details
    fi
fi
## Now we compress the snapshots themselves. This calls a python script, also 
##   included in the GIZMO source, which losslessly compresses the hdf5 files 
##   for snapshots. This assumes you are writing in hdf5, and assumes that 
##   your default output/snapshot name prefix is "snapshot" (and if its is 
##   writing multi-component snapshots, the directory prefix for them is 
##   "snapdir"). Change those if needed. Make sure the python code is linked 
##   and working properly. Again the compression is lossless, this is akin to 
##   gzipping the files, but allows you to continue using them like a 'normal'
##   hdf5 file, so its ideal, and typically gives a factor of ~2 compression 
##   (though it can be much higher if the data is highly regular). Make sure 
##   if you use this that you link the correct location of "compress_gizmosnap", 
##   here just assumed to be ${HOME}
#touch -a */*/*
for snapshotname in ./snapshot_*.hdf5
do
    python3 ${HOME}/compress_gizmosnap.py $snapshotname
done
for snapshotname in ./snapdir_*/snapshot_*.hdf5
do
    python3 ${HOME}/compress_gizmosnap.py $snapshotname
done
    cd ..
fi
## All done!
