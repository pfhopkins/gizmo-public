import os
import sys
import h5py
import argparse
import functools
import tempfile

## this is the core compression subroutine
def copy_level0(fs, fd, name, node):
    if isinstance(node, h5py.Dataset):
        # create the new dataset, defaulting to gzip level 4 compression,
        #  including chunking/shuffle options, which work well.
        dnew = fd.create_dataset(name, data=node, dtype=node.dtype, chunks=True,
            shuffle=True, compression="gzip", compression_opts=4, fletcher32=True);
        print('   ..dataset {} ({}) compressed'.format(name,node.dtype));
    elif isinstance(node, h5py.Group) and name == 'Header':
        # header entries don't get compressed, but are negligible for storage
        fs.copy(name, fd, name=name);
        print('  Header copied');
    else:
        print('  Group {}'.format(name));

## this is the main loop, called by parser, which does all the os-level operations
def main(filename):
    print('Losslessly compressing snapshot {}'.format(filename)); # start
    fs = h5py.File(filename, 'r'); # open file
    # check that the file has not already been compressed, in which case, exit
    if('CompactLevel' in fs['Header'].attrs):
        print(' .. this snapshot has already been compressed, done');
        fs.close();
        return;
    # create a temporary working file for intermediate steps
    tmpfilename = filename+'__tmp__';
    fd = h5py.File(tmpfilename, 'w');
    # recursively do the work here on the entries (call and compress)
    copy_datasets = functools.partial(copy_level0, fs, fd);
    fs.visititems(copy_datasets);
    # encode that the file is compressed, so it can be skipped in future
    fd['Header'].attrs['CompactLevel'] = 0;
    # close and rename temp file to replace original file
    fd.close(); fs.close(); os.rename(tmpfilename,filename);
    print('  Completed');
    return;

## this is the parser, meant to be called from the command-line
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Losslessly compress GIZMO hdf5 snapshots.'
                    '(c) R. Feldmann 2017, modified by PFH 2020')
    parser.add_argument('filename', help='hdf5 file to be compactified '
                        '(file will be over-written)')

    if len(sys.argv[1:]) == 0:
        print('Error: {} requires more parameters'.format(__file__));
        parser.print_help();
        parser.exit();

    args = parser.parse_args();
    main(args.filename)
