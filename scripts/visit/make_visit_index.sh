#!/bin/bash

function print_help() {
    echo "Usage: $1 [nchunks] [file list]"
}

if [ -z "$1" ]; then
    print_help $0
    exit 1
fi

if ! [[ $1 =~ [0-9]+$ ]]; then
    print_help $0
    exit 1
fi

echo "!NBLOCKS $1"
shift
for f in $*; do
    echo $f
done
