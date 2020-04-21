#!/bin/bash -e

if [ -z "$1" ] || [ -z "$2" ] ; then
    echo "usage: $(basename $0) git_ref_hash cache_dir_path"
    echo "generate a cache of IFO budget traces from a particular git commit"
    exit 1
fi

ref_hash="$1"
cache_dir="$2"

mkdir -p $cache_dir
cache_dir=$(cd $cache_dir && pwd)
gwinc_dir=$cache_dir/gwinc
mkdir -p $gwinc_dir

git archive $ref_hash | tar -x -C $gwinc_dir

cd $gwinc_dir

export LOG_LEVEL=INFO
for ifo in $(python3 -c "import gwinc; print(' '.join(gwinc.IFOS))") ; do
    python3 -m gwinc --save $cache_dir/${ifo}.h5 $ifo
done

echo $ref_hash > $cache_dir/ref_hash
