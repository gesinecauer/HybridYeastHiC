#!/bin/bash
# make_matrix.sh
# shell script for make_matrix.cpp
# Seungsoo Kim

dir="../references"
anns="../nobackup/annotations"
out="../nobackup"
masks="../masks"

samp=$1
ref=$2
bsize=$3
mask=$4
minavg=$5

override_homobin=false

if [ $stype = "single" ]; then
  if [ ! -r $anns/$ref.$bsize.allcombos.homology_noisolated.matrixbin ] || $override_homobin; then
    echo "Identify homologous bins"
    ./homology.sh $bsize
  fi
fi

./make_matrix $anns/$ref.chrom_lengths $bsize \
$anns/$ref.$bsize.allcombos.homology_noisolated.matrixbin \
$masks/$mask.bed <(grep -v invalid $out/assigned/$ref/$samp.assigned) $minavg \
$out/matrix/$ref/$bsize/$samp.$mask.rowsums $out/matrix/$ref/$bsize/$samp.$mask.matrix
