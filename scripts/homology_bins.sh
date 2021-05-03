#!/bin/bash
# homology_bins.sh
# calculates homologous bins
# Seungsoo Kim

refs="../references"
anns="../nobackup/annotations"
bsize=$1
ref=$2
thr=$3

# compile homology_matrices
if [ ! -x homology_matrices ]; then
  g++ homology_matrices.cpp -o homology_matrices -O3
fi

# ======== Seungsoo's method: allowing duplicates ref2 but not ref1 ========

# identify best homologous bins
join -a 1 <(join -a 1 <(awk -v t=$thr 'BEGIN{OFS="\t"; 
                                   a = ""; b = 0} 
                             $1 > t {if (a != $2 || b != $3) 
                                         print $2, $3, $4, $5; 
                                      a = $2; b = $3}' $anns/$ref.$bsize.homology | sort -k1,1 -k2,2n -k3,3 -k4,4n) \
                      <(cut -f 1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1 -k2,2n) | \
            awk '{OFS="\t"; print $3, $4, $2+$5}' | sort -k1,1 -k2,2n) \
          <(cut -f1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1) | \
  awk '{OFS="\t"; print $3, $2+$4}' | sort -k1,1n -k2,2n > $anns/$ref.$bsize.homology.updated

# exclude isolated homologous bins and find neighbors (to exclude from nonhomologs)
./homology_matrices $anns/$ref.$bsize.bin_annotations $anns/$ref.$bsize.homology.updated 2 $anns/$ref.$bsize.homology_noisolated $anns/$ref.$bsize.homology_neighbors

# ======== Seungsoo's method (modified): allowing duplicates ref2 but not ref1 ========

# allowing duplicates in ref1: identify best homologous bins
join -a 1 <(join -a 1 <(awk -v t=$thr 'BEGIN{OFS="\t"}
                             $1 > t {print $2, $3, $4, $5}' $anns/$ref.$bsize.homology | sort -k1,1 -k2,2n -k3,3 -k4,4n) \
                      <(cut -f 1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1 -k2,2n) | \
            awk '{OFS="\t"; print $3, $4, $1, $2, $2+$5}' | sort -k1,1 -k2,2n) \
          <(cut -f1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1) | \
  awk '{OFS="\t"; print $5, $2+$6, $3, $4, $1, $2}' | sort -k1,1n -k2,2n > $anns/$ref.$bsize.homology.matrixbin

# allowing duplicates in ref1: exclude isolated homologous bins and find neighbors (to exclude from nonhomologs)
./homology_matrices $anns/$ref.$bsize.bin_annotations <(cut -f1,2 $anns/$ref.$bsize.homology.matrixbin) 2 $anns/$ref.$bsize.homology_noisolated $anns/$ref.$bsize.homology_neighbors

# allowing duplicates in ref1: get full bin info for best homologous bins, excluding isolated homologous bins
join -a 2 -j 1 \
  <(awk '{OFS="\t"; print $1 ":" $2, $3, $4, $5, $6}' $anns/$ref.$bsize.homology.matrixbin | sort -k1,1) \
  <(sort $anns/$ref.$bsize.homology_noisolated | tr '\t' ':' | sort -k1,1) | tr  ': ' '\t' | sort -k1,1n -k2,2n > \
  $anns/$ref.$bsize.homology_noisolated.matrixbin


# ======== New: allowing duplicates in ref1 & ref2 ========

# allowing duplicates in ref1: identify best homologous bins
join -a 1 <(join -a 1 <(awk -v t=$thr 'BEGIN{OFS="\t"}
                             $1 > t {print $2, $3, $4, $5}' $anns/$ref.$bsize.homology | sort -k1,1 -k2,2n -k3,3 -k4,4n) \
                      <(cut -f 1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1 -k2,2n) | \
            awk '{OFS="\t"; print $3, $4, $1, $2, $2+$5}' | sort -k1,1 -k2,2n) \
          <(cut -f1,2 $anns/$ref.$bsize.chr_annotations | sort -k1,1) | \
  awk '{OFS="\t"; print $5, $2+$6, $3, $4, $1, $2}' | sort -k1,1n -k2,2n > $anns/$ref.$bsize.allcombos.homology.matrixbin

# allowing duplicates in ref1: exclude isolated homologous bins and find neighbors (to exclude from nonhomologs)
./homology_matrices $anns/$ref.$bsize.bin_annotations <(cut -f1,2 $anns/$ref.$bsize.allcombos.homology.matrixbin) 2 $anns/$ref.$bsize.allcombos.homology_noisolated $anns/$ref.$bsize.allcombos.homology_neighbors

# allowing duplicates in ref1: get full bin info for best homologous bins, excluding isolated homologous bins
join -a 2 -j 1 \
	<(awk '{OFS="\t"; print $1 ":" $2, $3, $4, $5, $6}' $anns/$ref.$bsize.allcombos.homology.matrixbin | sort -k1,1) \
	<(sort $anns/$ref.$bsize.allcombos.homology_noisolated | tr '\t' ':' | sort -k1,1) | tr  ': ' '\t' | sort -k1,1n -k2,2n > \
	$anns/$ref.$bsize.allcombos.homology_noisolated.matrixbin