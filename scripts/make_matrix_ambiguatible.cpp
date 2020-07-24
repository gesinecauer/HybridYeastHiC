/* make_matrix.cpp
 * takes assigned read pairs and creates matrix of interaction counts (long form)
 * Seungsoo Kim
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;

int main(int argc, char** argv) {

  // arg 1: chromosome lengths
  ifstream chromlengths(argv[1]);

  // arg 2: bin size
  stringstream bsizearg(argv[2]);
  long bsize;
  bsizearg >> bsize;

  // arg 3: homology file
  ifstream homologyfile(argv[3]);

  // arg 4: mask file
  ifstream maskfile(argv[4]);

  // arg 5: assigned reads file
  ifstream countfile(argv[5]);

  // arg 6: minimum average value across row
  stringstream minrowavg_arg(argv[6]);
  double minrowavg;
  minrowavg_arg >> minrowavg;

  // arg 7: output row sum file
  ofstream outrowsumf(argv[7]);

  // arg 8: output matrix
  ofstream outmatrf(argv[8]);

  // build homologous bins map
  map<string, vector<string > > homologous_map;
  string line;
  int bin;
  string chr1, chr2;
  string chrbin1, chrbin2;
  while (getline(homologyfile,line)) {
    // parse line of homologous bins
    stringstream ss(line);
    ss >> bin;
    ss >> bin;
    ss >> chr1;
    ss >> chrbin1;
    ss >> chr2;
    ss >> chrbin2;
    string homo1 = chr1 + "." + chrbin1;
    string homo2 = chr2 + "." + chrbin2;
    // if new homologous bin, add new vector
    if (homologous_map.find(homo1) == homologous_map.end()) { // FIXME
      homologous_map.insert(make_pair(homo1,vector<string >()));
    }
    if (homologous_map.find(homo2) == homologous_map.end()) {  // FIXME
      homologous_map.insert(make_pair(homo2,vector<string >()));
    }
    // add new element to vector for each homologous bin
    homologous_map[homo1].push_back(homo2);
    homologous_map[homo2].push_back(homo1);
  }
  homologyfile.close();

  // parse chromosome lengths
  string line;
  vector<string> chrs;
  vector<long> chrlens;
  string chr;
  long len;
  while (getline(chromlengths,line)) {
    stringstream ss(line);
    ss >> chr;
    ss >> len;
    chrs.push_back(chr);
    chrlens.push_back(len);
  }
  chromlengths.close();

  // calculate total number of bins, and cumulative number of bins for each chromosome
  long nbin=0;
  map<string,int> cumbins;
  for (int i = 0; i < chrs.size(); i++) {
    cumbins.insert(make_pair(chrs[i],nbin));
    nbin += chrlens[i]/bsize+1;
  }
  long minrowsum = (long) (minrowavg * nbin);

  // initialize matrices
  long* outmatr = new long[nbin*nbin];
  double* normmatr = new double[nbin*nbin];
  for (long i = 0; i < nbin*nbin; i++) {
    outmatr[i]=0;
    normmatr[i]=nan("");
  }

  // parse mask file
  vector<string> mask_chrs;
  vector<long> mask_starts;
  vector<long> mask_ends;
  long start;
  long end;
  while (getline(maskfile,line)) {
    stringstream ss(line);
    ss >> chr;
    ss >> start;
    ss >> end;
    mask_chrs.push_back(chr);
    mask_starts.push_back(start);
    mask_ends.push_back(end);
  }
  maskfile.close();

  string chr1, chr2;
  long re1, re2;
  long pos1, pos2;
  long pend1, pend2;
  long map1, map2;
  int dir1, dir2;
  int bin1, bin2;
  while (getline(countfile,line)) {
    // parse one read pair
    stringstream ss(line);
    ss >> chr1;
    ss >> pos1;
    ss >> pend1;
    ss >> map1;
    ss >> re1;
    ss >> dir1;
    ss >> chr2;
    ss >> pos2;
    ss >> pend2;
    ss >> map2;
    ss >> re2;
    ss >> dir2;
    // skip pair if either read maps in masked region
    bool skippair = false;
    for (int i = 0; i < mask_chrs.size(); i++) {
      if ((mask_chrs[i].compare(chr1)==0) && !((mask_starts[i] > pend1) || (mask_ends[i] < pos1))) {
        skippair = true;
        break;
      }
      if ((mask_chrs[i].compare(chr2)==0) && !((mask_starts[i] > pend2) || (mask_ends[i] < pos2))) {
        skippair = true;
        break;
      }
    }
    // add to matrix count (symmetrically) if both reads are from chromosomes included in the matrix and not in same bin
    if (!skippair && (cumbins.count(chr1) == 1) && (cumbins.count(chr2) == 1)) {
      bin1 = cumbins[chr1] + re1/bsize;
      bin2 = cumbins[chr2] + re2/bsize;
      if (bin1 != bin2) {
        outmatr[nbin*bin2 + bin1]++;
        outmatr[nbin*bin1 + bin2]++;
      }
    }
  }

  // calculate row sums
  long* rowsums = new long[nbin];
  for (int i = 0; i < nbin; i++) {
    rowsums[i] = 0;
    for (int j = 0; j < nbin; j++) {
      rowsums[i] += outmatr[nbin*i + j];
    }
    outrowsumf << i << '\t' << rowsums[i] << '\n';
  }

  // calculate total sum, excluding rows that don't meet minimum
  long allsum = 0;
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      allsum += rowsums[i];
    }
  }

  // calculate coverage-normalized matrix
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      for (int j = 0; j < nbin; j++) {
        if (rowsums[j] >= minrowsum) {
          normmatr[nbin*i + j] = ((double)outmatr[nbin*i + j]*allsum)/((double)rowsums[i]*rowsums[j]);
        }
      }
    }
  }

  // print matrix, one line per entry of matrix
  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      outmatrf << i << '\t' << j << '\t' << outmatr[nbin*i + j] << '\t' << normmatr[nbin*i + j] << '\n';
    }
  }

  delete[] outmatr;
  delete[] normmatr;
  delete[] rowsums;
  countfile.close();
  outrowsumf.close();
  outmatrf.close();
}
