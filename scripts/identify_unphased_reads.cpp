/* identify_unphased_reads.cpp
 * takes mapped reads, determines which multiple alignments are due to ambiguity in phasing
 * Gesine Cauer
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <utility>
#include <unordered_set>
#include <regex>


template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

using namespace std;

// function to find end position of read alignment from CIGAR string
int readend(string cig) {
  stringstream ss(cig);
  int pos = 0;
  int nts;
  char type;
  while (ss >> nts, ss >> type, !ss.eof()) {
    if (type=='M' || type=='D') {
      pos += nts;
    }
  }
  return pos;
}

// data type for read pairs
struct readpair {
  bool strand1;
  string chrom1;
  long coord1;
  bool strand2;
  string chrom2;
  long coord2;
  bool operator==(const readpair & rp) const {
    return (rp.chrom1.compare(chrom1) == 0 &&
      rp.chrom2.compare(chrom2) == 0 &&
      rp.strand1 == strand1 &&
      rp.strand2 == strand2 &&
      rp.coord1 == coord1 &&
      rp.coord2 == coord2);
  }
};

// hash function for deduplicating read pairs
namespace std {
  template <> struct hash<readpair> {
    inline size_t operator()(const readpair & a) const {
      size_t seed = 0;
      hash_combine(seed, a.chrom1);
      hash_combine(seed, a.chrom2);
      hash_combine(seed, a.strand1);
      hash_combine(seed, a.strand2);
      hash_combine(seed, a.coord1);
      hash_combine(seed, a.coord2);
      return seed;
    }
  };
}

int main(int argc, char** argv) {
  // parameters

  if (argc < 8) {
    cout << "USAGE: ./identify_unphased_reads RESTRICTION_SITE GENOME_DIGEST MAPQ BINSIZE HOMOLOGY_FILE READ_FILE READPAIR_OUT STAT_OUT\n";
    return 1;
  }

  // build map (binary search tree) of restriction fragments, with a vector for each chromosome
  map<string, vector<pair<long,long> > > restr_map;

  // build map of homologous bins
  map<string, vector<string > > homologous_map;

  // arg 1: restriction site
  string resite;
  stringstream ss(argv[1]);
  ss >> resite;

  // arg 2: restriction digest file
  ifstream genome_digest(argv[2]);

  // arg 3: mapq filter
  int mapq_filt;
  stringstream ss2(argv[3]);
  ss2 >> mapq_filt;

  // arg 4: bin size
  stringstream bsizearg(argv[4]);
  long bsize;
  bsizearg >> bsize;

  // arg 5: homology file
  ifstream homologyfile(argv[5]);

  // arg 6: read file
  ifstream readfile(argv[6]);
  // arg 7: output file for processed reads
  ofstream readout(argv[7]);
  // arg 8: output file for stats
  ofstream statout(argv[8]);

  // alignment score filter
  int as_filt = 10;

  string line;
  string chr;
  long pos1, pos2;

  // build restriction map
  //int counter = 0;
  while (getline(genome_digest,line)) {
    // parse line of restriction digest
    stringstream ss(line);
    ss >> chr;
    ss >> pos1;
    ss >> pos2;
    // if new chr, add new vector
    if (restr_map.find(chr) == restr_map.end()) {  // FIXME
      //cout << chr << " before new vector = " << int(restr_map[chr].size()) << ".\n" << flush;
      restr_map.insert(make_pair(chr,vector<pair<long,long> >()));
      //cout << chr << " after new vector = " << int(restr_map[chr].size()) << ".\n" << flush;
    }
    // add new element to vector for each restriction fragment
    //cout << chr << " before new element = " << int(restr_map[chr].size()) << ".\n" << flush;
    restr_map[chr].push_back(make_pair(pos1+1,pos2));
    //cout << chr << " after new element = " << int(restr_map[chr].size()) << ".\n\n" << flush;
    //counter++; if (counter > 20) {return 0;}
  }
  genome_digest.close();

  int bin;
  string chr1, chr2;
  string chrbin1, chrbin2;

  // build homologous bins map
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

  string rname1;
  int flag;
  long pos;
  int mapq;
  string cig;
  bool mult_align, secondary;
  int as, xs;
  size_t aspos, xspos;
  int phased;
  string homo;

  // data for primary alignment
  int primary_wins;
  int primary_flag;
  string primary_chr;
  long primary_pos;
  int primary_mapq;
  string primary_cig;
  string primary_homo;
  int primary_as;

  // for splitting species & chrom number
  regex species_chrom_re ("^([^_]+)_([^_]+)$");

  // counters for read pair categories
  long num_total = 0;
  long num_mult_align = 0;
  long num_phased = 0;
  long num_unphased = 0;
  long num_invalid = 0;
  long num_no_homo_bin_phased = 0;
  long num_primary_wins = 0;
  long num_no_homo_bin = 0;
  long num_nohomolog_loci = 0;
  long num_nohomolog_chr = 0;

  while (getline(readfile,line)) {
    // parse line from read 1 file
    stringstream ss(line);
    ss >> rname1;
    ss >> flag;
    if (flag >= 256) {
      flag = flag - 256;
      secondary = true;
    }
    else {
      secondary = false;
    }
    ss >> chr;
    ss >> pos;
    ss >> mapq;
    ss >> cig;

    // parse alignment scores for primary (and possibly secondary) alignment
    as = 0;
    aspos = line.find("AS:i:");
    if (aspos != string::npos) {
      stringstream ss(line.substr(aspos+5));
      ss >> as;
    }
    xs = -40;
    xspos = line.find("XS:i:");
    if (xspos != string::npos) {
      stringstream ss(line.substr(xspos+5));
      ss >> xs;
      mult_align = true;
    }
    else {
      mult_align = false;
    }

    if (cig == "*") {
      homo = "NA";
    }
    else {
      // find end positions of read alignments
      long posend = pos - 1 + readend(cig);

      // find restriction fragment ends
      long re1start, re1end;
      // binary search for restriction fragment 1 start
      int upper = restr_map[chr].size();
      int lower = 0;
      while (upper != lower) {
        int test = (upper+lower)/2;
        if (restr_map[chr][test].first <= pos) {
          lower = test;
        }
        if (restr_map[chr][test].second >= pos) {
          upper = test;
        }
      }
      re1start = restr_map[chr][lower].first;
      // search for restriction fragment 1 end starting from fragment 1 start
      while (restr_map[chr][upper].second + resite.length() < posend) {
        upper++;
      }
      re1end = restr_map[chr][upper].second + resite.length();
      // use restricton fragment location to determine which bin of this chromosome the read falls on
      long re = (flag == 0) ? re1end : re1start;
      long chrbin = re/bsize + 0.5;
      // define homologos bin
      homo = chr + "." + to_string(chrbin);

    }

    // if no multiple alignments exist, read is phased (maps to only one homolog)
    if (!mult_align) {
      if ((mapq < mapq_filt) || (as-xs < as_filt)) {
        phased = -1;
        num_invalid++;
      }
      else {
        // if homologous bin for alignment can't be identified
        if (homologous_map.find(homo) == homologous_map.end()) {
          // read is invalid: it can't be determined whether read maps to homologous loci
          phased = -1;
          num_invalid++;
          num_no_homo_bin_phased++;
        }
        else { // read is phased: read is valid and maps only to homolog of primary alignment
          phased = 1;
          num_phased++;
        }
      }
      num_total++;
    }
    // Process multiple alignments
    else {
      // if this is the primary alignment
      if (!secondary) {
        num_mult_align++;
        num_total++;
        // if primary MAPQ fails or the difference between primary vs secondary alignment scores fails threshold
        if ((mapq < mapq_filt) || (as-xs < as_filt)) {
          primary_wins = 0;
          // note primary alignment data for comparsion with secondary alignment data
          primary_flag = flag;
          primary_chr = chr;
          primary_pos = pos;
          primary_mapq = mapq;
          primary_cig = cig;
          primary_homo = homo;
          primary_as = as;
          // proceed to secondary alignment
          continue;
        }
        else { // read maps best to homolog of primary alignment
          // if homologous bin for primary alignment can't be identified
          if (homologous_map.find(homo) == homologous_map.end()) {
            // read is invalid: it can't be determined whether primary read maps to homologous loci
            primary_wins = -1;
            phased = -1;
            num_invalid++;
            num_no_homo_bin_phased++;
          }
          else { // read is phased: read is valid and maps best to homolog of primary alignment
            primary_wins = 1;
            phased = 1;
            num_phased++;
            num_primary_wins++;
          }
        }
      }
      else { // this is the secondary alignment
        // if primary MAPQ passes & the difference between primary vs secondary alignment scores passes threshold
        if (primary_wins != 0) {
          // discard data for secondary alignment (data for primary alignment was already recorded)
          continue;
        }
        // if primary vs secondary species is different and chromosome num is the same
        string species  = regex_replace (chr,species_chrom_re,"$1");
        string chromnum  = regex_replace (chr,species_chrom_re,"$2");
        string primary_species  = regex_replace (primary_chr,species_chrom_re,"$1");
        string primary_chromnum  = regex_replace (primary_chr,species_chrom_re,"$2");
        //if ((species != primary_species) && (chromnum == primary_chromnum)) {
        if ((species != primary_species)) {
          auto iter = homologous_map.find(homo);
          auto primary_iter = homologous_map.find(primary_homo);
          // if homologous bin for either alignment can't be identified
          if ((iter == homologous_map.end()) || (primary_iter == homologous_map.end()) ) {
            // read is invalid: it can't be determined whether reads map to homologous loci
            phased = -1;
            num_invalid++;
            num_no_homo_bin++;
          }
          else {
            // if primary & secondary alignment loci are in homologous bins
            if (find(iter->second.begin(), iter->second.end(), primary_homo) != iter->second.end() ) {
              // read is ambiguous: maps to 2 homologous loci
              phased = 0;
              num_unphased++;
              // record primary alignment
              flag = primary_flag;
              chr = primary_chr;
              pos = primary_pos;
              mapq = primary_mapq;
              cig = primary_cig;
            }
            else { // primary & secondary alignments aren't homologous loci
              // read is invalid: maps to 2 non-homologous loci
              phased = -1;
              num_invalid++;
              num_nohomolog_loci++;
            }
          }

        }
        else { // primary & secondary alignments aren't homologous loci (species/chrom mismatch)
          // read is invalid: maps to 2 non-homologous loci
          phased = -1;
          num_invalid++;
          num_nohomolog_chr++;
        }
      }
    }

    // print read pair info
    readout << rname1 << '\t' << flag << '\t' << chr << '\t' << pos << '\t' << mapq << '\t' << cig << '\t' << phased << '\n';
  }
  readout.close();

  // print read pair statistics
  statout << "num_total\t" << num_total << '\n';
  statout << "num_total_valid\t" << num_total - num_invalid << '\n';
  statout << "num_mult_align\t" << num_mult_align << '\n';
  statout << "num_phased\t" << num_phased << '\n';
  statout << "num_unphased\t" << num_unphased << '\n';
  statout << "num_invalid\t" << num_invalid << '\n';
  statout << "num_no_homo_bin_phased\t" << num_no_homo_bin_phased << '\n';
  statout << "num_primary_wins\t" << num_primary_wins << '\n';
  statout << "num_no_homo_bin\t" << num_no_homo_bin << '\n';
  statout << "num_nohomolog_loci\t" << num_nohomolog_loci << '\n';
  statout << "num_nohomolog_chr\t" << num_nohomolog_chr << '\n';
  statout.close();
}
