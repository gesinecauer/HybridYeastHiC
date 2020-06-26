This is supplmented by Dejun Lin documenting how to reproduce the homology
interaction loci between S.cerevisiae and S.uvarum. In the scripts used below,
Go over the "module load" commands to adjust to whatever available in your environment
because I found some of the modules used isn't available any more.
1) Download the sequencing data from SRA and put them in data. Modify the
   refs.txt, matrices.txt and rawfiles.txt to only include the downloaded data.
2) Follow the following steps in the README: 
        0. cd scripts
        1. ./compile.sh (for c+00x) # compile the c++ programs for processings
        2. ./prepare_genomes.sh # correct for structural rearangement in the reference genomes 
        3. ./pipeline.sh -r ../refs.txt -t ../rawfiles.txt  # index the reference genome and trim adaptors from reads
        5. ./centromere_annotations.sh # find centromere locations?
        6. ./annotate_references.sh # annotate  centromeres and restriction enzyme fragments?
        8. ./masks.sh # mask rDNA repeats 
3) In the "scripts" direction, run "./homology.sh"
    to generate the homologous interaction bins for 32000 bp bin size. 
4) Map the reads and process the read pairs:
         cd scripts && ./pipeline.sh -p ../samples.txt
5) generate the matrix file:
         cd scripts && ./pipeline.sh -m ../matrices.txt
The scripts hardcoded some module load commands in them which doesn't work any more. I used my own software from homebrew so I simply
commented the module commands out.
