# CircRNA Benchmarking

## nf-core/circrna Info
Version:    
nf-core/circrna v0.0.1dev 

Branch:   
https://github.com/nf-core/circrna/tree/new_benchmarking   


## Directory Structure
After all four runs of `nf-core/circrna` on both data types of both datasets, 
move the pipeline output of the BJS detection (`bsj_detection`) into eiter `polya/bed`, or `total/bed`, based on 
what datatype and dataset was used. Make sure the pipeline was executed using the `blacklist`.

Extract `bam` files located in `tools/unified/start/2nd_pass` into a shared directory `polya/bam` or `total/bam`
based on datatype and dataset.


The data directory structure should then be:

```
data
├── blacklist
│   └── hg38-blacklist.v2.bed
├── DEEP
│   ├── polya
│   │   ├── bam
│   │   │   ├── 1_polyA_41_Hf01_LiHe_L005.Aligned.out.bam
│   │   │   └── ...
│   │   ├── bed
│   │   │   ├── circexplorer2
│   │   │   │   ├── blacklist
│   │   │   │   │   ├── 1_polyA_41_Hf01_LiHe_L005.circexplorer2.blacklist.bed
│   │   │   │   │   └── ...
│   │   │   ├── ciriquant
│   │   │   │   ├── blacklist
│   │   │   │   │   ├── 1_polyA_41_Hf01_LiHe_L005.ciriquant.blacklist.bed
│   │   │   │   │   └── ...
│   │   │   ├── dcc
│   │   │   │   ├── blacklist
│   │   │   │   │   ├── 1_polyA_41_Hf01_LiHe_L005.dcc.blacklist.bed
│   │   │   │   │   └── ...
│   │   │   ├── find_circ
│   │   │   │   ├── blacklist
│   │   │   │   │   ├── 1_polyA_41_Hf01_LiHe_L005.find_circ.blacklist.bed
│   │   │   │   │   └── ...
│   │   │   └── segemehl
│   │   │       ├── blacklist
│   │   │       │   ├── ...
│   │   │       │   └── 9_polyA_51_Hf04_BlEM_L002.segemehl.unified.blacklist.bed
│   └── total
│       ├── bed
│       │   ├── circexplorer2
│       │   │   ├── blacklist
│       │   │   │   └── 9_polyA_51_Hf04_BlEM_L002.segemehl.unified.blacklist.bed
│       │   ├── ciriquant
│       │   ├── ...
│       │   └── segemehl
├── GSE138734
│   ├── polya
│   │   └── ...
│   └── total
│       └── ...
└── rrna_gtf
    └── ensembl_rRNA.gtf
```

## Additional Files

Blacklist: https://github.com/Boyle-Lab/Blacklist/     
Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

rRNA GTF: https://github.com/zxl124/rRNA_gtfs


## Run
After having set up everything, run
```
cd src/
bash run_analysis.sh
```
