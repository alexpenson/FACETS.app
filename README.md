# FACETS.app

## Version 0.9.5 (2015.08.30)

Single wrapper script `facets.py` which takes a tab-delimited file `tumor_normal_pairs.txt` listing input BAM files.
It must contain columns Tumor_Sample_Barcode, t_bamfile & n_bamfile. 

* counts the base coverage over SNPs
* creates a join tumor/normal counts file
* runs facets

The same wrapper script also does:
* checking of output
* gene level calling
* annotation of maf files with local copy number information

## FACETS results from bam files

usage::

    ./facets.py runlsf tumor_normal_pairs.txt --cval 50

remaining arguments are passed to doFacets.R


## Gene Level calling

usage:

    ./facets.py calls tumor_normal_pairs.txt

a text file listing, for each input cncf.txt file and each IMPACT341 gene, the integer copy number, TCGA-style copy number (-2, -1, 0, 1, 2) as well as more advanced copy number calling, including CNLOH for example.


* Waterfall supression from A. Penson

The tumor depth and allele counts are normalized to remove the dependence of the log-ratio on normal depth, using a lowess fit. This is designed to remove the "waterfall" effect produced when degraded tumor DNA with small fragment size means that SNPs adjacent to target exons are consistently covered less well in the tumor than in the normal.


Thanks to Zeng Zhang for GetBaseCounts https://github.com/zengzheng123/GetBaseCounts
