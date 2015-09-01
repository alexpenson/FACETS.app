# FACETS.app

## Version 0.9.5 (2015.08.30)

Wrapper script which takes a tab-delimited file `tumor_normal_pairs.txt` listing input BAM files.
It must contain columns Tumor_Sample_Barcode, t_bamfile & n_bamfile. 

* counts the base coverage over SNPs
* creates a join tumor/normal counts file
* runs facets
* gene level calling
* annotation of maf files with local copy number information


usage::
```bash
usage:: facets.py runlsf tumor_normal_pairs.txt
remaining arguments are passed to doFacets.R
```

## Updates

* Gene Level calling: get_gene_level_calls.R

usage:

    ./get_gene_level_calls.R output_file.txt *_cncf.txt

output:

a text file listing, for each input cncf.txt file and each IMPACT341 gene, the integer copy number, TCGA-style copy number (-2, -1, 0, 1, 2) as well as more advanced copy number calling, including CNLOH for example.

* dualFacets.R -  runs FACETS twice with two sets of input parameters dipLogR from the first iteration is used for the second

usage:

    ./dualFacets.R --purity_cval 300 --cval 100 countsMerged.dat.gz

output:

output files for both iterations are retained

* run.sh no longer recounts if counts files already present; can call it to re-run facets instead of calling doFacets.R directly. 

* ```--ndepth``` (```-n```) option added to ```doFacets.R``` added: reduction to 25 recommended by Venkat for high sensitivity

* Waterfall supression from A. Penson

The tumor depth and allele counts are normalized to remove the dependence of the log-ratio on normal depth, using a lowess fit. This is designed to remove the "waterfall" effect produced when degraded tumor DNA with small fragment size means that SNPs adjacent to target exons are consistently covered less well in the tumor than in the normal.


Thanks to Zeng Zhang for GetBaseCounts https://github.com/zengzheng123/GetBaseCounts
