# TrackCluster
Trackcluster is an isoform calling and quantification pipeline for Nanopore direct-RNA long reads.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Walkthrough](#walkthrough)


## <a name="overview"></a>Overview
A pipeline for reference-based identification of isoform calling using Nanopore direct-RNA long reads. 
This pipeline was designed to use **only** long and nosisy reads to make a valid transcriptome. 
An indicator for the intact 5' could be very helpful to the pipeline, i.e, the splicing leader in the mRNA of nematodes. 

It is recommended to combine all samples together to generate an new transcriptome reference. 
After this process, the expression of isoforms in each sample can be fetched by providing an "name:sample" table. 

The output format for this pipeline is ["bigGenePred"](https://github.com/Runsheng/trackcluster/blob/master/script/bigGenePred.as). 

## <a name="requirements"></a>Requirements

1. python 2.7.10+ (developed under python 2.7.12)
2. python modules: pysam, pandas, numpy
3. samtools V2.0+ , bedtools V2.24+
4. [minimap2](https://github.com/lh3/minimap2)

## Recommendations
1. UCSC Kent source tree (for generating binary track)
2. scipy for hclust function
3. matplotlib and pylab for plotting

## <a name="walkthrough"></a>Walkthrough

An walkthrough example can be found in the [ipython notebook file](https://github.com/Runsheng/trackcluster/blob/master/trackcluster_run_example.ipynb). 

## Updates:

- Add main work flow script, read all reads information from intersection result
- Disabled tempfiles and bedtools modules, fixed "No spaces left" issue
- Fixed a bug found when calculating distance matrix
- Disabled frefilter function
- Multiple mapped reads are now supported when concatenating long reads
- Performance improved especially for large genome and loci with very high long reads coverage 

## Main flow script

```
usage: trackcluster_mainflow.py [-h] [--inter INTERFILE]
                                [--annotation ANNOFILE] [--cpu THREADS]
                                [--out OUTFILE]

Main work flow of trackcluster

optional arguments:
  -h, --help            show this help message and exit
  --inter INTERFILE     minimap2 and annotation intersect file, named as
                        nano_inter_ref.bgp
  --annotation ANNOFILE
                        genome annotation, type column must be "isoform_anno"
  --cpu THREADS         Number of threads [4]
  --out OUTFILE         output file

```
