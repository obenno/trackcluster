#! /usr/bin/env python

import re, os, sys
from multiprocessing import Pool
from trackcluster.tracklist import read_bigg, write_bigg, bigglist_to_bedfile, read_bigg_fromDataFrame
from collections import OrderedDict
from trackcluster.batch import *
import pandas
import traceback
import argparse

parser = argparse.ArgumentParser(description='Main work flow of trackcluster')
parser.add_argument('--inter', dest='interFile',
                    help='minimap2 and annotation intersect file, named as nano_inter_ref.bgp')
parser.add_argument('--annotation', dest='annoFile',
                    help='genome annotation, type column must be "isoform_anno"')
parser.add_argument('--cpu', dest='threads', default = 4, type = int,
                    help='Number of threads [4]')
parser.add_argument('--out', dest='outFile',
                    help='output file')

if(len(sys.argv)==1):
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
## Give a better name of intersect list
nano_trans_BEDintersect = pandas.read_csv(args.interFile, sep = "\t", header = None)

def group_bigg_by_gene(bigglist):
    gene_bigg=OrderedDict()

    for bigg in bigglist:
        try:
            gene_bigg[bigg.geneName].append(bigg)
        except KeyError:
            gene_bigg[bigg.geneName]=[]
            gene_bigg[bigg.geneName].append(bigg)
    return gene_bigg

#gene_nano=group_bigg_by_gene(nano_bigg_new)
anno_bigg=read_bigg(args.annoFile)
gene_anno=group_bigg_by_gene(anno_bigg)
print("length of gene_anno.keys is ", len(gene_anno.keys()))
errors_ll=[]


def process_one_subsample_try(key, intronweight=0.5, by="ratio_all", full=False):

    print("Processing: " + key)

    nano_df = nano_trans_BEDintersect.loc[nano_trans_BEDintersect[37]==key, 0:19]
    bigg_nano_raw = read_bigg_fromDataFrame(nano_df)

    bigg_gff = gene_anno[key]

    print("Read bigg_gff and bigg_nano_raw.")

    ## bigg_nano = prefilter_smallexon(bigg_nano_raw, bigg_gff, cutoff=0)
    ## Set cutoff=0 to disable filter 2 of prefilter_smallexon
    ## prefilter_smallexon filtration will remove intergenic and anti-sense transcripts
    ## so it was disabled from main flow.
    bigg_nano = bigg_nano_raw

    if bigg_nano is None:
        print("bigg_nano is None")
        return 0

    try:
        D, bigg_nano_new = flow_cluster(bigg_nano, bigg_gff, by, intronweight=intronweight)
        bigg_nano_new = add_subread_bigg(bigg_nano_new)
        print("Finished flow_cluster.")
        ##print "len(bigg_nano_new) is ", len(bigg_nano_new)
        ### save nessary files
        ## for bigg in bigg_nano_new:
        ##     bigg.write_subread()
        ## print("Finished write_subread().")
        ##print len(bigg_nano_new)
        ##bigg_count_write(bigg_nano_new, out=biggout)
    except Exception as e:
        bigg_nano_new = []
        print("Exception occured when processing: " + key)
        errors_ll.append((key,e))

    return(bigg_nano_new)


## Test main flow without multiprocessing
bigg_out = []
## clustered_nano_parallel = []
## for gene in gene_anno.keys():
##     bigg_nano_new = process_one_subsample_try(gene)
##     print "type is", type(bigg_nano_new)
##     clustered_nano_parallel.append(bigg_nano_new)

clustered_nano_parallel = []
p = Pool(args.threads) ## Define how many threads to use
clustered_nano_parallel = p.map(process_one_subsample_try, gene_anno.keys())

print("clustered_nano_parallel length is "+str(len(clustered_nano_parallel)))
print("len(errors_ll) is ", len(errors_ll))

##print(clustered_nano_parallel)

for bigg_list in clustered_nano_parallel:
    ## i is list of bigglist
    bigg_count_write(bigg_list)
    for bigg_one in bigg_list:
        bigg_out.append(bigg_one.to_str())

with open(args.outFile, "w") as fw:
    fw.write("\n".join(bigg_out))

print("clustered_nano length is "+str(len(bigg_out)))

#bigg_count_write(clustered_nano, out= args.outFile, append=True)
##print("len(errors_ll) is ", len(errors_ll))