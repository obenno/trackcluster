#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:01 AM
# @Author  : Runsheng     
# @File    : cluster.py
"""
Generating similarity matrix using differ between gene models
make dendrogram for further plotting
The input is a list of
"""
# self import
from track import bigGenePred

# third part import
import scipy
import operator
from tracklist import wrapper_bedtools_intersect2, bigglist_to_bedfile, pandas_summary, add_subread_bigg, get_readall_bigg
## Mofidied function
from tracklist import wrapper_bedtools_intersect2_v2, bigglist_to_beddf, pandas_summary_v2, read_bigg_fromDataFrame, bigglist_pairwise_intersect
from utils import del_files
import pandas
from datetime import datetime # Add timestamp for message

def get_time():
    t=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return "[" + t + "]"

def flow_cluster(bigg_nano, bigg_gff, by="ratio_all", cutoff="auto", intronweight=0.5):

    bigg_nano.sort(key=operator.attrgetter("chromStart"))

    if by=="ratio_all":
        by1="ratio"
        by2="ratio_short"
    else:
        by1=by
        by2=by

    if cutoff=="auto":
        cutoff1=0.025
        cutoff2=0.001
    else: # expect cutoff as a tuple (0.05, 0.01)
        cutoff1, cutoff2= cutoff

    # hard code first filter of overalpping of 50 bp
    #bigg_l1=prefilter_smallexon(bigg_nano, bigg_gff, cutoff=50) # using default cutoff 0.95
    bigg_list=add_subread_bigg(bigg_gff+bigg_nano)
    print(len(bigg_list))
    ##bigg_df=pandas.DataFrame()
    bigg_tmp=[]
    for bigg in bigg_list:
        bigg=bigg.to_list()
        bigg_tmp.append(bigg)
    bigg_df=pandas.DataFrame(bigg_tmp)
    ## drop duplicated reads, as long as they
    ## have the same information on
    ## 0:chr, 1:start, 2:end, 5:strand, 9:blockCount, 10:blockSizes, 11:chromStarts
    bigg_df=bigg_df.drop_duplicates(subset=[0,1,2,5,9,10,11])
    ## Convert back to bigg class
    bigg_list = read_bigg_fromDataFrame(bigg_df)

    print("reducing bigg_list:")
    print("reduced bigg_list has " + str(len(bigg_list)) + " records")

    ## If locus has too many reads, will consume large mem and time
    ## subset original data for cal_distance and then merge the result
    if len(bigg_list) > 1000:
        print ("Locus has too many reads, preprocessing... " + get_time())
        a=range(0, len(bigg_list), len(bigg_list)/3) ## get index of subset
        ## subset
        bigg_list_tmp1 = bigg_list[a[0]:a[1]]
        bigg_list_tmp2 = bigg_list[a[1]:a[2]]
        if(len(a)==3):
            bigg_list_tmp3 = bigg_list[a[2]:]
            bigg_list_tmp4 = []
        else:
            bigg_list_tmp3 = bigg_list[a[2]:a[3]]
            bigg_list_tmp4 = bigg_list[a[3]:]
        ## Define subset flow function:
        def pre_sub_flow(bigg_list):
            bigg_list_tmp1 = bigg_list
            tmp_D1, tmp_bigg_list_by1 = cal_distance_v2(bigg_list_tmp1, intronweight=intronweight, by=by1)
            tmp_D, bigg_list_tmp1_l2 = filter_D(tmp_D1, tmp_bigg_list_by1, by=by1, cutoff=cutoff1)
            return bigg_list_tmp1_l2
        ## Process subset flow:
        bigg_list_sub1 = pre_sub_flow(bigg_list_tmp1)
        bigg_list_sub2 = pre_sub_flow(bigg_list_tmp2)
        bigg_list_sub3 = pre_sub_flow(bigg_list_tmp3)
        bigg_list_sub4 = pre_sub_flow(bigg_list_tmp4)
        ## Combine two subset:
        bigg_list = []
        bigg_list = add_subread_bigg(bigg_list_sub1+bigg_list_sub2+bigg_list_sub3+bigg_list_sub4)
        print("Finished preprocessing.")
    print("Start first time of cal_distance_v2:" + get_time())
    # can change filters
    D1, bigg_list_by1=cal_distance_v2(bigg_list, intronweight=intronweight, by=by1)
    ##print("len(bigg_list_by1)", len(bigg_list_by1))
    ##write_D(D1, bigg_list_by1, "new_d1.csv") # debug distance matrix
    _, bigg_l2=filter_D(D1, bigg_list_by1, by=by1, cutoff=cutoff1)
    print("Start second time of cal_distance_v2:" + get_time())
    D2, bigg_list_by2=cal_distance_v2(bigg_l2, intronweight=intronweight, by=by2)
    print("len(bigg_list_by2)", len(bigg_list_by2))
    D_remain, bigg_l3=filter_D(D2, bigg_list_by2, by=by2, cutoff=cutoff2)
    # add sanity check
    # the bigg_l3 subreads number together with read number+ bigg_l3=bigg_ll

    #missed_2=get_readall_bigg(bigg_list_by1)-get_readall_bigg(bigg_l2)
    #missed_3=get_readall_bigg(bigg_list_by2)-get_readall_bigg(bigg_l3)

    print "flow cluster", len(bigg_list),  len(bigg_l2), len(bigg_l3)

    ##D_remain=[]
    ##bigg_l3=[]
    return D_remain, bigg_l3


def getij(bigg_list):
    ij_list=[]
    for i in range(len(bigg_list)):
        for j in range(1, len(bigg_list)):
            ij_list.append((i,j))
    return ij_list


def get_pos_dic(bigg_list):
    pos_dic={}
    for n, bigg in enumerate(bigg_list):
        pos_dic[bigg.name]=n
    return pos_dic


def select_list(bigg_list, keep):
    # re_order D and bigg_list
    bigg_list_new=[]
    for i in keep:
        bigg_list_new.append(bigg_list[i])
    return bigg_list_new

def select_D(D, keep):
    D=D[keep,:]
    D=D[:,keep]

    return D

def prefilter_smallexon(bigg_list,bigg_list_gff, cutoff=50):
    """
    remove two kind of reads:

    1. not same strand as in current annotation
    2. with less than cutoff intersection with current annotation

    :param bigg_list:
    :param cutoff: at least 50bp intersection with current annotation
    :return: retained list
    """

    ## The anti-sense reads should not be removed
    ## in direct RNA sequencing, or novel anti-sence
    ## transcripts will be lost, thus filter 1 is disabled
    
    ## Intergenic transcripts will be lost due to filter 2,
    ## so the whole function will not be used in the main flow

    if len(bigg_list_gff)==0:
        return bigg_list

    ## # filter 1
    ## strand=bigg_list_gff[0].strand
    ## bigg_list_strand=[x for x in bigg_list if x.strand==strand]
    bigg_list_strand = bigg_list
    
    if len(bigg_list_strand)==0:
        return None

    # filter 2
    nano_exon, nano_intron=bigglist_to_bedfile(bigg_list_strand)
    gff_exon, gff_intron=bigglist_to_bedfile(bigg_list_gff)

    exonfile=wrapper_bedtools_intersect2(nano_exon, gff_exon)
    out_d=pandas_summary(exonfile)
    keep_name=set()
    for k, intersection in out_d.items():
        nano_name, gff_name=k
        if intersection > cutoff:
            keep_name.add(nano_name)

    bigg_list_new=[]
    for bigg in bigg_list:
        if bigg.name in keep_name:
            bigg_list_new.append(bigg)

    try:
        ### clean up
        del_files([exonfile, nano_intron, gff_intron])
        ### Also clean up exon fils
        del_files([nano_exon, gff_exon])
    except Exception as e:
        print("Cleanup in prefilter_smallexon is not successful: ", e)

    return bigg_list_new


def cal_distance(bigg_list, intronweight=0.5, by="ratio"):
    """
    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    #wkdir=set_tmp()
    #os.chdir(wkdir)

    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.get_exon()
        i.to_bedstr()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # get an pos combination and the name of bigg for each i
    # ij_list=getij(bigg_list)
    pos_dic=get_pos_dic(bigg_list)

    # flow begin
    file_exon, file_intron = bigglist_to_bedfile(bigg_list)

    exon_out=wrapper_bedtools_intersect2(file_exon, file_exon)
    exon_i=pandas_summary(exon_out)
    del_files([file_exon, exon_out])

    intron_out=wrapper_bedtools_intersect2(file_intron, file_intron)
    intron_i=pandas_summary(intron_out)
    del_files([file_intron, intron_out])

    for k, intersection in exon_i.items():
        name1, name2=k
        i=pos_dic[name1]
        j=pos_dic[name2]
        min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
        union = bigg_list[i].exonlen + bigg_list[j].exonlen - intersection
        # debug insanity
        if union <=0:
            print "exon", name1, name2, bigg_list[i].exonlen,  bigg_list[j].exonlen, union, intersection
        # debug over

        if by == "ratio":
            # exon could be 0?
            if min_length == 0:
                D_exon[i, j] = 1
            else:
                similar = float(intersection) / union
                D_exon[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 1
            else:
                D_exon[i, j] = 1 - float(intersection) / min_length

    for k, intersection in intron_i.items():
        name1, name2 = k
        i = pos_dic[name1]
        j = pos_dic[name2]
        min_length = min(bigg_list[i].intronlen, bigg_list[j].intronlen)
        union = bigg_list[i].intronlen + bigg_list[j].intronlen - intersection

        #### debug
        ## Intron union could equal to 0, for single exon transcripts
        if union <=0:
            print "intron",name1, name2, bigg_list[i].intronlen,  bigg_list[j].intronlen, union, intersection
        #### debug over

        if by == "ratio":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 1
            else:
                #print union
                similar = float(intersection) / union
                D_intron[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 1
            else:
                D_intron[i, j] = 1 - float(intersection) / min_length


    D=(D_exon+intronweight*D_intron)/float(1+intronweight)

    print("exon_out is ", exon_out)
    print("intron_out is ", intron_out)
    print("file_exon is ", file_exon)
    print("file_intron is ", file_intron)
    ## try:
    ##     # cleanup
    ##     del_files([exon_out, intron_out, file_exon, file_intron])
    ## except Exception as e:
    ##     print("Cleanup in flow_cluster is not successful: ", e)

    # debug:
    #print("D_exon",D_exon)
    #print("D_intron", D_intron)
    #print("D",D)

    #cleanup(remove_all=True)

    return D, bigg_list

def cal_distance_v2(bigg_list, intronweight=0.5, by="ratio"):
    """
    :param bigg_list:
    :param intronweight: if 0, do not cal the intron to save time
    :param by: used to cal the distance between two bigg object, can be "ratio", "ratio_short", "length", "length_short"
    :return: D: distance matrix
    """
    #wkdir=set_tmp()
    #os.chdir(wkdir)

    bigg_list.sort(key=operator.attrgetter("chromStart"))

    for i in bigg_list:
        i.get_exon()
        i.to_bedstr()

    length=len(bigg_list)
    D_exon=scipy.zeros([length, length])
    D_intron=scipy.zeros([length, length])

    # get an pos combination and the name of bigg for each i
    # ij_list=getij(bigg_list)
    pos_dic=get_pos_dic(bigg_list)
    ##print("pos_dic", pos_dic)
    # flow begin
    ## Generate bed6 dataframe for exon and intron
    ##exon_df, intron_df = bigglist_to_beddf(bigg_list)
    ##print(exon_df)
    ##exon_out=wrapper_bedtools_intersect2_v2(exon_df, exon_df)
    ##exon_i=pandas_summary_v2(exon_out)
    ## No need to del_files here, new function handels all
    ##del_files([file_exon, exon_out])

    ##intron_out=wrapper_bedtools_intersect2_v2(intron_df, intron_df)
    ##intron_i=pandas_summary_v2(intron_out)

    exon_i, intron_i = bigglist_pairwise_intersect(bigg_list)
    ##print(exon_i)
    print("len(exon_i) is", len(exon_i))
    for k, intersection in exon_i.items():
        ##print k, intersection
        name1, name2=k
        i=pos_dic[name1]
        j=pos_dic[name2]
        min_length = min(bigg_list[i].exonlen, bigg_list[j].exonlen)
        union = bigg_list[i].exonlen + bigg_list[j].exonlen - intersection
        # debug insanity
        if union <=0:
            print "exon", name1, name2, bigg_list[i].exonlen,  bigg_list[j].exonlen, union, intersection
        # debug over

        if by == "ratio":
            # exon could be 0?
            if min_length == 0:
                D_exon[i, j] = 1
            else:
                similar = float(intersection) / union
                D_exon[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_exon[i, j] = 1
            else:
                D_exon[i, j] = 1 - float(intersection) / min_length

    for k, intersection in intron_i.items():
        name1, name2 = k
        i = pos_dic[name1]
        j = pos_dic[name2]
        min_length = min(bigg_list[i].intronlen, bigg_list[j].intronlen)
        union = bigg_list[i].intronlen + bigg_list[j].intronlen - intersection

        #### debug
        ## Intron union could equal to 0, for single exon transcripts
        if union <0:
            print "intron",name1, name2, bigg_list[i].intronlen,  bigg_list[j].intronlen, union, intersection
        #### debug over

        if by == "ratio":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 1
            else:
                #print union
                similar = float(intersection) / union
                D_intron[i, j] = 1 - similar

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                D_intron[i, j] = 1
            else:
                D_intron[i, j] = 1 - float(intersection) / min_length


    D=(D_exon+intronweight*D_intron)/float(1+intronweight)
    print("len(D) is", len(D))
    # debug:
    #print("D_exon",D_exon)
    #print("D_intron", D_intron)
    #print("D",D)

    #cleanup(remove_all=True)

    return D, bigg_list



def write_D(D, bigg_list_new, outfile="./d.csv"):
    bigg_name=[x.name for x in bigg_list_new]

    if outfile is None:
        pass
    else:
        with open(outfile, "w") as fw:
            fw.write(",".join(bigg_name))
            fw.write("\n")
            for i in D:
                str_l=[str(x) for x in i]
                fw.write(",".join(str_l))
                fw.write("\n")


def filter_D(D, bigg_list, by="ratio", cutoff="auto", add_miss=False):

    """
    cutoff selection:
    learn from unc52, <0.025 in ratio_short
    return: index of the matrix that can be retained
    """
    if cutoff=="auto":
        if by=="ratio":
            cutoff=0.025
        elif by=="ratio_short":
            cutoff=0.001 # may need to add to 0.01
        elif by=="length" or by=="length_short":
            cutoff=100

    else: # expect two numbers for the cutoff
        cutoff=cutoff

    # hard code a cutoff for sw score of SL
    sw_score=11

    fullset=set(range(len(D)))
    drop=set()

    # first
    #print "first"
    #print len(get_readall_bigg(bigg_list))

    #for bigg in bigg_list:
    #    bigg.get_exon()


    # same list
    ij_list=getij(D)

    # add a filter function for D, if the distance between exon is 0, merge the small one with
    # unless need to parer the intronD and exonD separately, or else the filter should be outer function
    for i,j in ij_list:
        if D[i,j]<cutoff:
            if i==j:
                pass
            else:
                if by=="ratio":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen:
                        drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                    elif bigg_list[i].exonlen>bigg_list[j].exonlen:
                        drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)
                        # to add a subread add here

                if by=="ratio_short":
                    if bigg_list[i].exonlen<bigg_list[j].exonlen and bigg_list[i].score<sw_score:
                        drop.add(i)
                        bigg_list[j].subread.add(bigg_list[i].name)
                        bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)

                    elif bigg_list[i].exonlen==bigg_list[j].exonlen:
                        if i<j: # retain the later one
                            drop.add(i)
                            bigg_list[j].subread.add(bigg_list[i].name)
                            bigg_list[j].subread=bigg_list[j].subread.union(bigg_list[i].subread)
                        else:
                            drop.add(j)
                            bigg_list[i].subread.add(bigg_list[j].name)
                            bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

                    elif bigg_list[i].exonlen>bigg_list[j].exonlen and bigg_list[j].score<sw_score:
                        drop.add(j)
                        bigg_list[i].subread.add(bigg_list[j].name)
                        bigg_list[i].subread=bigg_list[i].subread.union(bigg_list[j].subread)

    keep=fullset-drop
    # change the default score of gene, no need to add
    for n, bigg in enumerate(bigg_list):
        if bigg.ttype=="isoform_anno":
            keep.add(n)


    # re_order D and bigg_list
    keepl = sorted(list(keep))
    bigg_list_new = select_list(bigg_list, keepl)


    #----------------------------------------#
    #### sanity check for missed ones
    ## collect the missed ones
    pos_dic=get_pos_dic(bigg_list)
    missed_name=get_readall_bigg(bigg_list)-get_readall_bigg(bigg_list_new)

    #print missed_name
    #print "inside"
    #print len(get_readall_bigg(bigg_list_new))
    if add_miss:
        if len(missed_name)>0:
            #print "{} missing bigg found, added back but may affect the isoforms".format(len(missed_name))
            missed_num=set()
            for k in missed_name:
                missed_num.add(pos_dic[k])

            keep=keep.union(missed_num)

            keepl=sorted(list(keep))
            bigg_list_new=select_list(bigg_list, keepl)

            #### end of sanity check


    D = select_D(D, keepl)

    return D, bigg_list_new


