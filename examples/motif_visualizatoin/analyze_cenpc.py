# -*- coding: utf-8 -*-
"""
This script just download all the sequences from our seqs.cvs file - the histone database.

"""
from __future__ import division
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import cPickle as pickle
from Bio import SeqIO

from Bio import motifs
from Bio.SeqUtils.CheckSum import seguid

import uuid
from Bio.SeqFeature import SeqFeature, FeatureLocation

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess
import pylab
import pandas as pd
import networkx as nx
from StringIO import StringIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
#from pylab import *
import itertools
import time
from multiprocessing import Pool, current_process
import cPickle as pickle

sys.path.insert(0,os.path.join(os.path.dirname(__file__),'MYSOFT/ALIGNMENT_TOOLS/'))
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'libs/ete/'))


# from hist_ss import get_hist_ss
from L_aln_tools import taxo_seq_architecture,features_via_hmm, get_prot_gbrec_by_gis, muscle_aln,taxo_msa,gen_fake_msa
from L_aln2html import aln2html
from L_taxonomy_tools import get_taxid_from_gbrec, gbrec_to_clickhtml

Entrez.email = "alexey.shaytan@nih.gov" 


def main():
#Loading sequences through gis - this is one option
    if(not os.path.isfile('int_data/fasta_dict.p')):
        gis_df=pd.read_csv('inp_data/cenpc_gis.csv',dtype=object)
        gis=list(gis_df['gi'])
        print gis
        gis=[i for i in gis if re.match('\d+',i)]
        fasta_dict=get_prot_gbrec_by_gis(gis)
        pickle.dump( fasta_dict, open("int_data/fasta_dict.p", "wb" ) )
    fasta_dict=pickle.load(open("int_data/fasta_dict.p", "r" ))

    #Uncomment below to get a manual list that we blasted.
    handle = open("inp_data/big_set.gb", "rU")
    records = SeqIO.parse(handle, "gb")
    non_red_rec=[]
    recidset=set()
    for record in records:
        if record.id in recidset:
            print "Ignoring %s" % record.id
            continue
        recidset.add(record.id)
        non_red_rec.append(record)
    handle.close()

    fasta_dict = SeqIO.to_dict(non_red_rec)

    print fasta_dict


    #Let's try to make alignments
    # msa=muscle_aln(fasta_dict.values())
    # print msa
    # taxids=map(get_taxid_from_gbrec,fasta_dict.values())
    # print taxids
    # aln2html(gen_fake_msa(fasta_dict.values()),'out_data/aln.html')
    # taxo_msa('out_data/aln.svg',taxids,'',fasta_dict.values())

    #Let's output list of species with clickable links in html.

    gbrec_to_clickhtml(fasta_dict.values(),'out_data/seqrecdata.html')
    #Let's start playing with feature annotations
    seqreclist=fasta_dict.values()
    taxids=map(get_taxid_from_gbrec,seqreclist)
    #now we will reset the feature annotation, although we could just append in future
    
    #We will be searching for at-hook also


    def at_hook(seq):
        hooks=["RGR",
        "PGR",
        "GGR",
        "MGR"]
        features=[]
        for h in hooks:
            p = re.compile(h)
            for m in p.finditer(str(seq)):
                features.append(SeqFeature(FeatureLocation(m.start()-2,m.start()+5), type="motif",qualifiers={'name':'AT-hook'}))
        return features

    def at_hook_canonical(seq):
        hooks=[r"(?=(\S\S[RPKST]GRP[RPKS]))"]
        # hooks=[r"PRGRP"]

        features=[]
        for h in hooks:
            p = re.compile(h)
            for m in p.finditer(str(seq)):
                features.append(SeqFeature(FeatureLocation(m.start(),m.start()+7), type="motif",qualifiers={'name':'Canonical-AT-hook'}))
        return features

    def cenpc_mot(seq):
        features=[]
        p = re.compile(r'R\S\S\S\SP\S\S[YFW]W')
        for m in p.finditer(str(seq)):
            features.append(SeqFeature(FeatureLocation(m.start(),m.start()+len(m.group())), type="motif",qualifiers={'name':'CENP-C-motif'}))

        return features

    for i in seqreclist:
        i.features=features_via_hmm(i.seq,'inp_data/comb.hmm',eval_thresh=0.1)
        i.features.extend(at_hook_canonical(i.seq))
        i.features.extend(at_hook(i.seq))
        i.features.extend(cenpc_mot(i.seq))

    # print seqreclist[0].features
    # features=features_via_hmm(fasta_dict.values()[0].seq,'inp_data/comb.hmm')
    # print features
    # print features[0].qualifiers['name']
    taxo_seq_architecture(seqreclist,outfile='out_data/taxo_arch.svg',taxids=taxids)

if __name__ == '__main__':
    main()
    
            