# -*- coding: utf-8 -*-
"""
Library to work with sets of histone sequences.
This was pre HistoneDB and will be depricated.

"""
import uuid
from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import cPickle as pickle
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess
import pylab
import pandas as pd
import networkx as nx
from StringIO import StringIO
import L_shade_hist_aln
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Emboss.Applications import NeedleCommandline
#from pylab import *
sys.path.append('../sec_str')

from hist_ss import get_hist_ss_in_aln
 
Entrez.email = "alexey.shaytan@nih.gov" 

PATH_to_NCBI_nodes_dmp="nodes.dmp"


os.environ['PATH']='/Users/alexeyshaytan/soft/x3dna-v2.1/bin:/Users/alexeyshaytan/soft/amber12/bin:/Users/alexeyshaytan/soft/sratoolkit/bin:/Users/alexeyshaytan/soft/bins/gromacs-4.6.3/bin:/opt/local/bin:/opt/local/sbin:/Users/alexeyshaytan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/ncbi/blast/bin:/usr/texbin'


def tax_filter_hist(hist_dataframe,clade_top_tax_id):
    """Should take a PANDAS list of histones and leave only those that are inside taxa"""
    
    p=pd.read_table(PATH_to_NCBI_nodes_dmp,sep='|',usecols=[0,1],header=None)
    G=nx.DiGraph()
    #Load all taxonomy as graph
    G.add_edges_from(zip(p.ix[:,1],p.ix[:,0]))
    clade=nx.dfs_tree(G,clade_top_tax_id)
    if(not clade.nodes()):
        clade.add_node(clade_top_tax_id)
    return(hist_dataframe[map(clade.has_node,hist_dataframe.TaxID)])


def hist_table_filter(h,hist_name='H2A',
    var_types=['all'],incl_putative='Yes',incl_predicted='Yes',
    ACCtype=['NP'],score_thresh=0.5,length_deviation=1.2,core_diff_thrs_min=0.95,
    core_diff_thrs_max=1.05,gi_blacklist=[]):
    """Filter histones by parameters"""

    hist_table=h[(h.FullSeq=='Yes')&((h.Putative=='No')|(h.Putative==incl_putative))&((h.Predicted=='No')|(h.Predicted==incl_predicted))&(h.Lendiff<length_deviation)&(h.Lendiff>(1/length_deviation))&(h.Score>score_thresh)&(h.Core_Lendiff<core_diff_thrs_max)&(h.Core_Lendiff>core_diff_thrs_min)]
    # hist_table=h[(h.FullSeq=='Yes')&((h.Putative=='No')|(h.Putative==incl_putative))&((h.Predicted=='No')|(h.Predicted==incl_predicted))&(h.Lendiff<length_deviation)&(h.Score>score_thresh)]
    hist_table=hist_table.reset_index(drop=True)
    # print hist_table.head()
    if(ACCtype[0]=='all'):
        hist_table=hist_table
    else:
        hist_table=hist_table[(hist_table['ACCtype'].apply(lambda x: x in ACCtype))&(hist_table['GI'].apply(lambda x: str(x) not in map(str,gi_blacklist)))]
    # print hist_table.head()
    # print hist_table
    h2=hist_table.reset_index(drop=True)

    if(var_types[0]=='all'):
        h2=h2[h2.Type==hist_name]
        title="%s , RefSeq, Predicted %s, Putative %s, Length Dev %.2f, Score thresh %.2f, Core lendiff thresh %.2f %.2f ACCtype %s"%(hist_name,incl_predicted,incl_putative,length_deviation,score_thresh,core_diff_thrs_min,core_diff_thrs_max,ACCtype)
    else:
        h2=h2[ (h2.Type==hist_name) & (h2['Variant'].apply(lambda x: bool(re.search('|'.join(var_types),x))))]
        title="%s, var %s , RefSeq, Predicted %s, Putative %s, Length Dev %.2f, Score thresh %.2f, Core lendiff thresh %.2f %.2f ACCtype %s"%(hist_name,'-'.join(var_types),incl_predicted,incl_putative,length_deviation,score_thresh,core_diff_thrs_min,core_diff_thrs_max,ACCtype)
    title=title.replace('$','eol')
    title=title.replace('^','bol')

    # print h2.head()
    return h2, title

def get_best_taxid(hist_table,tid,taxid_list):
    """Adding closest matching taxid from tax_dict, given this is a tid clade"""
    p=pd.read_table(PATH_to_NCBI_nodes_dmp,sep='|',usecols=[0,1],header=None)
    G=nx.DiGraph()
    G.add_edges_from(zip(p.ix[:,1],p.ix[:,0]))
    
    G.reverse(copy=False)

    hist_table['BestTaxID']=0
    for i in range(len(hist_table)):
        x=hist_table['TaxID'].iloc[i]
        dist=list()
        tidlist=list()
        for t in taxid_list:
            try:
                d=nx.shortest_path_length(G,x,t)
            except:
                continue
            dist.append(d)
            tidlist.append(t)
        hist_table['BestTaxID'].iloc[i]=tidlist[dist.index(min(dist))]
    return hist_table

def tax_split(hist_table,tax_dict,add_top_taxid=False,add_best_taxid=False):

    hist_tables=collections.defaultdict()

    for tname,tid in tax_dict.iteritems():
        print("#####Filtering by taxonomy %s###"%tname)
        hist_tables[tname]=tax_filter_hist(hist_table,tid)
        if(add_top_taxid):
            hist_tables[tname]['TopTaxID']=tid
        if(add_best_taxid):
            hist_tables[tname]=get_best_taxid(hist_tables[tname],tid,tax_dict.values())
            print hist_tables[tname].head()
    return hist_tables
 
def report_hist_cont(hist_table):
    """Reporting histone content"""
    print("Total histones %d"%len(hist_table))
    print("H3 histones %d"%len(hist_table[(hist_table['Type'].apply(lambda x: 'H3' in x))]))
    print("cenH3 histones %d"%len(hist_table[(hist_table['Variant'].apply(lambda x: 'cenH3' in x))]))
    print("H3.3 histones %d"%len(hist_table[(hist_table['Variant']=='H3.3')]))

    print("H2A histones %d"%len(hist_table[(hist_table['Type']=='H2A')]))
    print("H2A.Z histones %d"%len(hist_table[(hist_table['Variant'].apply(lambda x: 'H2A.Z' in x))]))
    print("H2A.X histones %d"%len(hist_table[(hist_table['Variant']=='H2A.X')]))


if __name__ == '__main__':
    p=pd.read_table(PATH_to_NCBI_nodes_dmp,sep='|',usecols=[0,1],header=None)
    G=nx.DiGraph()
    G.add_edges_from(zip(p.ix[:,1],p.ix[:,0]))
    # G=G.to_undirected()
    x=7227
    tax_dict_full={'human':9606,'mouse':10090,'elegans':6239,'xenopus':262014,'fly':7227,'gallus':9031,'mammals':40674,'metazoa':33208,'yeast':4932,'saccharomyces':716545,'vertebrates':7742,'fungi':4751,'gplants':33090,'eukaryota':2759,'archea':2157, 'bacteria':2}

    taxid_list=tax_dict_full.values()
    # print [nx.shortest_path_length(G,x,y) for y in taxid_list]
    print G.has_node(7227)
    print nx.shortest_path_length(G,33208,7227)
    
            