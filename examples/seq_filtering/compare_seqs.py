# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
Here we compare yeast and Xenopus (1kx5) histone sequences.
The gis are:
    1kx5        yeast
H3 27573749 6319482
H4 27573750 6319481
H2A 27573751 398366187
H2B 27573752 6319471

"""

import sys
import os
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS')
os.environ['PATH']='/Users/alexeyshaytan/soft/x3dna-v2.1/bin:/Users/alexeyshaytan/soft/amber12/bin:/Users/alexeyshaytan/soft/sratoolkit/bin:/Users/alexeyshaytan/soft/bins/gromacs-4.6.3/bin:/opt/local/bin:/opt/local/sbin:/Users/alexeyshaytan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/ncbi/blast/bin:/usr/texbin'

import pandas as pd
from L_aln_tools import *
from L_shade_hist_aln import *


# from Bio.Blast.Applications import NcbiblastpCommandline
# from Bio import ExPASy
# from Bio import SwissProt
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord


# from Bio import AlignIO
# from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.Polypeptide import PPBuilder
# import csv
# import collections
# from Bio import Entrez
# import cPickle as pickle
# from Bio import SeqIO

# import uuid

# from Bio.Align import MultipleSeqAlignment
# import re
# from Bio import AlignIO
# from Bio.Align.Applications import MuscleCommandline
# import subprocess
# import pylab
# import networkx as nx
# from StringIO import StringIO
# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML
# #from pylab import *
# import itertools
# import time
# from multiprocessing import Pool, current_process


# # from hist_ss import get_hist_ss
# from hist_ss import get_core_lendiff

Entrez.email = "alexey.shaytan@nih.gov" 



def main():


    hist_df=pd.read_csv('gis.csv',dtype={'xen_1kx5':object,'yeast_rs':object,'yeast_1id3':object},index_col='Type')
    fasta_records=get_prot_seq_by_gis(list(hist_df.xen_1kx5)+list(hist_df.yeast_rs)+list(hist_df.yeast_1id3))

    hist=['H3','H4','H2A','H2B']
    hist_aln=dict()
    for i in hist:
        hist_aln[i]=muscle_aln([SeqRecord(fasta_records[hist_df.loc[i,'xen_1kx5']],id="Xen|1kx5|%s"%i),SeqRecord(fasta_records[hist_df.loc[i,'yeast_rs']],id="Yeast|rs|%s"%i),SeqRecord(fasta_records[hist_df.loc[i,'yeast_1id3']],id="Yeast|1id3|%s"%i)])
        #Let's shade!
        hist_aln[i].sort()
        get_pdf(i,hist_aln[i],i)


    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    main()
    
            