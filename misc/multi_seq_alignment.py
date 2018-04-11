# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script works with sequnces

"""

DBFILEINIT="hist_seq.db" #intialdatabase

#DBFILE='init_struct.db'
#DBFILE='hist_seq.db'
DBFILE='h3_1_2_3_homo.db'



import pprint
from libs import markup
from libs import HTML
from Bio import Entrez
#from Bio import SeqIO
#import urllib2
#from Bio.PDB.PDBParser import PDBParser
#from xml.dom import minidom 
from Bio import ExPASy
#from Bio import SwissProt
from Bio import SeqIO
import re
import csv
#import sqlite3
import cPickle as pickle
import collections
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from seq_DB import Proteins, Prot_Fams, History
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO

Entrez.email = "alexey.shaytan@nih.gov" 

def main():
    
# do Loading    
    
    db = create_engine('sqlite:///data/'+DBFILE, echo=False)
    Session = sessionmaker(bind=db)
    session = Session()
    prot_fams_objs = session.query(Prot_Fams).all()
    UNP_recs = session.query(Proteins).order_by(Proteins.descr).all()
    histrec=session.query(History).all()
    
    
    UNP_to_seqrecord=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"UNP_seqrecords.dat", "rb" ))
    
    seqlist=[]
    ##Let's get the sequences
    for i in UNP_recs:
        seqlist.append(UNP_to_seqrecord[i.unp_id])
        print UNP_to_seqrecord[i.unp_id].format('fasta')
        
    #let's write to file
    output_handle = open("data/temp_mult_alignment.fasta", "w")
    SeqIO.write(seqlist, output_handle, "fasta")
    output_handle.close()
    
    p = re.compile('Full=(.+)$',re.IGNORECASE)
        
    muscle_cline = MuscleCommandline('bin/muscle',input="data/temp_mult_alignment.fasta")
    stdout, stderr = muscle_cline()
   # print stderr
    align = AlignIO.read(StringIO(stdout), "fasta")
    align.sort(key=lambda f: p.search(f.description).group(1))
  #  print align
    count = AlignIO.write(align, "data/"+DBFILE.replace('.db','_aln')+".fasta", "fasta")
    print "Saved %i alignments" % count
    
    
    session.close()
    os.remove("data/temp_mult_alignment.fasta")
    
if __name__ == '__main__':
    main()
    
    
            