# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script filters the database

"""

black_list1=['Q54F38','Q54Z07']
black_list2=['O15819','Q6C0C4','Q8SS77','P90543','P0CO05','P0CO04']

DBFILEINIT="hist_seq.db" #initial database
DBFILE="hist_seq.db" #database to filter
DBFILE2="h3_homo_seq.db" #filtered database
DBDESCR="Filtered database for h3 sequences of Homo sapience"
HISTREC="""
Protocol:
1) Filtered inital database by looking for *H3* in UNP description
2) Discarding *CENP*
3) Discarding *Fragment*
4) Discarding untypical and unclear entries: %s
5) Discarding *like*
6) Tuning discarding entries with insertions: %s
7) organism *Homo*
Should be full unambigous sequences of H3 variants
""" % (black_list1, black_list2)


OUTPUT_LINKS_TO_FILES=False #otherwise webpages

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

Entrez.email = "alexey.shaytan@nih.gov" 

def main():
    
# do Loading    
    
    try:
        os.remove('data/'+DBFILE2)
    except: pass
    db = create_engine('sqlite:///data/'+DBFILE, echo=False)
    Session = sessionmaker(bind=db)
    session = Session()
    
    prot_fams_objs = session.query(Prot_Fams).all()
    UNP_recs = session.query(Proteins).order_by(Proteins.descr).all()
    histrec=session.query(History).all()
    UNP_to_seqrecord=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"UNP_seqrecords.dat", "rb" ))
    
    
    #    
#    #Let prepare files for filtered DB
    db2 = create_engine('sqlite:///data/'+DBFILE2, echo=False)
    Proteins.metadata.create_all(db2)
    
    Session2 = sessionmaker(bind=db2)
    session2 = Session2()
    
    #Let's copy history records and add new one

    for rec in  histrec:
        r=History(rec.dbname,rec.dbdescr,rec.histrec,rec.time)
        session2.add(r)
    r=History(DBFILE2,DBDESCR,HISTREC)
    session2.add(r)
    
    for rec in  prot_fams_objs:
        r=Prot_Fams(rec.name)
        session2.add(r)
        session2.commit()
        
    
    p1 = re.compile('H3',re.IGNORECASE)
    p2 = re.compile('CEN',re.IGNORECASE)
    p3 = re.compile('Fragment',re.IGNORECASE)
    p4 = re.compile('like',re.IGNORECASE)
    q1= re.compile('homo',re.IGNORECASE)
    for rec in  UNP_recs:
        m1=p1.search(rec.descr)
        m2=p2.search(rec.descr)
        m3=p3.search(rec.descr)
        m4=p4.search(rec.descr)
        n1=q1.search(rec.organism)
        if(m1 and not m2 and not m3 and not m4 and (rec.unp_id not in (black_list1+black_list2)) and n1):
                r=Proteins(rec.unp_id,rec.descr,rec.organism)
                r.family=session2.query(Prot_Fams).filter(Prot_Fams.name==rec.family.name).first()
                session2.add(r)
    
    session2.commit()
    session2.close()
    session.close()
    

if __name__ == '__main__':
    main()
    
    
            