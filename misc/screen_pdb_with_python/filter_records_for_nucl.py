# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
filter database
"""
DBFILEINIT="init_nucl_struct.db" #initial database
DBFILE="init_nucl_struct.db" #database to filter
DBFILE2="nucl.db" #filtered database
DBDESCR="Filtered database of nucleosomes and octamers from initial nucleosome database, includes multinucleosomes"
HISTREC="""
Protocol:
1) Select PDBs with no less than 8 protein chains in assymetric unit
2) Manual review: get rid of 4DRA 4FIP 2X4X 3C9K (tubular crystall) 4FJC 3F9X 3F9Z 2YFW 2XQL 4DRB 4E45 2X4Y 3U5P 3U5O 2Z3F 3F9W
"""

import pprint
from libs import markup
from Bio import Entrez
import os
#from Bio import SeqIO
#import urllib2
#from Bio.PDB.PDBParser import PDBParser
#from xml.dom import minidom 
from Bio import ExPASy
from Bio import SwissProt
import re
import csv
#import sqlite3
import cPickle as pickle
import collections
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from PDB_UNP_DB import PDB_UNP,History

Entrez.email = "alexey.shaytan@nih.gov" 

def main():
    
# do Loading    
    try:
        os.remove('data/'+DBFILE2)
    except: pass
    db = create_engine('sqlite:///data/'+DBFILE, echo=False)
    Session = sessionmaker(bind=db)
    session = Session()
    res = session.query(PDB_UNP).all()
    histrec=session.query(History).all()
    mmdb_descr=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"MMDB_descr.dat", "rb" ))
    UNP_to_record=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"UNP_records.dat", "rb" ))
    
    UNPs=collections.Counter()
    PDBs=collections.Counter()
    TAXIDs=collections.Counter()
    ORGNs=collections.Counter()
    CHAINs_in_PDB=collections.Counter()
    NUMPDBs_with_NUMCHAINS=collections.Counter()
    PDBs_with_NUMCHAINS=collections.defaultdict(list)
    print "Read %d lines in table" % len(res)
    for rec in res:
        UNPs[rec.unp_id]+=1
        PDBs[rec.pdb_id]+=1
        TAXIDs[UNP_to_record[rec.unp_id].taxonomy_id[0]]+=1
        ORGNs[UNP_to_record[rec.unp_id].organism]+=1
        CHAINs_in_PDB[rec.pdb_id]+=rec.unp_num
    for pdb_id in CHAINs_in_PDB:
        NUMPDBs_with_NUMCHAINS[CHAINs_in_PDB[pdb_id]]+=1
        PDBs_with_NUMCHAINS[CHAINs_in_PDB[pdb_id]].append(pdb_id)
        


#    
#    #Let prepare files for filtered DB
    db2 = create_engine('sqlite:///data/'+DBFILE2, echo=False)
    PDB_UNP.metadata.create_all(db2)
    
    Session2 = sessionmaker(bind=db2)
    session2 = Session2()
    
    #Let's copy history records and add new one

    for rec in  histrec:
        r=History(rec.dbname,rec.dbdescr,rec.histrec,rec.time)
        session2.add(r)
    r=History(DBFILE2,DBDESCR,HISTREC)
    session2.add(r)
    

    blacklist=['4DRA', '4FIP', '2X4X', '3C9K', '4FJC', '3F9X', '3F9Z', '3F9W', '2YFW', '2XQL', '4DRB', '4E45', '2X4Y', '3U5P', '3U5O', '2Z3F']
#    for i in UNPs:
#        if(i not in blacklist):
#            print i, ' : ', UNPs[i], ' : ', UNP_to_record[i].description    
    
    ###Filtering and copy
    p = re.compile('3\.3',re.IGNORECASE)
    k=0
#    needed_pdbs=[]
#    for rec in res:
#        if(p.search(UNP_to_record[rec.unp_id].description)):
#            needed_pdbs.append(rec.pdb_id)
    
            
    for rec in res:
         if(rec.pdb_id not in blacklist):
            num_chains=0
            for rec2 in session.query(PDB_UNP).filter(PDB_UNP.pdb_id==rec.pdb_id).all():
                num_chains+=rec2.unp_num
            if(num_chains>7):
                r=PDB_UNP(rec.pdb_id,rec.unp_id,rec.unp_num)           
                session2.add(r)
                #  print rec.pdb_id, rec.unp_id
                k+=1
    session2.commit()
    session2.close()
    session.close()
    print "Output %d lines" % k

if __name__ == '__main__':
    main()
    
    
            