# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script searches the MMDB database for term histone, gets pdb codes of relevant structures,
then looks in SIFTS file for uniprot refs for all chains, goes to uniprot, looks if the title of the
record contains 'Histone', if yes remembers the pdb record and uniprot record.
Writes PDB to UNIPROT records table and saves uniprot records to file.
"""

#import pprint
from Bio import Entrez
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
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Entrez.email = "alexey.shaytan@nih.gov" 


Base = declarative_base()



class Hist_db(Base):
    """"""
    __tablename__ = "hist_db"
 
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String(4))
    unp_id = Column(String(4))
    unp_num = Column(Integer)
 
    def __init__(self,pdb_id, unp_id, unp_num):
        """"""
        self.pdb_id = pdb_id
        self.unp_id = unp_id
        self.unp_num = unp_num   


def main():
    
# do Loading    
    
    db = create_engine('sqlite:///hist_struct.db', echo=True)
    Session = sessionmaker(bind=db)
    session = Session()
    res = session.query(Hist_db).all()
  #  session.close()
    UNP_to_record=pickle.load(open( "UNP_records.dat", "rb" ))
    UNPs=collections.Counter()
    print "Read %d lines in table" % len(res)
    for rec in res:
        UNPs[rec.unp_id]+=1
       # print rec.pdb_id, rec.unp_id, UNP_to_record[rec.unp_id].organism    

####Manual filtering########  
# These are blacklisted for beening not histones   
    blacklist1=['Q8NB78', 'Q9V452', 'P19267', 'A8MT69', 'Q96L73', 'O15550', 'Q9V444', 'P48781', \
    'O41094', 'Q8N2Z9', 'Q9UPP1', 'O60341', 'Q04089', 'P29375', 'Q03330', 'P53165', 'Q9Y294',\
    'O74515', 'Q16576', 'Q24572', 'P40019', 'O15379', 'Q03164', 'P32447', 'Q92831', 'O93641',\
    'Q07794', 'Q8TEK3', 'Q09028', 'Q9NQR1', 'Q80TJ7', 'O75164', 'P40161', 'Q9H9B1', 'Q6ZMT4']
    
# These are blacklisted for beening not histones 

    
#    for i in UNPs:
#        if(i not in blacklist):
#            print i, ' : ', UNPs[i], ' : ', UNP_to_record[i].description
#    
#    #Let's make filtered DB
    db2 = create_engine('sqlite:///hist_struct_filtered.db', echo=True)
    Base.metadata.create_all(db2)
    
    Session2 = sessionmaker(bind=db2)
    session2 = Session2()
    k=0
    for rec in res:
        if(rec.unp_id not in blacklist1):
            r=Hist_db(rec.pdb_id,rec.unp_id,rec.unp_num)           
            session2.add(r)
          #  print rec.pdb_id, rec.unp_id
            k+=1
    session2.commit()
    session2.close()
    session.close()
    print "Output %d lines" % k

if __name__ == '__main__':
    main()
    
    
            