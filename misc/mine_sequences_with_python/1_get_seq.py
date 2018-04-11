# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script gets initial set of uniprot records for analysis

"""

DBFILE="hist_seq.db"
DBDESCR="Histone sequence structure database from UNIPROT"
HISTREC="""
Protocol:
get reviewed protein records from UNIPROT for different histone protein families
"""

#import pprint
#from Bio import Entrez
#from Bio import SeqIO
#import urllib2
#from Bio.PDB.PDBParser import PDBParser
#from xml.dom import minidom
from Bio import ExPASy
#from Bio import SwissProt
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
import re
import csv
#import sqlite3
import cPickle as pickle
import collections
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from seq_DB import Proteins, Prot_Fams, History

#Entrez.email = "alexey.shaytan@nih.gov" 


def main():
    
# do some initialization    
    try:
        os.remove('data/'+DBFILE)
    except: pass    
    
    db = create_engine('sqlite:///data/'+DBFILE, echo=True)
    Proteins.metadata.create_all(db)
    Session = sessionmaker(bind=db)
    session = Session()
    
    r=History(DBFILE,DBDESCR,HISTREC)
    session.add(r)
    
#    pdb_ids=[]
#    mmdb_descr={}
#    p = re.compile('Histone.+',re.IGNORECASE)
    
    print '#########Launching#######################\n\
    #########################################'*4
#    Define families to fetch
    prot_fams=['Histone H3','Histone H4', 'Histone H2A', 'Histone H2B', 'Histone H1 H5']
#    prot_fams=['test']
    print '''   Let's load saved UNIPROT ids  '''
    UNP_ids_by_fam = collections.defaultdict(list)
    for k in prot_fams:
        UNP_ids_file_name=k.replace(' ','_')
        UNP_ids_file=open('extDB/'+UNP_ids_file_name+'.list','r')
        csv_reader = csv.reader(UNP_ids_file, delimiter=' ') 
        for i in csv_reader:
            UNP_ids_by_fam[k].append(i[0])
        UNP_ids_file.close()
        print 'Read %s UNPs for %s' % (len(UNP_ids_by_fam[k]),k)
        
    print '''   Let's load UNIPROT records  and populate database'''
    
    prot_fams_obj=collections.defaultdict(Prot_Fams)
    UNP_to_record={} 
    UNP_to_seqrecord={} 
    UNP_rec_obj=collections.defaultdict(Proteins)
    for k in prot_fams:
        print '####   For %s' % k
        prot_fams_obj[k]=Prot_Fams(k)
        session.add(prot_fams_obj[k])
        for i in UNP_ids_by_fam[k]:
            handle = ExPASy.get_sprot_raw(i)   
          #  UNP_to_seqrecord[i] = SeqIO.read(handle, "swiss")
            try:
             #   UNP_to_record[i]=SwissProt.read(handle)
                UNP_to_seqrecord[i] = SeqIO.read(handle, "swiss")
            except ValueError:
                print 'Outdated UniProt identifier %s \n' % i
            else:
               # UNP_rec_obj[i]=Proteins(i,UNP_to_record[i].description,UNP_to_record[i].organism)
                UNP_rec_obj[i]=Proteins(i,UNP_to_seqrecord[i].description,UNP_to_seqrecord[i].annotations['organism'])
                UNP_rec_obj[i].family=prot_fams_obj[k]
                session.add(UNP_rec_obj[i])
                print i, 'added'
                
      
    pickle.dump(UNP_to_seqrecord, open( "data/"+DBFILE.replace('.db','_')+"UNP_seqrecords.dat", "wb" ) )   
    session.commit()
    session.close()
#
#    print 'After matching objects with SIFTS database', len(pdb_to_UNP), 'objects left'
#    print 'Following PDBs were not found in SIFTS databse'
#    for i in pdb_ids: 
#        if(i not in pdb_to_UNP): 
#            print i
#       
#    
#    for pdbid in pdb_ids:
#        counter=collections.Counter(pdb_to_UNP[pdbid])
#        print "Analyzing %s" % pdbid
#        records=[]
#        flag=0
#        unpids=[]
#        for unpid in counter:           
#            handle = ExPASy.get_sprot_raw(unpid)           
#            try:
#                records.append(SwissProt.read(handle))
#                unpids.append(unpid)
#            except ValueError:
#                print 'Outdated UniProt identifier %s for %s | Quantity: %s\n' %(unpid, pdbid, counter[unpid])
#            else:
#                m=p.search(records[-1].description)
#                if(m): flag=1
#        if(flag):
#            for rec,unpid in zip(records,unpids):
#                r=PDB_UNP(pdbid,unpid,counter[unpid])
#                session.add(r)
#                UNP_to_record[unpid]=rec
#                 #   f.write('%s | Organism %s | Quantity: %s\n' %(m.group(0), record.organism, counter[unpid]))
#                   # add_to_DB(pdbid,record)
#    session.commit()
#    session.close()
#    pickle.dump(UNP_to_record, open( "data/"+DBFILE.replace('.db','_')+"UNP_records.dat", "wb" ) )
#    
#    
#    #Now let's make relevant mmdb descr header database and save it
#    mmdb_descr2={}
#    for key in pdb_to_UNP:
#        mmdb_descr2[key]=mmdb_descr[key]
#    pickle.dump(mmdb_descr2, open( "data/"+DBFILE.replace('.db','_')+"MMDB_descr.dat", "wb" ) )

if __name__ == '__main__':
    main()
    
    
            