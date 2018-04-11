# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script get initial set of records for further filtering and analysis

"""

DBFILE="init_hist_struct.db"
DBDESCR="Initial histone structure database obtained from MMDB"
HISTREC="""
Protocol:
1) MMDB search for *histone*,
2) SIFTS uniprot refs obtained for all chains
3) Entries left where at least one chain had uniprot descr *histone* 
"""

#import pprint
from Bio import Entrez
#from Bio import SeqIO
#import urllib2
#from Bio.PDB.PDBParser import PDBParser
#from xml.dom import minidom
from Bio import ExPASy
from Bio import SwissProt
from Bio.PDB.PDBParser import PDBParser
import re
import csv
#import sqlite3
import cPickle as pickle
import collections
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from PDB_UNP_DB import PDB_UNP, History

Entrez.email = "alexey.shaytan@nih.gov" 


def main():
    
# do some initialization    
    try:
        os.remove('data/'+DBFILE)
    except: pass    
    
    db = create_engine('sqlite:///data/'+DBFILE, echo=True)
    PDB_UNP.metadata.create_all(db)
    Session = sessionmaker(bind=db)
    session = Session()
    
    r=History(DBFILE,DBDESCR,HISTREC)
    session.add(r)
    
    pdb_ids=[]
    mmdb_descr={}
    p = re.compile('Histone.+',re.IGNORECASE)
    
    print '#########Launching#######################\n\
    #########################################'*4
#    handle = Entrez.egquery(term="histone[TITL]")
#    record = Entrez.read(handle)
#    for row in record["eGQueryResult"]: print row["DbName"], row["Count"]    
    handle = Entrez.esearch(db="structure", term="histone", retmax="3000")
    record = Entrez.read(handle)
    print 'Number of records to analyze: ', record['Count']
    list_rec = record["IdList"]
    str = ",".join(list_rec)
    request=Entrez.epost(db="structure",id=str)
    result=Entrez.read(request)
    webEnv=result["WebEnv"]
    queryKey=result["QueryKey"]
    handle=Entrez.esummary(db="structure",webenv=webEnv, query_key=queryKey)
    record = Entrez.read(handle)
    for i in record:
        pdb_ids.append(i['PdbAcc'])
        mmdb_descr[i['PdbAcc']]=i
        
    pdb_to_UNP_file=open('extDB/pdb_chain_uniprot.csv','r')
    csv_reader = csv.reader(pdb_to_UNP_file, delimiter=',') 
    pdb_to_UNP = collections.defaultdict(list)
    for i in csv_reader:
        if(i[0].upper() in pdb_ids):
            pdb_to_UNP[i[0].upper()].append(i[2])            
    pdb_to_UNP_file.close()
        
    print 'After matching objects with SIFTS database', len(pdb_to_UNP), 'objects left'
    print 'Following PDBs were not found in SIFTS databse'
    for i in pdb_ids: 
        if(i not in pdb_to_UNP): 
            print i
    UNP_to_record={}    
    
    for pdbid in pdb_ids:
        counter=collections.Counter(pdb_to_UNP[pdbid])
        print "Analyzing %s" % pdbid
        records=[]
        flag=0
        unpids=[]
        for unpid in counter:           
            handle = ExPASy.get_sprot_raw(unpid)           
            try:
                records.append(SwissProt.read(handle))
                unpids.append(unpid)
            except ValueError:
                print 'Outdated UniProt identifier %s for %s | Quantity: %s\n' %(unpid, pdbid, counter[unpid])
            else:
                m=p.search(records[-1].description)
                if(m): flag=1
        if(flag):
            for rec,unpid in zip(records,unpids):
                r=PDB_UNP(pdbid,unpid,counter[unpid])
                session.add(r)
                UNP_to_record[unpid]=rec
                 #   f.write('%s | Organism %s | Quantity: %s\n' %(m.group(0), record.organism, counter[unpid]))
                   # add_to_DB(pdbid,record)
    session.commit()
    session.close()
    pickle.dump(UNP_to_record, open( "data/"+DBFILE.replace('.db','_')+"UNP_records.dat", "wb" ) )
    
    
    #Now let's make relevant mmdb descr header database and save it
    mmdb_descr2={}
    for key in pdb_to_UNP:
        mmdb_descr2[key]=mmdb_descr[key]
    pickle.dump(mmdb_descr2, open( "data/"+DBFILE.replace('.db','_')+"MMDB_descr.dat", "wb" ) )

if __name__ == '__main__':
    main()
    
    
            