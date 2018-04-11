# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script analyzes the contents of the database, whichever supplied
outputs useful text data and html file with pdb links
"""
#DBFILEINIT="init_nucl_struct.db" #intialdatabase
DBFILEINIT="init_hist_struct.db" #intialdatabase
#DBFILE='init_struct.db'
#DBFILE='nucl_homo.db'
DBFILE='h3_struct.db'

OUTPUT_LINKS_TO_FILES=False #otherwise webpages

import pprint
from libs import markup
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
from sqlalchemy.orm import sessionmaker
from PDB_UNP_DB import PDB_UNP,History

Entrez.email = "alexey.shaytan@nih.gov" 

def main():
    
# do Loading    
    
    db = create_engine('sqlite:///data/'+DBFILE, echo=False)
    Session = sessionmaker(bind=db)
    session = Session()
    res = session.query(PDB_UNP).all()
    histrec=session.query(History).all()
    session.close()
    UNP_to_record=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"UNP_records.dat", "rb" ))
    mmdb_descr=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"MMDB_descr.dat", "rb" ))
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
  #      print rec.pdb_id, rec.unp_id, UNP_to_record[rec.unp_id].taxonomy_id[0]   
  #  print vars(UNP_to_record[rec.unp_id])
    print """#############################################
#        Let's do the analysis              #
#############################################
"""
    gen_info= """********************General*******************
###There are %d unique pdb strucutres 
###There are %d unique protein ids
###There are %d unique taxonomic organisms present 
###Total protein chains in set %d""" % (len(PDBs),len(UNPs),len(TAXIDs),sum([a*b for a,b in zip(NUMPDBs_with_NUMCHAINS.values(),map(int,NUMPDBs_with_NUMCHAINS.keys()))]))   

       
    oligo_distr= "********************Oligomeric state distribution*******************\n"
    oligo_distr+= "###Olig state : Num PDBIDs : PDBIDs \n"
    for key in sorted(NUMPDBs_with_NUMCHAINS.iterkeys()): 
        oligo_distr += '### %s : %d ' %( key, NUMPDBs_with_NUMCHAINS[key])
        for item in PDBs_with_NUMCHAINS[key]: oligo_distr += ' %s ' % item
        oligo_distr+='\n'
    organism_distr='*****************Organismic distribution**************************\n'
    organism_distr+='### Num PDBs with at least one protein from organism | Organism name\n'
    for key in sorted(ORGNs, key=ORGNs.get, reverse=True): organism_distr+='### %8s     %s\n' % (ORGNs[key],key)
        
    prot_distr='*****************Protein distribution**************************\n'
    prot_distr+='### UniProtID | Num PDBs in database having it | Descr\n'
    for key in sorted(UNPs, key=UNPs.get, reverse=True):
        number_pdbs_with_unp=len(session.query(PDB_UNP).filter(PDB_UNP.unp_id==key).all())
        prot_distr+='### %8s %4d     %s\n' % (key,number_pdbs_with_unp,UNP_to_record[key].description)    
        
    
    print gen_info,oligo_distr,organism_distr,prot_distr

    
    
    htmlfile=open('data/'+DBFILE.replace('.db','.html'),'wt')    
    
    title = "PDBids from %s" % DBFILE
    header = "PDBids from %s" % DBFILE
    header+='''<script type="text/javascript" src="https://dl.dropboxusercontent.com/u/201202/js/tooltip.js"></script>
    <script type="text/javascript" src="https://dl.dropboxusercontent.com/u/201202/js/ajax.js"></script>
    <link rel="stylesheet" href="http://code.jquery.com/ui/1.10.3/themes/smoothness/jquery-ui.css" />
    <script src="http://code.jquery.com/jquery-1.9.1.js"></script>
    <script src="http://code.jquery.com/ui/1.10.3/jquery-ui.js"></script>
    <link rel="stylesheet" href="/resources/demos/style.css" />
    <script>
    $(function() {
    $( document ).tooltip();
    });
    </script>
    <style>
    label {
    display: inline-block;
    width: 5em;
    }
    </style>'''
  
    footer = ""
    styles = ( 'layout.css', 'alt.css', 'images.css' )
    page = markup.page( )
    page.init( css=styles, title=title, header=header, footer=footer )
    page.p("<h1> %s </h1>" % histrec[-1].dbdescr)
    page.p("<h1> Database analysis:</h1>")
    page.br( )
    page.p('<h2>History of how the database was created </h2>')
    for h in histrec:
        page.p(str(h.time.strftime("%Y-%m-%d %H:%M"))+' | '+ h.dbname+' | '+h.histrec)
    page.br()
    page.p('<h2>Now comes the analysis</h2>')
    page.p(gen_info.replace('\n','<br>'))
    page.br( )
    page.p("********************Oligomeric state distribution*******************")
    page.p("###Olig state : Num PDBIDs : PDBIDs ")
    for key in sorted(NUMPDBs_with_NUMCHAINS.iterkeys()):
        page.p('### %s : %d' % (key, NUMPDBs_with_NUMCHAINS[key])) # , ' : ', PDBs_with_NUMCHAINS[key]
        for i in PDBs_with_NUMCHAINS[key]:
            if(OUTPUT_LINKS_TO_FILES): page.a( i, class_='internal', href='http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s' % i, title=mmdb_descr[i]['PdbDescr'], onmouseover="showtrail(1,1,'http://www.rcsb.org/pdb/images/%s_asym_r_500.jpg');" % i, onmouseout="hidetrail();" )            
            else: page.a( i, class_='internal', href='http://www.rcsb.org/pdb/explore.do?structureId=%s' % i, title=mmdb_descr[i]['PdbDescr'], onmouseover="showtrail(1,1,'http://www.rcsb.org/pdb/images/%s_asym_r_500.jpg');" % i, onmouseout="hidetrail();" )
    page.br( )
    page.p(organism_distr.replace('\n','<br>').replace(' ','&nbsp;'))
    page.br()
    
    prot_distr='*****************Protein distribution**************************\n'
    prot_distr+='### UniProtID | Num PDBs in database having it | Descr\n'
    page.p(prot_distr.replace('\n','<br>'))
    for key in sorted(UNPs, key=UNPs.get, reverse=True):
        number_pdbs_with_unp=len(session.query(PDB_UNP).filter(PDB_UNP.unp_id==key).all())
        pdb_ids_for_unp=""
        for i in session.query(PDB_UNP).filter(PDB_UNP.unp_id==key).all(): pdb_ids_for_unp+=i.pdb_id+' '
        page.add('### &nbsp;&nbsp;&nbsp;&nbsp;')
        page.a('%8s' % key, href="http://www.uniprot.org/uniprot/%s" % key)
        page.add('&nbsp;&nbsp;')
        page.a('&nbsp;&nbsp; %4d&nbsp;&nbsp;' % number_pdbs_with_unp, href='#', title=pdb_ids_for_unp)
        page.add('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
        page.add('  %s' % UNP_to_record[key].description) 
        page.br()
   
        
        
    htmlfile.write(page.__call__())
    htmlfile.close()



if __name__ == '__main__':
    main()
    
    
            