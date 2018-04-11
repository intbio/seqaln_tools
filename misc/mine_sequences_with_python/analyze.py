# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
this script analyzes the database

"""

DBFILEINIT="hist_seq.db" #intialdatabase

#DBFILE='init_struct.db'
#DBFILE='hist_seq.db'

DBFILE='h3_1_2_3_homo.db'

SHORT_NAMES=True #otherwise webpages

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
    
    db = create_engine('sqlite:///data/'+DBFILE, echo=False)
    Session = sessionmaker(bind=db)
    session = Session()
    prot_fams_objs = session.query(Prot_Fams).all()
    UNP_rec_obj = session.query(Proteins).all()
    histrec=session.query(History).all()
    
    
    UNP_to_seqrecord=pickle.load(open( "data/"+DBFILEINIT.replace('.db','_')+"UNP_seqrecords.dat", "rb" ))
    
    
    print """#############################################
#        Let's do the analysis              #
#############################################
"""
    gen_info= """********************General*******************
###There are %d unique seq strucutres 
""" % (len(UNP_rec_obj))   

    p = re.compile('Full=(.+);',re.IGNORECASE)

    print '#####Summary table'  
    summary_data=[['Histone Family','Protein name','Organisms','UniProt IDs']]
    for k in prot_fams_objs:
        UNP_rec_objs=session.query(Proteins).filter(Proteins.family==k).order_by(Proteins.descr).all()
        summary_data.append(['######### Family: %s' % k.name])
        for i in UNP_rec_objs:
            m=p.search(i.descr)
            href='<a href="http://www.uniprot.org/uniprot/%s">%s</a>' % (i.unp_id,i.unp_id)
            if(SHORT_NAMES):
                prot_data=[k.name,m.group(1),i.organism,href]
            else:
                prot_data=[k.name,i.descr,i.organism,href]
            summary_data.append(prot_data)
        
     
       
    print gen_info, summary_data #,oligo_distr,organism_distr,prot_distr

    
    
    htmlfile=open('data/'+DBFILE.replace('.db','.html'),'wt')    
    
    title = "Seqs from %s" % DBFILE
    header = "Seqs from %s" % DBFILE
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
    page.p('<h2>Summary table </h2>')
    htmlcode = HTML.table(summary_data)
    page.add(htmlcode)
    
#    page.p("********************Oligomeric state distribution*******************")
#    page.p("###Olig state : Num PDBIDs : PDBIDs ")
#    for key in sorted(NUMPDBs_with_NUMCHAINS.iterkeys()):
#        page.p('### %s : %d' % (key, NUMPDBs_with_NUMCHAINS[key])) # , ' : ', PDBs_with_NUMCHAINS[key]
#        for i in PDBs_with_NUMCHAINS[key]:
#            if(OUTPUT_LINKS_TO_FILES): page.a( i, class_='internal', href='http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s' % i, title=mmdb_descr[i]['PdbDescr'], onmouseover="showtrail(1,1,'http://www.rcsb.org/pdb/images/%s_asym_r_500.jpg');" % i, onmouseout="hidetrail();" )            
#            else: page.a( i, class_='internal', href='http://www.rcsb.org/pdb/explore.do?structureId=%s' % i, title=mmdb_descr[i]['PdbDescr'], onmouseover="showtrail(1,1,'http://www.rcsb.org/pdb/images/%s_asym_r_500.jpg');" % i, onmouseout="hidetrail();" )
#    page.br( )
#    page.p(organism_distr.replace('\n','<br>').replace(' ','&nbsp;'))
#    page.br()
#    
#    prot_distr='*****************Protein distribution**************************\n'
#    prot_distr+='### UniProtID | Num PDBs in database having it | Descr\n'
#    page.p(prot_distr.replace('\n','<br>'))
#    for key in sorted(UNPs, key=UNPs.get, reverse=True):
#        number_pdbs_with_unp=len(session.query(PDB_UNP).filter(PDB_UNP.unp_id==key).all())
#        pdb_ids_for_unp=""
#        for i in session.query(PDB_UNP).filter(PDB_UNP.unp_id==key).all(): pdb_ids_for_unp+=i.pdb_id+' '
#        page.add('### &nbsp;&nbsp;&nbsp;&nbsp;')
#        page.a('%8s' % key, href="http://www.uniprot.org/uniprot/%s" % key)
#        page.add('&nbsp;&nbsp;')
#        page.a('&nbsp;&nbsp; %4d&nbsp;&nbsp;' % number_pdbs_with_unp, href='#', title=pdb_ids_for_unp)
#        page.add('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
#        page.add('  %s' % UNP_to_record[key].description) 
#        page.br()
#   
#        
        
    htmlfile.write(page.__call__())
    htmlfile.close()
    session.close()

if __name__ == '__main__':
    main()
    
    
            