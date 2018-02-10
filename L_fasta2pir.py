#!/usr/bin/env python2.7
"""
This library converts FASTA multiple alignment files to
PIR files as needed by MODELLER by adding additional information.
"""
__author__="Alexey Shaytan"



from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO



class aln():
    def __init__(self,BioPythonMSA=None,FASTAFile=None):
        if (BioPythonMSA!=None):
            self.MSA=BioPythonMSA
        else:
            self.MSA=AlignIO.read(FASTAFile, "fasta")
        self.pir_info=dict()
   
    def add_pir_info(self,id,type,pdb_file,res_beg=' ',chain_beg='.',res_end=' ',chain_end='.',prot_name='',src='',resolution='',Rfactor=''):
    	self.pir_info[id]={'type':type,'pdb_file':pdb_file,'res_beg':res_beg,'chain_beg':chain_beg,'res_end':res_end,'chain_end':chain_end,'prot_name':prot_name,'src':src,'resolution':resolution,'Rfactor':Rfactor}
        
    def write(self,file_name):
    	f = open(file_name, 'w')
    	f.write("C; this is a PIR alignment file for MODELLER\n\n")
    	for item in self.MSA:
    		f.write(">P1;%s\n"%item.id)
    		f.write("%s:%s:%-5s:%s:%-5s:%s:%s:%s:%s:%s\n"%(self.pir_info[item.id]['type'],\
    			self.pir_info[item.id]['pdb_file'],\
    			self.pir_info[item.id]['res_beg'],\
    			self.pir_info[item.id]['chain_beg'],\
    			self.pir_info[item.id]['res_end'],\
    			self.pir_info[item.id]['chain_end'],\
    			self.pir_info[item.id]['prot_name'],\
    			self.pir_info[item.id]['src'],\
    			self.pir_info[item.id]['resolution'],\
    			self.pir_info[item.id]['Rfactor']))
    		f.write("%s*\n\n"%item.seq)
    	f.close()

class aln_mult_chains():
    """
    This class deals with PIR fromat for multiple chains
    (a modeller hack in PIR format to handle homology modelling for multiple chain proteins)
    """
    def __init__(self,aln_list=[]):
        """
        Here we initialize with a list of several alignments in aln-class format.
        The ids of sequences in MSA should be the same
        """
        self.aln_list=aln_list
        #Impotant we sort by id, so that all rows in all alignements are in correspondence.
        #This is because Muscle often mixes the rows
        for i in range(len(self.aln_list)):
            self.aln_list[i].MSA.sort()


    def write(self,file_name):
        f = open(file_name, 'w')
        f.write("C; this is a PIR alignment file for MODELLER\n\n")
        n=0
        for item in self.aln_list[0].MSA:
            f.write(">P1;%s\n"%item.id)
            f.write("%s:%s:%-5s:%s:%-5s:%s:%s:%s:%s:%s\n"%(self.aln_list[0].pir_info[item.id]['type'],\
                self.aln_list[0].pir_info[item.id]['pdb_file'],\
                self.aln_list[0].pir_info[item.id]['res_beg'],\
                self.aln_list[0].pir_info[item.id]['chain_beg'],\
                self.aln_list[-1].pir_info[item.id]['res_end'],\
                self.aln_list[-1].pir_info[item.id]['chain_end'],\
                self.aln_list[0].pir_info[item.id]['prot_name'],\
                self.aln_list[0].pir_info[item.id]['src'],\
                self.aln_list[0].pir_info[item.id]['resolution'],\
                self.aln_list[0].pir_info[item.id]['Rfactor']))
            f.write("%s/\n"%item.seq)
            for i in self.aln_list[1:-1]:
                f.write("%s/\n"%i.MSA[n].seq)
            f.write("%s*\n\n"%self.aln_list[-1].MSA[n].seq)
            n=n+1
        f.close()


if __name__ == '__main__':
	MSA=MultipleSeqAlignment([SeqRecord(Seq('ABCDEFG'),id='h2a_xen'),SeqRecord(Seq('QRSTUVW'),id='h2a_yeast')])
	a=aln(MSA)
	a.add_pir_info('h2a_xen','structureX','h2a_x','1')
	a.add_pir_info('h2a_yeast','structureX','h2a_y','1')

	a.write('testaln.pir')




