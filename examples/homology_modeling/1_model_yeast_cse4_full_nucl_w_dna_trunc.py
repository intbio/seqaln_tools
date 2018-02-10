#This script homology models yeast cenp-a nucleosome....

#Try1
#based on 1kx5 nucleosome structure
#we take 
#yeast HTA1 gene (gi 398366187) - as H2A,
# HTB1 (gi 398366183) as H2B,
#HHT1 (gi 6319482) - this should be H3
# but we will take CSE4 gi 27808712
# HHF1 (gi 6319481) - as H4
# we put the 601 DNA sequence
#we truncate the tails according to our NCP model simulations (see Shaytan et al. NAR(?) 2015)
#

#DNA is following
#aaGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAG T TTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAaaga

#Alexey Shaytan 2015

import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from chimera import openModels
import sys
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS')
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/X3DNA_analysis_modeling')


sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/structure_analysis')
sys.path.append('/Library/modeller-9.14/modlib')
sys.path.append('/Library/modeller-9.14/lib/mac10v4')
# sys.path.append(';'.join(['', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/setuptools-0.9.6-py2.7.egg', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/PROPKA-3.1-py2.7.egg', '/Library/modeller-9.14/examples/atom_files', '/Applications/Bioinf/biana', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/readline', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PIL', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyObjC', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/gtk-2.0']))
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import L_fasta2pir
from L_aln_tools import muscle_aln

# from dna_tools_simple import change_dna_seq_in_pdb



os.environ['LD_LIBRARY_PATH']='/Library/modeller-9.14/lib/mac10v4'

from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

#####
#Extract chains with Chimera
nucl_xen=openModels.open('1kx5.pdb',type='PDB')
# nucl_yeast=openModels.open('1id3.pdb',type='PDB')
# h2a_xen=Seq(str(nucl[0].sequence('C')))
# h2a_yeast=Seq(str(nuclZ[0].sequence('C')))
rc('select :.A :.B :.C :.D :.E :.F :.G :.H')
rc('write selected format pdb #0 xen_nucl.pdb')
# rc('write selected format pdb #1 h2a_yeast_xray.pdb')

#generate alignments
#Biopython extracts seqs from pdb
p = PDBParser(PERMISSIVE=1)
s = p.get_structure('1kx5', '1kx5.pdb')
ppb=PPBuilder()
seqs_xen=dict()
for i in ['A','B','C','D','E','F','G','H']:
	seqs_xen[i]=ppb.build_peptides(s[0][i])[0].get_sequence()

#Here we manually input the trunctaed versions of yeast sequences
seqs_yeast=dict()
#CSE4 - H3
#SSKQQWVSSAIQSDSSGRSLSNVNRLAGDQQSINDRALSLLQRTRATKNLFPRREERRRYESSKSDLDIETDYEDQAGNLEIETENEEEAEMETEVPAPVRTHSYALDRYVRQKRREKQRKQSLKR
#VEKKYTPSELALYEIRKYQRSTDLLISKIPFARLVKEVTDEFTTKDQDLRWQSMAIMALQEASEAYLVGLLEHTNLLALHAKRITIMKKDMQLARRIRGQFI
#VEKKYTPSELALYEIRKYQRSTDLLISKIPFARLVKEVTDEFTTKDQDLRWQSMAIMALQEASEAYLVGLLEHTNLLALHAKRITIMKKDMQLARRIRGQFI
seqs_yeast['A']=Seq('VEKKYTPSELALYEIRKYQRSTDLLISKIPFARLVKEVTDEFTTKDQDLRWQSMAIMALQEASEAYLVGLLEHTNLLALHAKRITIMKKDMQLARRIRGQFI')
seqs_yeast['E']=Seq('VEKKYTPSELALYEIRKYQRSTDLLISKIPFARLVKEVTDEFTTKDQDLRWQSMAIMALQEASEAYLVGLLEHTNLLALHAKRITIMKKDMQLARRIRGQFI')
#H4
seqs_yeast['B']=Seq('KRHRKILRDNIQGITKPAIRRLARRGGVKRISGLIYEEVRAVLKSFLESVIRDSVTYTEHAKRKTVTSLDVVYALKRQGRTLYGFGG')
seqs_yeast['F']=Seq('KRHRKILRDNIQGITKPAIRRLARRGGVKRISGLIYEEVRAVLKSFLESVIRDSVTYTEHAKRKTVTSLDVVYALKRQGRTLYGFGG')
#H2A
#MSGGKGGKAGSAA KASQ SRSAKAGLTFPVGRVHRLLRRGNYAQRIGSGAPVYLTAVLEYLAAEILELAGNAARDNKKTRIIPRHLQLAIRNDDELNKLLGNVTIAQGGVLPNIHQNLLPK KSAKATKASQEL
# SGRGKQGGKTR  AKAK TRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPK KTESSKSKSK
seqs_yeast['C']=Seq('KASQSRSAKAGLTFPVGRVHRLLRRGNYAQRIGSGAPVYLTAVLEYLAAEILELAGNAARDNKKTRIIPRHLQLAIRNDDELNKLLGNVTIAQGGVLPNIHQNLLPK')
seqs_yeast['G']=Seq('KASQSRSAKAGLTFPVGRVHRLLRRGNYAQRIGSGAPVYLTAVLEYLAAEILELAGNAARDNKKTRIIPRHLQLAIRNDDELNKLLGNVTIAQGGVLPNIHQNLLPK')

#H2B
#       AKSAPAPKKGSKKAVTKTQKKDGKKRRKT RKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK
#MSAKAEKKPASKAPAEKKPAAKKTSTSTDGKKRSKA RKETYSSYIYKVLKQTHPDTGISQKSMSILNSFVNDIFERIATEASKLAAYNKKSTISAREIQTAVRLILPGELAKHAVSEGTRAVTKYSSSTQA
seqs_yeast['D']=Seq('RKETYSSYIYKVLKQTHPDTGISQKSMSILNSFVNDIFERIATEASKLAAYNKKSTISAREIQTAVRLILPGELAKHAVSEGTRAVTKYSSSTQA')
seqs_yeast['H']=Seq('RKETYSSYIYKVLKQTHPDTGISQKSMSILNSFVNDIFERIATEASKLAAYNKKSTISAREIQTAVRLILPGELAKHAVSEGTRAVTKYSSSTQA')

# s = p.get_structure('1id3', '1id3.pdb')
# ppb=PPBuilder()
# seqs_yeast=dict()
# for i in ['A','B','C','D','E','F','G','H']:
	# seqs_yeast[i]=ppb.build_peptides(s[0][i])[0].get_sequence()
# print nucl_yeast

#Biopython aligns them and prepares for PIR format
#Force gaps to be with high penalty.
seq_aln=dict()
for i in ['A','B','C','D','E','F','G','H']:
	seqlist=[SeqRecord(seqs_xen[i],id='xen'),SeqRecord(seqs_yeast[i],id='yeast')]
	aln_pir=L_fasta2pir.aln(muscle_aln(seqlist))
	aln_pir.add_pir_info('xen','structureX','xen_nucl', 'FIRST', i,'LAST',i)
	aln_pir.add_pir_info('yeast','sequence','nucl_yeast')
	seq_aln[i]=aln_pir

mult_aln=L_fasta2pir.aln_mult_chains([seq_aln[i] for i in ['A','B','C','D','E','F','G','H']])

mult_aln.write('aln.pir')

#####
#Now let's do MODELLER
env = environ()  # create a new MODELLER environment to build this model in
env.io.atom_files_directory = ['.','data']
# change to folder with data files
log.verbose()    # request verbose output
# directories for input atom files
# env.io.atom_files_directory = ['data']
a = automodel(env,
              alnfile  = 'aln.pir',
              knowns   = 'xen',
              sequence = 'yeast')
# a.library_schedule = autosched.slow
# a.max_var_iterations = 300
# automodel.write_intermediates=True
a.starting_model= 1
a.ending_model  = 1
a.rand_method=None
a.max_sc_mc_distance=10
a.max_sc_sc_distance=10
a.md_level=None # We will do MD anyway, and the structurea are similar.

a.make()

print a.outputs[0]['name']


rc('close session')


#let's generate a new DNA sequence
#aaGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAG
#T
#TTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAaaga
#len 162 - need to add to 147, i.e. 12 to each side
# cen3_seq='AAGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAGTTTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAAAGA'
# cen3_seq=map(lambda x: x, cen3_seq)
# cen3_seq_flanked=['A']*12+cen3_seq+['A']*12

# change_dna_seq_in_pdb('1kx5.pdb','1kx5_cen3.pdb',cen3_seq_flanked)



rc('open %s'%a.outputs[0]['name'])
rc('open 1kx5.pdb')
rc('open 1kx5_cen3.pdb')

# rc('open h2a_yeast_xray.pdb')
rc('matchmaker #1 #0')

rc('match #2:@N1,N9 #1:@N1,N9')
# rc('delete ~#2:.A & ~#2:.B & ~#0')
rc('changechains A,B I,J #2')

rc('combine #0,2 newchainids false close false')
rc('delete :MN')
rc('delete :HOH')

rc('resrenumber 127 #3:.A')
rc('resrenumber 127 #3:.E')

rc('resrenumber 16 #3:.B')
rc('resrenumber 16 #3:.F')

rc('resrenumber 13 #3:.C')
rc('resrenumber 13 #3:.G')

rc('resrenumber 36 #3:.D')
rc('resrenumber 36 #3:.H')

rc('resrenumber -73 #3:.I')
rc('resrenumber -73 #3:.J')


rc('write format pdb #3 temp.pdb')

rc('close session')

rc('open nucl_aln.pdb')
rc('open temp.pdb')
rc('matchmaker #0 #1')
rc('write format pdb #1 yeast_nucl_cse4_trunc_aln_v1.pdb')

os.system('cp yeast_nucl_cse4_trunc_aln_v1.pdb ../2_MD_relax/yeast_cse4/prep/')



# print("KUKU")
# quit()
# #Position of the last base on I strand that is curled in DNA
# for str_post in range(-15,70):
# #str_post=-1



# 	str_seq=seq_i[73-18:str_post+73+1]
# 	os.system('./gen_dna.x %s dna.pdb'%str_seq)

# 	rc('open 1kx5.pdb')
# 	rc('delete #0:HOH')
# 	rc('delete #0:CL')
# 	rc('delete #0:-73-%d.I,%d-73.J'%(str_post-1,(str_post-1)*(-1)))
# 	rc('delete #0:MN')


# 	rc('open dna.pdb')
# 	rc('resrenumber -18 #1:.I')
# 	rc('resrenumber %d #1:.J'%(str_post*(-1)))
# 	rc('write format pdb #1 dna.pdb')

# 	rc('open 2o5i_renum.pdb')
# 	rc('delete #2:HOH')
# 	rc('delete #2:MN')
# 	rc('delete #2:ZN.J')


# 	rc('changechains A,B,C,D,E,H K,L,M,N,O,P #2')

# 	##Now let's do superpositions
# 	rc('match #1:-18.I,18.J@N1,N2,N3,N4,N7,N9 #2:-18.I,18.J@N1,N2,N3,N4,N7,N9 ')

# 	rc('match #0:%d.I,%d.J@N1,N2,N3,N4,N7,N9 #1:%d.I,%d.J@N1,N2,N3,N4,N7,N9 '%(str_post,str_post*(-1),str_post,(-1)*str_post))





# 	rc('delete #1:-18.I,18.J')
# 	rc('delete #1:%d.I,%d.J'%(str_post,(-1)*str_post))


# 	rc('combine #0,1,2 newchainids false close false')

# 	#Let's add bonds
# 	rc('bond #3:-18.I@O3\':-17.I@P')
# 	rc('bond #3:%d.I@O3\':%d.I@P'%(str_post-1,str_post))
# 	rc('bond #3:18.J@P:17.J@O3\'')
# 	rc('bond #3:%d.J@P:%d.J@O3\''%((-1)*str_post+1,(-1)*str_post))


# 	rc('write format pdb #3 complex.pdb')

# 	os.system('cp complex.pdb models/complex_%d.pdb'%str_post)

# 	rc('close session')
# # # gather the names of .pdb files in the folder
# # file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

# # # loop through the files, opening, processing, and closing each in turn
# # for fn in file_names:
# # 	replyobj.status("Processing " + fn) # show what file we're working on
# # 	rc("open " + fn)
# # 	rc("align ligand ~ligand") # put ligand in front of remainder of molecule
# # 	rc("focus ligand") # center/zoom ligand
# # 	rc("surf") # surface receptor
# # 	rc("preset apply publication 1") # make everything look nice
# # 	rc("surftransp 15") # make the surface a little bit see-through
# # 	# save image to a file that ends in .png rather than .pdb
# # 	png_name = fn[:-3] + "png"
# # 	rc("copy file " + png_name + " supersample 3")
# # 	rc("close all")
# # uncommenting the line below will cause Chimera to exit when the script is done
# #rc("stop now")
# # note that indentation is significant in Python; the fact that
# # the above command is exdented means that it is executed after
# # the loop completes, whereas the indented commands that 
# # preceded it are executed as part of the loop.
