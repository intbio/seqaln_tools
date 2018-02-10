# -*- coding: utf-8 -*-
"""
Visualizing histone MSAs via TEXSHADE.
A very powerful example.

"""

from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import cPickle as pickle
from Bio import SeqIO

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

from StringIO import StringIO

from hist_ss import get_hist_ss
from hist_ss import get_hist_ss_in_aln
from Bio.Align.AlignInfo import SummaryInfo
import L_aln_tools
#from pylab import *
 
Entrez.email = "alexey.shaytan@nih.gov" 

PATH_to_NCBI_nodes_dmp="nodes.dmp"

TEMP_DIR='/Users/alexeyshaytan/junk/temptemp'


def write_texshade(file_handle,aln_fname,features,res_per_line=120,showlegend=True,shading_modes=['similar'],logo=False,hideseqs=False,setends=[],ruler=False,numbering_seq='consensus',hide_ns=False):

    for shading in shading_modes:
        shading=str(shading)
        file_handle.write("""
    \\begin{texshade}{%s}
    \\residuesperline*{%d}
    """%(aln_fname,res_per_line))
       

        file_handle.write(r"""
    \seqtype{P}
    \defconsensus{{}}{*}{upper}
    """)
        #a very dirty hack
        if(setends):
            if(numbering_seq=='consensus'):
                file_handle.write("""
    \\setends[%d]{%s}{%d..%d}
    """%(setends[0]+1,numbering_seq,setends[0]+setends[0],setends[1]+setends[0]+setends[0]-1))
            else:
                    file_handle.write("""
    \\setends{%s}{%d..%d}
    """%(numbering_seq,setends[0],setends[1]))
           
        if(ruler):
            file_handle.write("""
    \\showruler{top}{%s}
    """%numbering_seq)
        if(hide_ns):
            file_handle.write("""
    \\hideseq{%s}
    """%numbering_seq)

        file_handle.write(r"""
    \seqtype{P}
    """)
        if((shading=='similar')|(shading=='0')):
            file_handle.write(r"""
    \shadingmode{similar}
    \threshold[80]{50}
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo{top} \showlogoscale{leftright}
    \namesequencelogo{logo}
    """)

        if((shading=='hydropathy_functional')|(shading=='1')):
            file_handle.write(r"""
    \shadingmode[hydropathy]{functional}
    \shadeallresidues
    \threshold[80]{50}
    
    """)     
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[hydropathy]{top} \showlogoscale{leftright}

    """)
        if((shading=='chemical_functional')|(shading=='2')):
            file_handle.write(r"""

    \shadingmode[chemical]{functional}
    \shadeallresidues
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[chemical]{top} \showlogoscale{leftright}

    """)
        if((shading=='structure_functional')|(shading=='3')):
            file_handle.write(r"""
    \shadingmode[structure]{functional}
    \shadeallresidues
    """)

        if((shading=='charge_functional')|(shading=='4')):
            file_handle.write(r"""
    \shadingmode[charge]{functional}
    \shadeallresidues
    """)

        if((shading=='diverse')|(shading=='5')):
            file_handle.write(r"""
    \shadingmode{diverse}

    """)
       
        if(hideseqs):
            file_handle.write(r"""
    \hideseqs

    """)
        file_handle.write(r"""

    %\setends{consensus}{1..160}
    %\setends{consensus}{1..160}


    %\feature{ttop}{1}{1..160}{bar:conservation}{}
    %\showfeaturestylename{ttop}{conserv.}
    \ttopspace{-\baselineskip}

    %\feature{top}{1}{1..160}{color:charge}{}
    %\showfeaturestylename{top}{charge}

    %\feature{bottom}{1}{1..160}{color:molweight[ColdHot]}{}
    %\showfeaturestylename{bottom}{weight}

    %\bbottomspace{-\baselineskip}
    %\feature{bbottom}{2}{1..160}{bar:hydrophobicity[Red,Gray10]}{}
    %\showfeaturestylename{bbottom}{hydrophob.}

    %\bargraphstretch{3}
    %\featurestylenamescolor{Red}
    %\featurestylenamesrm  \featurestylenamesit

    %\showsequencelogo{top}


    %\showconsensus[ColdHot]{bottom}
    \showconsensus[black]{top}

    %\defconsensus{.}{lower}{upper}
    %\defconsensus{{}}{lower}{upper}
    %\defconsensus{{}}{*}{upper}

    """)
        file_handle.write(features)
        file_handle.write(r"""
    \featurerule{1mm}""")
        if(showlegend):
            file_handle.write(r"""
    \showlegend""")

        align = AlignIO.read(open(aln_fname,'r'), "fasta")
        for a,i in zip(align,range(len(align))):
            # print a.id.replace('|',' | ')
            file_handle.write("""
    \\nameseq{%d}{%s}"""%(i+1,a.id.replace('|',' | ')))

        file_handle.write(r"""

    \end{texshade}
    """)





def get_pdf(hist_name,align,title,shading_modes=['similar'],logo=False,hideseqs=False,splitN=20,setends=[],ruler=False):

    # align.sort()
    # a=SummaryInfo(align)
    ns='consensus'
    a_len=len(align)
    num=int(a_len/splitN)+1
    if((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    if((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    if((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    print a_len, splitN
    for i in range(num):
        #if(setends or ruler):
        #    ns='techconsensus'
        #    t_aln=L_aln_tools.add_consensus(align[(i*splitN):((i+1)*splitN)],threshold=0.01, ambiguous='X',name=ns)
        #else:
        t_aln=align[(i*splitN):((i+1)*splitN)]
        AlignIO.write(t_aln,open(TEMP_DIR+'/alignment%d.fasta'%i,'w'), 'fasta')
    # align.append(SeqRecord(a.dumb_consensus(threshold=0.1, ambiguous='X'),id='Zcons'))
    # AlignIO.write(align, open('data/alignment.fasta','w'), 'fasta')
    res_per_line=len(align[0])

    #Let's make some drawing
    hv,ss=get_hist_ss_in_aln(align,debug=1)
    # print hv,ss
    #prepare feature section
    features=''
    # features='\\feature{top}{1}{1..%d}{---}{loop}'%len(align[0])
    for i in ss:
        if(re.search('alpha',i)):
            features=features+"\\feature{tttop}{consensus}{%d..%d}{helix}{%s}"%(ss[i][0]+1,ss[i][1]+1,i)
        if(re.search('beta',i)):
            features=features+"\\feature{tttop}{consensus}{%d..%d}{-->}{%s}"%(ss[i][0]+1,ss[i][1]+1,i)
        if(re.search('loop',i)):
            features=features+"\\feature{ttttop}{consensus}{%d..%d}{loop}{%s}"%(ss[i][0]+1,ss[i][1]+1,i)
        if(re.search('domain',i)):
            features=features+"\\feature{ttttop}{consensus}{%d..%d}{loop}{%s}"%(ss[i][0]+1,ss[i][1]+1,i)
        if(re.search('tail',i)):
            features=features+"\\feature{ttttop}{consensus}{%d..%d}{loop}{%s}"%(ss[i][0]+1,ss[i][1]+1,i)
        if(re.search('mgarg',i)):
            features=features+"\\frameblock{consensus}{%d..%d}{Red[1.5pt]}"%(ss[i][0]+1,ss[i][1]+1)
    a=open(TEMP_DIR+'/align.tex','w')

    a.write(r"""\documentclass[11pt,landscape]{article}
%\documentclass{standalone}
%\usepackage[a0paper]{geometry}
\usepackage{hyperref}
""")

    a.write("""
\\usepackage[paperwidth=%fin, paperheight=18in]{geometry}
        """%(22/200.*res_per_line+2.5))

    a.write(r"""
\usepackage{texshade}

\begin{document}""")
    a.write("""
\\Huge{%s}"""%title)

    for i in range(num-1):
        write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%i , features, res_per_line,False,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False)
    write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%(num-1) , features, res_per_line,True,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False)

    a.write(r"""
\end{document} """)
    a.close()

    command='/usr/texbin/pdflatex --file-line-error --synctex=1 -output-directory=%s --save-size=10000  %s/align.tex > /dev/null'%(TEMP_DIR,TEMP_DIR)

    print('Launcning command:')
    print(command)
    os.system(command)
    os.system('mv '+TEMP_DIR+'/align.pdf %s.pdf'%title.replace(' ','_'))



    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')

    get_pdf('H2A',MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2A.Z',name='H2A.Z')]),'H2AvsH2A.Z',[0,5,1],True,True)
    # get_pdf('H2A',MultipleSeqAlignment([SeqRecord(human_h2a_z_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='1H2A.Z',name='H2A.Z')]),'H2AvsH2A.Z',[0,5,1])



    
            