# -*- coding: utf-8 -*-
"""
This library outputs MSA to HTML and allows annotation.
This is my take on mview.

"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
import re
from hist_ss import get_hist_ss_in_aln_for_html

restypedict={'A':'hphob','C':'hphob','D':'neg','E':'neg','F':'hphob','G':'pol','H':'pol','I':'hphob','K':'pos','L':'hphob','M':'hphob','N':'pol','P':'hphob','Q':'pol','R':'pos','S':'pol','T':'pol','V':'hphob','W':'hphob','Y':'pol'}

def aln2html(msa,filename,features=None,title=None,description=True,field1w=20,field2w=35):
    """ 
    This function outputs HTML from msa and annotates features.
    msa - Biopython MSA,
    filename - html file to output the result.
    features - a dictionary of features, organized as follows:
    {(begin,end):{'level':0(default),'symbol':'H','description':'desc'}}
    if features overlap and not levels, they will be split to different levels.
    Only three levels (0,1,2) are available.
    """


    style="""
pre,td{margin: 0px;padding: 0px;border: 0px;}
.pos{color:blue;}
.neg{color:red;}
.pol{color:green;}
.hphob{color:grey;}
.def{color:black;}
.conserved{background:lightblue;}
.nonconserved{background:white;}
"""
    sinfo=SummaryInfo(msa)
    cons=sinfo.gap_consensus(threshold=0.9, ambiguous='X')
    
    #Let's work on features
    f_description=''
    msatext=''
    annot_line=[0,1,2]
    if(features):
        annot_line[0]=list(' '*len(cons))
        annot_line[1]=list(' '*len(cons))
        annot_line[2]=list(' '*len(cons))
        keys=sorted(list(features.keys()),key=lambda x: x[0])
        for k in keys:
            if(features[k].get('description',0)):
                f_description+='{0}-{1};'.format(features[k]['symbol'],features[k]['description'])
            lev=features[k].get('level',0)
            if(re.match('^\s+$',''.join(annot_line[lev][k[0]:k[1]+1]))):
                annot_line[lev][k[0]:k[1]+1]=features[k]['symbol']*(k[1]-k[0]+1)
            else:
                lev+=1
                if(re.match('^\s+$',''.join(annot_line[lev][k[0]:k[1]+1]))):
                    annot_line[lev][k[0]:k[1]+1]=features[k]['symbol']*(k[1]-k[0]+1)
                else:
                    lev+=1
                    if(re.match('^\s+$',''.join(annot_line[lev][k[0]:k[1]+1]))):
                        annot_line[lev][k[0]:k[1]+1]=features[k]['symbol']*(k[1]-k[0]+1)
        if(not re.match('^\s+$',''.join(annot_line[2]))):
            msatext='<TR><TD><PRE>{0:<{field1w}}</PRE></TD>'.format('annotation',field1w=field1w+2)
            if(description):
                msatext+='<TD><PRE>{0:<{field2w}}</PRE></TD>'.format('level 2',field2w=field2w+2)
            for c in annot_line[2]:
                msatext+='<TD><PRE>{0}</PRE></TD>'.format(c)
            msatext+='</TR>'
        if(not re.match('^\s+$',''.join(annot_line[1]))):
            msatext='<TR><TD><PRE>{0:<{field1w}}</PRE></TD>'.format('annotation',field1w=field1w+2)
            if(description):
                msatext+='<TD><PRE>{0:<{field2w}}</PRE></TD>'.format('level 1',field2w=field2w+2)
            for c in annot_line[1]:
                msatext+='<TD><PRE>{0}</PRE></TD>'.format(c)
            msatext+='</TR>'
        if(not re.match('^\s+$',''.join(annot_line[0]))):
            msatext='<TR><TD><PRE>{0:<{field1w}}</PRE></TD>'.format('annotation',field1w=field1w+2)
            if(description):
                msatext+='<TD><PRE>{0:<{field2w}}</PRE></TD>'.format('level 0',field2w=field2w+20)
            for c in annot_line[0]:
                msatext+='<TD><PRE>{0}</PRE></TD>'.format(c)
            msatext+='</TR>'
        f_description+='<BR><BR>'



    for s in msa:
        if(re.search(r'\d\d+',s.id)):
            gi=re.search(r'(\d\d+)',s.id).group(1)
            line='<TR><TD><PRE><a href="http://www.ncbi.nlm.nih.gov/protein/?term={0}">{1:<{field1w}}</a></PRE></TD>'.format(gi,s.id[:field1w],field1w=field1w+2)
            if(description):
                line+='<TD><PRE>{0:<{field2w}}</PRE></TD>'.format(s.description[:field2w],field2w=field2w+2)
        else:
            line='<TR><TD><PRE>{0:<{field1w}}</PRE></TD>'.format(s.id[:field1w],field1w=field1w+2)
            if(description):
                line+='<TD><PRE>{0:<{field2w}}</PRE></TD>'.format(s.description[:field2w],field1w=field1w+2)
        for c,i in zip(s.seq,range(len(s.seq))):
            line+='<TD><PRE class="{0} {1}">{2}</PRE></TD>'.format(restypedict.get(c,'def'),'conserved' if c==cons[i] and c!='-' else 'nonconserved',c)
        line+='</TR>'
        msatext=msatext+line


    a=open(filename,'w')
    a.write("""
<!DOCTYPE html>
<HTML>
<HEAD>
<META http-equiv="Content-Type" content="text/html; charset=utf-8"/>
<TITLE>MultipleSequenceAlignment</TITLE>
<style>
{style}
</style>
</HEAD>
<BODY style="background-color:white; color:black; a:link:blue; a:active:red; a:visited:purple">
{title}<BR><BR>
{features}
<TABLE style="border:0px; border-spacing:0px; background-color:white; color:black; a:link:blue; a:active:red; a:visited:purple;">

{msatext}

</TABLE>
</BODY>
</HTML>
""".format(\
title=title,\
msatext=msatext,\
style=style,\
features=f_description
))

    a.close()

    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    test_h2a_core=Seq('TRSTRAHLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2A.Z',name='H2A.Z'),SeqRecord(test_h2a_core,id='test',name='test')])
    f=get_hist_ss_in_aln_for_html(msa)
    print f
    aln2html(msa,"test.html",title="MyMSA",features=f)
    # get_pdf('H2A',MultipleSeqAlignment([SeqRecord(human_h2a_z_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='1H2A.Z',name='H2A.Z')]),'H2AvsH2A.Z',[0,5,1])



    
            