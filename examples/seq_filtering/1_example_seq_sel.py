# -*- coding: utf-8 -*-
"""
This is an exmple script that interactively explores H2A
alignments and conservation using tools and examples developed in libraries.
First of all it employs an elaborate sequence selection mechansm.
Makes alignments in fasta, html and SVG format with ete2, including
taxonomic tree.
"""
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "libs"))
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS/')

from L_aln_tools import muscle_aln, trim_aln_gaps, cons_prof, add_consensus, cluster_seq_support
from L_hist_aln import *
from L_seq_subset import *
from L_shade_aln import *
from L_plot4seq import *
from L_taxonomy_tools import *
from hist_ss import get_hist_ss_in_aln_for_html, get_hist_ss_in_aln_as_string
from L_aln2html import aln2html
from ete2 import NCBITaxa
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node,AttrFace,TextFace
from Bio.Align import MultipleSeqAlignment
from math import log
from itertools import chain


ncbi = NCBITaxa()

def main():
    title=''
    #1. Getting data
    ########################################################
    ########################################################
    # df=pd.read_csv('int_data/seqs_rs_redef.csv') #Histone types info #Does not really seem that we need to redefine variants based on best score.
    df=pd.read_csv('int_data/seqs_rs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
    
    #2. Filtering - filter initial dataset by type, variant and other parameters
    ########################################################
    ########################################################

    #2.1. Narrow by variant/type
    ########################################################
    title+='H2A'
    # f_df=df[(df['hist_var']=='canonical_H4')]
    # f_df['hist_var']='canonical_H4'
    f_df=df[((df['hist_var']=='canonical_H2A')|(df['hist_var']=='H2A.X'))&(df['partial']==False)&(df['non_st_aa']==False)]
    # f_df=df[((df['hist_var']=='H2A.Z'))&(df['partial']==False)&(df['non_st_aa']==False)]

    # f_df=df[(df['hist_type']=='H2A')]

    print "Number of seqs after narrowing by hist type/var:", len(f_df)
    

    #2.2. Filter by list of taxonomy clades - restrict sequences to certain taxonomic clades
    #########################################################
    title+=' across cellular organisms'
    # parent_nodes=[9443] #131567 - cellular organisms, 7215 4930 Drosophila and yeast, 9443 - primates
    parent_nodes=[131567] #131567 - cellular organisms, 7215 4930 Drosophila and yeast, 9443 - primates
    #33682 - euglenozoa
    #6656 - arthropods
    # 4751 - fungi
    #5782 - dictostelium
    #This is akin manual removal of bad species
    del_nodes=[5782,5690]

    print "Selecting taxonomic subset for taxids: ",parent_nodes
    print "while removing taxonomic subset for taxids: ",del_nodes

    taxids=set(parent_nodes)
    for i in parent_nodes:
        taxids.update(ncbi.get_descendant_taxa(i,intermediate_nodes=True))
    for i in del_nodes:
        taxids=taxids.difference(set([i]))
        taxids=taxids.difference(set(ncbi.get_descendant_taxa(i,intermediate_nodes=True)))

    f_df=f_df[f_df['taxid'].isin(taxids)]
    print "Number of seq after taxonomic subset: ",len(f_df)
    


    #2.3.0 Marking number of identical sequence within each species and subspecies.
    #This will simplify further analysis of sequence filtering on similarity
    #We know that all refseqs are duplicated for instance.
    ################################################
    ident=dict()
    new_gis=list()
    tids=set(list(f_df['taxid']))
    for i in tids:
        # print i.name, i.sci_name
        temp_df=f_df[(f_df['taxid']==i)]
        gis=list(temp_df['gi']) #this is to limit exec time
        # print gis
        if(len(gis)>1):
            res=cluster_seq_support({gi:fasta_dict[str(gi)] for gi in gis},ident_thresh=1.00)
            ident.update(res)
        else:
            ident.update({gis[0]:1})

    f_df['ident']=[ident.get(k,1) for k in f_df['gi']]
    #where ident - number of identical sequnces for current sepecies/subspecies.
    print "Identity of sequence inside each taxid determined"

    #2.3.1. Calculate number of similar seqs for every seq in tax group
    #########################################################
    # Use powerful method, to get rid of random errors is to identify identical sequences
    # if a sequence is supported by two or more entires - this is good.
    # Here we add a degen column to our data set - showing how many similar sequences are found
    # for a given sequence in its taxonomic clade (genus currently) 

    #We will traverse the species tree by species, genus or family, and determine degeneracy level
    degen=dict()
    new_gis=list()
    tids=list(f_df['taxid']) 
    t = ncbi.get_topology(tids,intermediate_nodes=True)
    for i in t.search_nodes(rank='family'):
        # print i.name, i.sci_name
        nodeset=list()
        for k in i.traverse():
            nodeset.append(int(k.name))
        temp_df=f_df[(f_df['taxid'].isin(nodeset))]
        gis=list(temp_df['gi']) #this is to limit exec time
        # print gis
        res=cluster_seq_support({gi:fasta_dict[str(gi)] for gi in gis},ident_thresh=1.00)
        degen.update(res)

    # print degen
    f_df['degen']=[degen.get(k,1) for k in f_df['gi']]

    #2.3.2. Remove seqs that do not have support outside their species
    # if they are not curated or RefSeq NP.
    ###########################################################

    f_df=f_df.sort(['RefSeq','degen'],ascending=False) # so that RefSeq record get priority on removing duplicates
    f_df=f_df[(f_df['degen']>f_df['ident'])|(f_df['curated']==True)|(f_df['RefSeq']==2)]
    print "After removing mined seqs with no support in neighboring species: ",len(f_df)

    #2.3.3. Shuffle sequnces, so that upon further selection, RefSeq and high degeneracy get priority
    ###########################################################
    #RefSeq and degenerate sequence get priority
    # title+=' 1ptax'
    f_df=f_df.sort(['RefSeq','degen'],ascending=False) # so that RefSeq record get priority on removing duplicates
    # print f_df[0:10]
    # f_df=f_df.drop_duplicates(['taxid','hist_var'])


    #2.4 Take one best representative per specific taxonomic rank (e.g. genus)
    ############################################################
    pruningrank='genus'
    print "Pruning taxonomy by ", pruningrank
    
    title+=' , one seq. per %s'%pruningrank
    #Common ranks: superorder-order-suborder-infraorder-parvorder-superfamily-family-subfamily-genus-species-subspecies
    seqtaxids=list(f_df['taxid']) #old list
    grouped_taxids=group_taxids(seqtaxids,rank=pruningrank)
    # print seqtaxids
    # print grouped_taxids
    #Now we need to take best representative
    #refseq NP, curated, or the one with largest degeneracy
    new_gis=list()
    for tids in grouped_taxids:
        t_df=f_df[f_df['taxid'].isin(tids)]
        #try take curated first
        if(len(t_df[t_df['curated']==True])>0):
            new_gis.append(t_df.loc[t_df.curated==True,'gi'].values[0])
            continue
        #try take NP records nest
        #RefSeq 2 means NP, 1 means XP
        if(len(t_df[t_df['RefSeq']==2])>0):
            new_gis.append(t_df.loc[t_df.RefSeq==2,'gi'].values[0])
            continue
        # take best degenerate otherwise
        else:
            t_df=t_df.sort(['degen','RefSeq'],ascending=False) 
            new_gis.append(t_df['gi'].iloc[0])

    f_df=f_df[f_df['gi'].isin(new_gis)]

    print "After pruning taxonomy we have: ",len(f_df)


    #2.5. Check seq for sanity - needs to be checked!
    ##############################################
    # title+=' seqQC '

    # print "Checkig sequence quality"
    # newgis=list()
    # for i,row in f_df.iterrows():
    #     gi=row['gi']
    #     seq=fasta_dict[str(gi)].seq
    #     hist_type=row['hist_type']
    #     hist_var=row['hist_var']
    #     if(check_hist_length(seq,hist_type,hist_var,5)&check_hist_core_length(seq,hist_type,5)):
    #         newgis.append(gi)
    # f_df=f_df[f_df['gi'].isin(newgis)] #remake the dataframe
    # print len(f_df)

    #3. Make a list of seq with good ids and descriptions
    ##############################################

    f_fasta_dict={key: value for (key,value) in fasta_dict.iteritems() if int(key) in list(f_df['gi'])}
    print len(f_fasta_dict)
    taxid2name = ncbi.get_taxid_translator(list(f_df['taxid']))
    #Relabel sequences gi=> type and organism
    f_fasta_dict={key: SeqRecord(id=key, description=f_df.loc[f_df.gi==int(key),'hist_var'].values[0]+' '+taxid2name[f_df.loc[f_df.gi==int(key),'taxid'].values[0]],seq=value.seq) for (key,value) in f_fasta_dict.iteritems() }
    #with arbitrary index
    # f_fasta_dict_rel={key: SeqRecord(id=str(index), description=f_hist_df.loc[f_hist_df.gi==key,'hist_var'].values[0]+' '+taxid2names[f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]],seq=f_fasta_dict[key].seq) for (index,key) in enumerate(f_fasta_dict) }
    # exit()

    #4. Make MSA
    #################
    #Here we construct MSA
    msa=muscle_aln(f_fasta_dict.values(),gapopen=float(-20))
    AlignIO.write(msa, "int_data/example_msa.fasta", "fasta")

    msa_annot=MultipleSeqAlignment([SeqRecord(Seq(''.join(get_hist_ss_in_aln_as_string(msa)).replace(' ','-')),id='annotation',name='')])
    msa_annot.extend(msa)
    AlignIO.write(msa_annot, "int_data/example_msa_annot.fasta", "fasta")

    for i in range(len(msa)):
        gi=msa[i].id
        msa[i].description=f_fasta_dict[gi].description.replace('canonical','ca')
    msa.sort(key=lambda x: x.description)


    #5. Visualize MSA############
    aln2html(msa,'example_h2a.html',features=get_hist_ss_in_aln_for_html(msa,'H2A',0),title="canonical H2A alignment",description=True,field1w=10,field2w=35)

    #6. Trim alignment - this is optional
    #6.1. Trim gaps
    # title+=' gaptrim'
    # msa_tr=trim_aln_gaps(msa,threshold=0.8)

    #6.2. Trim to histone core sequence
    msa_tr=trim_hist_aln_to_core(msa)
    # msa_tr=msa
    # print get_hist_ss_in_aln_for_shade(msa_tr,below=True)

    # exit()

    #7. Vizualize MSA with ete2.##########
    taxid2gi={f_df.loc[f_df.gi==int(gi),'taxid'].values[0]:gi for gi in list(f_df['gi'])}
    gi2variant={gi:f_df.loc[f_df.gi==int(gi),'hist_var'].values[0] for gi in list(f_df['gi'])}

    msa_dict={i.id:i.seq for i in msa_tr}
    t = ncbi.get_topology(list(f_df['taxid']),intermediate_nodes=False)
    a=t.add_child(name='annotation')
    a.add_feature('sci_name','annotation')
    t.sort_descendants(attr='sci_name')
    ts = TreeStyle()
    def layout(node):
        # print node.rank
        # print node.sci_name
        if getattr(node, "rank", None):
            if(node.rank in ['order','class','phylum','kingdom']):   
                rank_face = AttrFace("sci_name", fsize=7, fgcolor="indianred")
                node.add_face(rank_face, column=0, position="branch-top")
        if node.is_leaf():
            sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
            node.add_face(sciname_face, column=0, position="branch-right")
        if node.is_leaf() and not node.name=='annotation':
            s=str(msa_dict[str(taxid2gi[int(node.name)])])
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            gi=taxid2gi[int(node.name)]
            add_face_to_node(TextFace(' '+str(gi)+' '),node,column=1, position = "aligned")
            add_face_to_node(TextFace('      '+str(int(node.name))+' '),node,column=2, position = "aligned")
            add_face_to_node(TextFace('      '+str(gi2variant[gi])+' '),node,column=3, position = "aligned")

        if node.is_leaf() and node.name=='annotation':
            s=get_hist_ss_in_aln_as_string(msa_tr)
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            add_face_to_node(TextFace(' '+'NCBI_GI'+' '),node,column=1, position = "aligned")
            add_face_to_node(TextFace('       '+'NCBI_TAXID'+' '),node,column=2, position = "aligned")
            add_face_to_node(TextFace('       '+'Variant'+'       '),node,column=3, position = "aligned")



    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(title, fsize=20), column=0)
    t.render("example_motifs_H2A.svg", w=6000, dpi=300, tree_style=ts)

    #10. Conservation############
    #############################
    features=get_hist_ss_in_aln_for_shade(msa_tr,below=True)
    cn=add_consensus(msa_tr,threshold=0.5)[-2:-1]
    # Below are three methods that we find useful.
    # plot_prof4seq('cons_sofp_psic',map(float,cons_prof(msa_tr,f=2,c=2)),cn,features,axis='conservation')
    plot_prof4seq('example_cons_ent_unw',map(lambda x:log(20)+x,map(float,cons_prof(msa_tr,f=0,c=0))),cn,features,axis='conservation')
    plot_prof4seq('example_cons_ent_unw_norm',map(lambda x:log(20)+x,map(float,cons_prof(msa_tr,f=0,c=0,norm="T"))),cn,features,axis='conservation')
    
    # plot_prof4seq('cons_sofp_unw',map(float,cons_prof(msa_tr,f=0,c=2)),cn,features,axis='conservation')
    plot_prof4seq('example_cons_sofp_unw_renorm1',map(float,cons_prof(msa_tr,f=0,c=2,m=1)),cn,features,axis='conservation')
    plot_prof4seq('example_cons_sofp_unw',map(float,cons_prof(msa_tr,f=0,c=2,m=0)),cn,features,axis='conservation')
    plot_prof4seq('example_cons_sofp_psic_renorm1',map(float,cons_prof(msa_tr,f=2,c=2,m=1)),cn,features,axis='conservation')

    
    # plot_prof4seq('cons_ent_psic',map(lambda x:log(20)+x,map(float,cons_prof(msa_tr,f=2,c=0))),cn,features,axis='conservation')

   
#we get an alignment - we need to get conservation profile and visualize it on annotated histone sequence
    
    #11. Subfamily specific sites





    #12.Phylogenetic trees




if __name__ == '__main__':
    main()
    
            