# -*- coding: utf-8 -*-
"""
Taxonomy tools.
Functions of version 1 - were implemented using networkx
A better implementation using ETE2 is now available - see below.

"""
import os
import sys
import cPickle as pickle
import pandas as pd
import networkx as nx
from ete2 import NCBITaxa

# Entrez.email = "alexey.shaytan@nih.gov" 

PATH_to_NCBI_nodes_dmp="nodes.dmp"


os.environ['PATH']+=os.path.sep+'/Users/alexeyshaytan/soft/x3dna-v2.1/bin:/Users/alexeyshaytan/soft/amber12/bin:/Users/alexeyshaytan/soft/sratoolkit/bin:/Users/alexeyshaytan/soft/bins/gromacs-4.6.3/bin:/opt/local/bin:/opt/local/sbin:/Users/alexeyshaytan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/ncbi/blast/bin:/usr/texbin'

#This function is DEPRECATED
def check_tax_id_clade(clade_top_tax_id,check_tax_id):
    """Checks if given tax_id is inside a clade formed by taxa described by its top taxid"""

    p=pd.read_table(PATH_to_NCBI_nodes_dmp,sep='|',usecols=[0,1],header=None)
    G=nx.DiGraph()
    #Load all taxonomy as graph
    G.add_edges_from(zip(p.ix[:,1],p.ix[:,0]))
    clade=nx.dfs_tree(G,clade_top_tax_id)
    if(not clade.nodes()):
        clade.add_node(clade_top_tax_id)
    return(clade.has_node(check_tax_id))


ncbi = NCBITaxa()


def subsample_taxids(taxids,rank='species'):
    """
    For a given set of taxids leaves only one representative per selected rank
    Eg. for a set of subspecies - leave only species.
    """
    rank_dict={'superkingdom':0,'kingdom':1,'phylum':2,'class':3,'superorder':4,'order':5,'suborder':6,'infraorder':7,'parvorder':8,'superfamily':9,'family':10,'subfamily':11,'genus':12,'subgenus':13,'species':14,'subspecies':15}

    tree = ncbi.get_topology(taxids,intermediate_nodes=True)
    #We have now a phylogenetic tree with all annotations for our taxids.
    subsampled_taxids=set()
    #We are iterating through the taxids and for every we are determining only one representative from this group.
    #These representatives will be the same for the taxids in one group - and hence subsampling will happen.
    for t in taxids:
        #From a taxid we need to go up the tree till we reach the desired rank.
        node=tree.search_nodes(name=str(t))[0]
        # while node.rank!=rank:
        #     node=node.up
        while rank_dict.get(node.rank,100)>rank_dict.get(rank):
            node=node.up
        #And now we go down taking the first child taxid.
        while str(node.name) not in map(str,taxids):
            node=node.children[0]
        subsampled_taxids.add(int(node.name))
    return list(subsampled_taxids)


# def group_taxids(taxids,rank='species'):
#     """
#     For a given set of taxids group them, so that each group covers a single rank clade
#     Eg. for a set of subspecies - group by species.
#     """
#     rank_dict={'superkingdom':0,'kingdom':1,'phylum':2,'class':3,'superorder':4,'order':5,'suborder':6,'infraorder':7,'parvorder':8,'superfamily':9,'family':10,'subfamily':11,'genus':12,'subgenus':13,'species':14,'subspecies':15}

#     tree = ncbi.get_topology(taxids,intermediate_nodes=True)
#     #We have now a phylogenetic tree with all annotations for our taxids.
#     grouped_taxids=list()
#     print tree.get_ascii(attributes=["name", "sci_name", "taxid","rank"])
#     last_known_rank=0
#     for t in tree.traverse():
#         #sometimes the taxonomy lacks certain level of classification, here is a fix
#         #needs to be tested
#         flag = 0
#         if(last_known_rank < rank_dict.get(rank)) & (rank_dict.get(t.rank,-100) > rank_dict.get(rank)):
#             flag=1     
 
#         if t.rank==rank or flag:
#             taxa=set([int(t.name)])
#             dec=set([int(i.name) for i in t.get_descendants()])
#             taxa.update(dec)
#             grouped_taxids.append(list(taxa.intersection(map(int,taxids))))

#         if(rank_dict.get(t.rank,-1)>-1):
#             last_known_rank=rank_dict.get(t.rank)
#         # print "rank set",last_known_rank
#     return grouped_taxids

def group_taxids(taxids,rank='species'):
    """
    For a given set of taxids group them, so that each group covers a single rank clade
    Eg. for a set of subspecies - group by species.
    """
    rank_dict={'superkingdom':0,'kingdom':1,'phylum':2,'class':3,'superorder':4,'order':5,'suborder':6,'infraorder':7,'parvorder':8,'superfamily':9,'family':10,'subfamily':11,'genus':12,'subgenus':13,'species':14,'subspecies':15}
    inv_rank_dict={v:k for k,v in rank_dict.iteritems()}
    tree = ncbi.get_topology(taxids,intermediate_nodes=True)
    #We have now a phylogenetic tree with all annotations for our taxids.
    grouped_taxids=list()
    class_taxids=set()
    # print tree
    rnodes=tree.search_nodes(rank=str(rank))
    # print rnodes[0]
    for n in rnodes:
        taxa=set([int(n.name)])
        dec=set([int(i.name) for i in n.get_descendants()])
        taxa.update(dec)
        grouped_taxids.append(list(taxa.intersection(map(int,taxids))))
        class_taxids.update(grouped_taxids[-1])
    #The problem is some taxids are still unclassified due to lack of respective levels in tree
    non_class=set(taxids).difference(class_taxids)
    if(len(non_class)>0):
        print "non classified taxids left:",non_class

        #try to do it with lower rank
        new_rank=inv_rank_dict[rank_dict[rank]+1]
        print "reclassifying them at level: ",new_rank
        grouped_taxids.extend(group_taxids(list(non_class),new_rank))        
    return grouped_taxids

def get_taxid_from_gbrec(seqrec):
    """
    Get taxid from GB seqrecrod in SeqRec Format
    """
    sources=[]
    for f in seqrec.features:
        if f.type=='source':
            if f.qualifiers['db_xref'][0][0:6]=='taxon:':
                sources.append(f.qualifiers['db_xref'][0][6:])
    if len(sources)>0:
        return sources[0]
    else:
        return None

def gbrec_to_clickhtml(seqreclist,filename='seqrecdata.html',title=''):
    """
    Outputs an html file with seqrecord information and with clickable links.
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
    text=''
    for s in seqreclist:
        line='<TR><TD><PRE><a href="http://www.ncbi.nlm.nih.gov/protein/?term={0}">{1:<{field1w}}</a></PRE></TD>'.format(s.id,s.id,field1w=17)
        line+='<TD><PRE>{0}</PRE></TD>'.format(s.annotations['organism'])
        line+='<TD><PRE>{0}</PRE></TD>'.format(s.description)

        line+='</TR>'
        text=text+line


    a=open(filename,'w')
    a.write("""
<!DOCTYPE html>
<HTML>
<HEAD>
<META http-equiv="Content-Type" content="text/html; charset=utf-8"/>
<TITLE>SeqRecData</TITLE>
<style>
{style}
</style>
</HEAD>
<BODY style="background-color:white; color:black; a:link:blue; a:active:red; a:visited:purple">
{title}<BR><BR>

<TABLE style="border:0px; border-spacing:0px; background-color:white; color:black; a:link:blue; a:active:red; a:visited:purple;">

{text}

</TABLE>
</BODY>
</HTML>
""".format(\
title=title,\
text=text,\
style=style
))

    a.close()

def sort_taxids(taxids):
    """
    ideas is to sort taxids so that they hierarchically cluster as in tree.
    Good to group sequences of similar species together.
    """
    newtaxids=[]
    tree = ncbi.get_topology(taxids,intermediate_nodes=False)
    #Start traversing tree
    for node in tree.traverse("preorder"):
        if(node.name in taxids):
            for i in range(taxids.count(node.name)):
                newtaxids.append(node.name)
    return newtaxids





# def tax_filter_hist(hist_dataframe,clade_top_tax_id):
#     """Should take a PANDAS list of histones and leave only those that are inside taxa"""
    
#     p=pd.read_table(PATH_to_NCBI_nodes_dmp,sep='|',usecols=[0,1],header=None)
#     G=nx.DiGraph()
#     #Load all taxonomy as graph
#     G.add_edges_from(zip(p.ix[:,1],p.ix[:,0]))
#     clade=nx.dfs_tree(G,clade_top_tax_id)
#     if(not clade.nodes()):
#         clade.add_node(clade_top_tax_id)
#     return(hist_dataframe[map(clade.has_node,hist_dataframe.TaxID)])
#     #pylab.plot(prof)
if __name__ == '__main__':
    #some examples

    tids=[9606,639525]

    print group_taxids(tids,rank='class')
    # print ncbi.get_descendant_taxa(tree.name)

    exit()
    #1. Getting data
    #################
    hist_df=pd.read_csv('inp_data/seqs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences

    #2. Filter dataframe by histone variant
    #################
    f_hist_df=hist_df[(hist_df['hist_var']=='canonical_H2B')]

    #3. Select one variant per taxid
    #################
    f_hist_df=f_hist_df.drop_duplicates(['taxid','hist_var'])
    
    #4. Filter by list of taxonomy clades   
    ################
    parent_nodes=[9443] #131567 - cellular organisms
    taxids=list()
    for i in parent_nodes:
        taxids.extend(ncbi.get_descendant_taxa(i))

    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(taxids)]

    #5. Take one representative per species or specific rank.
    ################
    #Common ranks: superorder-order-suborder-infraorder-parvorder-superfamily-family-subfamily-genus-species-subspecies
    seqtaxids=list(f_hist_df['taxid']) #old list
    new_seqtaxids=subsample_taxids(seqtaxids,rank='family') #new subsampled list
    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(new_seqtaxids)] #remake the dataframe
    
    #---------------
    #Output tree before subsampline
    tree = ncbi.get_topology(seqtaxids,intermediate_nodes=False)
    print tree.get_ascii(attributes=["sci_name", "rank","taxid"])
    #Output after subsampling
    tree = ncbi.get_topology(new_seqtaxids,intermediate_nodes=False)
    print tree.get_ascii(attributes=["sci_name", "rank","taxid"])
    #Output new and old sets
    print seqtaxids
    print new_seqtaxids

    #TODO:#6. Filter sequences by their quality: align to known histones and remove those that have strange indels.
    ################ We can use the features.json file ideally.
    f_fasta_dict={key: value for (key,value) in fasta_dict.iteritems() if key in list(f_hist_df['gi'])} # get fasta dict


    
            