# seqaln_tools
Different python libraries that automate downloading, filtering and analysis of sequences and sequence alignments.
Dig here if you want to know how:
- Download seqeunces from NCBI databases and filter them.
- Display a multiple sequence alignment using TexShade and python wrapper.
- Example scripts of doing homology modeling with Modeller including converting MSA fromats to PIR format.
- Map sequences or there features to a taxonomic tree through ete2
- Calculate mutual information between two sequences.
- Make dataplots on top of seqeunces or sequence alignments.
- Find and display secondary structure elements in sequences.

## L_aln2html.py
This library outputs MSA to HTML and allows some annotation.
Functions: aln2html
Python dependencies: biopython

## L_aln_tools.py
This is a complex and powerful library to retrieve sequences from NCBI and build, edit, analyze and display MSAs.
Function: get_prot_seq_by_gis, get_prot_seqrec_by_gis, get_prot_gbrec_by_gis, aln_undup, cons_prof, trim_aln_gaps, trim_aln_to_key_seq, trim_aln_to_seq, trim_aln_to_seq_length, muscle_aln, add_consensus, cluster_seq_support_nw, cluster_seq_support, taxo_msa, gen_fake_msa, features_via_hmm, taxo_seq_architecture, debump_features   
Python dependencies: biopython, networkx, PyQT, modified ete2 https://github.com/molsim/ete/tree/2.3
External programs: AL2CO, Muscle, EMBOSS
Installation: TBA

## L_fasta2pir.py
This library converts FASTA multiple alignment files to
PIR files as needed by MODELLER by adding additional information.

## L_mut_info.py
Library caclulate mutual information of two sequences
via R.
Sequences should be simple strigns of capital letters.

## L_plot4seq.py

This library plots profiles on top of sequences.
Using combination of R and python.
It can generate visual representation of sequences by itself.

## L_seq_table_tools.py
Deprecated

## L_shade_aln.py
This is a library that makes good images of shaded alignments 
through TeXShade.

## L_shade_hist_aln.py
Visualizing histone MSAs via TEXSHADE.
A very powerful example.

## L_taxonomy_tools.py

Makes use of ETE2 library to filter sequences using there taxonomy.
Functions: check_tax_id_clade, subsample_taxids, group_taxids, get_taxid_from_gbrec, gbrec_to_clickhtml
