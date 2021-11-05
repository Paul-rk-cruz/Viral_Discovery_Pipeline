import pandas as pd
import numpy as np
import Bio
from Bio import Entrez
import os.path
import argparse
from argparse import ArgumentParser

'''
input an output from rapsearch2 or diamond alignment and produce excel sheet of top 100 hits
'''

Entrez.email = 'katej16@uw.edu'
paramEutils = { 'usehistory':'Y' }


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')                                                                                       
    parser.add_argument('nr_alignment_output', help='nr alignment file from rapsearch (.m8) or diamond (.tsv)')                            
    args = parser.parse_args()                                                                                                             

    protein_output_file = args.nr_alignment_output

def get_lineage(protein_output):
    rapsearch_file = pd.read_csv(protein_output, sep='\t')
    rapsearch_file.columns = ['1','2','3','4','5','6','7','8','9','10','11','12']
    
    #entrez direct is computationally expensive. Downsample data 
    drop_duplicates = rapsearch_file.groupby('1').head(2)
    subsampled = drop_duplicates.head(100)
    subsampled.columns = ['1','2','3','4','5','6','7','8','9','10','11','12']

    #take proteins from subsampled df
    protein_v_lst = subsampled['2'].to_list()
    prot_lst2 = []
    locus_lst = []
    protein_lst = []
    taxonomy_lst = []

    #get genbankformatted efetch of query
    #if query accession is invalid, ignore that shit
    for i in protein_v_lst:
        try:
            handle = Entrez.efetch(db="protein", id=i, retmode="xml")
            records = Entrez.parse(handle)
            #parse handle object for lineage, Accession, and protein description 
            for record in records:
                taxonomy = record["GBSeq_taxonomy"]
                locus = record["GBSeq_locus"]
                protein = record["GBSeq_definition"]
                locus_lst.append(locus)
                protein_lst.append(protein)
                taxonomy_lst.append(taxonomy)
        except:
            pass

        #make df to export as csv and merge with other csvs if needed
        df = pd.DataFrame(zip(locus_lst,taxonomy_lst))
        df.columns = ['Protein Accession','Lineage']
        df['Description'] = protein_lst
    return df

def find_viruses_3(protein_output):
    lineagedf = get_lineage(protein_output)
    lineagedf.columns = ['Protein Accession','Lineage','Description']
    diamond_df = pd.read_csv(protein_output, sep='\t')
    diamond_df.columns = ['Query', 'Subject', '% Identical Match', 'Alignment Len', '# Mismatch', '# Gap Openings', 'Start in Query', 'End in Query', 'Start in Subject', 'End in Subject','evalue','bitscore']

    #need to re-drop duplicates and subsample upon merge
    diamond_df = diamond_df.groupby('Query').head(2)
    diamond_df = diamond_df.head(100)
    acc_col = diamond_df['Subject'].tolist()
    new_acc = []
    for acc in acc_col:
        #get rid of protein version
        acc_new = acc[:-2]
        new_acc.append(acc_new)
    diamond_df['Subject'] = new_acc

    #merge the dataframes
    merged_df = diamond_df.merge(lineagedf, left_on = 'Subject', right_on = 'Protein Accession')
    #pick what you want from dataframes
    select_merged_df = merged_df[['Query','Subject', 'evalue', '% Identical Match', 'Alignment Len', 'Lineage','Description']]
    #get node len for sorting
    query_list = select_merged_df['Query'].to_list()
    node_length = []
    for query in query_list:
        split_query = query.split("_")
        lengths = split_query[1:2]
        for value in lengths:
            ints = int(value)
            node_length.append(ints)
    select_merged_df.insert(loc=0, column = 'Node Number', value=node_length, allow_duplicates=False)
    #sort according to node len and evalue and take first match (there will be multiple hits of the same protein)
    select_sorted_df = select_merged_df.sort_values(["Node Number", "evalue"], ascending = (True, True))
    final_df = select_sorted_df.drop_duplicates()

    ##for python script return .xlsx
    #name file with file input

    name = protein_output
    filename = "%s.xlsx" % name
    excel = final_df.to_excel(filename)
    return excel

find_viruses_3(protein_output_file)
