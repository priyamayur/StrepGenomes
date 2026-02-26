import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime
import pandas as pd


directory = "embl_short"
output_directory = "output/"
Entrez.email = "pb11@iastate.edu"


with open(output_directory +'all_neighbors_updated_new_strep.json') as f:
    all_neighbours = json.load(f)

print("all neighbours length: ", len(all_neighbours.keys()))

with open(output_directory +'RGG_gene_found_dict_new_strep.json') as f:
    rgg_genes_only = json.load(f)

with open(output_directory +'source_information_new_strep.json') as f:
    source = json.load(f)

files_gene_information =[]

for gca_key in all_neighbours.keys():
    count_neighbour = 0
    count_rgg = 0
    species= "Unknown"
    new = False
    if gca_key.count("_") > 1:
        new = True
    # print("------------------------------------------------------------------")
    if all_neighbours[gca_key] != None:
        # print("GCA key: ", gca_key, " neighbours: ", all_neighbours[gca_key] )
        # print("count of neighbours: ", len(all_neighbours[gca_key]))
        count_neighbour = len(all_neighbours[gca_key])

    if rgg_genes_only[gca_key] != None and len(rgg_genes_only[gca_key])> 0:
        # print("GCA key with RGG: ", gca_key, " RGG genes: ", rgg_genes_only[gca_key] )
        # print("count of RGG genes: ", len(rgg_genes_only[gca_key]['hit_ids']))
        count_rgg = len(rgg_genes_only[gca_key]['hit_ids'])
    if gca_key in source.keys():
        for chro in source[gca_key]:
            plasmid = chro['plasmid']
            if plasmid == False:
                species = chro['description']  
    print("File: ", gca_key, " Species: ", species, " Neighbour count: ", count_neighbour, " RGG gene count: ", count_rgg)  
    files_gene_information.append([gca_key, new, species, count_neighbour, count_rgg])
    
files_gene_information_df = pd.DataFrame(files_gene_information, columns=['file_name','new addition', 'species','neighbour_count', 'rgg_gene_count'])
files_gene_information_df.to_csv(output_directory + "gene_count_information_including_new_strep.tsv", sep='\t', index=False)
 
        
