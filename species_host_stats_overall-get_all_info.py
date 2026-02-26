import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime
import pandas as pd


pmid_neighbours_updated = pd.read_csv("only_neighbours_new_strep.tsv", sep='\t')
chatgpt_file = pd.read_excel( "abstract_improved_merged_new_strep_sept_12.xlsx")
previous_chatgpt_file = pd.read_excel("final_chatgpt_aug_combined_results_21.xlsx")

with open('filename_pmid_species_new_strep.json') as f:
    species = json.load(f)

host_information = pd.read_excel( "host_information_new_strep.xlsx")
gene_count_info = pd.read_csv("gene_count_information_including_new_strep.tsv", sep='\t')

pathogenecity=pd.read_excel("Streptococcus_Pathogenicity.xlsx")
pathogenecity_dict = pathogenecity.set_index("Species")["Pathogenicity"].to_dict()

host_information_chromosome = host_information[host_information["plasmid"] == False]
print("chromosome records=", len(host_information_chromosome))
file_species_host = []
species_pathogenicity = "Unknown"
for file_name in species.keys():
    print("Processing file:", file_name)
    new_file = False
    if len(file_name.split("_")) >2 :
        new_file = True
    found_pmid = False
    species_name = species[file_name]['species']
    pmids = species[file_name]['pmids']

    host_info_row = host_information_chromosome[(host_information_chromosome['filename'] == file_name)]
    host = host_info_row['host'].values[0] if not host_info_row.empty else "Unknown"
    isolation_source = host_info_row['isolation_source'].values[0] if not host_info_row.empty else "Unknown"

    neighbours_found = pmid_neighbours_updated[pmid_neighbours_updated['file_name'] == file_name]
    
    
    gene_count_info_neighbour_count = gene_count_info[(gene_count_info['file_name'] == file_name)]['neighbour_count'].values[0]
    gene_count_info_rgg_count = gene_count_info[(gene_count_info['file_name'] == file_name)]['rgg_gene_count'].values[0]
    
    
    for path_sp in pathogenecity_dict:
        if path_sp in species_name or species_name in path_sp:
            species_pathogenicity = pathogenecity_dict[path_sp]
            break
    if species_pathogenicity == "Unknown":    
        print("No pathogenicity info for species:", species_name)    
        

    if len(neighbours_found) > 0:
        neghbours = 1
        # continue 
    else:
        neghbours = 0 

    for pmid in pmids:
        found_pmid = True
        # print("Processing file:", file_name, "with pmid:", pmid)
        previous_chatgpt_file_row = previous_chatgpt_file[(previous_chatgpt_file['pmid'] == pmid)]
        # print("previous_chatgpt_file_row: ", len(previous_chatgpt_file_row), pmid, file_name)
        chatgpt_file_row = chatgpt_file[(chatgpt_file['pmid'] == int(pmid))]
        # chatgpt_host = chatgpt_file_row['host'].values[0] if not chatgpt_file_row.empty else "Unknown"
        chatgpt_host_host_site = chatgpt_file_row['host and host-site'].values[0] if not chatgpt_file_row.empty else "Unknown"
        chatgpt_pathogenicity_rating = chatgpt_file_row['Pathogenicity Rating'].values[0] if not chatgpt_file_row.empty else "Unknown"
        host_normalized = previous_chatgpt_file_row['Host normalized'].values[0] if not previous_chatgpt_file_row.empty else "unknown_pmid"

        file_species_host.append([file_name, pmid, species_name, neghbours,gene_count_info_neighbour_count, gene_count_info_rgg_count,
                                host, isolation_source,
                                chatgpt_host_host_site, host_normalized,
                                species_pathogenicity, chatgpt_pathogenicity_rating])
        # print("file_species_host: ", file_species_host)
        
    if not found_pmid:
        pmid = "no pmid"
        chatgpt_isolation_process = "Unknown"
        chatgpt_location = "Unknown"
        chatgpt_host_site = "Unknown"
        chatgpt_non_host_source = "Unknown"
        chatgpt_disease = "Unknown"
        chatgpt_pathogenicity_rating = "Unknown"
        chatgpt_host_host_site = "Unknown"
        previous_chatgpt_file_row = previous_chatgpt_file[(previous_chatgpt_file['file_name'] == file_name)]

        host_normalized = previous_chatgpt_file_row['Host normalized'].values[0] if not previous_chatgpt_file_row.empty else "unknown_filename"
        file_species_host.append([file_name, pmid, species_name, neghbours,gene_count_info_neighbour_count, gene_count_info_rgg_count,
                                host, isolation_source, 
                                chatgpt_host_host_site,host_normalized,
                                species_pathogenicity, chatgpt_pathogenicity_rating
                                ])
    

file_species_host_df = pd.DataFrame(file_species_host, 
                                    columns=['file_name', 'pmid', 'species', 'neighbours_found', 'gene_neighbour_count', 'rgg_gene_count',
                                             'host', 'isolation_source', 
                                              'chatgpt_host_host_site', 'host_normalized',
                                             'pathogenecity','chatgpt_pathogenicity_rating'])
print("length final==>",len(file_species_host_df))
print("file_species_host_df shape: ", file_species_host_df.head())


file_species_host_df.to_excel('chatgpt_plus_original_sep_12_new_strep.xlsx', index=False)


