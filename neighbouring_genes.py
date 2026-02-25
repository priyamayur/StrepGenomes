import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime
import pandas as pd

# directory = "/work/idoerg/britta_strep/embl"
# output_directory = "/work/idoerg/britta_strep/output/"
directory = "embl_short"
output_directory = "output/"
output_filename = "neighbours_and_pmids_updated_with_more_annotated.json"
Entrez.email = "pb11@iastate.edu"



# with open(output_directory +'all_neighbors_updated.json') as f:
#     all_neighbours = json.load(f)

with open(output_directory +'all_neighbors_updated_new_strep.json') as f:
    all_neighbours = json.load(f)

    # print(all_neighbours.keys())
# print("all neighbours length: ", len(all_neighbours.keys()))

# with open(output_directory +'pubmedIds_updated2_chromosomes_plasmids_corrected.json') as f:
#     pubmedIds = json.load(f)
    
with open(output_directory +'pubmedIds_chromosomes_plasmids_plus_new.json') as f:
    pubmedIds = json.load(f)

total_count_overlap = 0
count_overlap_neighbours = 0
neighbours_pmid = {}
britta_txt = []
txt_only_neighbours = []
only_neighbours = 0
all_pmids= []

for recordname in pubmedIds.keys():
    
    recordname_clean = recordname.split("_")[0] + "_" + recordname.split("_")[1] 
    record_type = recordname.split("_")[-1]
    
    if record_type == 'p':
        continue
    
    for pmid_this in pubmedIds[recordname]['pmids']:
        all_pmids.append([recordname_clean, pmid_this])
    print("Processing recordname: ", recordname_clean)
    if recordname_clean in all_neighbours.keys():
        total_count_overlap += 1
        # print("recordname_clean before: ", recordname_clean)
        if all_neighbours[recordname_clean] == None:
            continue
        neighbours_count = len(all_neighbours[recordname_clean])
        pmids = (pubmedIds[recordname]['pmids'])
        if neighbours_count > 0:
            only_neighbours += 1
            for pmid in pmids:
                txt_only_neighbours.append([recordname_clean, pmid])
            if len(pmids) == 0:    
                txt_only_neighbours.append([recordname_clean, "no pmid"])
        if neighbours_count > 0 and len(pmids) > 0:
            if len(pmids) > 1:
                print("recordname for more than one pmid: ", recordname_clean)
            # print("recordname_clean: ", recordname_clean)
            count_overlap_neighbours += 1
            neighbours_pmid[recordname_clean] = pubmedIds[recordname] 
            neighbours_pmid[recordname_clean]['neighbours'] = all_neighbours[recordname_clean]
            for pmid in pmids:
                britta_txt.append([recordname_clean, pmid])


britta_txt_df = pd.DataFrame(britta_txt, columns=['file_name', 'pmid'])
txt_only_neighbours_df = pd.DataFrame(txt_only_neighbours, columns=['file_name', 'pmid'])
all_pmids_df = pd.DataFrame(all_pmids, columns=['file_name', 'pmid'])
all_pmids_df = all_pmids_df.drop_duplicates()
# all_pmids_df.to_csv(output_directory + "all_pmids_new_strep.tsv", sep='\t', index=False)
txt_only_neighbours_df.to_csv(output_directory + "only_neighbours_new_strep.tsv", sep='\t', index=False)
# britta_txt_df.to_csv(output_directory + "pmid_neighbours_updated_new_strep.tsv", sep='\t',index=False)
# with open(output_directory + output_filename , "w") as json_file:
#     json.dump(neighbours_pmid, json_file, indent=4)
    
print("total_count_overlap: ", total_count_overlap)
print("count_overlap_neighbours: ", count_overlap_neighbours)
print("only_neighbours: ", only_neighbours)