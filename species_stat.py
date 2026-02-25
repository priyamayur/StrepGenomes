import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime
import pandas as pd

output_directory = "output/"
# with open(output_directory +'pubmedIds_updated2_chromosomes_plasmids_corrected.json') as f:
#     pubmedIds = json.load(f)
    
with open(output_directory +'pubmedIds_chromosomes_plasmids_plus_new.json') as f:
    pubmedIds = json.load(f)

species_to_pmid = {}
filename_pmids = []
filename_pmids_all = {}

for recordname in pubmedIds.keys():
  
  if len(recordname.split("_")) <= 4:
    
    recordname_clean = recordname.split("_")[0] + "_" + recordname.split("_")[1] 
  else:
    print("Processing record:", recordname)
    recordname_clean = recordname.split("_")[0] + "_" + recordname.split("_")[1] + "_" + recordname.split("_")[2]  
  record_type = recordname.split("_")[-1]
  if record_type == 'p':
      continue
  pmids = (pubmedIds[recordname]['pmids'])
  description = pubmedIds[recordname]['description']
  # print(recordname_clean)
  # print(description)
  if (description) == None:
    print("No description found for record: ", recordname_clean)
    continue
  species_name = description.split(" ")[0] + " " + description.split(" ")[1]
  filename_pmids_all[recordname_clean] = {'pmids':pmids,'species':species_name}
  


# with open(output_directory + 'filename_pmid_species.json' , "w") as json_file:
#     json.dump(filename_pmids_all, json_file, indent=4)

with open(output_directory + 'filename_pmid_species_new_strep.json' , "w") as json_file:
    json.dump(filename_pmids_all, json_file, indent=4)