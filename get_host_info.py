import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime
import pandas as pd


# out_direc = "c:\\Users\\pb11\\Documents\\Projects\\ras_ripp\\combined_bakta_original_annotations\\"
out_direc = "c:\\Users\\pb11\\Documents\\Projects\\ras_ripp\\combined_bakta_original_annotations_plus_new\\"
file_name_to_source = {}
file_name_to_host=[]
count = 0
for filename in os.listdir(out_direc):
    count += 1
    file_path = os.path.join(out_direc, filename)
    # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    
    filename_clean = filename.split(".")[0]
    print("file name =", filename_clean, "count=",count)
    # file_path = out_direc + file + ".embl.gz"
    try:
        with gzip.open(file_path, "rt") as handle:
            new_records = []
            records = list(SeqIO.parse(handle, "embl"))  # Use SeqIO.parse() for multiple records
            # if len(list(records)) == 0:
            #     print("No records found in", file_path)
            # if len(records) > 1:
            file_name_to_source[filename_clean] = []
            
            for record in records:
                plasmid = False
                print("ID:", record.id)
                host =''
                isolation_source=''
                for feature in record.features:
                    if feature.type == "source":
                        # print("---- SOURCE FEATURE ----")
                        
                        for key, value in feature.qualifiers.items():
                            if "host" in key.lower():
                                host = value[0] if value else ''
                            if "isolation_source" in key.lower():
                                isolation_source = value[0] if value else ''    
                            if "plasmid" in key.lower():
                                plasmid =True
                           # print(f"{key}: {value[0]}")
                file_name_to_source[filename_clean].append({
                    "filename": filename_clean,
                    "id": record.id,
                    "description": record.description,
                    "length": len(record.seq),
                    "host": host,
                    "isolation_source": isolation_source,
                    "plasmid": plasmid
                })    
                file_name_to_host.append(
                   [filename_clean,
                    record.id,
                    record.description,
                    len(record.seq),
                    host,
                    isolation_source,
                    plasmid]
                )    




    except Exception as e:
        print("Error processing file:", filename, "Error:", e)
        logging.error(f"Error processing file {filename}: {e}")

print("all files=",len(file_name_to_source))        

with open("c:\\Users\\pb11\\Documents\\Projects\\ras_ripp\\output\\" +  "host_information_new_strep.json", "w") as json_file:
    json.dump(file_name_to_source, json_file, indent=4)

file_name_to_source_df = pd.DataFrame(file_name_to_host, columns=["filename", "id", "description", "length", "host", "isolation_source", "plasmid"])
file_name_to_source_df.to_excel("c:\\Users\\pb11\\Documents\\Projects\\ras_ripp\\output\\" +  "host_information_new_strep.xlsx", index=False)