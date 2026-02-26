import pandas as pd
from Bio import Entrez
import json


pmid_neighbours_old = pd.read_csv( 'pmid_filenames_all.tsv', sep='\t')
pmid_neighbours_new_strep = pd.read_csv('pmid_filenames_all_new_strep.tsv', sep='\t')

pmid_neighbours = pmid_neighbours_new_strep[~pmid_neighbours_new_strep['pmid'].isin(pmid_neighbours_old['pmid'])]

pmid_neighbours = pmid_neighbours.drop_duplicates(subset=['pmid'])
print("Number of pmids: ", len(pmid_neighbours))
print(len(set(pmid_neighbours['pmid'])))
# print(pmid_neighbours.head())

Entrez.email = "pb11@iastate.edu"
def get_pubmed_abstract(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
    abstract = (handle.read())
    cleaned_abstract = abstract.replace('\n', ' ')

    handle.close()
    return cleaned_abstract

file_abstract_dict = {}
for index, row in pmid_neighbours.iterrows():
    # if index == 5:
    #     break
    pmid = int(row['pmid'])
    file_name = row['file_name']
    print("pmid: ", pmid, " file_name: ", file_name)
    abstract = get_pubmed_abstract(pmid)
    # print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    # print(abstract)
    if file_name in file_abstract_dict.keys():
        file_abstract_dict[file_name].append({'pmid':pmid,'abstract':abstract})
    else:
        file_abstract_dict[file_name] = [{'pmid':pmid,'abstract':abstract}]        

# json_output_path = source_folder + "pmid_abstracts.json"

json_output_path = "pmid_abstracts_new_strep.json"

with open(json_output_path, 'w') as json_file:
    json.dump(file_abstract_dict, json_file)
    

    
# I will give you a json file, can you extract the following information in a tabular form with columns "Filename", "Species", "Isolation host", "Where and how it was isolated".   The "Filename" is the key startting with "GCA". The other columns are extracted from the abstract as follows 1. "Species" is the species mentioned in the abstract, 2. "Isolation host" is the host it was isolated from as mentioned in the abstract 3. "Where and how it was isolated" is where and how it was isolated as mentioned in the abstract.
# some keys starting with "GCA" might have more than one abstract, please extract both in separate rows wth same "Filename"



