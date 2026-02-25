import json
import re
from Bio import Entrez
import os
import gzip
from Bio import SeqIO
import logging
from datetime import datetime

# directory = "/work/idoerg/britta_strep/embl"
# output_directory = "/work/idoerg/britta_strep/output/"
directory = "embl_short"
# directory = "combined_bakta_original_annotations_plus_new"
output_directory = "output/"
output_filename = "pubmedIds_short_sept_10.json"
output_filename = "pubmedIds_chromosomes_plasmids_plus_new.json"
Entrez.email = "pb11@iastate.edu"

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log_filename = f"log_{timestamp}.log"
logging.basicConfig(filename=output_directory+log_filename, level=logging.INFO, format="%(asctime)s - %(message)s")

def jaccard_similarity(str1, str2):
    """
    Compute Jaccard similarity between two strings.
    """
    # Convert strings to sets of words or characters
    set1 = set(str1.lower().split())  # Tokenize by words
    set2 = set(str2.lower().split())  # Tokenize by words

    # Compute intersection and union
    intersection = set1.intersection(set2)
    union = set1.union(set2)

    # Compute Jaccard similarity
    return len(intersection) / len(union) if union else 0.0

def clean_list(lst):
    """
    Remove None, empty strings, and whitespace strings from a list.
    """
    clean_lst =  [item for item in lst if item not in [None, "", " "]]
    clean_lst_unique = list(set(clean_lst))
    
    return clean_lst_unique

def remove_special_chars(text):
    """
    Remove special characters from a string.
    """
    return re.sub(r'[^a-zA-Z0-9\s]', '', text)  # Keeps letters, numbers, and spaces

def verify_pubmed_id(pmid, title, first_author):
  """
  Verify that the PubMed ID (PMID) corresponds to the correct article by comparing the title and first author.

  """
  jaccard_similarity_threshold = 0.7
  same_article_name = False
  same_first_author = False
  handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
  records = Entrez.read(handle)
  handle.close()

  ##check
  if ("PubmedArticle" not in records.keys()) :
     return same_first_author and same_article_name
  if (len(records["PubmedArticle"]) == 0) :
     return same_first_author and same_article_name
  if "MedlineCitation" not in records["PubmedArticle"][0].keys():
     return same_first_author and same_article_name
  if "Article" not in records["PubmedArticle"][0]["MedlineCitation"].keys():
     return same_first_author and same_article_name
  if "ArticleTitle" not in records["PubmedArticle"][0]["MedlineCitation"]["Article"].keys():
     return same_first_author and same_article_name
  if "AuthorList" not in records["PubmedArticle"][0]["MedlineCitation"]["Article"].keys():
     return same_first_author and same_article_name
  if len(records["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"]) == 0:
     return same_first_author and same_article_name

  # verify the article title
  found_article_title = records["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"]
  clean_found_article_title = remove_special_chars(found_article_title)
  similarity = jaccard_similarity(clean_found_article_title, title)
  if similarity > jaccard_similarity_threshold:
    same_article_name =  True

  # verify the first author
  found_first_authors = records["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"][0]
  if "LastName" in found_first_authors.keys():
      found_first_author_lastname = found_first_authors['LastName']
      first_author_lastname = first_author.split(" ")[0]
      if found_first_author_lastname == first_author_lastname:
            same_first_author = True

  return same_first_author and same_article_name

def find_pubmed_id(title, first_author):
  """
  Find the PubMed ID (PMID) of a paper by searching PubMed with different search criterias.
  """
  
  clean_title = remove_special_chars(title)
  # Search PubMed by title and first author
  handle = Entrez.esearch(db="pubmed", term= clean_title+"[title]"+" AND "+ first_author + "[Author - First]", retmax=40)
  record = Entrez.read(handle)
  handle.close()

  # Get PubMed ID (PMID) of the matching papers
  pmid_list = record["IdList"]

  if len(pmid_list) > 0:
    return pmid_list[0]
  else:
    # Search PubMed by first n words in the title
    first_n_words_list = list(range(5, len(clean_title.split(" ")), 5))
    max_size_pmid_list = []
    for n_word in first_n_words_list:
      first_n_words = n_word
      if len(clean_title.split(" ")) > first_n_words:
        shorter_title_list = title.split(" ")[0:first_n_words]
        shorter_title = " ".join(shorter_title_list)
        handle = Entrez.esearch(db="pubmed", term= shorter_title+"[title]", retmax=40)
        record = Entrez.read(handle)
        handle.close()
        pmid_list = record["IdList"]
        if len(pmid_list) > 0:
          max_size_pmid_list.extend(pmid_list)

    pmid_list_final = list(set(max_size_pmid_list))

  ## verify the pmids
  for pmid in pmid_list_final:
    is_pubmed_id = verify_pubmed_id(pmid, title, first_author)
    if is_pubmed_id:
      return pmid

def get_pubmedid_from_database(reference):
  """
  Get the PubMed ID (PMID) of a paper from the PubMed database.
  """
  pubmed_id = None
  notes = ""

  if len(reference.title) > 2: # title is not empty
    title = reference.title
    if len(reference.authors) > 2: # authors is not empty
      first_author = reference.authors.split(",")[0]
      pubmed_id = find_pubmed_id(title, first_author)
      if pubmed_id != None:
        notes = pubmed_id+": found in pubmed database;"  
        print("FOUND FROM PUBMED DB")
        logging.info("FOUND FROM PUBMED DB")
      else:
        notes = title+": not found in pubmed database;"  

  return pubmed_id, notes
   
def get_complete_pubmedid_info(record):
  """
  Get the complete information of the PubMed ID (PMID).
  """
  pubmed_ids = []  
  ## get first accession and description    
  accession = record.annotations.get("accessions", ["Unknown"])[0]
  description = record.description
  record_id = record.id
  record_seq_length = len(record.seq)

  ## find the pubmed ids
  notes = ''
  for ref in record.annotations.get("references", []):
      pubmed_id = ''
      ## pubmed id already found in the file
      if len(ref.pubmed_id) > 0:
        pubmed_id =str(ref.pubmed_id)
        notes +=  pubmed_id+": found in file;"
      ## search for pubmed id from the database by the title  
      else:
        pubmed_id, notes_db = get_pubmedid_from_database(ref)
        notes += notes_db     
      pubmed_ids.append(pubmed_id)
  pubmed_ids_clean = clean_list(pubmed_ids)
  
  if notes == '':
      notes="no paper"

  complete_info = {'pmids': pubmed_ids_clean,'accesion':accession, 'description':description, 'notes':notes }  
  extra_info = {'record_id':record_id, 'record_seq_length':record_seq_length}  

  return complete_info, extra_info

def get_chromosome_and_plasmids(complete_info_all_records, complete_info_metadata, filename_clean):
  """
  Get the chromosome and plasmids' information separated from all the records in the embl file.
  """
  chromosomes_and_plasmids = {}

  # find the longest sequence as chromosome
  max_chromosome_length = 0
  chromosome_record_id = list(complete_info_metadata.keys())[0]
  for record_id in complete_info_metadata.keys():
    record_seq_length = complete_info_metadata[record_id]['record_seq_length']
    if record_seq_length > max_chromosome_length:
      max_chromosome_length = record_seq_length
      chromosome_record_id = record_id
  
  chromosome_file_name = filename_clean + "_" + chromosome_record_id + "_c"
  chromosomes_and_plasmids[chromosome_file_name] = complete_info_all_records[chromosome_record_id]

  # everything else checked for a plasmid
  for record_id in complete_info_metadata.keys():
    if record_id != chromosome_record_id:
      record_description = complete_info_metadata[record_id]['record_description']
      if "plasmid" in record_description.lower() or "extrachromosomal" in record_description.lower():
        plasmid_file_name = filename_clean + "_" + record_id + "_p"
        chromosomes_and_plasmids[plasmid_file_name] = complete_info_all_records[record_id]
      else:
         print("cannot classify into chromosome or plasmid: ", filename_clean, record_id, record_description)  
         logging.info("cannot classify into chromosome or plasmid: ", filename_clean, record_id, record_description)

  return chromosomes_and_plasmids


results ={}
for filename in os.listdir(directory):
    file_path = os.path.join(directory, filename)
    if os.path.isfile(file_path):
      print(filename)
      filename_clean = filename.split(".")[0]
      try :
        with gzip.open(file_path, "rt") as handle: # read as text mode
          
          complete_info_metadata = {}
          complete_info_all_records = {}
          for record in list(SeqIO.parse(handle, "embl")):
            complete_info, extra_info = get_complete_pubmedid_info(record)
            record_id = extra_info['record_id']
            record_seq_length = extra_info['record_seq_length']
            record_description = complete_info['description']
            complete_info_metadata[record_id] = {'record_seq_length':record_seq_length, 'record_description':record_description}
            complete_info_all_records[record_id] = complete_info
            
          chromosomes_and_plasmids = get_chromosome_and_plasmids(complete_info_all_records, complete_info_metadata, filename_clean)
          results |= chromosomes_and_plasmids
      except Exception as e:
        pubmed_ids_clean = []
        notes = "error"
        accession = None
        description = None
        print("exception: ",filename)
        print(f"An error occurred: {e}")
        logging.info("exception: ",filename)
        logging.info(f"An error occurred: {e}")
       
        results[filename_clean] = {'pmids': pubmed_ids_clean,'accesion':accession, 'description':description, 'notes':notes }


with open(output_directory + output_filename , "w") as json_file:
    json.dump(results, json_file, indent=4)

print("done")

