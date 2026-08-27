# StrepGenomes

### Genome files 

All the genome files used can be found here: [strep genomes](https://drive.google.com/drive/folders/1F-sI_ro7On3kW2aRG4iMaWyJcQLgv1it?usp=drive_link)

## finding Ras-RiPP/RGG neighbors

**script: find_rgg_neighbors.py**

This script identifies neighboring genes in annotated bacterial genomes based on Pfam HMM hits. It extracts proteomes from EMBL genome files ([strep genomes](https://drive.google.com/drive/folders/1F-sI_ro7On3kW2aRG4iMaWyJcQLgv1it?usp=drive_link)), performs HMM searches using HMMER, and detects gene pairs that occur within a defined genomic distance. By default, it will find neighbouring genes between RaS-RiPPs (PF04055) and Rggs (PF21259). In addition the script finds single hit results for Rggs (PF21259).

Python dependencies:
Biopython

usage: 

python find_rgg_neighbors.py \
  --embldir path/to/embl_files \
  --pfam1 PF04055 \
  --pfam2 PF21259 \
  --outfolder output_directory

| Argument      | Description                                  | Default                      |
| ------------- | -------------------------------------------- | ---------------------------- |
| `--embldir`   | Directory containing `.embl.gz` genome files | current path         |
| `--pfam1`     | First Pfam ID                                | `PF04055`                    |
| `--pfam2`     | Second Pfam ID                               | `PF21259`                    |
| `--outfolder` | Directory for output files                   | `output_ras_ripps_neighbors` |

### Output:
- {pfam1}_{pfam2}_neighbors.txt
Human-readable neighbor relationships

- {pfam1}_{pfam2}_neighbors.json
Structured JSON of neighbor results

- {pfam2}_hits.json
All single-gene hits for Pfam2 (Rggs (PF21259)) across genomes

## finding pubmedids for missing ids

**script: findpubmedid.py**

This script parses compressed EMBL genome files, extracts or retrieves verified PubMed IDs for each sequence record, and classifies sequences into chromosome (longest record) and plasmids based on description keywords. It outputs a JSON file containing accession details, descriptions, associated PMIDs, and notes on how each ID was obtained, with logging for errors and processing steps.

## find host information from embl files

**script: get_host_info.py**

This script parses compressed EMBL files to extract source metadata for each sequence record, including host, isolation source, sequence length, and whether the record represents a plasmid. It compiles this information per file and record, then outputs the results both as a structured JSON file and as an Excel spreadsheet for easier downstream analysis and review.

## extract species name

**script: species_stat.py**

This script reads a JSON file containing chromosome and plasmid PubMed ID mappings, filters to retain only chromosome records, and extracts the species name from each record’s description. It then creates a simplified mapping of filename to associated PMIDs and species, exporting the results as a new JSON file for downstream species-level analysis.

## extract neighbourhood information

**script: gene_count_and_analysis.py**

This script integrates multiple precomputed datasets to summarize gene neighborhood and RGG gene information per genome. For each genome key, it counts the number of neighboring genes and identified RGG hits, retrieves species information from source metadata (chromosome only), flags newly added entries, and compiles the results into a tab-delimited file for downstream comparative analysis.

## combine all data with chatGPT derived data

**script: species_host_stats_overall-get_all_info.py**

This script consolidates genomic, publication, host, gene neighborhood, and pathogenicity data into a single integrated dataset. For each genome file, it combines species and PMID information, chromosome-level host metadata, neighbor and RGG gene counts, curated pathogenicity data, and ChatGPT-derived annotations (e.g., host-site and pathogenicity ratings). It then compiles all matched records—including genomes without PMIDs—into a final Excel file for comprehensive downstream analysis and interpretation.
