#!/usr/bin/env python
import gzip
import re
import json
import sys
import logging
import subprocess
import glob
from pathlib import Path
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


BASE_PATH = Path('/mnt/c/Users/pb11/Documents/Projects/ras_ripp')

def neighboring_genes(json_list):
    plus_list = []
    minus_list = []
    all_list = []
    # + strand dictionaries
    locus_to_num_plus = {} # Key: locus tag; value: index of gene on strand
    num_to_locus_plus = {} # Key: index of gene on strand; value: locus tag 

    # - strand dictionaries
    locus_to_num_minus = {}
    num_to_locus_minus = {}

    # both strand dictionaries
    locus_to_num_all = {}
    num_to_locus_all = {}

    # build and sort neighbor vectors
    for gene in json_list:
        if gene['strand'] == 1:
            plus_list.append((gene['start'], gene['end'], gene['locus_tag'], gene['strand']))
        else:
            minus_list.append((gene['start'], gene['end'], gene['locus_tag'], gene['strand']))
        all_list.append((gene['start'], gene['end'], gene['locus_tag'], gene['strand']))
    plus_list.sort()
    minus_list.sort()
    all_list.sort()

    for (i,plusdata) in enumerate(plus_list):
        locus_tag = plusdata[2]
        num_to_locus_plus[i] = locus_tag
        locus_to_num_plus[locus_tag] = i

    for (i,minusdata) in enumerate(minus_list):
        locus_tag = minusdata[2]
        num_to_locus_minus[i] = locus_tag
        locus_to_num_minus[locus_tag] = i

    for (i,alldata) in enumerate(all_list):
        locus_tag = alldata[2]
        num_to_locus_all[i] = locus_tag
        locus_to_num_all[locus_tag] = i
    return(locus_to_num_plus, locus_to_num_minus, locus_to_num_all)

def find_chromosome(dna_molecules):
    """
    Find the largest chromosome in a list of DNA molecules
    """
    maxlen = (0, None)
    for i, dna in enumerate(dna_molecules):
        if len(dna.seq) > maxlen[0]:
            maxlen = (len(dna.seq), i)
    return(dna_molecules[maxlen[1]])
             
    
def create_proteome_files(infile, proteome_path, infmt="embl", outfmt="fasta"):
    """
    Create two proteome files,from the genomic EMBL file. FASTA and JSON
    """
    dna_molecules = []
    outfile_stem = proteome_path /  f"{Path(infile).with_suffix('').stem}_aa"
    json_outpath = f'{outfile_stem}.json'
    seqfile_outpath = f'{outfile_stem}.{outfmt}'
    # Read embl file. May contain more than one DNA molecule for its organism
    # E.g. chromosome and more than one plasmid
    try:
        with gzip.open(infile, 'rt') as f:
            for dna_mol in SeqIO.parse(f,infmt):
                dna_molecules.append(dna_mol)
    except:
        logging.error(f'{infile} Corrupt file?')
        return(0)
    # If there is more than one DNA molecule in the EMBL
    # file, find the longest one (the chromosome)
    if len(dna_molecules) == 0:
        logging.error(f'No molecules found in {infile}. Corrupt file?')
        return(0)
    elif len(dna_molecules) > 1:
        genome = find_chromosome(dna_molecules)
        logging.info(f'More than one DNA mol in file. Chromosome ID={genome.id}')
    elif len(dna_molecules) == 1:
        genome = dna_molecules[0]
        logging.info(f'One DNA mol in file. Chromosome ID={genome.id}')
    json_list = []
    outseq_list = []
    if len(genome.features) < 5:
        logging.error(f'{infile} probably unannotated. Stopping.')
        return(0)
    logging.info(f'infile: {infile}')
#    logging.info(genome.features[2])
    for feature in genome.features:
        if feature.type == 'CDS':
            # don't do pseudogenes
            if ('pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers):
                logging.info(f'{infile} pseudogene skipped {feature.location.start+1}..{feature.location.end}')
                continue
            if 'translation' not in feature.qualifiers:
                logging.info(f'{infile} no translation {feature.location.start+1}..{feature.location.end}')
                continue
            protseq = feature.qualifiers['translation'][0]
            protid = feature.qualifiers['protein_id'][0]
            start = feature.location.start
            end = feature.location.end
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
            else:
                locus_tag = f'ltag_{start+1}_{end}'
                logging.info(f'{infile} no locus tag {start+1}..{end}')

            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0]
            else:
                gene = ''
            strand = feature.location.strand
            # note = feature.qualifiers['note']
            # create a fasta record
            outrec = SeqRecord(
                Seq(protseq),
                id=locus_tag,
                name=gene,
                description=f'({start+1}:{end}),{strand},{protid}'
                )
            outseq_list.append(outrec)

            # Create a dictionary for later json output           
            outrecdict = {'sequence':protseq,
                          'id': protid,
                          'locus_tag': locus_tag,
                          'name': gene,
                          # 'description': note,
                          'start': start+1,
                          'end': end,
                          'strand': strand}
            json_list.append(outrecdict)
    # write the FASTA proteome file 
    logging.info(f'writing {seqfile_outpath} with {len(outseq_list)} sequences')
    SeqIO.write(outseq_list, seqfile_outpath, outfmt)

    # map loci to their strand, coordinates, and relative position on their strands
    locus_to_num_plus, locus_to_num_minus, locus_to_num_all = \
        neighboring_genes(json_list)

    # Write proteome files with only + strand proteins, - strand proteins,
    # and both strands
    logging.info(f'writing {json_outpath} with {len(json_list)} sequences')
    with open(json_outpath,"w") as f:
        json.dump(json_list,f)
    l2n_plus_path = f'{outfile_stem}_l2n+.json' 
    l2n_minus_path = f'{outfile_stem}_l2n-.json' 
    l2n_all_path = f'{outfile_stem}_l2na.json' 

    with open(l2n_plus_path,"w") as f:
        json.dump(locus_to_num_plus,f)
    with open(l2n_minus_path,"w") as f:
        json.dump(locus_to_num_minus,f)
    with open(l2n_all_path,"w") as f:
        json.dump(locus_to_num_all,f)

def hmmsearch(hmm_path, proteome_path):
    """
    hmmsearch using the hmm in hmmpath (e.g. a Pfam file)
    for hits in the protein multifasta file proteome_path
    Creates a tab-formatted Hmmr search results file
    """
    hmmoutpath = BASE_PATH / 'hmms'
    hmmoutpath.mkdir(parents=True, exist_ok=True)
    proteome_stem = Path(proteome_path).with_suffix('').stem
    hmm_stem = Path(hmm_path).with_suffix('').stem
    subprocess.run(["hmmsearch",  "--tblout", f"{hmmoutpath}/{hmm_stem}_{proteome_stem}.hmmsrch", "-o","/dev/null", hmm_path, proteome_path])

def all_hmmsearch(hmm_path, proteome_dir):
    """
    Create hmmsrch files with a single Pfam query for all investigated proteomes.
    """
    for proteome_file in Path(proteome_dir).glob('*.fasta'):
        hmmsearch(hmm_path, proteome_file)

def get_hmmr_hits(hmmsrch_path,eval_threshold=0.001):
    """
    Get all the hits for one hmmsrch file, which is all the
    Pfam hits in that proteome for one pfam hmm
    """
    # create a regexp pattern to parse the Pfam ID and the genome sequence ID from the hmmrsrch filename
    # e.g. PF04055_GCA_965111745_aa.hmmsrch
    pattern = r"^(PF\d+)_(GCA_\d+)"
    hmmr_hits = {}
    hit_ids = []
    hmmsrch_stem = Path(hmmsrch_path).with_suffix('').stem
    
    match = re.match(pattern, hmmsrch_stem)
    # print(hmmsrch_stem, "====>",match)
    pfam_id = match.group(1)
    
    genome_id = match.group(2)
    for qresult in SearchIO.read(BASE_PATH / 'hmms' / hmmsrch_path, 'hmmer3-tab'):
        if qresult.evalue < eval_threshold:
            hit_ids.append(qresult.id)
    hmmr_hits['pfam'] = pfam_id
    hmmr_hits['genome_id'] = genome_id
    hmmr_hits['hit_ids'] = hit_ids  
    return(hmmr_hits)

def get_hits_positions(hmmr_hits):
    """
    Get the relative positions in the genome (gene number on strand, or gene number in both strands)
    from a hmmr hit list
    """
    proteome_path = BASE_PATH / 'proteome'
    json_file_plus = proteome_path / f'{hmmr_hits["genome_id"]}_aa_l2n+.json'
    json_file_minus = proteome_path / f'{hmmr_hits["genome_id"]}_aa_l2n-.json'
    json_file_all = proteome_path / f'{hmmr_hits["genome_id"]}_aa_l2na.json'
    hmmr_hits_all = {}
    hmmr_hits_all['hit_ids'] = {} 
    hmmr_hits_all['genome_id'] = hmmr_hits['genome_id']
    hmmr_hits_all['pfam'] = hmmr_hits['pfam']

    with open(json_file_all,'r') as f:
        all_strand_genes = json.load(f)
    for hit_id in hmmr_hits['hit_ids']:
        if hit_id in all_strand_genes:
            hit_pos = all_strand_genes[hit_id]
            hmmr_hits_all['hit_ids'][hit_id] =  hit_pos 
    return hmmr_hits_all

def get_hit_distances(hmmr_hits_1, hmmr_hits_2,max_distance = 4):
    """
    Do and all-vs-all search between the two hit lists to see if any hits are within 
    max_distance distance in the same proteome
    """
    neighbors = []
    found_pairs = {}
    for hit_id_1 in hmmr_hits_1['hit_ids']:
        for hit_id_2 in hmmr_hits_2['hit_ids']:
            if hit_id_1 == hit_id_2:
                print("what the fuck?")
                raise
            hit_pos_1 = hmmr_hits_1['hit_ids'][hit_id_1]
            hit_pos_2 = hmmr_hits_2['hit_ids'][hit_id_2]
            gene_distance = abs(hit_pos_1 - hit_pos_2)
            hit_pair = [hit_id_1, hit_id_2]
            hit_pair.sort()
            hit_pair = tuple(hit_pair)
            if hit_pair in found_pairs:
                continue
            else:
                found_pairs[hit_pair] = 1
                if gene_distance <= max_distance:
                    neighbors.append((hit_pair, gene_distance))
    return(neighbors)

def get_neighboring_hmm_hits(pfam1, pfam2, proteome):
    # proteome_path = BASE_PATH / 'proteome'
    hmmrsrch_path_1 = f'{pfam1}_{proteome}_aa.hmmsrch'
    hmmrsrch_path_2 = f'{pfam2}_{proteome}_aa.hmmsrch'
    try:
        hmmr_hits_1 = get_hmmr_hits(hmmrsrch_path_1)
    except ValueError:
        print (f'{hmmrsrch_path_1} empty')
        return ([])
    try:
        hmmr_hits_2 = get_hmmr_hits(hmmrsrch_path_2)
    except ValueError:
        print (f'{hmmrsrch_path_2} empty')
        return ([])
    hmmr_hits_pos_1 = get_hits_positions(hmmr_hits_1)
    hmmr_hits_pos_2 = get_hits_positions(hmmr_hits_2)
    neighbors = get_hit_distances(hmmr_hits_pos_1, hmmr_hits_pos_2)
    return(neighbors)

def get_all_neighbors(pfam1, pfam2,outfolder,embl_dir):
    outfile = f'{outfolder}/{pfam1}_{pfam2}_neighbors.txt'
    fout = open(outfile,"w")
    # proteome_dir = BASE_PATH / 'embl'
    proteome_files = glob.glob(f'{embl_dir}/*.embl.gz')
    neighbors_dict = {}
    for pfile in proteome_files:
        gca_id = f'{Path(pfile).stem.split(".")[0]}'
        try:
            neighbors = get_neighboring_hmm_hits(pfam1, pfam2, gca_id)
        except FileNotFoundError:
            print(f'problem with {gca_id}: not annotated?')
            neighbors_dict[gca_id] = None
            continue
        fout.write(f'{gca_id}\n')
        neighbors_dict[gca_id] = neighbors
        if not neighbors:
            fout.write('NN\n')
        else:
            for i in neighbors:
                fout.write(f'{i[0][0]}\t{i[0][1]}\t{i[1]}\n')
        fout.write('//\n')
    fout.close()
    z = 0
    n = 0
    for i in neighbors_dict:
        if neighbors_dict[i] is  None:
            n += 1
        elif len(neighbors_dict[i]) > 0:
            z += 1
    print(f'Pairs with neighbors: {z}, {n} not analyzed,  out of {len(neighbors_dict)}')

    with open(f'{outfolder}/{pfam1}_{pfam2}_neighbors.json','w') as f:
        json.dump(neighbors_dict,f)
        
    print(f'All neighbors written to {outfile} and {outfolder}/{pfam1}_{pfam2}_neighbors.json')
        

def create_search_files(hmm_path, fasta_paths):
    for fasta_path in fasta_paths: 
        hmmsearch(hmm_path, fasta_path)

def create_proteome_db(embldir, proteome_path):
    for embl_file in embldir.glob('*.embl.gz'):
        create_proteome_files(embl_file, proteome_path)


def get_single_gene_hit(pfam,outfile,embl_dir):
    proteome_files = glob.glob(f'{embl_dir}/*.embl.gz')
    
    gene_found_dict = {}
    for pfile in proteome_files:
        gca_id = f'{Path(pfile).stem.split(".")[0]}'
        try:
            hmmrsrch_path_1 = f'{pfam}_{gca_id}_aa.hmmsrch'
            # print(hmmrsrch_path_1)
            try:
                hmmr_hits_1 = get_hmmr_hits(hmmrsrch_path_1)
            except ValueError:
                print (f'{hmmrsrch_path_1} empty')
                hmmr_hits_1 = ([])
            # neighbors = get_hmmr_hits(pfam1, pfam2, gca_id)
        except FileNotFoundError:
            print(f'problem with {gca_id}: not annotated?')
            gene_found_dict[gca_id] = None
            continue
        gene_found_dict[gca_id] = hmmr_hits_1

    with open(outfile,'w') as f:
        json.dump(gene_found_dict,f,indent=2)       
        
        
if __name__ == '__main__':

    logging.basicConfig(
        filename='ras_rips_hmm.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(
        description="Run HMM neighbor detection pipeline"
    )

    parser.add_argument(
        "--embldir",
        default=str(BASE_PATH / "combined_bakta_original_annotations_plus_new"),
        help="Path to directory containing .embl.gz files"
    )

    parser.add_argument(
        "--pfam1",
        default="PF04055",
        help="First Pfam ID (default: PF04055)"## rasripps
    )

    parser.add_argument(
        "--pfam2",
        default="PF21259",
        help="Second Pfam ID (default: PF21259)"### RGG
    )

    parser.add_argument(
        "--outfolder",
        default="output_ras_ripps_neighbors",
        help="Output filename (default: fooneighbors)"
    )

    args = parser.parse_args()

    proteome_path = BASE_PATH / 'proteome'
    proteome_path.mkdir(parents=True, exist_ok=True)

    embldir = Path(args.embldir)

    print(f"Using EMBL directory: {embldir}")
    print(f"Using Pfam IDs: {args.pfam1}, {args.pfam2}")
    print(f"Output will be written to: {args.outfile}")
    
    print("Creating proteome database...")
    create_proteome_db(embldir, proteome_path)
    print("Running HMM searches...")
    all_hmmsearch(f"{args.pfam1}.hmm", proteome_path)
    print(f"Completed HMM search for {args.pfam1}") 
    all_hmmsearch(f"{args.pfam2}.hmm", proteome_path)
    print(f"Completed HMM search for {args.pfam2}") 
    
    print("Finding neighboring hits...")
    
    get_all_neighbors(args.pfam1, args.pfam2, args.outfolder, embldir)
    print("Neighbor detection complete.")
    
    print("Getting single gene hits for Pfam1...")
    
    get_single_gene_hit(args.pfam1, f'{args.outfolder}/{args.pfam1}_hits.json', embldir)
    print(f"Single gene hit detection for {args.pfam1} complete. Saved to {args.pfam1}_hits.json")
    

