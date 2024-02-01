#!/usr/bin/env python

"""
File: gbgenerecord.py
Author: Jake Leyhr
GitHub: https://github.com/jakeleyhr/GetVISTA/
Date: January 2024
Description: Query the GenBank database with a species and gene name to obtain a list of different \
    records containing the sequence to inform use of gbvistagene.py
"""

# Import dependencies
import argparse
from Bio import Entrez
from getvista.version_check import check_for_updates

# Function #1 - get gene record
def search_gene_info(species, gene_name):
    # Set your email address
    Entrez.email = "dummy@gmail.com"

    # Build the query
    query = f"{species}[ORGN] AND {gene_name}[Gene Name]" # Strict check on species name and gene name (also searches gene name synonyms)

    # Search Entrez Gene
    handle = Entrez.esearch(db="gene", term=query, retmode="xml")
    record = Entrez.read(handle)
    # print(f' Entrez record: {record}')

    # Fetch gene information from each record
    if "IdList" in record and record["IdList"]:
        gene_id = record["IdList"][0]
        #print(gene_id)
        # Fetch gene information
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        gene_record = Entrez.read(handle)
        return gene_record

    return None


def gbrecord(species, gene_name):
    # Get target gene record from Entrez
    gene_info = search_gene_info(species, gene_name)
    
    # Extract the relevant information
    if gene_info:
        #print(gene_info[0])
        gene_ref_name = gene_info[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
        print(f"Query gene: {gene_ref_name}")
        try:
            gene_ref_desc = gene_info[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
            print(f'Description: {gene_ref_desc}') 
        except KeyError:
            print('Description: None available')
        try:
            synonyms = gene_info[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
            print(f'Synonyms: {synonyms}')
        except KeyError:
            print('Synonyms: None available')
        try:
            locus = gene_info[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_maploc']
            print(f'Locus: {locus}') 
        except KeyError:
            print('Locus: None available')
        
        # Print each record information
        for i in range(0, 999):
            try:
                assembly = gene_info[0]['Entrezgene_locus'][i]['Gene-commentary_label']
                accession = gene_info[0]['Entrezgene_locus'][i]['Gene-commentary_accession']
                start = int(gene_info[0]['Entrezgene_locus'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from'])
                end = int(gene_info[0]['Entrezgene_locus'][i]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to'])
                length = end - start + 1

                print("") 
                print(f"RECORD {i}")
                print(f"Assembly: {assembly}")
                print(f'Accession: {accession}')
                print(f'Location: {start}:{end}')
                print(f'Length: {length}')
            except IndexError:
                continue
    else:
        print("No records found. Check species and gene names.")

def main():
# Create an ArgumentParser
    parser = argparse.ArgumentParser(description="Query the GenBank database with a species and gene name \
                                     to obtain a list of different records containing the sequence to inform \
                                     use of the gbgene module.")

    # Add arguments for species and gene_name
    parser.add_argument("-s", "--species", help="Species name (e.g., 'Homo_sapiens' or 'Human')", required=True)
    parser.add_argument("-g", "--gene_name", help="Gene name (e.g. BRCA1 or brca1)", required=True)

    # Parse the command-line arguments
    args = parser.parse_args()

    gbrecord(
        args.species, 
        args.gene_name,
    )

if __name__ == "__main__":
    check_for_updates()
    main()

    