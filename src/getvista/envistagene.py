#!/usr/bin/env python

"""
File: envistagene.py
Author: Jake Leyhr
GitHub: https://github.com/jakeleyhr/GetVISTA/
Date: February 2024
Description: Query the Ensembl database with species and gene name to obtain FASTA file and gene feature coordinates in pipmaker format
"""

 # Import dependencies
import sys
import json
import time
import argparse
import requests
from shutil import get_terminal_size
from Bio import SeqIO
from Bio.Seq import Seq
from getvista.version_check import check_for_updates

# Python 2/3 adaptability
try:
    from urllib.parse import urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('ERROR: Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
                print('Check species and/or gene name')
                sys.exit()
        return data

    def get_genes_in_genomic_coordinates(self, species, genomic_coordinates):
        genes = self.perform_rest_action(
            endpoint='/overlap/region/{0}/{1}'.format(species, genomic_coordinates),
            params={'feature': 'gene'}
        )
        genes_sorted = sorted(genes, key=lambda x: (x['end'] + x['start']) / 2) # Sort list by middle gene coordinate
        return genes_sorted


# Function 1 - get genomic coordinates and FASTA sequence
def download_dna_sequence(species, gene_name, start_adjust, end_adjust, gene_oriented_adjustment):
    # Ensembl REST API endpoint for gene lookup
    endpoint = f"http://rest.ensembl.org/lookup/symbol/{species}/{gene_name}"

    # Make the request to the Ensembl REST API (times out after 60 seconds)
    response = requests.get(endpoint, headers={"Content-Type": "application/json"}, timeout = 60)
    
    # Get start and end coordinates from gene lookup
    if response.status_code == 200:
        gene_info = response.json()
        if 'start' in gene_info and 'end' in gene_info:
            gene_start_coordinate = gene_info['start']
            gene_end_coordinate = gene_info['end']
            gene_ref_name = gene_info['display_name']
        else:
            print(f"Coordinates not found for gene: {gene_name}")
    else:
        print(f"ERROR: Failed to retrieve gene information. Status code: {response.status_code}")
        print("Check species and gene name are correct.")
        sys.exit()
    
    # Get strand information:
    strandint = gene_info['strand']
    if strandint == 1:
        strand = 'forward'
    elif strandint == -1:
        strand = 'reverse'

    if gene_oriented_adjustment and strand == 'reverse':
        # Adjust genomic coordinates according to base adjust inputs oriented with the gene direction
        start = gene_start_coordinate-end_adjust
        if start < 1:
            print(f"WARNING: {start} is not a valid start coordinate, changing to 1.")
            start = 1
        end = gene_end_coordinate+start_adjust
    else:
        # Adjust genomic coordinates according to gene and base adjust inputs
        start = gene_start_coordinate-start_adjust
        if start < 1:
            print(f"WARNING: {start} is not a valid start coordinate, changing to 1.")
            start = 1
        end = gene_end_coordinate+end_adjust

    genomic_coordinates = f"{gene_info['seq_region_name']}:{start}-{end}" 

    if (end - start + 1) > 5000000:
        print("ERROR: Query sequence must be under 5Mb")
        sys.exit()
    
    # Ensembl REST API endpoint for DNA sequences in FASTA format
    url = f"https://rest.ensembl.org/sequence/region/{species}/{genomic_coordinates}?format=fasta" 
 
    # Specify the headers with the required Content-Type
    headers = {"Content-Type": "text/x-fasta"} 
    # Make the request to the Ensembl REST API with the headers (times out after 60 seconds)
    response = requests.get(url, headers=headers, timeout = 60) 
    
    # Print details
    print(f"Assembly name: {gene_info['assembly_name']}")
    print(f"{species} {gene_name} coordinates: {gene_info['seq_region_name']}:{gene_start_coordinate}-{gene_end_coordinate}")
    print(f"{species} {gene_name} is on {strand} strand")
    print(f"{species} {gene_name} sequence length: {gene_end_coordinate-gene_start_coordinate+1}bp")
    print(f"Specified coordinates: {genomic_coordinates}")

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the DNA sequence from the response
        fasta_lines = response.text.strip().split('\n')
    else:
        print(f"Error: Unable to retrieve DNA sequence. Status code: {response.status_code}")
        print(f"Response content: {response.text}")
        sys.exit()
    
    return genomic_coordinates, fasta_lines, strand, gene_ref_name


# Function 2 - get gene feature coordinates in pipmaker format
def pipmaker(genes, genomic_coordinates, apply_reverse_complement, nocut, all_transcripts):
    client = EnsemblRestClient()
    
    input_region_start = int(genomic_coordinates.split(":")[1].split("-")[0])
    input_region_end = int(genomic_coordinates.split(":")[1].split("-")[1])
    sequence_length = (input_region_end - input_region_start) + 1

    print(f"Specified sequence length: {sequence_length}bp")
    print("\nTranscripts included in region:")

    coordinates_content = ""
    for gene in genes:
        gene_id = gene['id']
        gene_info = client.perform_rest_action(
            '/lookup/id/{0}'.format(gene_id),
            params={'expand': '1', 'utr': '1'}
        )
        if gene_info and 'Transcript' in gene_info:
            transcripts = gene_info['Transcript']
            new_start = 1

            for transcript in transcripts:
                if all_transcripts:
                    filter_type = 'all'
                else:
                    filter_type = 'canonical'

                if filter_type == 'all' or ('is_canonical' in transcript and transcript['is_canonical'] == 1):
                    strand_indicator = ">" if gene_info['strand'] == 1 else "<"
                    start_position = transcript.get('start') - input_region_start + new_start
                    end_position = transcript.get('end') - input_region_start + new_start
                    transcript_name = transcript.get('display_name', transcript['id'])
                    
                    if start_position < 0 and end_position < 0:
                        print(f"({transcript_name} transcript out of 5' range:{start_position}:{end_position})")
                    elif start_position > sequence_length and end_position > sequence_length:
                        print(f"({transcript_name} transcript out of 3' range:{start_position}:{end_position})")
                    else:
                        print(transcript_name)

                    if nocut == False:
                        if not ((start_position < 0 and end_position < 0) or (start_position > sequence_length and end_position > sequence_length)):
                            if start_position < 0:
                                if apply_reverse_complement:
                                    transcript_name += f"-cut3':{1 - start_position}bp"
                                else:
                                    transcript_name += f"-cut5':{1 - start_position}bp"
                                start_position = 1
                            if end_position > sequence_length:
                                if apply_reverse_complement:
                                    transcript_name += f"-cut5':{end_position - sequence_length}bp"
                                else:
                                    transcript_name += f"-cut3':{end_position - sequence_length}bp"
                                end_position = sequence_length
                            coordinates_content += f"{strand_indicator} {start_position} {end_position} {transcript_name}\n" # Assemble header line
                    else: 
                        coordinates_content += f"{strand_indicator} {start_position} {end_position} {transcript_name}\n" # Assemble header line

                    coordinates = []

                    if 'Exon' in transcript:
                        exons = transcript['Exon']
                        for exon in exons:
                            start = exon.get('start') - input_region_start + new_start
                            end = exon.get('end') - input_region_start + new_start
                            
                            if nocut == False:
                                if start < 0 and end < 0:
                                    continue
                                if start > sequence_length and end > sequence_length:
                                    continue
                                if start < 0:
                                    start = 1
                                if end > sequence_length:
                                    end = sequence_length

                            utr_start = 0
                            utr_end = 0

                            if 'UTR' in transcript:
                                utrs = transcript['UTR']
                                for utr in utrs:
                                    utr_start = utr.get('start') - input_region_start + new_start
                                    utr_end = utr.get('end') - input_region_start + new_start
                                    
                                    if nocut == False:
                                        if utr_start < sequence_length and utr_end > sequence_length:
                                            utr_end = sequence_length
                                        #print(f"exon start: {start}")
                                        #print(f"exon end: {end}")
                                        #print(f"utr start: {utr_start}")
                                        #print(f"utr end: {utr_end}")
                                        # if CDS and UTR exon is 5' cut off, make sure UTR takes priority
                                        if end == utr_end and utr_start<start:
                                            print('true')
                                            start = 0
                                            end = 0
                                            utr_start=1

                                    if start == utr_start and end == utr_end:
                                        start = 0
                                        end = 0
                                    elif utr_start <= start <= utr_end:
                                        start = utr_end + 1
                                    elif utr_start <= end <= utr_end:
                                        end = utr_start - 1

                                    if nocut == False:
                                        if start > sequence_length and end > sequence_length:
                                            start = 0
                                            end = 0
                                            continue
                                        if end > sequence_length:
                                            end = sequence_length
                                            continue
                                        if utr_start < 0 and utr_end < 0:
                                            continue
                                        if utr_start < 0:
                                            utr_start = 1
                                        

                            coordinates.append((f"{start} {end} exon", start))

                        if 'UTR' in transcript:
                            utrs = transcript['UTR']
                            for utr in utrs:
                                utr_start_abs = utr.get('start')
                                utr_end_abs = utr.get('end')
                                utr_start = utr_start_abs - input_region_start + new_start
                                utr_end = utr_end_abs - input_region_start + new_start

                                if nocut == False:
                                    if utr_start < 0 and utr_end <= 0:
                                        continue
                                    if utr_start > sequence_length and utr_end > sequence_length:
                                        continue
                                    if utr_start < 0 and utr_end > 0:
                                        utr_start = 1
                                    if utr_start < sequence_length and utr_end > sequence_length:
                                        utr_end = sequence_length

                                coordinates.append((f"{utr_start} {utr_end} UTR", utr_start))

                    else:
                        print("No exons found in the transcript.")

                    coordinates.sort(key=lambda x: x[1])
                    for coord, _ in coordinates:
                        if coord.strip() != "0 0 exon" and "0 0 UTR":
                            coordinates_content += f"{coord}\n"
                    coordinates_content += "\n\n"

    return coordinates_content, sequence_length

# Function A - reverse pipmaker coordinates
def reverse_coordinates(coordinates, sequence_length):
    # Prepare to collect reversed coordianates
    reversed_coordinates = []

    # Prepare to collect lines for each transcript
    transcript_lines = []
    for line in coordinates.split('\n'):
        if line:
            fields = line.split()
            # For lines with 3 fields (without strand information i.e. not header lines):
            if len(fields) == 3: 
                start = int(fields[0])
                end = int(fields[1])
                feature = fields[2]

                # Swap start and end values, subtract from sequence_length to flip coordinates
                new_start = sequence_length - end + 1
                new_end = sequence_length - start + 1

                transcript_lines.append((new_start, f"{new_start} {new_end} {feature}")) # Append new coordinates to object
            
            # For lines with 4 fields (including strand information i.e. header lines)
            elif len(fields) == 4:
                strand = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                feature = fields[3]

                # Swap start and end values, subtract from sequence_length to flip coordinates
                new_start = sequence_length - end + 1
                new_end = sequence_length - start + 1

                # Swap < for > and vice versa
                new_strand = '>' if strand == '<' else '<'

                transcript_lines.append((new_start, f"{new_strand} {new_start} {new_end} {feature}")) # Append new coordinates to object

        elif transcript_lines:
            sorted_lines = sorted(transcript_lines, key=lambda x: x[0]) # Sort lines with 3 fields based on the new start coordinate
            reversed_coordinates.extend(line for _, line in sorted_lines)
            reversed_coordinates.append("")  # Add an empty line between transcripts
            reversed_coordinates.append("")  # Add an empty line between transcripts
            transcript_lines = []  # Reset for the next transcript

    return '\n'.join(reversed_coordinates)

#Function B - use flanking genes
def useflanks(species, genes, genomic_coordinates, flank):
    if flank not in ["in", "ex"]:
        print("ERROR: invalid flank argument; must be 'in' or 'ex'")
        sys.exit()

    chrom = genomic_coordinates.split(":")[0]
    input_region_start = int(genomic_coordinates.split(":")[1].split("-")[0])
    input_region_end = int(genomic_coordinates.split(":")[1].split("-")[1])
    sequence_length = (input_region_end - input_region_start) + 1

    print(f"Specified sequence length: {sequence_length}bp")
    print("\nGenes included in region:")
    genes_data = []  # Initialize an empty list to store gene data
    for gene in genes:
        gene_id = gene['id']
        try:
            gene_name = gene['external_name']
        except KeyError:
            gene_name=gene_id

        print(gene_name)

        # Store gene data as a dictionary
        gene_data = {
            'gene_name': gene_name,
            'start_position': gene['start'],
            'end_position': gene['end'],
        }
    
        # Append the dictionary to the list
        genes_data.append(gene_data)

    # Initialize an empty dictionary to map gene names to start positions
    gene_start_positions = {} 
    gene_end_positions = {} 
    for gene_data in genes_data:
        gene_name = gene_data['gene_name']
        start_position = min(gene_data['start_position'], gene_data['end_position'])
        end_position = max(gene_data['start_position'], gene_data['end_position'])

        if flank == "in":
            # Store gene start positions in the dictionary
            gene_start_positions[gene_name] = start_position
            gene_end_positions[gene_name] = end_position
        elif flank == "ex":
            gene_start_positions[gene_name] = end_position+1
            gene_end_positions[gene_name] = start_position-1
    while True:      
        # Prompt the user for input
        gene_name_input1 = input("Please enter the first gene name (case sensitive): ")
        # Check if the gene name is present in the dictionary
        if gene_name_input1 in gene_start_positions:
            start_coordinate = gene_start_positions[gene_name_input1]
            if flank == "in":
                print(f"The start coordinate for {gene_name_input1} is: {start_coordinate}")
            if flank == "ex":
                print(f"The start coordinate after {gene_name_input1} is: {start_coordinate}")
            break
        else:
            choice = input("Invalid input. Do you want to try again? (yes/no): ")
            if choice.lower() != "yes":
                print("Terminating the script.")
                sys.exit()  # Terminate the loop if the user chooses not to try again

    while True:      
        # Prompt the user for input
        gene_name_input2 = input("Please enter the second gene name (case sensitive): ")
        # Check if the gene name is present in the dictionary
        if gene_name_input2 in gene_end_positions:
            end_coordinate = gene_end_positions[gene_name_input2]
            if flank == "in":
                print(f"The end coordinate for {gene_name_input2} is: {end_coordinate}")
            if flank == "ex":
                print(f"The end coordinate before {gene_name_input2} is: {end_coordinate}")
            break
        else:
            choice = input("Invalid input. Do you want to try again? (yes/no): ")
            if choice.lower() != "yes":
                print("Terminating the script.")
                sys.exit()  # Terminate the loop if the user chooses not to try again

    new_genomic_coordinates = (f'{chrom}:{start_coordinate}-{end_coordinate}')      
    print(f'\nNew genomic coordinates: {new_genomic_coordinates}')

    # Download the FASTA file again 
    # Ensembl REST API endpoint for DNA sequences in FASTA format
    url = f"https://rest.ensembl.org/sequence/region/{species}/{new_genomic_coordinates}?format=fasta" 
 
    # Specify the headers with the required Content-Type
    headers = {"Content-Type": "text/x-fasta"} 
    # Make the request to the Ensembl REST API with the headers (times out after 60 seconds)
    response = requests.get(url, headers=headers, timeout = 60) 
    
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the DNA sequence from the response
        fasta_lines = response.text.strip().split('\n')
    else:
        print(f"Error: Unable to retrieve DNA sequence. Status code: {response.status_code}")
        print(f"Response content: {response.text}")
        sys.exit()

    return new_genomic_coordinates, gene_name_input1, gene_name_input2, fasta_lines


#Function C - add graphical display
def graphic(coordinates_content, sequence_length, all_transcripts):
    graphic_coordinates = []
    for line in coordinates_content.split('\n'):
        if line.startswith('>') or line.startswith('<'):
            parts = line.split()
            gene_start = int(parts[1])
            gene_end = int(parts[2])
            gene_name = parts[3]
            direction = parts[0]
            graphic_coordinates.append((gene_start, gene_end, gene_name, direction))
    if all_transcripts:
        bp_per_dash = 30 
    else:
        bp_per_dash = 70
    num_dashes = sequence_length // bp_per_dash
    genomic_map = ['='] * num_dashes

    # Sort gene coordinates by start and end coordinates
    graphic_coordinates.sort(key=lambda x: (x[0], x[1]))
    #print(graphic_coordinates)
    adjusted = False
    # Adjust start coordinates that are under bp_per_dash
    for i, (gene_start, _, _, _) in enumerate(graphic_coordinates):
        if gene_start <= bp_per_dash and not adjusted:
            adjusted = True
            continue
        if gene_start <= bp_per_dash and adjusted:
            graphic_coordinates[i] = (graphic_coordinates[i][0] + bp_per_dash*i, *graphic_coordinates[i][1:])
            adjusted = True

    # Proceed to generate graphic
    for i in range(len(graphic_coordinates)):
        gene_start, gene_end, gene_name, direction = graphic_coordinates[i]
        gene_start_dash = (gene_start - 1) // bp_per_dash
        gene_end_dash = (gene_end - 1) // bp_per_dash

        if gene_start_dash == gene_end_dash:
            # Gene fits within one dash
            genomic_map[gene_start_dash] = f"{direction}{gene_name}{direction}"
        else:
            # Gene spans multiple dashes
            genomic_map[gene_start_dash] = f"{direction}{gene_name}{direction}"
            for j in range(gene_start_dash + 1, min(gene_end_dash + 1, len(genomic_map))):
                genomic_map[j] = ''

        # Adjust start coordinate if there are overlapping genes with same start
        if i < len(graphic_coordinates) - 1:
            next_gene_start, next_gene_end, _, _ = graphic_coordinates[i + 1]
            if gene_start == next_gene_start:
                if gene_end > next_gene_end:
                    graphic_coordinates[i] = list(graphic_coordinates[i])
                    graphic_coordinates[i][0] += 1
                    graphic_coordinates[i] = tuple(graphic_coordinates[i])
                else:
                    graphic_coordinates[i + 1] = list(graphic_coordinates[i + 1])
                    graphic_coordinates[i + 1][0] += 1
                    graphic_coordinates[i + 1] = tuple(graphic_coordinates[i + 1])

    # Combine adjacent directionality markers
    genomic_map_combined = ''.join(genomic_map).replace("<>", "#").replace("><", "#").replace("<<", "#").replace(">>", "#")
    print(f"\nGraphical representation of specified sequence region:")
    print(genomic_map_combined)


def engene(
    species, 
    gene_name, 
    start_adjust, 
    end_adjust, 
    gene_oriented_adjustment=False,
    fasta_output_file=None, 
    coordinates_output_file=None, 
    all_transcripts=None, 
    nocut=None, 
    apply_reverse_complement=False, 
    autoname=False, 
    fw=False, 
    flank=None,
    vis=False
):
    # Print run header
    terminalwidth = get_terminal_size()[0]
    nameswidth = len(f" {species} {gene_name} ")
    leftindent = ((terminalwidth-nameswidth)//2)
    print("▒"*terminalwidth+
          "▒"*leftindent+ 
          f" {species} {gene_name} "+ 
          "▒"*(terminalwidth-leftindent-nameswidth)+
          "▒"*terminalwidth)
    
     # Get genomic coordinates from gene record and adjustments
    genomic_coordinates, fasta_lines, strand, gene_ref_name = download_dna_sequence(species, gene_name, start_adjust, end_adjust, gene_oriented_adjustment)
 

    # Get genes info
    client = EnsemblRestClient()
    genes = client.get_genes_in_genomic_coordinates(species, genomic_coordinates)

    # If no genes in genes list, terminate script
    if len(genes) == 0:
        print("\nNo genes found in region. Terminating script\n")
        sys.exit()

    if flank:  
        genomic_coordinates, fgene1, fgene2, fasta_lines = useflanks(species, genes, genomic_coordinates, flank)
        genes = client.get_genes_in_genomic_coordinates(species, genomic_coordinates)

    # If no genes in genes list
    if len(genes) == 0:
        print("No genes found in region\n")
        #sys.exit()

    # Get coordinates and sequence length
    coordinates_content, sequence_length = pipmaker(genes, genomic_coordinates, apply_reverse_complement, nocut, all_transcripts)
    #print(coordinates_content)

    # If -fw argument used and the gene is on the reverse strand, force reverse complement
    if fw and strand == 'reverse':
        print(f'\n{gene_name} is on the reverse strand, flipped automatically.')
        apply_reverse_complement = True

    genomic_coordinates_fixed = genomic_coordinates.replace(":", ".")
        
    if autoname and (start_adjust or end_adjust) and not flank:
        if not fasta_output_file:
            if not apply_reverse_complement:
                fasta_output_file = f"{species}_{gene_ref_name}_{genomic_coordinates_fixed}.fasta.txt"
            else:
                fasta_output_file = (f"{species}_{gene_ref_name}_{genomic_coordinates_fixed}_revcomp.fasta.txt")
        if not coordinates_output_file:
            if not apply_reverse_complement:
                coordinates_output_file = (f"{species}_{gene_ref_name}_{genomic_coordinates_fixed}.annotation.txt")
            else:
                coordinates_output_file = (f"{species}_{gene_ref_name}_{genomic_coordinates_fixed}_revcomp.annotation.txt")
    if autoname and (flank == "in"):
        if not fasta_output_file:
            if not apply_reverse_complement:
                fasta_output_file = f"{species}_{gene_ref_name}__{fgene1}-{fgene2}.fasta.txt"
            else:
                fasta_output_file = (f"{species}_{gene_ref_name}__{fgene2}-{fgene1}_revcomp.fasta.txt")
        if not coordinates_output_file:
            if not apply_reverse_complement:
                coordinates_output_file = (f"{species}_{gene_ref_name}__{fgene1}-{fgene2}.annotation.txt")
            else:
                coordinates_output_file = (f"{species}_{gene_ref_name}__{fgene2}-{fgene1}_revcomp.annotation.txt")    
    if autoname and (flank == "ex"):
        if not fasta_output_file:
            if not apply_reverse_complement:
                fasta_output_file = f"{species}_{gene_ref_name}__{fgene1}-{fgene2}.ex.fasta.txt"
            else:
                fasta_output_file = (f"{species}_{gene_ref_name}__{fgene2}-{fgene1}.ex_revcomp.fasta.txt")
        if not coordinates_output_file:
            if not apply_reverse_complement:
                coordinates_output_file = (f"{species}_{gene_ref_name}__{fgene1}-{fgene2}.ex.annotation.txt")
            else:
                coordinates_output_file = (f"{species}_{gene_ref_name}__{fgene2}-{fgene1}.ex_revcomp.annotation.txt")       
    if autoname and not (start_adjust or end_adjust or flank):
        if not fasta_output_file:
            if not apply_reverse_complement:
                fasta_output_file = f"{species}_{gene_ref_name}.fasta.txt"
            else:
                fasta_output_file = (f"{species}_{gene_ref_name}_revcomp.fasta.txt")
        if not coordinates_output_file:
            if not apply_reverse_complement:
                coordinates_output_file = (f"{species}_{gene_name}.annotation.txt")
            else:
                coordinates_output_file = (f"{species}_{gene_name}_revcomp.annotation.txt")       

    # Check if ".txt" is already at the end of the coordinates_output_file argument
    if coordinates_output_file and not coordinates_output_file.endswith(".txt"):
        coordinates_output_file += ".txt"

    if apply_reverse_complement:
        coordinates_content = reverse_coordinates(coordinates_content, sequence_length) # Reverse the coordinates
        message=f"\nReversed coordinates saved to "
        if vis:
            graphic(coordinates_content, sequence_length, all_transcripts)
    else:
        message=f"\nCoordinates saved to "
        if vis:
            graphic(coordinates_content, sequence_length, all_transcripts)
            
    # If run with -anno (save annotation coordinates):
    if coordinates_output_file:
        with open(coordinates_output_file, 'w') as coordinates_file:
                coordinates_file.write(coordinates_content)
        print(f"{message}{coordinates_output_file}")
    else:
        print("\nNo coordinates output file specified")

    # Check if ".txt" is already at the end of the fasta_output_file argument
    if fasta_output_file and not fasta_output_file.endswith(".txt"):
        fasta_output_file += ".txt"

    # Split FASTA header line from DNA sequence
    header_line = fasta_lines[0][1:] # first line is the header, remove the >
    dna_sequence = Seq(''.join(fasta_lines[1:]))
    # Save the DNA sequence to the specified output file
    if fasta_output_file:
        if not apply_reverse_complement:
            # Create a SeqRecord object and save it in FASTA format
            fasta_record = SeqIO.SeqRecord(
                dna_sequence, 
                id=str(header_line), 
                description="",
                )
            print(f"DNA sequence saved to {fasta_output_file}\n")
        else:
            # Reverse complement the sequence
            reverse_complement_sequence = str(dna_sequence.reverse_complement())
            last_colon_1_index = header_line.rfind(":1") # Find the strand indicator in FASTA header
            reversed_header_line = header_line[:last_colon_1_index] + ":-1" + header_line[last_colon_1_index+2:] # Replace ":1" with ":-1" in the header line
            # Create a SeqRecord object and save it in FASTA format
            fasta_record = SeqIO.SeqRecord(
                Seq(reverse_complement_sequence),
                id=reversed_header_line,
                description="",
            )
            print(f"Reverse complement DNA sequence saved to {fasta_output_file}\n")
        # Write to file
        SeqIO.write(fasta_record, f"{fasta_output_file}", "fasta")   
    else:
        print("No FASTA output file specified\n") # Print to notify user neither sequence specified as output


def main():
    #Check for updates
    check_for_updates()
    
    # Create an ArgumentParser
    parser = argparse.ArgumentParser(description="Query the Ensembl database with a species and gene name to obtain DNA sequences in FASTA format and gene feature coordinates in pipmaker format.")
    
    # Add arguments for species, gene_name, etc
    parser.add_argument("-s", "--species", nargs='+', required=True, help="Species name(s) (e.g., 'Homo_sapiens' or 'Human')")
    parser.add_argument("-g", "--gene_name", required=True, help="Gene name (e.g. BRCA1 or brca1)")
    parser.add_argument("-sa", "--start_adjust", type=int, default=0, help="Number to subtract from the start coordinate (default: 0)")
    parser.add_argument("-ea", "--end_adjust", type=int, default=0, help="Number to add to the end coordinate (default: 0)")
    parser.add_argument("-goa", "--gene_oriented_adjustment", action="store_true", default=False, help="Make start/end adjustments follow gene orientation, not assembly orientation")
    parser.add_argument("-fasta", "--fasta_output_file", default=None, help="Output file name for the DNA sequence in FASTA format")
    parser.add_argument("-anno", "--coordinates_output_file", default=None, help="Output file name for the gene coordinates in pipmaker format")
    parser.add_argument("-all", "--all_transcripts", action="store_true", default=False, help="Include all transcripts (instead of canonical transcript only)")
    parser.add_argument("-nocut", action="store_true", default=False, help="Keep feature annotations not included in sequence")
    parser.add_argument("-rev", action="store_true", default=False, help="Reverse complement DNA sequence and coordinates")
    parser.add_argument("-autoname", action="store_true", default=False, help="Automatically generate output file names based on species and gene name")
    parser.add_argument("-fw", action="store_true", default=False, help="Automatically orient the gene in the forward strand by reverse complementing if needed")
    parser.add_argument("-flank", default=None, help="Select 2 genes to specify new range. 'in' to include the flanking genes, 'ex' to exclude them")
    parser.add_argument("-vis", action="store_true", default=False, help="Display graphical representation of the sequence in the terminal")


    # Parse the command-line arguments
    args = parser.parse_args()

    # Loop through multi-species inputs
    for species in args.species:
        # Pass arguments in the correct order
        engene(
            species, 
            args.gene_name, 
            args.start_adjust, 
            args.end_adjust, 
            args.gene_oriented_adjustment,
            args.fasta_output_file, 
            args.coordinates_output_file, 
            args.all_transcripts, 
            args.nocut, 
            args.rev,
            args.autoname,
            args.fw,
            args.flank,
            args.vis,
        )
        

if __name__ == '__main__':
    main()


