#!/usr/bin/env python

import sys
import json
import time
import requests
import argparse

# Python 2/3 adaptability
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
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
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_genes_in_genomic_coordinates(self, species, genomic_coordinates):
        genes = self.perform_rest_action(
            endpoint='/overlap/region/{0}/{1}'.format(species, genomic_coordinates),
            params={'feature': 'gene'}
        )
        return genes

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    cleaned_sequence = sequence.replace('\n', '').upper()
    valid_sequence = ''.join(base for base in cleaned_sequence if base in complement_dict)
    return ''.join(complement_dict[base] for base in reversed(valid_sequence))


def reverse_coordinates(coordinates, sequence_length):
    reversed_coordinates = []
    transcript_lines = []  # Collect lines for each transcript
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


def download_dna_sequence(species, gene_name, start_adjust, end_adjust, output_filename=None, apply_reverse_complement=False, verbose=True):
    # Ensembl REST API endpoint for gene lookup
    endpoint = f"http://rest.ensembl.org/lookup/symbol/{species}/{gene_name}"

    # Make the request to the Ensembl REST API
    response = requests.get(endpoint, headers={"Content-Type": "application/json"})
    
    # Get start and end coordinates from gene lookup
    if response.status_code == 200:
        gene_info = response.json()
        if 'start' in gene_info and 'end' in gene_info:
            gene_start_coordinate = gene_info['start']
            gene_end_coordinate = gene_info['end']
        else:
            print(f"Coordinates not found for gene: {gene_name}")
    else:
        print(f"Failed to retrieve gene information. Status code: {response.status_code}")

    genomic_coordinates = f"{gene_info['seq_region_name']}:{gene_start_coordinate-start_adjust}-{gene_end_coordinate+end_adjust}" # Adjust genomic coordinates according to gene and base adjust inputs

    url = f"https://rest.ensembl.org/sequence/region/{species}/{genomic_coordinates}?format=fasta" # Ensembl REST API endpoint for DNA sequences in FASTA format
    headers = {"Content-Type": "text/x-fasta"} # Specify the headers with the required Content-Type
    response = requests.get(url, headers=headers) # Make the request to the Ensembl REST API with the headers
    if verbose:
        print(f"Assembly name: {gene_info['assembly_name']}")
        print(f"{species} {gene_name} coordinates: {gene_info['seq_region_name']}:{gene_start_coordinate}-{gene_end_coordinate}")
        print(f"{species} {gene_name} is on {gene_info['strand']} strand")
        print(f"{species} {gene_name} sequence length: {gene_end_coordinate-gene_start_coordinate+1}bp")
        print(f"Extracting coordinates: {genomic_coordinates}")

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the DNA sequence from the response
        fasta_lines = response.text.strip().split('\n')
        header_line = fasta_lines[0]
        dna_sequence = ''.join(fasta_lines[1:])

        # Save the DNA sequence to the specified output file
        if output_filename:
            with open(output_filename, 'w') as output_file:
                # If run with -rev option, reverse complement DNA and change strand indicator in header
                if apply_reverse_complement:
                    reverse_sequence = reverse_complement(dna_sequence) # Apply reverse_complement function to DNA sequence
                    last_colon_1_index = header_line.rfind(":1") # Find the strand indicator in FASTA header
                    reversed_header_line = header_line[:last_colon_1_index] + ":-1" + header_line[last_colon_1_index+2:] # Replace ":1" with ":-1" in the header line
                    output_file.write(reversed_header_line + '\n' + reverse_sequence + '\n') # Write the header and sequence to the output file 
                    print(f"Reverse complement DNA sequence saved to {output_filename}")
                
                # Otherwise, write original FASTA to file
                else:
                    output_file.write(header_line + '\n' + dna_sequence + '\n') # Write the header and sequence to the output file 
                    if verbose:
                        print(f"DNA sequence saved to {output_filename}")
        else:
            print("No sequence output file specified")
    else:
        print(f"Error: Unable to retrieve DNA sequence. Status code: {response.status_code}")
        print(f"Response content: {response.text}")
    
    return genomic_coordinates


def run_gene_coordinates(species, gene_name, start_adjust, end_adjust, fasta_output_file=None, coordinates_output_file=None, all_transcripts=None, nocut=None, apply_reverse_complement=False):
    genomic_coordinates = download_dna_sequence(species, gene_name, start_adjust, end_adjust, fasta_output_file, apply_reverse_complement)
    client = EnsemblRestClient()
    gene_info = client.get_genes_in_genomic_coordinates(species, genomic_coordinates)
    if gene_info:
        if coordinates_output_file:
            with open(coordinates_output_file, 'w') as coordinates_file:
                input_region_start = int(genomic_coordinates.split(":")[1].split("-")[0])  # Extracting start position from the input region
                input_region_end = int(genomic_coordinates.split(":")[1].split("-")[1])  # Extracting end position from the input region
                sequence_length = (input_region_end - input_region_start) + 1 # Calculating total sequence length
                print(f"Total sequence length: {sequence_length}bp")

                for gene in gene_info:
                    gene_id = gene['id']
                    gene_info = client.perform_rest_action(
                        '/lookup/id/{0}'.format(gene_id),
                        params={'expand': '1', 'utr': '1'} # Essential for obtaining UTR information
                    )
                    #print("Gene Info:", gene_info)  # Add this line for debugging
                    if gene_info and 'Transcript' in gene_info:
                        transcripts = gene_info['Transcript']
                        gene_start = gene_info['start']
                        new_start = 1

                        for transcript in transcripts:
                            #print("Transcript:", transcript)  # Add this line for debugging

                            # Choose to return annotations for all transcripts or just the canonical ones
                            if args.all_transcripts:
                                filter_type = 'all'
                            else:
                                filter_type = 'canonical'

                            if filter_type == 'all' or ('is_canonical' in transcript and transcript['is_canonical'] == 1):
                                strand_indicator = ">" if gene_info['strand'] == 1 else "<" # Determines strand direction
                                start_position = transcript.get('start') - input_region_start + new_start # Calculate transcript start coordinate relative to region
                                #print(f'transcript start: {start_position}')
                                end_position = transcript.get('end') - input_region_start + new_start # Calculate transcript end coordinate relative to region
                                #print(f'transcript end: {end_position}')
                                transcript_name = transcript.get('display_name', transcript['id']) # Get transcript name

                                # If run with -nocut option:
                                if nocut:
                                    if start_position < 0 and end_position < 0: # If entire transcript is out of region range, ignore it. (Shouldn't occur.)
                                        print("Some transcripts out of range.")
                                        continue

                                    if start_position < 0: # If only the start of the transcript is out of range, set the start position to 1 and add cut flag to transcript name
                                        if apply_reverse_complement:
                                            transcript_name += f"-cut3':{1-start_position}bp"
                                        else:
                                            transcript_name += f"-cut5':{1-start_position}bp"
                                        start_position = 1

                                    if end_position > sequence_length: # If only the end of the transcript is out of range, set the end position to sequence_length (maximum of range) and add cut flag to transcript name
                                        if apply_reverse_complement:
                                            transcript_name += f"-cut5':{end_position-sequence_length}bp"
                                        else:
                                            transcript_name += f"-cut3':{end_position-sequence_length}bp"
                                        end_position = sequence_length
                                        
                                coordinates_file.write(f"{strand_indicator} {start_position} {end_position} {transcript_name}\n") # Write the transcript header line to coordinates file

                                coordinates = []

                                ## Gather exon information
                                if 'Exon' in transcript:
                                    exons = transcript['Exon']
                                    for exon in exons:
                                        start = exon.get('start') - input_region_start + new_start # Calculate exon start coordinate relative to region
                                        #print(f'exon start: {start}')
                                        end = exon.get('end') - input_region_start + new_start # Calculate exon end coordinate relative to region
                                        #print(f'exon end: {end}')

                                        # If run with -nocut option:
                                        if nocut:
                                            if start < 1 and end < 1: # If whole exon is 5' out of region, ignore it
                                                continue
                                            if start > sequence_length and end > sequence_length: # If whole exon is 3' out of region, ignore it
                                                continue
                                            if start < 1: # If only the start of the exon is 5' out of region, set the start coordinate to 1
                                                start = 1
                                            if end > sequence_length: # If only the end of the exon is 3' out of region, set the end coordinate to sequence_length (maximum of range)
                                                end = sequence_length

                                        # Initialize separate UTR start and end
                                        utr_start = 0  
                                        utr_end = 0

                                        ## Examine UTRs to modify exon extents
                                        if 'UTR' in transcript:
                                            utrs = transcript['UTR']
                                            for utr in utrs:
                                                utr_start = utr.get('start') - input_region_start + new_start # Calculate UTR start coordinate relative to region
                                                #print(f'utr start: {utr_start}')
                                                utr_end = utr.get('end') - input_region_start + new_start # Calculate UTR end coordinate relative to region
                                                #print(f'utr end: {utr_end}')
                                                
                                                # If run with -nocut option:
                                                if nocut:
                                                    if utr_start < sequence_length and utr_end > sequence_length: # If UTR end is 3' out of region, set the end coordinate to sequence_length (maximum of range)
                                                        utr_end = sequence_length
                                                    
                                                if start == utr_start and end == utr_end: # If 'exon' overlaps entirely with annotated UTR, set 'exon' coordinates to 0 0 (to be removed later) - gives UTR priority
                                                    start = 0
                                                    end = 0
                                                elif utr_start <= start <= utr_end: # If exon start lies inside UTR, give UTR priority and set exon start coordinate to begin after UTR ends
                                                    start = utr_end + 1
                                                    #print(f'start_final {start}')
                                                elif utr_start <= end <= utr_end: # If exon end lies inside UTR, give UTR priority and set exon end coordinate to begin before UTR starts
                                                    end = utr_start - 1
                                                    #print(f'end_final {end}')
                                                
                                                # If run with -nocut option:
                                                if nocut:
                                                    if start > sequence_length and end > sequence_length: # If whole exon is 3' out of region, set exon coordinates to 0 0 (to be removed later)
                                                        start = 0
                                                        end = 0
                                                        continue
                                                    if end > sequence_length: # If only the end of the UTR is 3' out of region, set the end coordinate to sequence_length (maximum of range)
                                                        end = sequence_length
                                                        continue
                                                    if utr_start < 0 and utr_end < 0: # If whole UTR is 5' out of region, ignore it (DELETE?)
                                                        continue
                                                    if utr_start < 0: # If only the start of the UTR is 5' out of region, set the start coordinate to 1 (DELETE?)
                                                        utr_start = 1

                                        coordinates.append((f"{start} {end} exon", start)) # Append the exon lines to coordinates object
                                    
                                    ## Gather UTR information
                                    if 'UTR' in transcript:
                                        utrs = transcript['UTR']
                                        for utr in utrs:
                                            utr_start_abs = utr.get('start') # Get absolute coordinate of UTR start
                                            utr_end_abs = utr.get('end') # Get absolute coordinate of UTR end
                                            utr_start = utr_start_abs - input_region_start + new_start # Calculate UTR start coordinate relative to region
                                            #print(f'utr_start_abs {utr_start_abs}')
                                            #print(f'utr_start {utr_start}')
                                            utr_end = utr_end_abs - input_region_start + new_start # Calculate UTR end coordinate relative to region
                                            #print(f'utr_end_abs {utr_end_abs}')
                                            #print(f'utr_end {utr_end}')

                                            # If run with -nocut option:
                                            if nocut:
                                                if utr_start < 1 and utr_end < 1: # If whole UTR is 5' out of region, ignore it
                                                    continue
                                                if utr_start > sequence_length and utr_end > sequence_length: # If whole UTR is 3' out of region, ignore it
                                                    continue
                                                if utr_start < 1 and utr_end > 1: # If UTR start is 5' out of region, set the start coordinate to 1
                                                    utr_start = 1
                                                if utr_start < sequence_length and utr_end > sequence_length: # If UTR end is 3' out of region, set the end coordinate to sequence_length (maximum of range)
                                                    utr_end = sequence_length

                                            coordinates.append((f"{utr_start} {utr_end} UTR", utr_start)) # Append the UTR lines to coordinates object
                                
                                else:
                                    print("No exons found in the transcript.")

                                coordinates.sort(key=lambda x: x[1]) # Sort coordinate lines by start coordinate

                                # Write coordinate lines to coordinates file, except lines with 0 0 coordinates (flagged for removal earlier)
                                for coord, _ in coordinates:
                                    if coord.strip() != "0 0 exon" and "0 0 UTR":
                                        coordinates_file.write(f"{coord}\n")
                                coordinates_file.write("\n\n") # Add 2 blank lines after each transcript block

            # If run with -anno (save annotation coordinates):
            if coordinates_output_file:
                with open(coordinates_output_file, 'r') as coordinates_file:
                    coordinates_content = coordinates_file.read()

                # If run with -rev option, reverse complement coordinates
                if apply_reverse_complement:
                    reversed_coordinates = reverse_coordinates(coordinates_content, sequence_length) # Reverse the coordinates
                    # Save the reversed coordinates to the specified output file
                    with open(coordinates_output_file, 'w') as reversed_coordinates_file:
                        reversed_coordinates_file.write(reversed_coordinates)
                    print(f"Reversed coordinates saved to {coordinates_output_file}")

                # Otherwise, write original coordinates to file
                else:
                    with open(coordinates_output_file, 'w') as coordinates_file:
                        coordinates_file.write(coordinates_content)
                    print(f"Coordinates saved to {coordinates_output_file}")  

        else:
            print("No coordinates output file specified.") # Print to notify user neither annotation nor fasta file specified as outputs

    else:
        print("No genes found in the specified region.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download DNA sequences in FASTA format and gene annotation coordinates in VISTA format from Ensembl.")
    parser.add_argument("-s", "--species", help="Species name (e.g., 'Homo_sapiens' or 'Human')")
    parser.add_argument("-gene", "--gene_name", help="Gene name")
    parser.add_argument("-sa", "--start_adjust", type=int, default=0, help="Number to subtract from the start coordinate (default: 0)")
    parser.add_argument("-ea", "--end_adjust", type=int, default=0, help="Number to add to the end coordinate (default: 0)")
    parser.add_argument("-fasta", "--fasta_output_file", default=None, help="Output file name for the DNA sequence in VISTA format")
    parser.add_argument("-anno", "--coordinates_output_file", default=None, help="Output file name for the gene coordinates")
    parser.add_argument("-all", "--all_transcripts", action="store_true", help="Include all transcripts (instead of canonical transcript only)")
    parser.add_argument("-nocut", action="store_false", default=True, help="Delete annotations not included in sequence")
    parser.add_argument("-rev", action="store_true", help="Reverse complement DNA sequence and coordinates")
    parser.add_argument("-autoname", action="store_true", help="Automatically generate output file names based on species and gene name")
    args = parser.parse_args()

    # Automatically generate output file names if -autoname is provided
    if args.autoname:
        if not args.fasta_output_file:
            args.fasta_output_file = f"{args.species}_{args.gene_name}.fasta.txt"
        if not args.coordinates_output_file:
            args.coordinates_output_file = f"{args.species}_{args.gene_name}.annotation.txt"

    # Pass arguments in the correct order
    run_gene_coordinates(
        args.species, 
        args.gene_name, 
        args.start_adjust, 
        args.end_adjust, 
        args.fasta_output_file, 
        args.coordinates_output_file, 
        args.all_transcripts, 
        args.nocut, 
        args.rev
    )