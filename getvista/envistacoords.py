#!/usr/bin/env python

"""
File: envistacoords.py
Author: Jake Leyhr
GitHub: https://github.com/jakeleyhr/GetVISTA/
Date: January 2024
Description: Query the Ensembl database with species and genomic coordinates to obtain FASTA file and gene feature coordinates in pipmaker format
"""

 # Import dependencies
import sys
import json
import time
import argparse
import requests
import re
from Bio import SeqIO
from Bio.Seq import Seq
from getvista.version_check import check_for_updates

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

    def get_genes_in_region(self, species, region):
        genes = self.perform_rest_action(
            endpoint='/overlap/region/{0}/{1}'.format(species, region),
            params={'feature': 'gene'}
        )
        return genes

# Function 1 - get genomic coordinates and FASTA sequence
def download_dna_sequence(species, genomic_coordinates):
    # Ensembl REST API endpoint for DNA sequences in FASTA format
    url = f"https://rest.ensembl.org/sequence/region/{species}/{genomic_coordinates}?format=fasta" 
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
    
    return fasta_lines

# Function 2 - get gene feature coordinates in pipmaker format
def pipmaker(genes, genomic_coordinates, apply_reverse_complement, nocut, all_transcripts):
    client = EnsemblRestClient()
    
    input_region_start = int(genomic_coordinates.split(":")[1].split("-")[0])
    input_region_end = int(genomic_coordinates.split(":")[1].split("-")[1])
    sequence_length = (input_region_end - input_region_start) + 1

    print(f"Specified sequence length: {sequence_length}bp")
    print("")
    print(f"Transcripts included in region:")

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
                    print(transcript_name)

                    if start_position < 0 and end_position < 0:
                        print(f"({transcript_name} transcript out of 5' range:{start_position}:{end_position})")
                    if start_position > sequence_length and end_position > sequence_length:
                        print(f"({transcript_name} transcript out of 3' range:{start_position}:{end_position})")
                        

                    if nocut == False:
                        if not (start_position < 0 and end_position < 0) or (start_position > sequence_length and end_position > sequence_length):
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
                                    start = 0
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
                                            utr_start = 0

                            coordinates.append((f"{start} {end} exon", start))

                        if 'UTR' in transcript:
                            utrs = transcript['UTR']
                            for utr in utrs:
                                utr_start_abs = utr.get('start')
                                utr_end_abs = utr.get('end')
                                utr_start = utr_start_abs - input_region_start + new_start
                                utr_end = utr_end_abs - input_region_start + new_start

                                if nocut == False:
                                    if utr_start < 0 and utr_end < 0:
                                        continue
                                    if utr_start > sequence_length and utr_end > sequence_length:
                                        continue
                                    if utr_start < 0 and utr_end > 0:
                                        utr_start = 0
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


def encoords(species, genomic_coordinates, fasta_output_file=None, coordinates_output_file=None, all_transcripts=None, nocut=None, apply_reverse_complement=False, autoname=False):
    #Check if query sequence is >5Mb
    coordssplit = re.match(r"(\d+):(\d+)-(\d+)", genomic_coordinates)
    start = int(coordssplit.group(2))
    end = int(coordssplit.group(3))
    if (end - start + 1) > 5000000:
        print("ERROR: Query sequence must be under 5Mb")
        sys.exit()

    # Get genes info
    client = EnsemblRestClient()
    genes = client.get_genes_in_region(species, genomic_coordinates)

    print(f"Assembly name: {genes[0]['assembly_name']}")

    # Get coordinates and sequence length
    coordinates_content, sequence_length = pipmaker(genes, genomic_coordinates, apply_reverse_complement, nocut, all_transcripts)

    # Automatically generate output file names if -autoname is provided
    chrom = coordssplit.group(1)
    if autoname:
        if not fasta_output_file:
            if not apply_reverse_complement:
                fasta_output_file = f"{species}_{chrom}_{start}-{end}.fasta.txt"
            else:
                fasta_output_file = (f"{species}_{chrom}_{start}-{end}_revcomp.fasta.txt")
        if not coordinates_output_file:
            if not apply_reverse_complement:
                coordinates_output_file = (f"{species}_{chrom}_{start}-{end}.annotation.txt")
            else:
                coordinates_output_file = (f"{species}_{chrom}_{start}-{end}_revcomp.annotation.txt")

    # Check if ".txt" is already at the end of the coordinates_output_file argument
    if coordinates_output_file and not coordinates_output_file.endswith(".txt"):
        coordinates_output_file += ".txt"

    # Check if ".txt" is already at the end of the fasta_output_file argument
    if fasta_output_file and not fasta_output_file.endswith(".txt"):
        fasta_output_file += ".txt"

    print("")
    # If run with -anno (save annotation coordinates):
    if coordinates_output_file:
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
        print("No coordinates output file specified.") # Print to notify user no annotation file specified as output


    fasta_lines = download_dna_sequence(species, genomic_coordinates)
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
            print(f"DNA sequence saved to {fasta_output_file}")
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
            print(f"Reverse complement DNA sequence saved to {fasta_output_file}")
        # Write to file
        SeqIO.write(fasta_record, f"{fasta_output_file}", "fasta")   
    else:
        print("No sequence output file specified") # Print to notify user neither sequence specified as output



def main():
    # Create an ArgumentParser
    parser = argparse.ArgumentParser(description="Query the Ensembl database with a species name and genomic coordinates to obtain DNA sequences in FASTA format and gene feature coordinates in pipmaker format.")
    
    # Add arguments for species, gene_symbol etc
    parser.add_argument("-s", "--species", required=True, help="Species name (e.g., 'Homo_sapiens' or 'Human')")
    parser.add_argument("-c", "--gencoordinates", required=True, help="Genomic coordinates (e.g., 1:1000-2000)")
    parser.add_argument("-fasta", "--fasta_output_file", default=None, help="Output file name for the DNA sequence in FASTA format")
    parser.add_argument("-anno", "--coordinates_output_file", default=None, help="Output file name for the gene coordinates in pipmaker format")
    parser.add_argument("-all", "--all_transcripts", action="store_true", default=False, help="Include all transcripts (instead of canonical transcript only)")
    parser.add_argument("-nocut", action="store_true", default=False, help="Don't delete annotations not included in sequence")
    parser.add_argument("-rev", action="store_true", default=False, help="Reverse complement DNA sequence and coordinates")
    parser.add_argument("-autoname", action="store_true", default=False, help="Automatically generate output file names based on species and gene name")
    
    # Parse the command-line arguments
    args = parser.parse_args() 

    # Pass arguments in the correct order
    encoords(
        args.species,
        args.gencoordinates,
        args.fasta_output_file,
        args.coordinates_output_file,
        args.all_transcripts,
        args.nocut,
        args.rev,
        args.autoname
    )

    
if __name__ == '__main__':
    check_for_updates()
    main()


