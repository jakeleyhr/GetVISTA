#!/usr/bin/env python

"""
File: gbvistagene.py
Author: Jake Leyhr
GitHub: https://github.com/jakeleyhr/GetVISTA/
Date: January 2024
Description: Query the GenBank database with a species and gene name \
    to obtain FASTA file and gene feature coordinates in pipmaker format
"""

# Import dependencies
import re
import argparse
import sys
from collections import defaultdict
import http.client
from Bio import Entrez, SeqIO
from Bio.Seq import Seq


http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = "HTTP/1.0"

# functions for debugging
def record_directories(gene_info):
    for key, value in gene_info[0].items():
        print(f"Key at top level: {key}")
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                print(f"  Sub-key: {sub_key}")
                if isinstance(sub_value, dict):
                    for sub_sub_key, sub_sub_value in sub_value.items():
                        print(f"    Sub-sub-key: {sub_sub_key}")
                        print(f"    Sub-sub-value: {sub_sub_value}")
                else:
                    print(f"  Sub-value: {sub_value}")
def search_key_value(data, target_key, indent=""):
    for key, value in data.items():
        print(key)
        if key == target_key:
            print(f"{indent}{key}: {value}")
            continue
        if isinstance(value, dict):
            search_key_value(value, target_key, indent + "  ")
def explore_structure(data, indent=""):
    if isinstance(data, dict):
        for key, value in data.items():
            print(f"{indent}{key}:")
            explore_structure(value, indent + "  ")
    elif isinstance(data, list):
        for i, item in enumerate(data):
            print(f"{indent}[{i}]:")
            explore_structure(item, indent + "  ")
    else:
        print(f"{indent}{data}")


# Function #1 - get gene record
def search_gene_info(species, gene_name):
    # Set your email address
    Entrez.email = "dummy@gmail.com"

    # Build the query - strict check on species name and gene name (also searches gene name synonyms)
    query = f"{species}[ORGN] AND {gene_name}[Gene Name] OR {gene_name}[Accession]"

    # Search Entrez Gene
    handle = Entrez.esearch(db="gene", term=query, retmode="xml")
    record = Entrez.read(handle)
    # print(f' Entrez record: {record}')

    # Fetch gene information from each record
    if "IdList" in record and record["IdList"]:
        gene_id = record["IdList"][0]
        # print(gene_id)
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        gene_record = Entrez.read(handle)
        # print(f'Entrez gene record: {gene_record}')
        return gene_record

    return None


# Function #2 - process the gene record and extract relevant information
def process_gene_info(gene_info, record_id, start_adjust, end_adjust, species, gene_symbol):
    # Extract the relevant information
    if gene_info:
        gene_ref_name = gene_info[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
        print(f"Query gene: {gene_ref_name}")
        try:
            gene_ref_desc = gene_info[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_desc"]
            print(f"Description: {gene_ref_desc}")
        except KeyError:
            print("Description: None available")
        try:
            synonyms = gene_info[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_syn"]
            print(f"Synonyms: {synonyms}")
        except KeyError:
            print("Synonyms: None available")
        try:
            locus = gene_info[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_maploc"]
            print(f"Locus: {locus}")
        except KeyError:
            print("Locus: None available")

        try:
            assembly = gene_info[0]["Entrezgene_locus"][record_id]["Gene-commentary_label"]
            accession_number = gene_info[0]["Entrezgene_locus"][record_id]["Gene-commentary_accession"]
            start = (int(
                gene_info[0]["Entrezgene_locus"][record_id]["Gene-commentary_seqs"][0]
                ["Seq-loc_int"]["Seq-interval"]["Seq-interval_from"]
            ) + 1)
            end = (int(
                    gene_info[0]["Entrezgene_locus"][record_id]["Gene-commentary_seqs"][0]
                    ["Seq-loc_int"]["Seq-interval"]["Seq-interval_to"]
            ) + 1)  # "Gene-commentary_seqs" gives merged gene seqeunce, "Gene-commentary-products" gives transcripts
            length = end - start + 1

            print("")
            print(f"Using record {record_id}:")
            print(f"Assembly: {assembly}")
            print(f"Accession: {accession_number}")
            print(f"Location: {start}:{end}")
            print(f"Length: {length}bp")
            print("")
        except IndexError:
            print(f"Record_ID #{record_id} not found. Try a different value")
            sys.exit()  # Exit the script

        # For debugging - explore the file format:
        # Assuming gene_info is a dictionary, explore the format
        #record_directories(gene_info)
        # Assuming gene_info is a list with a single dictionary, explore the format
        #search_key_value(gene_info[0], 'Seq-interval_from')
        # Print the structure of gene_info[0]
        #explore_structure(gene_info[0]['Entrezgene_locus'])

    else:
        print(f"No gene information found for {gene_symbol} in {species}.")

    # Calculate region start and end based on the gene start and end coordinates +/- the user-provided adjustment values
    requested_start_position = start - start_adjust
    requested_end_position = end + end_adjust


    return accession_number, requested_start_position, requested_end_position


# Function #3 - get list of genes and features in specified region
def get_genes_in_region(accession, requested_start_position, requested_end_position):
    # Retrieve GenBank record
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Check if start and end coordinates are within range of sequence record
    for feature in record.features:
        if feature.type == "source":
            source_str = str(feature.location)
            match = re.search(r"\[([<>]?\d+):([<>]?\d+[<>]?)\]\([+-]\)", source_str) # Look for location coordinates in particular format. '[<>]?' allows for < or >
            if match:
                source_start = int(match.group(1).lstrip('<').lstrip('>')) + 1 # +1 to start only because of 0-based indexing in Entrez record (not present on NCBI website)
                source_end = int(match.group(2).lstrip('<').lstrip('>'))
    
    start = requested_start_position
    end = requested_end_position
    if requested_start_position < source_start:
        start = source_start
        print("Input start coordinate is out of bounds, trimming to closest value:")
    if requested_end_position > source_end:
        end = source_end
        print("Input end coordinate is out of bounds, trimming to closest value:")

    # Calculate total user-specified sequence length in bp (need to add 1)
    sequence_length = end - start + 1
    
    start_position = start
    end_position = end

    print(f"Specified region: {accession}:{start}-{end}")
    print(f"Specified region length: {sequence_length}bp")
    print("")
    
    print("Finding features in region...")
    print("")


    # Prepare to collect genes and features
    genes_in_region, collected_features = [], []
    collect_features = False

    # Parse GenBank features to identify genes that overlap with specified sequence region
    for feature in record.features:
        if feature.type == "gene":
            location_str = str(feature.location)
            match = re.search(r"\[([<>]?\d+):([<>]?\d+[<>]?)\]\([+-]\)", location_str) # Look for location coordinates in particular format. '[<>]?' allows for < or >
            if match:
                gene_start = int(match.group(1).lstrip('<').lstrip('>')) + 1 # +1 to start only because of 0-based indexing in Entrez record (not present on NCBI website)
                gene_end = int(match.group(2).lstrip('<').lstrip('>'))
                if (
                    start <= gene_start <= end # Gene start inside region?
                    or start <= gene_end <= end # Gene end inside region?
                    or gene_start <= start <= end <= gene_end # Gene middle inside region?
                ):
                    # print(feature)
                    # print(f'gene start: {gene_start}, gene end: {gene_end}')
                    # print(f'region start: {start}, region end: {end}')
                    genes_in_region.append(feature.qualifiers["gene"][0]) # Collect feature
                    collect_features = True # Start collecting subsequent features from list
                    continue
            else:
                print("Location coordinates in unexpected format")

        if collect_features:
            if (
                (feature.type == "ncRNA"
                and "N" in feature.qualifiers.get("transcript_id", [""])[0] # N designates curated records as opposed to X (e.g. NR vs XR)
                )
                or (
                    feature.type == "mRNA"
                    and "N" in feature.qualifiers.get("transcript_id", [""])[0] # N designates curated records as opposed to X (e.g. NM vs XM)
                )
                or (
                    feature.type == "CDS"
                    and "N" in feature.qualifiers.get("protein_id", [""])[0] # N designates curated records as opposed to X (e.g. NP vs XP)
                )  
            ):
                gene_value = feature.qualifiers.get("gene", [""])[0]  # extract the gene name associated with the feature
                print(f"gene: {gene_value}")
                if feature.type == "mRNA":
                    transcript = feature.qualifiers.get("transcript_id", [""])[0] # extract the transcript_id associated with the feature
                    print(f"mRNA found: {transcript}")
                if feature.type == "CDS":
                    transcript = feature.qualifiers.get("protein_id", [""])[0]
                    print(f"CDS found: {transcript}")
                if feature.type == "ncRNA":
                    transcript = feature.qualifiers.get("transcript_id", [""])[0]
                    print(f"ncRNA found: {transcript}")

                # Remove extraneous characters from location and reorder the coordinates
                simplified_location = reorder_location(re.sub(r"[^0-9,:]", "", str(feature.location)))

                # Get the feature's strand direction
                strand = re.search(r"[+-]", str(feature.location)).group()  

                # Create dictionary format
                feature_dict = {
                    "gene": gene_value,
                    "type": feature.type,
                    "transcriptid": transcript,
                    "location": simplified_location,
                    "strand": strand,
                }
                # Add features to dictionary
                collected_features.append(feature_dict)

            # Stop collecting features when the next "gene" feature is found
            if feature.type == "gene":
                collect_features = False

        # print("Collected features:")
        # print(type(collected_features))
        # for collected_feature in collected_features:
        #    print(collected_feature)
        #    print("")

    return genes_in_region, collected_features, start_position, end_position, sequence_length


# Function #3.5 - simplify and reorder gene feature coordinates
def reorder_location(location_str):
    # Split the location string into features
    features = location_str.split(",")

    # Split each feature into pair of numbers (start and end coordinates) and convert to integers
    pairs = [tuple(map(int, feature.split(":"))) for feature in features]

    # Sort the pairs based on the leftmost value (start coordinate)
    ordered_pairs = sorted(pairs, key=lambda x: x[0])

    # Format the ordered pairs back into the desired string format e.g. "100:200, 300:400"
    ordered_location_str = ",".join([f"{int(left)+1}:{right}" for left, right in ordered_pairs]) # +1 to left(start) to account for 0-based numbering from Entrez record

    return ordered_location_str


# Function 4 - reformat feature information into lists groups by gene
def reformat(collected_features):
    # Initialize empty lists
    mrna_list, cds_list, ncrna_list = [], [], []
    # Loop through features and extract corresponding lists
    for feature in collected_features:
        if feature["type"] == "mRNA":
            mrna_list.extend(extract_features(feature))
        elif feature["type"] == "CDS":
            cds_list.extend(extract_features(feature))
        elif feature["type"] == "ncRNA":
            ncrna_list.extend(extract_features(feature))

    # Generate ID number for mRNA transcript based on gene and transcriptID
    gene_mapping = {}
    transcriptid_mapping = {}
    current_id = 1
    for mrna in mrna_list:
        gene = mrna['gene']
        transcriptid = mrna['transcriptid']
        key = (gene, transcriptid) # Use both gene name and transcriptid as the key
        firstkey = (gene)
        if firstkey not in transcriptid_mapping and key not in gene_mapping:
            # If both the gene and the combiantion of gene and transcriptid haven't appeared earlier, reset ID count to 1
            current_id = 1
            gene_mapping[key] = current_id
            transcriptid_mapping[firstkey] = current_id
        elif firstkey in transcriptid_mapping and key not in gene_mapping:
            # If the gene has appeared before but not the combination of gene and transcriptid, add 1 to ID count
            current_id = current_id + 1
            gene_mapping[key] = current_id
        mrna['ID']=gene_mapping[key]

    # Generate ID number for CDS transcript based on gene and transcriptID
    gene_mapping = {}
    transcriptid_mapping = {}
    current_id = 1
    for cds in cds_list:
        gene = cds['gene']
        transcriptid = cds['transcriptid']
        key = (gene, transcriptid) # Use both gene name and transcriptid as the key
        firstkey = (gene)
        if firstkey not in transcriptid_mapping and key not in gene_mapping:
            # If both the gene and the combiantion of gene and transcriptid haven't appeared earlier, reset ID count to 1
            current_id = 1
            gene_mapping[key] = current_id
            transcriptid_mapping[firstkey] = current_id
        elif firstkey in transcriptid_mapping and key not in gene_mapping:
            # If the gene has appeared before but not the combination of gene and transcriptid, add 1 to ID count
            current_id = current_id + 1
            gene_mapping[key] = current_id
        cds['ID']=gene_mapping[key]

    #Now make the gene name + ID as the "transcript key" for pairing mRNAs and corresponding CDSs

    # Create a dictionary to organize entries by transcript and type
    organized_dict = defaultdict(lambda: {"mRNA": [], "CDS": []})
    # Iterate through the mRNA list and organize entries based on transcriptid
    for mrna_entry in mrna_list:
        transcript_key = f"{mrna_entry['gene']}transcript{mrna_entry['ID']}"
        organized_dict[transcript_key]["mRNA"].append(mrna_entry)
    # Iterate through the CDS list and pair entries with the same transcriptid
    for cds_entry in cds_list:
        transcript_key = f"{cds_entry['gene']}transcript{cds_entry['ID']}"
        organized_dict[transcript_key]["CDS"].append(cds_entry)
    # Convert the values of the organized_dict to a list - the final output with mRNAs paired to CDSs
    paired_list = list(organized_dict.values())


    # Create a mapping of unique ncRNA transcriptids to numbers
    transcriptid_mapping_transcript = defaultdict(lambda: len(transcriptid_mapping_transcript) + 1)
    
    # Create a dictionary to organize ncRNA entries by transcript and type
    organized_ncrnas = defaultdict(lambda: {"ncRNA": []})

    # Iterate through the ncRNA list and pair entries with the same transcriptid
    for ncrna_entry in ncrna_list:
        transcriptid = ncrna_entry["transcriptid"]
        transcript_key = f"{ncrna_entry['gene']}transcript{transcriptid_mapping_transcript[transcriptid]}"
        organized_ncrnas[transcript_key]["ncRNA"].append(ncrna_entry)

    # Convert the values of the organized_dict to a list - the final output of ncRNAs
    ncrna_list = list(organized_ncrnas.values())

    return ncrna_list, paired_list


# Function 4.5 - extract feature coordinates 
def extract_features(feature):
    # Extract the 'location' string from the feature dictionary
    location_str = feature["location"]
    # Split the string into pairs
    pairs = location_str.split(",")
    # Loop through the pairs and extract the leftmost and rightmost values (start and end coords)
    feature_list = [
        {
            "gene": feature["gene"],
            "type": feature["type"],
            "strand": feature["strand"],
            "transcriptid": feature["transcriptid"],
            "start": int(pair.split(":")[0]),
            "end": int(pair.split(":")[1]),
        }
        for pair in pairs
    ]
    return feature_list


# Function 5 - reformat gene feature information into pipmaker format
def pipmaker(paired_list, ncrna_list, start_position):
    # Prepare for writing results:
    result_text = []

    # First, deal with protein-coding genes
    for transcript in paired_list:
        if "mRNA" in transcript:
            mrnas = transcript["mRNA"]
            if mrnas:
                first_mrna = mrnas[0]  # Assuming the first mRNA entry represents the transcript
                gene_name = first_mrna["gene"]
                transcript_id = first_mrna["transcriptid"]
                strand_indicator = (">" if first_mrna["strand"] == "+" else "<")  # Get strand direction < or >
                start = (min(mrna["start"] for mrna in mrnas) - start_position + 1)  # Add 1 to avoid 0 values
                end = (max(mrna["end"] for mrna in mrnas) - start_position + 1)  # Add 1 to avoid 0 values

                result_text.append(f"{strand_indicator} {start} {end} {gene_name}:{transcript_id}") # Assemble header line

                # Prepare to write feature lines
                feature_lines = []

                # Make UTR feature lines
                if "mRNA" in transcript:
                    utrs = transcript["mRNA"]
                    for utr in utrs:
                        utr_start = (utr["start"] - start_position + 1)  # Add 1 to avoid 0 values
                        utr_end = (utr["end"] - start_position + 1)  # Add 1 to avoid 0 values
                        feature_lines.append(f"{utr_start} {utr_end} UTR")

                # Make exon feature lines
                if "CDS" in transcript:
                    exons = transcript["CDS"]
                    for exon in exons:
                        exon_start = (exon["start"] - start_position + 1)  # Add 1 to avoid 0 values
                        exon_end = (exon["end"] - start_position + 1)  # Add 1 to avoid 0 values
                        feature_lines.append(f"{exon_start} {exon_end} exon")

                # Sort feature lines under each header
                feature_lines.sort(key=lambda x: (int(x.split()[0]), int(x.split()[1])))

                # Add feature lines
                result_text.extend(feature_lines)

    # Second, deal with ncRNA genes
    for transcript in ncrna_list:
        if "ncRNA" in transcript:
            ncrnas = transcript["ncRNA"]
            if ncrnas:
                first_ncrna = ncrnas[0]  # Assuming the first mRNA entry represents the transcript
                gene_name = first_ncrna["gene"]
                transcript_id = first_ncrna["transcriptid"]
                strand_indicator = (">" if first_ncrna["strand"] == "+" else "<")  # Get strand direction < or >
                start = (min(ncrna["start"] for ncrna in ncrnas) - start_position + 1)  # Add 1 to avoid 0 values
                end = (max(ncrna["end"] for ncrna in ncrnas) - start_position + 1)  # Add 1 to avoid 0 values

                result_text.append(f"{strand_indicator} {start} {end} {gene_name}:{transcript_id}") # Assemble header line

                # Prepare to write feature lines
                feature_lines = []

                # Make UTR feature lines
                if "ncRNA" in transcript:
                    ncrnas = transcript["ncRNA"]
                    for ncrna in ncrnas:
                        ncrna_start = (ncrna["start"] - start_position + 1)  # Add 1 to avoid 0 values
                        ncrna_end = (ncrna["end"] - start_position + 1)  # Add 1 to avoid 0 values
                        feature_lines.append(f"{ncrna_start} {ncrna_end} UTR")
                
                # Sort feature lines under each header
                feature_lines.sort(key=lambda x: (int(x.split()[0]), int(x.split()[1])))

                # Add feature lines
                result_text.extend(feature_lines)

    # Separate exons and UTRs by iterating through lines
    i = 0
    while i < len(result_text) - 1:
        current_line = result_text[i]
        next_line = result_text[i + 1]

        # Skip lines with '<' or '>'
        if (
            "<" in current_line
            or ">" in current_line
            or "<" in next_line
            or ">" in next_line
        ):
            i += 1
            continue

        # Split the lines into values (columns)
        current_values = current_line.split()
        next_values = next_line.split()

        # Check conditions and modify the lines
        # If the UTR and exon coordinates are the same:
        if current_values[0] == next_values[0] and current_values[1] == next_values[1]:
            result_text.pop(i) # Remove the top line from result_text
            #print("Redundant UTR removed")
        else:
            # If UTR overlaps with 3' exon, cut back the UTR
            if (
                int(current_values[1]) > int(next_values[0])
                and current_values[2] == "UTR"
                and next_values[2] == "exon"
            ):
                current_values[1] = str(int(next_values[0]) - 1)
                result_text[i] = " ".join(current_values)
            # If exon overlaps with 3' UTR, cut forward the UTR
            if (
                int(current_values[1]) > int(next_values[0])
                and current_values[2] == "exon"
                and next_values[2] == "UTR"
            ):
                next_values[0] = str(int(current_values[1]) + 1)
                result_text[i + 1] = " ".join(next_values)

            #print(f"start: {int(current_values[0])}, end: {int(current_values[1])}")
            #print({int(current_values[0]) == int(current_values[1]) + 1})

            # Last bodge fix
            if int(current_values[0]) == int(current_values[1]) + 1:
                result_text.pop(i) # Remove the line from result_text
                print("error UTR removed")

            # Move to the next pair of lines
            i += 1

    # print(result_text)
    coordinates = result_text

    return coordinates


# Function A - cut/remove out-of-range sequence features from the pipmaker file (unless -nocut argument included)
def cut(coordinates, sequence_length):
    processed_coordinates = []

    for line in coordinates:
        if line.startswith(">"):
            header_values = list(map(int, line.split()[1:3]))
            if header_values[0] < 1:
                line_parts = line.split()
                line_parts[3] = f"{line_parts[3]}:cut5':{1 - header_values[0]}bp"  # Update name with cut flag
                header_values[0] = 1  # Set to 1 if less than 1
                line = " ".join(map(str, line_parts))
            if header_values[1] > sequence_length:
                line_parts = line.split()
                line_parts[3] = f"{line_parts[3]}:cut3':{header_values[1]-sequence_length}bp"  # Update name with cut flag
                header_values[1] = min(sequence_length, header_values[1])  # Set to sequence_length if greater than sequence_length
                line = " ".join(map(str, line_parts))
            processed_line = f"> {header_values[0]} {header_values[1]} {line.split()[3]}"
            processed_coordinates.append(processed_line)
        elif line.startswith("<"):
            header_values = list(map(int, line.split()[1:3]))
            if header_values[0] < 1:
                line_parts = line.split()
                line_parts[3] = f"{line_parts[3]}:cut5':{1 - header_values[0]}bp"  # Update name with cut flag
                header_values[0] = 1  # Set to 1 if less than 1
                line = " ".join(map(str, line_parts))
            if header_values[1] > sequence_length:
                line_parts = line.split()
                line_parts[3] = f"{line_parts[3]}:cut3':{header_values[1]-sequence_length}bp"  # Update name with cut flag
                header_values[1] = min(sequence_length, header_values[1])  # Set to sequence_length if greater than sequence_length
                line = " ".join(map(str, line_parts))
            processed_line = (f"< {header_values[0]} {header_values[1]} {line.split()[3]}")
            processed_coordinates.append(processed_line)
        else:
            header_values = list(map(int, line.split()[0:2]))
            header_values[0] = max(1, header_values[0])  # Set to 1 if less than 1
            header_values[0] = min(sequence_length, header_values[0])  # Set to sequence_length if greater than sequence_length
            header_values[1] = max(1, header_values[1])  # Set to 1 if less than 1
            header_values[1] = min(sequence_length, header_values[1])  # Set to sequence_length if greater than sequence_length
            processed_line = f"{header_values[0]} {header_values[1]} {line.split()[2]}"
            if (header_values[0] == 1 and header_values[1] == 1) or (
                header_values[0] == sequence_length
                and header_values[1] == sequence_length
            ):
                continue
            else:
                processed_coordinates.append(processed_line)

    return processed_coordinates


# Function B - reverse the coordinates in the pipmaker file (if -rev argument included)
def reverse_coordinates(coordinates, sequence_length):
    # Prepare objects to append data
    reversed_coordinates, transcript_lines = [], [] 

    # Go though line by line and reverse the coordinates
    for line in coordinates.split('\n'):
        if line:
            fields = line.split()
            # For lines with 3 fields (without strand information i.e. not header lines):
            if len(fields) == 3: 
                start = int(fields[0])
                end = int(fields[1])
                feature = fields[2]

                # Swap start and end values, subtract from sequence_length to flip coordinates
                new_start = sequence_length - end + 1 #+1 makes sense
                new_end = sequence_length - start + 1 #+1 makes sense

                # Append new coordinates to transcript_lines object
                transcript_lines.append((new_start, f"{new_start} {new_end} {feature}"))
            
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

                # Append new coordinates to transcript_lines object
                transcript_lines.append((new_start, f"{new_strand} {new_start} {new_end} {feature}")) 

        elif transcript_lines:
            sorted_lines = sorted(transcript_lines, key=lambda x: x[0]) # Sort lines with 3 fields based on the new start coordinate
            reversed_coordinates.extend(line for _, line in sorted_lines)
            reversed_coordinates.append("")  # Add an empty line between transcripts
            reversed_coordinates.append("")  # Add an empty line between transcripts
            transcript_lines = []  # Reset for the next transcript

    return '\n'.join(reversed_coordinates)


# Function C - download DNA sequence region in FASTA format (if -fasta argument included)
def download_fasta(species, accession, start, end, fasta_output_file, apply_reverse_complement):
    print("Getting FASTA sequence...")
    # Retrieve GenBank record
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Extract sequence from GenBank record based on coordinates
    sequence = record.seq[start:end]

    if not apply_reverse_complement:
        # Create a SeqRecord object and save it in FASTA format
        fasta_record = SeqIO.SeqRecord(sequence, id=f"{species}:{accession}_{start}:{end}:1", description="")
    else:
        # Reverse complement the sequence
        reverse_complement_sequence = str(sequence.reverse_complement())
        # Create a SeqRecord object and save it in FASTA format
        fasta_record = SeqIO.SeqRecord(
            Seq(reverse_complement_sequence),
            id=f"{species}:{accession}_{start}:{end}:-1",
            description="",
        )

    # Write to file
    SeqIO.write(fasta_record, f"{fasta_output_file}", "fasta")


def gbgene(
    species,
    gene_symbol,
    record_id,
    start_adjust,
    end_adjust,
    fasta_output_file,
    coordinates_output_file,
    nocut=None,
    apply_reverse_complement=False,
    autoname=False,
):
    # Get target gene record from Entrez
    gene_info = search_gene_info(species, gene_symbol)

    # Extract key information about the sequence region to be processed
    accession_number, requested_start_position, requested_end_position = process_gene_info(
        gene_info, record_id, start_adjust, end_adjust, species, gene_symbol)

    print("Finding features in region...")

    # Get a list of genes and their features included in the sequence region
    genes, collected_features, start_position, end_position, sequence_length = get_genes_in_region(accession_number, requested_start_position, requested_end_position)

    print("")
    print("Genes in the specified region:", genes)
    print("")

    # If -anno argument is included (or -autoname), collect feature coordinates and write to .txt file:
    if coordinates_output_file or autoname:
        # Get 2 reformatted lists of genes and features: 1 for ncRNAs, and 1 for paired mRNA/CDS features
        ncrna_list, paired_list = reformat(collected_features)

        # print("Features in the specified region:")
        # print("List of Protein-coding gene features:", paired_list)
        # print("List of ncRNAs:", ncrna_list)

        # Convert the lists of features into the final pipmaker format
        coordinates = pipmaker(paired_list, ncrna_list, start_position)

        # If -nocut argument is included, continue. If not, cut the ends to remove coordinates out of range.
        if nocut is False:
            coordinates = cut(coordinates, sequence_length)

        # Format the coordinates into string with lines
        formatted_coordinates = ""
        for i, line in enumerate(coordinates):
            if i > 0 and (line.startswith("<") or line.startswith(">")):
                formatted_coordinates += "\n\n"  # Add 2 blank lines before the header line
            parts = line.split()
            if (line.startswith("<") or line.startswith(">")):
                formatted_coordinates += f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}\n"
            else:
                formatted_coordinates += f"{parts[0]} {parts[1]} {parts[2]}\n"

        # Automatically generate output file names if -autoname is provided
        if autoname:
            if not fasta_output_file:
                if not apply_reverse_complement:
                    fasta_output_file = f"{species}_{gene_symbol}_{start_position}-{end_position}.fasta.txt"
                else:
                    fasta_output_file = (f"{species}_{gene_symbol}_{start_position}-{end_position}_revcomp.fasta.txt")
            if not coordinates_output_file:
                if not apply_reverse_complement:
                    coordinates_output_file = (f"{species}_{gene_symbol}_{start_position}-{end_position}.annotation.txt")
                else:
                    coordinates_output_file = (f"{species}_{gene_symbol}_{start_position}-{end_position}_revcomp.annotation.txt")

        # Check if ".txt" is already at the end of the coordinates_output_file argument
        if coordinates_output_file and not coordinates_output_file.endswith(".txt"):
            coordinates_output_file += ".txt"

        # Check if ".txt" is already at the end of the fasta_output_file argument
        if fasta_output_file and not fasta_output_file.endswith(".txt"):
            fasta_output_file += ".txt"

        # Check if need to reverse complement and write to txt file
        if not apply_reverse_complement:
            with open(coordinates_output_file, 'w') as coordinates_file:
                coordinates_file.write(formatted_coordinates)
            print(f"Coordinates saved to {coordinates_output_file}") 
            print("")
        else:
            reversed_coordinates = reverse_coordinates(formatted_coordinates, sequence_length) # Reverse the coordinates
            with open(coordinates_output_file, 'w') as coordinates_file:
                coordinates_file.write(reversed_coordinates)
            print(f"Reversed feature coordinates saved to {coordinates_output_file}")
            print("")
                
    else:
        print("No pipmaker output file specified.")

    # If -fasta argument is included (or -autoname), write the DNA sequence to a txt file:
    if fasta_output_file or autoname:
        # First download the sequence
        download_fasta(species, accession_number, start_position, end_position, fasta_output_file, apply_reverse_complement) 
        # Check if need to reverse complement
        if not apply_reverse_complement:
            print(f"DNA sequence saved to {fasta_output_file}")
        else:
            print(f"Reverse complement DNA sequence saved to {fasta_output_file}")
    else:
        print("No FASTA output file specified.")


def main():
    # Create an ArgumentParser
    parser = argparse.ArgumentParser(description="Query the GenBank database with a species and gene name \
                                     to obtain FASTA file and gene feature coordinates in pipmaker format.")

    # Add arguments for species and gene_symbol
    parser.add_argument("-s", "--species", help="Species name", required=True)
    parser.add_argument("-g", "--gene_symbol", help="Gene symbol", required=True)
    parser.add_argument("-r", "--record_id", type=int, default=0, help="Record ID number (default=0, the top match)")
    parser.add_argument("-sa", "--start_adjust", type=int, default=0, help="Number to subtract from the start coordinate (default: 0)")
    parser.add_argument("-ea", "--end_adjust", type=int, default=0, help="Number to add to the end coordinate (default: 0)")
    parser.add_argument("-fasta", "--fasta_output_file", default=None, help="Output file name for the DNA sequence in VISTA format")
    parser.add_argument("-anno", "--coordinates_output_file", default=None, help="Output file name for the gene coordinates")
    parser.add_argument("-nocut", action="store_true", default=False, help="Delete annotations not included in sequence")
    parser.add_argument("-rev", action="store_true", default=False, help="Reverse complement DNA sequence and coordinates")
    parser.add_argument("-autoname", action="store_true", default=False, help="Automatically generate output file names based on accession and gene name")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Provide arguments to run
    gbgene(
        args.species,
        args.gene_symbol,
        args.record_id,
        args.start_adjust,
        args.end_adjust,
        args.fasta_output_file,
        args.coordinates_output_file,
        args.nocut,
        args.rev,
        args.autoname
    )


if __name__ == "__main__":
    main()

