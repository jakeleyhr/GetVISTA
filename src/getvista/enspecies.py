#!/usr/bin/env python

"""
File: enspecies.py
Author: Jake Leyhr
GitHub: https://github.com/jakeleyhr/GetVISTA/
Date: February 2024
Description: query Ensembl database with complete or partial species common name, \
    binomial name, or taxon to return binomial names (to use in engene or encoords modules)
"""
# Import dependencies
import sys
import argparse
import requests
from fuzzywuzzy import process
from getvista.version_check import check_for_updates

# Dictionary to store alternative search terms for common names
alternative_search_terms = {
    "C.elegans": "nematode",
    "Paramormyrops kingsleyae": "elephantfish",
    "Periophthalmus magnuspinnatus": "mudskipper",
    "ass": "donkey",
    "cattle": "cow",
    "hybrid cattle": "cow",
    "mummichog": "killifish"
}

def get_species_list():
    # Ensembl REST API endpoint for species list
    url = "https://rest.ensembl.org/info/species?content-type=application/json"
    try:
        # Sending GET request to the Ensembl API
        response = requests.get(url)

        # Checking if the request was successful (status code 200)
        if response.status_code == 200:
            species_data = response.json()

            # Extracting species information
            species_list = []
            for species_info in species_data["species"]:
                species_common_name = species_info["common_name"]
                species_latin_name = species_info["name"]
                species_list.append((species_common_name, species_latin_name))

            return species_list

        else:
            # If the request was not successful, print error message with status code
            print(f"Error: Unable to retrieve species list. Status code: {response.status_code}")

    except requests.exceptions.RequestException as e:
        # Handling request exceptions
        print("Error:", e)

    # If there's an issue with the response, return None
    return None


def common_name_search(common_name, species_list):
    search_keyword = common_name.lower()
    matched_species = []
    # Check if search keyword matches any common names directly
    matched_species += [(common_name, latin_name) for common_name, latin_name in species_list if search_keyword in common_name.lower()]
    # Check if search keyword matches any alternative search terms
    for common_name, alt_term in alternative_search_terms.items():
        if search_keyword.lower() in alt_term.lower():
            # Find the Latin name corresponding to the matching common name
            latin_name = next((latin_name for common, latin_name in species_list if common.lower() == common_name.lower()), None)
            # Append only the matching species to the matched_species list
            if latin_name:
                matched_species.append((common_name, latin_name))
    # Sort and print matched species
    if matched_species:
        matched_species.sort(key=lambda x: x[0].lower())
        max_common_name_length = max(len(common_name) for common_name, _ in matched_species)
        print("\nList of Matching Species in Ensembl:")
        for common_name, latin_name in matched_species:
            dashes = '-' * (max_common_name_length - len(common_name))
            print(f"{common_name} {dashes}-> {latin_name}")
    else:
        # If no exact match is found, perform a spellcheck-like search
        # Get a list of all common names
        all_common_names = [name for name, _ in species_list] + list(alternative_search_terms.values())
        # Use fuzzywuzzy to find the closest match
        closest_matches = process.extract(search_keyword, all_common_names, limit=3)
        for closest_match in closest_matches:
            closest_common_name = closest_match[0]
            # If the closest match is an alternative search term, find the corresponding common name
            if closest_common_name in alternative_search_terms.values():
                closest_common_name = next((common_name for common_name, alt_term in alternative_search_terms.items() if alt_term.lower() == closest_common_name), None)
            closest_latin_name = next((latin_name for common_name, latin_name in species_list if common_name == closest_common_name), None)
            matched_species.append((closest_common_name, closest_latin_name))
        max_common_name_length = max(len(common_name) for common_name, _ in matched_species)
        print("\nFinding up to 3 closest matches:")
        for common_name, latin_name in matched_species:
            dashes = '-' * (max_common_name_length - len(common_name))
            print(f"{common_name} {dashes}-> {latin_name}")


def latin_name_search(latin_name, species_list):
    search_keyword = latin_name.lower()
    # Filter species list based on search keyword if provided
    matched_species = [(common_name, latin_name) for common_name, latin_name in species_list if search_keyword in latin_name.lower()]
    # Sort and print matched species
    if matched_species:
        matched_species.sort(key=lambda x: x[0].lower())
        max_common_name_length = max(len(common_name) for common_name, _ in matched_species)
        print("\nList of Matching Species in Ensembl:")
        for common_name, latin_name in matched_species:
            dashes = '-' * (max_common_name_length - len(common_name))
            print(f"{common_name} {dashes}-> {latin_name}")
    else:  # For Latin name search
        closest_matches = process.extract(search_keyword, [ln for _, ln in species_list], limit=3)
        for closest_match in closest_matches:
            closest_latin_name = closest_match[0]
            closest_common_name = next((common_name for common_name, latin_name in species_list if latin_name == closest_latin_name), None)
            matched_species.append((closest_common_name, closest_latin_name))
            
        max_common_name_length = max(len(common_name) for common_name, _ in matched_species)
        print("\nFinding up to 3 closest matches:")
        for common_name, latin_name in matched_species:
            dashes = '-' * (max_common_name_length - len(common_name))
            print(f"{common_name} {dashes}-> {latin_name}")


def taxonomysearch(taxon, species_list):
    url = f"https://rest.ensembl.org/info/genomes/taxonomy/{taxon}?content-type=application/json"

    species_set = set()  # Set to store unique species names and Latin names
    try:
        # Make 10 requests to the API and combine species into set - because result can vary each time for some reason
        for _ in range(10):
            # Sending GET request to the Ensembl API
            response = requests.get(url)

            # Checking if the request was successful (status code 200)
            if response.status_code == 200:
                species_data = response.json()

                # Extracting species information
                for species_info in species_data:
                    species_common_name = species_info.get("display_name")
                    species_latin_name = species_info.get("name")
                    if species_common_name and (species_common_name, species_latin_name) not in species_set:
                        # Add species names to set
                        species_set.add((species_common_name, species_latin_name))

            #else:
                # If the request was not successful, print error message with status code
                #print(f"Error: Unable to retrieve taxonomic information. Status code: {response.status_code}")
                

    except requests.exceptions.RequestException as e:
        # Handling request exceptions
        print("Error:", e)
        
    # Print the unique species names collected
    if species_set:
        common_species=[]
        # Iterate over each pair (common_name, latin_name) in species_set
        for common_name1, latin_name1 in species_set:
            # Check if the Latin name is present in species_list
            for common_name2, latin_name2 in species_list:
                # If the Latin name is found in species_list, add the pair to common_species
                if latin_name1 == latin_name2:
                    formatted_string = common_name1.replace(" ", "_")
                    if latin_name1 == formatted_string.lower():
                        common_species.append((common_name2, latin_name1))
                    else:
                        common_species.append((common_name1, latin_name1))    

        max_common_name_length = max(len(common_name) for common_name, _ in common_species)
        print("\nList of Species in Ensembl:")
        for common_name, latin_name in sorted(common_species):
            # Calculate the number of dashes needed for alignment
            dashes = '-' * (max_common_name_length - len(common_name))
            # Print common name with aligned Latin name
            print(f"{common_name} {dashes}-> {latin_name}")
    else:
        print("\nNo species records retrieved. Check taxon name is correct\n")
        sys.exit()
    print(f"\n {len(common_species)} species records in {taxon}")


def getspecies(common_name, latin_name, taxon):
    if common_name:
        species_list = get_species_list()
        common_name_search(common_name, species_list)
    elif latin_name:
        species_list = get_species_list()
        latin_name_search(latin_name, species_list)
    elif taxon:
        species_list = get_species_list()
        taxonomysearch(taxon, species_list)
    else: 
        species_list = get_species_list()
        if species_list:
            # Find the maximum length of common names
            max_common_name_length = max(len(common_name) for common_name, _ in species_list)

            # Printing the list of species with common names and Latin species names
            print("List of Species in Ensembl:")
            for common_name, latin_name in species_list:
                # Calculate the number of dashes needed for alignment
                dashes = '-' * (max_common_name_length - len(common_name))
                # Print common name with aligned Latin name
                print(f"{common_name} {dashes}-> {latin_name}")
        else:
            print("No species records retrieved.")
    print("")


def main():
    #Check for updates
    check_for_updates()
    
    # Argument parser setup
    parser = argparse.ArgumentParser(description="List species names from the Ensembl genome browser by common name, bionomial latin name, or taxon.")

    parser.add_argument("-cn", "--common_name", help="Search for species by their common name or part of it (e.g. fly)")
    parser.add_argument("-ln", "--latin_name", help="Search for species by their Latin names or part of it (e.g. melano)")
    parser.add_argument("-tax", "--taxon", help="Search for species by taxon name (e.g. Carnivora)")


    # Parse command-line arguments
    args = parser.parse_args()

    getspecies(args.common_name,
               args.latin_name,
               args.taxon)


if __name__ == "__main__":
    main()
    
    