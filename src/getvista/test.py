# Example object
region_info = ['NC_000020', 'REGION:', '34383347..36454749']

# Find the element containing the ".." separator
separator_element = next((elem for elem in region_info if ".." in elem), None)

if separator_element:
    # Extract the numbers on either side of ".."
    start_number, end_number = map(int, separator_element.split('..'))

    print("Start number:", start_number)
    print("End number:", end_number)
else:
    print("Separator '..' not found in the list.")