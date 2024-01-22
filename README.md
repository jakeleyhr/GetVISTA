# GetVISTA
Query Ensembl to obtain genomic information in VISTA format.
* **vistacoords.py**: query with species and _genomic coordinates_
* **vistagene.py**: query with species and _gene name_

## Author
Jake Leyhr (@jakeleyhr)

## Dependencies
* Python 3.11
* requests
* argparse

## Quick start guide
* Install and open [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
* Create an environment with python 3.11 e.g:
```
conda create -n getvistaenv python=3.11
```
* Activate (enter) the environment:
```
conda activate getvistaenv
```
* install the package dependencies:
```
pip install requests argparse
```
* Download the vistacoords.py and/or vistagene.py files
* In the Miniconda terminal, navigate to the folder containing the .py files
* Then you're ready to start!

# vistacoords.py usage

```
$ python vistacoords.py -h
usage: vistacoords.py [-h] [-s SPECIES] [-r REGION] [-fasta FASTA_OUTPUT_FILE] [-anno COORDINATES_OUTPUT_FILE] [-all] [-nocut] [-rev]
                      [-autoname]

Download DNA sequences in FASTA format and gene annotation coordinates in pipmaker format from Ensembl.

options:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name (e.g., 'Homo_sapiens' or 'Human')
  -r REGION, --region REGION
                        Genomic coordinates (e.g., 1:1000-2000)
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in FASTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -nocut                Don't delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and gene name
```
## vistacoords.py inputs and outputs:
The simplest inputs are the species name (**-s**) and region coordinates (**-r**), along with the -autoname flag:
```
$ python vistacoords.py -s human -r 1:100000-200000 -autoname
```
This produces the following output in the terminal:
```
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to human_1:100000-200000.annotation.txt
DNA sequence saved to human_1:100000-200000.fasta.txt
Total sequence length 100001
```
Along with two text files - the first contains the coordinates of the exons and UTRs of all genes contained within the genomic region selected in pipmaker format, and the second contains the DNA sequence of the selected region in FASTA format. By using the **-autoname** flag, the names of these files were automatically generated from the species and region inputs.

Alternatively, the output file names can be specified manually using the **-anno** and **-fasta** arguments, e.g:
```
$ python vistacoords.py -s human -r 1:100000-200000 -anno annotationoutput.txt -fasta fastaoutput.txt
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to annotationoutput.txt
DNA sequence saved to fastaoutput.txt
Total sequence length 100001
```
Without **-anno**, **-fasta**, or **-autoname** arguments, the terminal output will be provided but no output .txt files. If, for example, only **-anno** is provided, **-autoname** can also be provided to generate the remaining (fasta) filename:
```
$ python vistacoords.py -s human -r 1:100000-200000 -anno annotationoutput.txt -autoname             
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to annotationoutput.txt
DNA sequence saved to human_1:100000-200000.fasta.txt
Total sequence length 100001
```
## vistacoords.py specific arguments:
By default, only the exon and UTR coordinates of the canonical gene transcripts are included in the annotation .txt file, e.g:
```
$ python vistacoords.py -s human -r 1:950000-1000000 -autoname
```
```
...

< 25199 32094 PERM1-202
25199 26172 UTR
26173 26270 exon
26500 26625 exon
28882 31030 exon
31031 31174 UTR
32066 32094 UTR
```
However, by including the **-all** flag, all transcripts are included:
```
$ python vistacoords.py -s human -r 1:950000-1000000 -autoname -all
```
```
...

< 25199 32094 PERM1-202
25199 26172 UTR
26173 26270 exon
26500 26625 exon
28882 31030 exon
31031 31174 UTR
32066 32094 UTR


< 25205 32094 PERM1-201
25205 26172 UTR
26173 26270 exon
26500 26625 exon
28882 30658 exon
31138 31167 exon
31168 31174 UTR
32066 32094 UTR


< 25206 26642 PERM1-203
25206 26270 exon
26500 26642 exon


< 26122 32118 PERM1-204
26122 26172 UTR
26173 26270 exon
26500 26625 exon
28882 31030 exon
31031 31048 UTR
32066 32118 UTR
```
Also by default, the script carefully trims the transcript coordinates to ensure that the reported coordinates fit inside the specified region. For example, the mouse Cenpa gene is located on chromosome 5:30824121-30832175. If those coordinates are input, the resulting annotation file appears like this:
```
$ python vistacoords.py -s mouse -r 5:30824121-30832175 -autoname
```
```
> 1 8054 Cenpa-205
1 252 UTR
253 337 exon
5688 5794 exon
6206 6283 exon
6517 6651 exon
6652 6667 UTR
7262 8054 UTR
```
If the region 5:30824621-30832074 is specified instead, which cuts off 500bp from the 5' end and 100bp from the 3' end, The resulting annotation file is adjusted to include only bases inside the region. In this case, the 5' UTR and 1st exon have been deleted entirely, and the 3' UTR has been truncated. Information about the truncation has been added to the transcript name (Cenpa-205-cut5':500bp-cut3':100bp) to make it clear to the user that the selection has cut off part of the gene.
```
$ python vistacoords.py -s mouse -r 5:30824621-30832074 -autoname
```
```
> 1 7454 Cenpa-205-cut5':500bp-cut3':100bp
5188 5294 exon
5706 5783 exon
6017 6151 exon
6152 6167 UTR
6762 7454 UTR
```
This option can be turned off by including the **-nocut** flag, such that cut off parts of the gene are still included in the annotation file, with negative coordinates or coordinates that extend beyond the end the sequence:
```
$ python vistacoords.py -s mouse -r 5:30824621-30832074 -autoname -nocut
```
```
> -499 7554 Cenpa-205
-499 -248 UTR
-247 -163 exon
5188 5294 exon
5706 5783 exon
6017 6151 exon
6152 6167 UTR
6762 7554 UTR
```
By default, the specified genomic region is read on the forward strand, but for some purposes a gene on the reverse strand may want to be collected in the 5'>3' direction. In such cases, the **-rev** flag can be included. This reverse complements the DNA sequence returned in the fasta file (in addition to modifying the header to reflect this by changing :1 to :-1). It also flips the annotation coordinates. Returning to the mouse Cenpa gene as an example, this is the output when extracting the whole gene with **-rev**:
```
$ python vistacoords.py -s mouse -r 5:30824121-30832174 -autoname -rev
```
```
< 1 8054 Cenpa-205
1 793 UTR
1388 1403 UTR
1404 1538 exon
1772 1849 exon
2261 2367 exon
7718 7802 exon
7803 8054 UTR
```
Note that the strand direction indicator has changed (> to <), and the 252bp 5' UTR is now at the bottom (3' end) of the file rather than the top, with the rest of the annotations following suit.


# vistagene.py usage
```
$ python vistagene.py -h
usage: vistagene.py [-h] [-s SPECIES] [-gene GENE_NAME] [-sa START_ADJUST] [-ea END_ADJUST] [-fasta FASTA_OUTPUT_FILE] [-anno COORDINATES_OUTPUT_FILE] [-all] [-nocut] [-rev] [-autoname]

Download DNA sequences in FASTA format and gene annotation coordinates in pipmaker format from Ensembl.

options:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name (e.g., 'Homo_sapiens' or 'Human')
  -gene GENE_NAME, --gene_name GENE_NAME
                        Gene name
  -sa START_ADJUST, --start_adjust START_ADJUST
                        Number to subtract from the start coordinate (default: 0)
  -ea END_ADJUST, --end_adjust END_ADJUST
                        Number to add to the end coordinate (default: 0)
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in FASTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -nocut                Delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and gene name
```
## vistagene.py inputs:
The output arguments, in addition to the **-all**, **-nocut**, **-rev** arguments are identical to vistacoords.py described above, but the inputs are quite different. Rather than defining a species and genomic region, a species and gene name are input. For example, mouse and the gdf5 gene. This script outputs a detailed log of the gene information and the sequence region extracted:
```
$ python vistagene.py -s mouse -gene gdf5 -autoname 
Assembly name: GRCm39
mouse gdf5 coordinates: 2:155782943-155787287
mouse gdf5 is on -1 strand
mouse gdf5 sequence length: 4345bp
Extracting coordinates: 2:155782943-155787287
DNA sequence saved to mouse_gdf5.fasta.txt
Total sequence length: 4345bp
Coordinates saved to mouse_gdf5.annotation.txt
```
Two additional arguments can be used to adjust the start (**-sa**) and end (**-ea**) coordinates beyond the gene start and end. For example, to extract the sequence and annotations for the gdf5 gene plus an additional 50,000bp from the 5' flank and an additional 20,000bp from the 3' flank (direction relative to the assembly forward strand):
```
$ python vistagene.py -s mouse -gene gdf5 -autoname -sa 50000 -ea 20000 
Assembly name: GRCm39
mouse gdf5 coordinates: 2:155782943-155787287
mouse gdf5 is on -1 strand
mouse gdf5 sequence length: 4345bp
Extracting coordinates: 2:155732943-155807287
DNA sequence saved to mouse_gdf5.fasta.txt
Total sequence length: 74345bp
Coordinates saved to mouse_gdf5.annotation.txt
```
In the output, note that the gene length is 4,345bp, but the total sequence length extracted is 74,345bp as a result of the 70,000bp flanking regions also being included. The annotation file also reflects these additional sequences, including the genes in the expanded region:
```
< 1 39288 Uqcc1-204-cut5':44102bp
10276 10333 exon
18331 18403 exon
19335 19442 exon
20719 20814 exon
30613 30705 exon
38991 39014 exon
39015 39288 UTR


> 49706 51728 Gm15557-201
49706 49816 UTR
50460 50770 UTR
50771 51475 exon
51476 51728 UTR


< 50001 54345 Gdf5-201
50001 50520 UTR
50521 51395 exon
53421 54033 exon
54034 54345 UTR


> 65536 74345 Cep250-204-cut3':33533bp
65536 65635 UTR
70850 70949 UTR
70950 71132 exon
71878 71934 exon
72691 72773 exon
73088 73253 exon
73967 74073 exon
74291 74345 exon
```
## Notes
* Per the [Ensembl REST API documentation](https://rest.ensembl.org/documentation/info/overlap_region), the maximum sequence length that can be queried is 5Mb. Requests above this limit will fail (Status code: 400 Reason: Bad Request).
* For species with common names more than one word long (e.g. Alpine marmot or Spotted gar, as opposed to human or mouse), the full species name according to Ensembl must be used with underscores separating the words. For the Alpine marmot: marmota_marmota_marmota, and for the Spotted Gar: lepisosteus_oculatus

## Bugs

Please submit via the [GitHub issues page](https://github.com/jakeleyhr/GetVISTA/issues).  

## Software Licence

[GPLv3](https://github.com/jakeleyhr/GetVISTA/blob/main/LICENSE)
