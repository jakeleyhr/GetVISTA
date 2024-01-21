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

## vistacoords.py usage

```
$ vistacoords.py -h
usage: vistacoords.py [-h] [-s SPECIES] [-r REGION] [-fasta FASTA_OUTPUT_FILE] [-anno COORDINATES_OUTPUT_FILE] [-all] [-start_value START_VALUE] [-nocut] [-rev]
                      [-autoname]

Download DNA sequences in FASTA format and gene annotation coordinates in VISTA format from Ensembl.

options:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name (e.g., 'Homo_sapiens' or 'Human')
  -r REGION, --region REGION
                        Genomic coordinates (e.g., 1:1000-2000)
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in VISTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -start_value START_VALUE
                        Start value for coordinates, 1 by default
  -nocut                Don't delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and gene name
```
### vistacoords.py inputs and outputs:
The simplest inputs are the species name (**-s**) and region coordinates (**-r**), along with the -autoname flag:
```
$ vistacoords.py -s human -r 1:100000-200000 -autoname
```
This produces the following output in the terminal:
```
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to human_1:100000-200000.annotation.txt
DNA sequence saved to human_1:100000-200000.fasta.txt
Total sequence length 100001
```
Along with two text files - the first contains the coordinates of the exons and UTRs of all genes contained within the genomic region selected, and the second contains the DNA sequence of the selected region in FASTA format. By using the **-autoname** flag, the names of these files were automatically generated from the species and region inputs.

Alternatively, the output file names can be specified manually using the **-anno** and **-fasta** arguments, e.g:
```
$ vistacoords.py -s human -r 1:100000-200000 -anno annotationoutput.txt -fasta fastaoutput.txt
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to annotationoutput.txt
DNA sequence saved to fastaoutput.txt
Total sequence length 100001
```
Without **-anno**, **-fasta**, or **-autoname** arguments, the terminal output will be provided but no output .txt files. If, for example, only **-anno** is provided, **-autoname** can also be provided to generate the remaining (fasta) filename:
```
$ vistacoords.py -s human -r 1:100000-200000 -anno annotationoutput.txt -autoname             
Extracting human coordinates: 1:100000-200000
Assembly name: GRCh38
Coordinates saved to annotationoutput.txt
DNA sequence saved to human_1:100000-200000.fasta.txt
Total sequence length 100001
```
### vistacoords.py specific arguments:
By default, only the exon and UTR coordinates of the canonical gene transcripts are included in the annotation .txt file, e.g:
```
python vistacoords.py -s human -r 1:950000-1000000 -anno -autoname
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
python vistacoords.py -s human -r 1:950000-1000000 -anno -autoname -all
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
vistacoords.py -s mouse -r 5:30824121-30832175 -anno annotationoutput2.txt
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
If the region 5:30824621-30832074 is specified instead, which cuts off 500bp from the 5' end and 100bp from the 3' end, The resulting annotation file is adjusted to include only bp inside the region. In this case, the 5' UTR and 1st exon have been deleted entirely, and the 3' UTR has been truncated. Information about the truncation have been added to the transcript name (Cenpa-205-cut5':500bp-cut3':100bp) to make it clear to the user that the selection has cut off part of the gene.
```
$ vistacoords.py -s mouse -r 5:30824621-30832074 -anno annotationoutput2.txt
```
```
> 1 7454 Cenpa-205-cut5':500bp-cut3':100bp
5188 5294 exon
5706 5783 exon
6017 6151 exon
6152 6167 UTR
6762 7454 UTR
```
This option can be turned off by including the -nocut flag, such that cut off parts of the gene are still included in the annotation file, with negative coordinates or coordinates that extend beyond the end the sequence:
```
vistacoords.py -s mouse -r 5:30824621-30832074 -anno annotationoutput2.txt -nocut
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
start_value, rev

## vistagene.py usage

```
$ vistagene.py -h
usage: vistagene.py [-h] [-s SPECIES] [-gene GENE_NAME] [-sa START_ADJUST] [-ea END_ADJUST] [-fasta FASTA_OUTPUT_FILE] [-anno COORDINATES_OUTPUT_FILE] [-all]
                    [-start_value START_VALUE] [-nocut] [-rev] [-autoname]

Download DNA sequences in FASTA format and gene annotation coordinates in VISTA format from Ensembl.

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
                        Output file name for the DNA sequence in VISTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -start_value START_VALUE
                        Start value for coordinates, 1 by default
  -nocut                Delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and gene name
```



## Bugs

Please submit via the [GitHub issues page](https://github.com/jakeleyhr/GetVISTA/issues).  

## Software Licence

[GPLv3](https://github.com/jakeleyhr/GetVISTA/blob/main/LICENSE)
