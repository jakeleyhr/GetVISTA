# GetVISTA
Query Ensembl or GenBank to obtain genomic information in VISTA format. Useful for collecting sequences to run in [mVISTA](https://genome.lbl.gov/vista/mvista/submit.shtml) multi-species sequence conservation analyses.
* **envistagene.py**: query _Ensembl_ database with species and _gene name_
* **gbvistagene.py**: query _GenBank_ database with species and _gene name_
* **envistacoords.py**: query _Ensembl_ database with species and _genomic coordinates_
* **gbvistacoords.py**: query _GenBank_ database with species and _genomic coordinates_
* **gbgenerecord.py**: query _GenBank_ with database species and _gene name_, get list of records to select (to choose in gbgene module)
* **emailaddress.py**: check and update the email address used to make GenBank Entrez queries
* **version_check.py**: check package version is up to date

<img src="https://github.com/jakeleyhr/GetVISTA/assets/154226340/41d2c750-5241-4962-9c00-8a7b973bb1ef" width="800">



## Author
Jake Leyhr ([@jakeleyhr](https://twitter.com/JakeLeyhr)) \
https://github.com/jakeleyhr/GetVISTA

## Dependencies
* Python 3.11
* requests
* packaging
* biopython
* configparser

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
* Install the package:
```
pip install getvista
```
* Navigate to the folder you want to deposit the output file in:
```
cd \path\to\working\directory
```
* Then you're ready to start!









&nbsp;
&nbsp;
# engene usage
```
$ engene -h
usage: engene [-h] -s SPECIES [SPECIES ...] -g GENE_NAME [-sa START_ADJUST] [-ea END_ADJUST] 
                            [-fasta FASTA_OUTPUT_FILE] [-anno COORDINATES_OUTPUT_FILE] 
                            [-goa] [-all] [-nocut] [-rev] [-autoname] [-fw] [-flank FLANK] 
                            [-vis]

Query the Ensembl database with a species and gene name to obtain DNA sequences in FASTA 
format and gene feature coordinates in pipmaker format.

options:
  -h, --help            show this help message and exit
  -s SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]
                        Species name(s) (e.g., 'Homo_sapiens' or 'Human')
  -g GENE_NAME, --gene_name GENE_NAME
                        Gene name (e.g. BRCA1 or brca1)
  -sa START_ADJUST, --start_adjust START_ADJUST
                        Number to subtract from the start coordinate (default: 0)
  -ea END_ADJUST, --end_adjust END_ADJUST
                        Number to add to the end coordinate (default: 0)
  -goa, --gene_oriented_adjustment
                        Make start/end adjustments follow gene orientation, not assembly
                        orientation
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in FASTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates in pipmaker format
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -nocut                Delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and 
                        gene name
  -fw                   Automatically orient the gene in the forward strand by reverse 
                        complementing if needed
  -flank FLANK          Select 2 genes to specify new range. 'in' to include the flanking 
                        genes, 'ex' to exclude them
  -vis                  Display graphical representation of the sequence in the terminal
```
## engene inputs:
## -s, -g, and -autoname
The simplest inputs are the species name (**-s**) and gene name (**-g**), along with the **-autoname** flag:
```
$ engene -s human -g gdf5 -autoname
```
This produces the following output in the terminal:
```
▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35433347-35454746
Specified sequence length: 21400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201

Coordinates saved to human_gdf5.annotation.txt
DNA sequence saved to human_gdf5.fasta.txt
```
The species name has to be entered in a form recognised by Ensembl, which includes binomial and one-word common names. For species with multi-word species names such as Spotted Gar, the binomial name must be used with an underscore i.e. Lepisosteus_oculatus. Multiple species names can be entered, not just one, as will be explained later in this document.

Two text files are generated in the working directory - the first contains the coordinates of the exons and UTRs of all genes contained within the genomic region selected in pipmaker format, and the second contains the DNA sequence of the selected region in fasta format. By using the **-autoname** flag, the names of these files were automatically generated from the species and gene name inputs.

## -anno and -fasta
Alternatively, the output file names can be specified manually using the **-anno** and **-fasta** arguments, e.g:
```
$ engene -s human -g gdf5 -fasta fastafilename -anno annotationfilename

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35433347-35454746
Specified sequence length: 21400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201

Coordinates saved to annotationfilename.txt
DNA sequence saved to fastafilename.txt
```

Without **-anno**, **-fasta**, or **-autoname** arguments, the terminal output will be provided but no output .txt files will be saved to the working directory. This can be useful for quickly testing a set of arguments without cluttering up your working directory with test output files.
```
$ engene -s human -g gdf5

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35433347-35454746
Specified sequence length: 21400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201

No coordinates output file specified
No sequence output file specified
```

If only **-anno** is provided, **-autoname** can also be provided to generate the remaining (fasta) filename:
```
$ engene -s human -g gdf5 -anno annotationfilename -autoname 

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35433347-35454746
Specified sequence length: 21400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201

Coordinates saved to annotationfilename.txt
DNA sequence saved to human_gdf5.fasta.txt
```
## -sa and -ea
Two additional arguments can be used to adjust the start (**-sa**) and end (**-ea**) coordinates beyond the gene start and end. For example, to extract the sequence and annotations for the human gdf5 gene plus an additional 50,000bp from the 5' flank and an additional 20,000bp from the 3' flank (strand direction relative to the assembly, NOT the gene):
```
$ engene -s human -g gdf5 -autoname -sa 50000 -ea 20000 

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35383347-35474746
Specified sequence length: 91400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201
CEP250-202
UQCC1-205

Coordinates saved to human_gdf5_20.35383347-35474746.annotation.txt
DNA sequence saved to human_gdf5_20.35383347-35474746.fasta.txt
```
In the output, note that the gene length is 21,400bp, but the total sequence length extracted is 91,400bp as a result of the 70,000bp flanking regions also being included. The automatically generated output file names also now include the new genomic coordinates. The annotation file also reflects these additional sequences, including the gene transcripts in the expanded region:
```
< 1 28627 UQCC1-205-cut5':80769bp
692 787 exon
10746 10850 exon
28594 28617 exon
28618 28627 UTR

> 49683 52104 GDF5-AS1-201
49683 49817 exon
50506 52104 exon

< 50001 54882 GDF5-201
50001 50562 UTR
50563 51437 exon
53952 54582 exon
54583 54882 UTR

< 70608 70751 MIR1289-1-201
70608 70751 exon

> 72252 91400 CEP250-202-cut3':44534bp
72252 72405 UTR
74958 75028 UTR
76637 76759 UTR
78919 79021 UTR
79022 79207 exon
80229 80285 exon
82397 82479 exon
82693 82858 exon
83620 83726 exon
83958 84209 exon
86544 86640 exon
88704 88805 exon
89327 89485 exon
90028 90206 exon
90524 90706 exon
```
## -goa
By default, when using **-sa** or **-ea**, these adjustments are done according to the assembly direction, regardless of the gene orientation. For example **-sa 1000** adds 1000bp from the 'left' of the specified gene, even if this is 'downstream' because when the gene is on the reverse strand. However, if you want to collect 1000bp from 'in front' of the gene (e.g. for collecting promoter sequences), you can use **-sa 1000** in combination with **-goa**. By adding the **-goa** argument to the command, the module ensures the 1000bp are always added 'upstream' of the gene (and vice-versa when using **-ea**).

## -nocut
By default, the module carefully trims the transcript coordinates to ensure that the reported coordinates fit inside the specified region and only features inside the region are included. In this case, the specified coordinate range cut off 80,769bp from the UQCC1 gene-205 transcript, and 44,534bp has been cut off the CEP250-205, and this information about the truncation has been added to the transcript names to make it clear to the user that the selection has cut off part of these transcripts.
This option can be turned off by including the **-nocut** flag, such that cut-off parts of the transcripts are still included in the annotation file, with negative coordinates or coordinates that extend beyond the end the sequence:
```
$ engene -s human -g gdf5 -autoname -sa 50000 -ea 20000 -nocut
```
produces the following text in the annotation file:
```
< -80768 28627 UQCC1-205
-80768 -79412 UTR
-79411 -79277 exon
-76680 -76567 exon
-68658 -68581 exon
-36182 -36074 exon
-16789 -16732 exon
-9162 -9090 exon
-1428 -1321 exon
692 787 exon
10746 10850 exon
28594 28617 exon
28618 28627 UTR

> 49683 52104 GDF5-AS1-201
49683 49817 exon
50506 52104 exon

< 50001 54882 GDF5-201
50001 50562 UTR
50563 51437 exon
53952 54582 exon
54583 54882 UTR

< 70608 70751 MIR1289-1-201
70608 70751 exon

> 72252 135934 CEP250-202
72252 72405 UTR
74958 75028 UTR
76637 76759 UTR
78919 79021 UTR
79022 79207 exon
80229 80285 exon
82397 82479 exon
82693 82858 exon
83620 83726 exon
83958 84209 exon
86544 86640 exon
88704 88805 exon
89327 89485 exon
90028 90206 exon
90524 90706 exon
92156 92300 exon
93103 93249 exon
94525 94755 exon
95885 96078 exon
96300 96427 exon
96630 96799 exon
107291 107458 exon
107866 108000 exon
110083 110226 exon
111178 111311 exon
113231 113369 exon
114373 114721 exon
115249 115370 exon
116703 116823 exon
118499 118620 exon
119044 121659 exon
124392 124505 exon
124689 124844 exon
125597 125698 exon
126652 126708 exon
128017 128280 exon
128281 135934 UTR
```
Bear in mind that this file won't work as an input to mVISTA, because the annotation coordinates have to match the fasta file and having negative coordinates or coordinates greater than the total sequence length will result in an error. However, this option can be useful to see the number of features that are cut off in the default output. It's worth noting that the results using the **-nocut** option match the output you get from using the Ensembl website's export data tool, which then often has to be manually curated to fit the sequence in order for the mVISTA analysis to run successfully. In this regard, the default mode of this module (without **-nocut**) represents another incentive to use the GetVISTA package as opposed to Ensembl website interface.

## -vis
When the **-vis** flag is included in the command, a simple graphical representation of the specified sequence region is printed in the terminal:
```
$ engene -s human -g gdf5 -sa 50000 -ea 20000 -vis

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35383347-35474746
Specified sequence length: 91400bp

Transcripts included in region:
UQCC1-205
GDF5-AS1-201
GDF5-201
MIR1289-1-201
CEP250-202

Graphical representation of specified sequence region:
<UQCC1-205-cut5':80769bp<========================================================
=================================================================================
=================================================================================
=================================================================================
=>GDF5-AS1-201#GDF5-201<=========================================================
=================================================================================
=================================================================================
====<MIR1289-1-201<=====================>CEP250-202-cut3':44534bp>

No coordinates output file specified
No FASTA output file specified
```
Here, the transcripts in the region are printed with their strand directions indicated by the flanking '>' or '<' symbols while the '=' symbols represent intergenic sequences. Partially overlapping transcripts, such as GDF5-AS-1 and GDF5-201, are shown by the presence of the '#' symbol between the two names. The lengths of the displayed intergenic sequences are approximately accurate, while genes are only represented by their printed names, regardless of size.

## -all
By default, only the exon and UTR coordinates of the canonical gene transcripts are included in the annotation .txt file, e.g:
```
< 1 4882 GDF5-201
1 562 UTR
563 1437 exon
3952 4582 exon
4583 4882 UTR
```
However, by including the **-all** flag in the command, all transcripts are included, including the non-canonical ones, e.g:
```
< 1 4882 GDF5-201
1 562 UTR
563 1437 exon
3952 4582 exon
4583 4882 UTR

< 1 21400 GDF5-202
1 562 UTR
563 1437 exon
3952 4582 exon
4583 4823 UTR
7886 8041 UTR
21294 21400 UTR
```
## -rev
Also by default, the specified genomic region is read on the genome assembly's forward strand, but for some purposes, a gene on the reverse strand may want to be collected in the 5'>3' direction (or indeed vice versa). In such cases, the **-rev** flag can be included. This reverse complements the DNA sequence returned in the fasta file (in addition to modifying the header to reflect this by changing :1 to :-1). It also flips the annotation coordinates. We have already seen that the human GDF5 gene is on the reverse strand (indicated by left-facing '<' in the annotation file), so using the **-rev** argument:
```
$ engene -s human -g gdf5 -autoname -rev              
```
Produces this output in the terminal:
```
▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35433347-35454746
Specified sequence length: 21400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201

Reversed coordinates saved to human_gdf5_revcomp.annotation.txt
Reverse complement DNA sequence saved to human_gdf5_revcomp.fasta.txt
```
The automatically generated output file names reflect the fact that the sequence and annotation coordinates have been reversed, the coordinates looking like this:
```
< 19297 21400 GDF5-AS1-201-cut3':318bp
19297 20895 exon

> 16519 21400 GDF5-201
16519 16818 UTR
16819 17449 exon
19964 20838 exon
20839 21400 UTR

> 650 793 MIR1289-1-201
650 793 exon
```
## -fw and multi-species entry
Instead of manually reversing of the strand direction, the **-fw** argument forces the use of the forward strand for the specified gene. For example, if you know you always want your target gene to be on the forward strand, you can add the **-fw** argument every time, because if the target gene is already on the forward strand, it won't reverse the direction, but if it's on the reverse strand, it will.

Crucially, multiple species names can be included as arguments with **-s**, for example:
```
$ engene -s human mouse chicken -g gdf5 -autoname -sa 50000 -ea 20000 -fw
```
In this case, the module behaves as previously described, just iterating through the different species names (separated by spaces). As a result, in this example the 6 output files from all 3 species will be saved in the working directory, and the output in the terminal looks like this:
```
▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ human gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCh38
human gdf5 coordinates: 20:35433347-35454746
human gdf5 is on reverse strand
human gdf5 sequence length: 21400bp
Specified coordinates: 20:35383347-35474746
Specified sequence length: 91400bp

Transcripts included in region:
GDF5-AS1-201
GDF5-201
MIR1289-1-201
CEP250-202
UQCC1-205

gdf5 is on the reverse strand, flipped automatically.

Reversed coordinates saved to human_gdf5_20.35383347-35474746_revcomp.annotation.txt
Reverse complement DNA sequence saved to human_gdf5_20.35383347-35474746_revcomp.fasta.txt

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ mouse gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: GRCm39
mouse gdf5 coordinates: 2:155782943-155787287
mouse gdf5 is on reverse strand
mouse gdf5 sequence length: 4345bp
Specified coordinates: 2:155732943-155807287
Specified sequence length: 74345bp

Transcripts included in region:
Uqcc1-204
Gm15557-201
Gdf5-201
Cep250-204

gdf5 is on the reverse strand, flipped automatically.

Reversed coordinates saved to mouse_gdf5_2.155732943-155807287_revcomp.annotation.txt
Reverse complement DNA sequence saved to mouse_gdf5_2.155732943-155807287_revcomp.fasta.txt

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒ chicken gdf5 ▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒

Assembly name: bGalGal1.mat.broiler.GRCg7b
chicken gdf5 coordinates: 20:1563813-1568758
chicken gdf5 is on forward strand
chicken gdf5 sequence length: 4946bp
Specified coordinates: 20:1513813-1588758
Specified sequence length: 74946bp

Transcripts included in region:
CEP250-203
GDF5-201
UQCC1-201
ERGIC3-201
ENSGALT00010061117
ENSGALT00010061118
ENSGALT00010061120

Coordinates saved to chicken_gdf5_20.1513813-1588758.annotation.txt
DNA sequence saved to chicken_gdf5_20.1513813-1588758.fasta.txt
```
The inclusion of the **-fw** argument caused the human and mouse, but not chicken, strand directions to be automatically flipped such that in all 3 species, Gdf5 is output in the same (forward) strand. This is particularly useful in conjunction with the input of multiple species, as the same gene will often be annotated on different strands in assemblies from different species, and for sequence alignment a common orientation is essential.
## -flank
Finally, the last option for specifying the sequence region is **-flank**. Adding this argument causes the module to pause after printing the list of genes in the specified sequence region, and prompts the user to input two gene names from this list. These two genes are then used to define the boundaries of a new sequence region, whereupon the module reruns to extract the feature coordinates and fasta sequence corresponding to this new region. **-flank** can be used in two modes: **-flank in** includes the two genes within the boundaries, such that the new sequence begins at the first base of the first gene and ends at the last base of the second gene (read left to right along the forward strand), while **-flank ex** excludes these genes, such that the new sequence begins immediately after the first gene, and ends immediately before the second gene. For example, imagine engene is run with Gene C as an input and **-sa 100000 -ea 100000** so the region includes Gene A, Gene B, Gene C, Gene D, and Gene E. When the **-flank in** argument is included and Gene B and Gene D are entered at the prompt, the new sequence includes just Genes B, C, and D, and the two intergenic sequences between them. On the other hand, when the **-flank ex** argument is used and Gene B and Gene D are entered at the prompt, the new sequence only includes Gene C and the intergenic sequences on either side of it.

&nbsp;
&nbsp;
# gbgene usage
```
$ gbgene -h
usage: gbgene [-h] -s SPECIES [SPECIES ...] -g GENE_NAME [-r RECORD_ID] [-sa START_ADJUST] 
                            [-ea END_ADJUST] [-goa] [-fasta FASTA_OUTPUT_FILE]
                            [-anno COORDINATES_OUTPUT_FILE] [-x] [-nocut] [-rev] 
                            [-autoname] [-fw] [-flank FLANK] [-vis]

Query the GenBank database with a species and gene name to obtain FASTA file and gene 
feature coordinates in pipmaker format.

options:
  -h, --help            show this help message and exit
  -s SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]
                        Species name(s) (e.g., 'Homo_sapiens' or 'Human')
  -g GENE_NAME, --gene_name GENE_NAME
                        Gene name (e.g. BRCA1 or brca1)
  -r RECORD_ID, --record_id RECORD_ID
                        Record ID number (default=0, the top match)
  -sa START_ADJUST, --start_adjust START_ADJUST
                        Number to subtract from the start coordinate 
                        (default: 0)
  -ea END_ADJUST, --end_adjust END_ADJUST
                        Number to add to the end coordinate (default: 0)
  -goa, --gene_oriented_adjustment
                        Make start/end adjustments follow gene orientation, not assembly
                        orientation
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in VISTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates
  -x                    Include predicted (not manually curated) transcripts in results
  -nocut                Delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on accession and 
                        gene name
  -fw                   Automatically orient the gene in the forward strand by reverse 
                        complementing if needed
  -flank FLANK          Select 2 genes to specify new range. 'in' to include the flanking 
                        genes, 'ex' to exclude them
  -vis                  Display graphical representation of the sequence in the terminal
```
This command functions almost identically to engene, except that it queries the GenBank nucleotide database rather than Ensembl. There is no **-all** option, as all transcripts are automatically included in the annotation file. However, this only includes the manually curated transcripts. To get all transcripts including the predicted ones, add the **-x** flag. This may be particularly relevant when exploring new genomes with few manually curated genes/transcripts. There is also an extra option **-r**, to specify the sequence record. By default it is 0 (the default record according to GenBank), but a different record may be desired in some cases (e.g. to use the human T2T assembly CHM13v2.0 instead of the GRCh38.14 assembly - see the gbrecords module description below). Unlike in engene, the species name can be entered in any form with underscores separating the words (e.g. carcharodon_carcharias or great_white_shark).









&nbsp;
&nbsp;
# encoords usage
```
$ encoords -h
usage: encoords.py [-h] -s SPECIES -c GENCOORDINATES [-fasta FASTA_OUTPUT_FILE]
                              [-anno COORDINATES_OUTPUT_FILE] [-all] [-nocut] 
                              [-rev] [-autoname] [-flank FLANK] [-vis]

Query the Ensembl database with a species name and genomic coordinates to obtain DNA 
sequences in FASTA format and gene feature coordinates in pipmaker format.

options:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name (e.g., 'Homo_sapiens' or 'Human')
  -c GENCOORDINATES, --gencoordinates GENCOORDINATES
                        Genomic coordinates (e.g., 1:1000-2000)
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in FASTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene coordinates in pipmaker format
  -all, --all_transcripts
                        Include all transcripts (instead of canonical transcript only)
  -nocut                Don't delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on species and 
                        genomic coordinates
  -flank FLANK          Select 2 genes to specify new range. 'in' to include the flanking 
                        genes, 'ex' to exclude them
  -vis                  Display graphical representation of the sequence in the terminal
```
The encoords module works similarly to engene, but with some critical differences as the intention is to query a set of genomic coordinates instead of a gene name. Therefore, the **-g** argument is replaced by **-c**, and the coordinates are entered in the format 'chromosome:start-end' e.g:
```
$ encoords -s mouse -c 5:30810000-30890000 -autoname

Assembly name: GRCm39
Specified sequence length: 80001bp

Transcripts included in region:
Slc35f6-201
Cenpa-205
Dpysl5-203

Coordinates saved to mouse_5_30810000-30890000.annotation.txt
DNA sequence saved to mouse_5_30810000-30890000.fasta.txt
```
As the coordinates are specified in place of the gene name, the **-sa**, **-ea**, and **-fw** arguments in engene are absent from encoords. Also, as the chromosomes and coordinates are not homologous between species, there is no option to enter multiple species names as in engene.


&nbsp;
&nbsp;
# gbcoords usage
```
$ gbcoords -h
usage: genecoords [-h] -a ACCESSION -c GENCOORDINATES [-fasta FASTA_OUTPUT_FILE]
                              [-anno COORDINATES_OUTPUT_FILE] [-x] [-nocut] [-rev] 
                              [-autoname] [-flank FLANK] [-vis]

Query the GenBank database with an accession and range of coordinates to obtain FASTA file 
and gene feature coordinates in pipmaker format.

options:
  -h, --help            show this help message and exit
  -a ACCESSION, --accession ACCESSION
                        accession code
  -c GENCOORDINATES, --gencoordinates GENCOORDINATES
                        Genomic coordinates
  -fasta FASTA_OUTPUT_FILE, --fasta_output_file FASTA_OUTPUT_FILE
                        Output file name for the DNA sequence in VISTA format
  -anno COORDINATES_OUTPUT_FILE, --coordinates_output_file COORDINATES_OUTPUT_FILE
                        Output file name for the gene annotation coordinates
  -x                    Include predicted (not manually curated) transcripts in results
  -nocut                Delete annotations not included in sequence
  -rev                  Reverse complement DNA sequence and coordinates
  -autoname             Automatically generate output file names based on accession and 
                        genomic coordinates
  -flank FLANK          Select 2 genes to specify new range. 'in' to include the flanking 
                        genes, 'ex' to exclude them
  -vis                  Display graphical representation of the sequence in the terminal
```
This command functions almost identically to encoords, except that it queries the GenBank nucleotide database rather than Ensembl so has aspects of gbcoords (the absence of the **-all** argument, replaced in a sense by **-x**). The other key difference with encoords is that an accession code (e.g. NC_000020 for human chromosome 20 in GRCh38.p14) must be specified with **-a** instead of a species name with **-s**, and the genomic coordinates therefore just require the base region, not the chromosome (e.g. 500000-600000 instead of 20:500000-600000).
```
$ gbcoords -a NC_000020 -c 500000-600000
```









&nbsp;
&nbsp;
# gbrecords usage
```
usage: gbrecord [-h] -s SPECIES -g GENE_NAME

Query the GenBank database with a species and gene name to obtain a list of different 
records containing the sequence to inform the use of the gbgene module.

options:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name (e.g., 'Homo_sapiens' or 'Human')
  -g GENE_NAME, --gene_name GENE_NAME
                        Gene name (e.g. BRCA1 or brca1)
```
This accessory module is intended to be used to see which genomic sequence records are available in GenBank for the given species and gene name. For example:
```
$ gbrecord -s human -g gdf5
```
gives this output in the terminal:
```
RECORD 0
Assembly: Chromosome 20 Reference GRCh38.p14 Primary Assembly
Accession: NC_000020
Location: 35433346:35454748
Length: 21403

RECORD 1
Assembly: RefSeqGene
Accession: NG_008076
Location: 21518:26399
Length: 4882

RECORD 2
Assembly: Chromosome 20 Alternate T2T-CHM13v2.0
Accession: NC_060944
Location: 37154207:37175639
Length: 21433
```
This shows that there are 3 different genomic records readily available in GenBank that contain the human GDF5 gene. Record 0 and Record 2 are different genomic assemblies, while Record 1 is a smaller RefSeqGene sequence (28kb). The transcript of GDF5 contained in this sequence is one of the shorter isoforms, only 4882bp as opposed to the longer 21403bp isoform in GRCh38.p14 or 21433bp isoform in CHM13v2.0. If you were only interested in the isolated region around the smaller core sequence of GDF5, you may want to use **-r 1** when running the gbgene command, as this would significantly speed up the request compared to using the default (**-r 0**) GRCh38.p14 record. For example, the command 
```
$ gbgene -s human -g gdf5 -autoname -r 0
```
takes ~68 seconds to complete, while the command:
```
$ gbgene -s human -g gdf5 -autoname -r 1
``` 
takes only ~8 seconds.



&nbsp;
&nbsp;
# gbemail usage
```
$ gbemail -h
usage: gbemail [-h] [-check] [-update]

Manage email address for GenBank Entrez queries. Only stored locally and sent with 
queries to NCBI, nowhere else.

options:
  -h, --help  show this help message and exit
  -check      Check the current saved email address
  -update     Update the email address
```
The first time you try to run gbgene, gbcoords, or gbrecord, you will be prompted for an email address as NCBI requires this to send queries via Entrez. After you enter this email address the first time, you won't need to enter it again, as it saved locally to a config file in the package directory. It is never shared with anyone except NCBI when you send queries via Entrez. This accessory module provides the option to check or update (change) this email address if you decide to use a different address at a later date. Any dummy email address or word (e.g. 'dummy@gmail.com' or 'dummy') can also be used if preferred. 


&nbsp;
&nbsp;
## Notes
* Per the [Ensembl REST API documentation](https://rest.ensembl.org/documentation/info/overlap_region), the maximum sequence length that can be queried with engene or encoords is 5Mb. Requests above this limit will fail (Status code: 400 Reason: Bad Request).
* Requests to GenBank with gbgene, gbcoords, or gbrecord sometimes fail for reasons unknown. If you get an "HTTP Error 400: Bad Request" when running gbcoords, gbgene, or gbrecord, try running the command once or twice again, and the query should go through.
* When running either gbcoords or gbgene on a large sequence record (e.g. a whole chromosome), it may take several tens of seconds to run, compared to the almost instant response from encoords and engene. This is because the gb scripts always search through the entire chromosomal record for gene features in the specified range, while the en scripts are able to narrow their search range to this range from the beginning.
* In the graphical representation, if two or more gene transcripts have identical start and end coordinates, only one of the transcripts will be visualised.

## Bugs
Please submit via the [GitHub issues page](https://github.com/jakeleyhr/GetVISTA/issues).  

## Software Licence
[GPLv3](https://github.com/jakeleyhr/GetVISTA/blob/main/LICENSE)
