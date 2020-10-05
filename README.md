# in-vivo-mRNA-Display
scripts for processing of in vivo mRNA display sequencing data

"In vivo mRNA display anables large-scale proteomics by next generation sequencing"
Panos Oikonomou, Roberto Salatino, Saeed Tavazoie
doi: 10.1073/pnas.2002650117
PNAS, October 2020

- command: python run_IMVD_alignment_yORF.py
- requirements: python 2.7, cutadapt v1.8.3, bowtie2 v2.2.6, samtools v1.2, bedtools v2.17.0
- sequencing files *.fastq.gz should be placed in the folder data/ in the same directory as the script.
- exemplar fastq.gz are provided with this repo
- if new sequencing files are to be analyzed they should be appropriately listed in run_IMVD_alignment_yORF.py, e.g.:  in_prefixes = ["ind2_CGATGTAT", "ind6_GCCAATAT", "ind8_CTTGTAAT"]; 
- sequencing files should follow the naming scheme R1_index_ind*_NNNNNNNN.fastq.gz, where NNNNNNNN is the illumina index used that should be present in the first column of data/barcodes/barcode_map.tsv
- all internal indexes listed in data/barcodes/myBarcodes.fa should be present in the second column of data/barcodes/barcode_map.tsv
- the output of this script is a count matrix file: new_strict_counts.tsv
