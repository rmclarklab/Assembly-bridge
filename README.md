# Assembly-bridge

The assembler.py program links together genomic scaffolds from a draft genome assembly into superscaffolds using BLASR-aligned sequence contigs from one or more other assemblies (typically the latter are more fragmented assemblies, as might be generated by Illumina short-read *de novo* assemblies, but the program is agnostic with respect to assembly source).

Specifically, assembler.py was developed to order and orient Sanger scaffolds in the draft genome sequence of the two-spotted spider mite (*Tetranychus urticae*) using Illumina *de novo* assemblies from multiple *T. urticae* strains. The command line used to do this is given below, and the output was used to construct Table S3 in the publication that describes the study and method (data files used as input for the program can be downloaded from the supporting material). The manuscript is:

- Wybouw, N., Kosterlitz, O., Kurlovs, A. H., Bajda, S., Greenhalgh, R., Snoeck, S., Bui, H., Bryon, A., Dermauw, W., Van Leeuwen, T., and Clark, R. M. 2019. Long-term population studies uncover the genome structure and genetic basis of xenobiotic and host plant adaptation in the herbivore *Tetranychus urticae*. Under review (for a preprint, see bioRxiv 474064; doi: https://doi.org/10.1101/474064)

For the data sets needed to test installation, and to replicate Table S3 in Wybouw et al. (2019), please contact Richard Clark (richard.m.clark@utah.edu). The data sets used as input will be available for public download as soon as the manuscript is accepted and appears online.

Although the program was written for and tested on draft genome sequences from *T. urticae*, it was designed to work for related genome projects for which comparable data sets and input files are available.

The assembler.py program was written by Robert Greenhalgh (robert.greenhalgh@utah.edu). Please contact either Robert Greenhalgh or Richard Clark with questions about the software.

## Citation

For publications resulting from use of the software, we request that users cite Wybouw et al. 2019 (see above).

## Requirements

The assembler.py program is written in Python and requires pysam as a dependency. The progam was developed and tested with [Python](https://www.python.org/download/releases/2.7/) 2.7 and [pysam](https://pysam.readthedocs.io/en/latest/index.html) 0.14.1.

## Usage

### The script requires the following arguments to run:

- -b/--bams: One or more BAM files of BLASR-aligned de novo assemblies.

- -l/--length: An integer value specifying the minimum length a contig must be mapped to two scaffold ends in order to join them.

- -d/--distance: An integer value specifying the maximum distance to look for contigs around scaffold ends.

### The following optional arguments may be specified:

- -f/--fraction: The maximum fractional length of a scaffold to look for contigs at that scaffold's ends. If this value is less than that specified for -d/--distance for a scaffold, it will be used instead. Example: if -d/--distance is set to 50,000, -f/--fraction is set to 0.25, and the length of the scaffold being analyzed is 100,000 bp long, then the script will only analyze the 25,000 bp (0.25*100,000) around the scaffold ends.
    
- -r/--breakpoints: A tab-delimited text file used to break scaffolds into subscaffolds. Column 1 contains the name of the scaffold to split, while columns 2+ contain the basepair coordinates at which to split the scaffold. Each scaffold to be split should only be listed once in this file. If a scaffold is listed on multiple lines, only the coordinates in the last line will be used.
    
- -c/--cutoff: The name of the last scaffold to consider when building superscaffolds. All sequences after this will be ignored.
    
- -k/--blacklist: One or more scaffolds that should be excluded from the analysis. Subscaffold IDs may specified as well.

- -t/--table: Generate a detailed table of the contigs and strains supporting each join in addition to the standard output of superscaffolds and orphan scaffolds.

### Example command used to generate the information for Table S3:

```./assembler.py -b File_S9_Heber.bam File_S10_Lon_Inb.bam File_S11_Parrott.bam File_S12_RS.bam File_S13_ShCo.bam File_S14_SR-VP.bam -l 7500 -d 75000 -f 0.25 -r Breakpoints.txt -c scaffold_44 -k scaffold_15.1 scaffold_42 -t```
