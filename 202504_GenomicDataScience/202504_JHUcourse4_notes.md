**Course 4 -> Command Line Tools for Genomic Data Science**

Course link: <https://www.coursera.org/learn/genomic-tools?specialization=genomic-data-science>

---

- [Basic Unix Commands](#basic-unix-commands)
  - [Related to directories](#related-to-directories)
  - [Viewing file content](#viewing-file-content)
  - [Redirecting](#redirecting)
  - [Querying and processing content](#querying-and-processing-content)
  - [Comparing content](#comparing-content)
  - [Archiving](#archiving)
- [Sequences and Genomic Features](#sequences-and-genomic-features)
  - [Genomic feature representations: BED/GTF](#genomic-feature-representations-bedgtf)
    - [BED format: Browser Extensible Data.](#bed-format-browser-extensible-data)
    - [GTF format](#gtf-format)
  - [Alignments: SAM and paired-end sequencing](#alignments-sam-and-paired-end-sequencing)
    - [SAM format: Sequence Alignment/Map](#sam-format-sequence-alignmentmap)
  - [Retrieving sequences and features](#retrieving-sequences-and-features)
  - [SAMtools and BEDtools](#samtools-and-bedtools)
    - [SAMtools](#samtools)
    - [Sorting and indexing](#sorting-and-indexing)
    - [BEDtools](#bedtools)
- [Alignment \& Sequence Variation](#alignment--sequence-variation)
  - [Genomic tools](#genomic-tools)
    - [Alignment tools: Bowtie and BWA](#alignment-tools-bowtie-and-bwa)
    - [Sequence variation tools: SAMtools mpileup and BCFtools](#sequence-variation-tools-samtools-mpileup-and-bcftools)
    - [Simple application: variant calling](#simple-application-variant-calling)
- [Tools for Transcriptomics](#tools-for-transcriptomics)
    - [Running a shell script](#running-a-shell-script)
    - [Alignment: TopHat](#alignment-tophat)
    - [Transcript assembly: Cufflinks](#transcript-assembly-cufflinks)
    - [Differntial analysis for expression and splicing: Cuffdiff](#differntial-analysis-for-expression-and-splicing-cuffdiff)
    - [Visualization: IGV](#visualization-igv)


---

## Basic Unix Commands

> Some background info to get things *rlly* clear... OvO

- **Unix** is the pioneering OS (Operating System, acts as a bridge between users/programs and physical hardware). Unix-like OS include macOS, Linux (Ubuntu), and Android (Linux-based).
- A **shell** is a command-line interpreter that lets users interact with the OS by typing commands. It’s a layer above the OS kernel. 
  - User types a command (eg. `ls`) -> shell interprets it and asks the OS kernel to execute -> kernel accesses hardware (eg. reads disk to list files) -> results are returned to the shell -> displayed to the user.
  - Shell types include Bash (default on most Linux systems) and Zsh (default on macOS, Zsh can run most Bash scripts).
- The **terminal** isn't the same thing as the shell! It provides a text-based interface (an input/output *window*) for users and **starts a shell**, but it's *the shell* that interprets the user's command by requesting the OS kernel to run it and retuning the output to the terminal. 
  - A shell can run without a terminal, eg. when running .sh scripts, *we have demos (and practices) in the last module of this course that demonstrates the use of scripts for RNA-seq tasks to make it run in the background* (so do stay tuned~~).



### Related to directories

- `pwd` prints working directory; `cd` changes directory, **relative path has no slash before the subdirectory name**; `cd .` denotes the current directory; `cd ../dir` locates to the **parent directory, and then** to its (the parent one's) subdirectory; `ls` lists directory content; `man command_name` opens manual page for that command.
  - Wildcards (like regex) match directories/files. (`*`, `?`, `[]`, `{}`)
  - `ls -l` or `ls -lt` also displays user operations (the former sorts in alphabetical order (default), the latter sorts by modification time since that's what `-t` denotes. **Note that when multiple options in the same location of the command are present we only write one dash**, like `-lt` rather than `-l -t`)

- `mkdir new_dir_name(s)` creates a new directory under current directory. 
- `cp file(s) dir` to copy the files to the target directory. `mv file(s) dir` to move. 
- `rm` to remove files. add option `-i` to ask for confirmation before removal. 
  - `rmdir` to remove *empty* directories. **add `-r` to remove recursively** (from the subdirectories and then the directory iteself). 

### Viewing file content

- `more file_name` views content of the file. Press Enter to view the next line and space to view the next page. Press q to exit. 
  - `less file_name` is more advanced (named as a joke 'less is more'). Allows for bidirectional scrolling via arrows, and press g (start) or G (end) to jump to start/end.
  - `/` slash + pattern to search for a specific line. Eg. for a multi-FASTA-formatted file in which it contains multiple DNA sequences each starting with *the header*: a `>` and then a description. Thus, we can search for each sequence by entering `/>` consecutively.

- `head/tail file_name` to show the first/last 10 lines of the file. `head -20 file_name` to view the first 20 (in this case) lines. 
- `cat file_name(s)` views *all* contents in the file(s). If multiple files, then the output is just concatenating the files. 
  - Add `-n` to return content of each line **preceded by the line number**, an option also applicable to many other commands. 
  - Use wildcards like `cat */*.md` to display the content of all matched files in a directory. 
- `wc file_name(s)` word count, returns 1. number of lines; 2. number of words; 3. number of characters; 4. the file name. Add `-l` to only return the number of lines. 

### Redirecting

- We can redirect the **standard output** ('stdout', shown in the **terminal**, eg. the output of `wc file_name(s)`) **into a file (standard input) under your current working directory** using `>`. The file created after this operation is called 'md_wordcounts' (located in the working directory), which is used in the following examples. 
  - Similarly, we can use `<` to redirect files into a standard output. Eg. `wc < *.md` returns an output of '1035   12033   79738'.
  - `2>` redirects standard errors ('stderr', the error messages). *Modern* shells use `&>` to combine both redirections (stdout and stderr) to a file (eg. output.log) or another destination *so that nothing is printed to the terminal*.  

- Another redirect operation is the **pipeline**, with a vertical bar `|`. It **redirects the output from the prior command into the input of the subsequent command**. Eg. `cat file_name(s) | more` to browse the file(s) one page at a time instead of having all the content displayed at once; `ls | wc -l` to give us the number of items in the current working directory (since it counts the number of lines in the output of 'ls').

### Querying and processing content
- `sort text_file_name(s)` displays the fields sorted alphabetically. Add `-r` to sort in reversed order. Add `-k3` to sort according to the 3rd (in this case) column.
  -  For a column with numerical values, it still sorts alphabetically by default, eg. '10' goes before '2'. So we use `-k 3n` (i.e. add `-n`) to specify that it's a numerical column, and `-k 3nr` for descending order rather than ascending. 
  -  Sorting by multiple columns (and prioritizing the ones to the left): eg. `sort -k 3n -k 4 file_name(s)`.
  -  Example: `sort -k 4r md_wordcounts "md_wordcounts copy"` (since the duplicated file **contains a space in its filename we have to have it quoted**), returns: 
```
    1035   12034   79738 total
    1035   12034   79738 total
      19     323    1931 202504_JHUcourse4_notes.md
      19     323    1931 202504_JHUcourse4_notes.md
     714    6020   39067 202504_GoogleDA_notes.md
     714    6020   39067 202504_GoogleDA_notes.md
     302    5691   38740 202504_DLbook_part1_math.md
     302    5691   38740 202504_DLbook_part1_math.md

```
    
  - Add `-u` to only return the unique fields. `sort -k 4r -u md_wordcounts "md_wordcounts copy"` returns: 

```
    1035   12034   79738 total
      19     323    1931 202504_JHUcourse4_notes.md
     714    6020   39067 202504_GoogleDA_notes.md
     302    5691   38740 202504_DLbook_part1_math.md
```

- `cut` accesses a specific **range in the fields**. By default, it uses tabs as delimiters. Eg. to *cut columns one to three*, use the `cut -f1-3 file_name` command (don't use a colon to specify the range!! it's a dash!!) 
  - Since it looks for tabs by default, add `-d 'custom_delimiter'`. Eg. `cut -d ' ' -f1 file_name` cuts column 1 by space.
- `uniq file_name` does almost the same thing as `sort -u`, but it only *removes duplicated fields that are **consecutive***. 
  - Add `-c` to count the number of *consecutive occurrences* of fields. 
- `grep 'pattern' file_name(s)` returns the **all the specific lines** that contain the pattern. Eg. `grep -n ' that contain ' *.md` returns (with the matched lines preceded by the line number via `-n`):
```
202504_GoogleDA_notes.md:295:Now in the results, we can see the entries that *violated* length = 2. But we still want to select those that contain a specific substring (even though the lengths are mismatched). 
202504_JHUcourse4_notes.md:63:- `grep 'pattern' file_name(s)` returns the **specific lines** that contain the pattern. Eg. `grep -n ' that contain ' *.md` returns (with the matched lines preceded by the line number via `-n`):
```
- `grep -c 'pattern' file_name(s)` counts the number of occurrences of the pattern, and does not return the lines that contain it (so this is equivalent to using the original `grep` first and then using pipeline `|` to `wc -l`). 
  - Eg. still `grep -nc ' that contain ' *.md` but with a 'c' in options. 
```
202504_DLbook_part1_math.md:0
202504_GoogleDA_notes.md:1
202504_JHUcourse4_notes.md:4
```
- `grep -v 'pattern' file_name(s)` says **ignore the pattern and return everything else**. 

**Additionally:**

- `awk 'pattern {action}' file_name`: super versatile, basically **scans the pattern and processes the action for lines in the input file that match the pattern**. 
  - It reads the input file **line by line**, and splits each line into **fields**. Default field separator is whitespace; use `-F` to set it. (eg. for a *comma-separated CS*V, `awk -F,'pattern {action}' file_name`)
  - Patterns can be: lines containing the text in `/text/`; relational expressions eg. `$1 > 100` (`$0` means the entire line, `$1` means the first field, ..., `$NF` means the last field) 
  - Actions: print, arithmetic operations, etc. 
  - `BEGIN {action}` and `END {action}` executed before and after reading/processing all input lines.

### Comparing content 
- `diff file1 file2` returns the locations and the types of differences of lines that differed in the two files, and the content they differed (`<` for the first file and `>` for the second).
- `comm file1 file21` returns a table-like thingy, if the line appears only in the first/second file, the content of it will be listed in the first/second column of the table. If it appears in both of them it's going to be listed in the third column, *but the lines in the input files must be sorted in the same order* to return the correct common lines.
  - If we only want to focus on one column (say, the 3rd one since it contains the common lines), add `-1 -2` to exclude the other columns.

### Archiving 
- `gzip file_name` to get a .gz file with the original one gone. 
  - `gunzip file_name.gz` to retrieve it. There's no loss of information to the compression and then the decompression. 
- `bzip2 file_name` has more compressing power. `bunzip2 file_name.bz2`.
- `zcat zipped_file` to **view the file without unzipping it**. 

**Compressing multiple files together**: we first create an archive that chains together the files, and then compress that resulting archive.
- `tar -cvf archive_name.tar file1 file2 ...`
  - Options: 'c' for 'create' to create an archive, 'v' for verbose mode to print the progress to the terminal, and 'f' means 'from the following file', which is archive_name.tar in our command.
  - After this step, the **.tar simply archives the separate files but does not compress them**. So we `gzip` the archive_name.tar file to get archive_name.tar.gz.
  - After we have the .tar.gz file, first `gunzip`, and then `untar -xvf archive_name.tar`, 'x' means to extract files. 


## Sequences and Genomic Features

### Genomic feature representations: BED/GTF
Genomic features: determined by gene annotations.
- Precise location: the start and end sites on the strand.
- Structure: eg. clustering of exons.

#### BED format: Browser Extensible Data.
- A single interval:
  - The basic format only needs 3 columns: the *genomic substrate* (which chromosome or scaffold that feature came from), the start and the end. 
  - These coordinates are 0-based!!! But the weird thing is that though the start is indeed 1 less than an 1-based coordinate, the **end** indices for both 0-based and 1-based are equal...
  - We can add more fields/columns to the basic format (i.e. only 3 columns) to show more information, such as the name of the feature (eg. Exon3), +/- strand, how the feature should be interpreted for visualizations, etc. 
- **To group related intervals (like a cluster of exons for a gene) together**: put additional fields (columns) at the end of **a BED record that represents the gene/the whole cluster**. 
  - **BED12** format, containing a total of 12 columns; features three additional columns named 'block count', 'block sizes' and 'block starts' for each block (each single interval). 

#### GTF format
- Columns 1-9 separated by tab `\t`; fields **within column 9** (this column provides more than one piece of identifier info, eg. 'gene_id "genA"; transcript_id "genA.1"', and other additional pieces of information) **separated by space** (this is very important since in `cut` we have to specify that the separator is ' '). 
  - Since the column 9 is so informative, one useful way to carry out analysis using a GTF file is via the `cut` command to **cut the 9th column** and pipelining it to other subsequent operations (like `sort -u`, `wc -l`, etc.) 
  - Say, we only want to view the output transcripts (eg. "genA.1" as denoted above). We can count that it is the 4th 'column' (again, separated by spaces) *within the 9th column*. So we use the pipeline: `cut -f9 file_name.gtf | cut -d ' ' -f4`. 
- **Each interval feature (eg. each exon) takes one line/row/record**. 
  - Intervals in the same cluster are grouped together by **stacking the rows in consecutive order**, and we can see from *column 9* that the 'gene_id' attribute is identical (eg. all "genA") for these consecutive rows. 
  - Distinguish from how **BED groups the clustered intervals by adding more columns to a single record about the whole cluster**, so BED has a more compact presentation. 
- Unlike BED, coordinates are 1-based. 
- Conceivably, regarding the operations via BEDtools, **the results on GTF format (separate exonic features) are basically just similar to using `split` on BED12 format**.

**GFF3** format: 
- Just like the GTF format, also one line for one single exon in the gene. Also 1-based.
- The last column gives some information about the single feature and the cluster that it belongs to (i.e., the upper-level entity). 
  - Eg. for an exon (for illustration purpose) with 'ID=exon00001', it also contains the info 'Parent=mrna002'.

### Alignments: SAM and paired-end sequencing

- Sliced alignments are for mRNA sequencing due to the presence of introns in reference genomes. 
- For longer effective coverage and better alignment accuracy, Most NGS (except for small RNAs) use **paired-end (PE) sequencing, which generates two reads** (R1 on the forward strand of the DNA fragment and R2 on the reverse strand) **from the same DNA fragment but its both ends and thus both strands**. The SAM/BAM format stores information about their relationship in fields like FLAG.
  - The term 'mate' just refers to the partner read in a paired-end sequencing experiment. 
  - In PE sequencing, **two separate FASTQ files (R1.fq and R2.fq)** are produced, since NGS sequencers (eg. Illumina) can't interleave R1/R2 during sequencing and *each read's data is stored separately in real time*. Mappers (BWA and Bowtie2, elaborated in a following section) require R1/R2 to be ordered identically in their files, eg. *Read 1 (say, represented by Line 1 of the file) in R1.fq always pairs with Read 1 in R2.fq*.
- When the partner reads are mapped on the genome in opposite directions and facing each other, and the distance between them is within the boundaries that are specified for the original set of fragments (the fragment size of reads should follow a normal distribution with a given mean and SD), we call that particular mapping as being **properly paired** or **concordant**. If any of those conditions is being violated, we call the mapping non-concordant.

#### SAM format: Sequence Alignment/Map

Representations of **alignments** (BAM (Binary Alignment/Map) is another format but compressed and *represented in binary*). 
- The first portion is the header, which contains the metadata. Every line in the header will start with an `@`, and then by a very short code that tells us what the line does. Eg. 'PG' tells us the alignment program.
- Following the header, we have information about the alignments of the reads to the reference genome. The second field is **FLAG**, and the fifth field is **CIGAR**. 
  - FLAG gives us a combination of multiple pieces of info (represented in binary) about this particular match for **paired-end reads**, eg. whether it is properly paired with the mate (i.e. partner read).
  - CIGAR gives information about **how that read aligns to the reference genome** (recall the section on 'edit distance' covered in Course 3). It can be represented as **a sequence of short strings that each contains a number followed by a letter**. The letter marks edit operations (eg. 'M' for matched bases between the read and the reference genome, 'N' for skipping (introns) and 'D' for deletion), and the number tell us how many such edit operations (i.e., the length of corresponding bases) we have.

### Retrieving sequences and features
**For sequences:** 
- Go to a database (NCBI GenBank in the example), select the specific database (eg. 'Nucleotide') and use keywords to search for a sequence. Download the item in FASTA format. 
- In other cases (like accessing FASTQ sequences in the example; we select SRA for RNA-seq data), we can use command line tools to download the SRA data. Use `wget` and the Internet address that it's going to fetch and download from. It deposits the data precisely in ur intended directory.
  - After it is downloaded via `wget` and saved, we have to convert it to the appropriate format to use. The command is `fastq-dump`, and we add `nohup` (no hangup) before the command to **execute it in the background**. This means even after closing the environment on my laptop, this is still going to run on the remote terminal; add a `&` after the `fastq-dump` command so that we can have access to the shell for other commands at the same time. 
    - While the command is being executed in the background, any errors reported will be automatically written to the file 'nohup.out'.
    - The `fastq-dump` command converts the SRA data to a .fastq file which contains the 4 lines for each read and we can view them  normally.

**For features**:
UCSC Table Browser (used in the example) or Ensembl.

We can use the program to retrieve data about a target genomic feature in a selected format (eg. BED or GTF), and then we save it to a .txt file. 

### SAMtools and BEDtools

#### SAMtools
 SAMtools is a suite of utilities to manipulate SAM and BAM formatted files. Allows us to process, sort, index, filter, and analyze alignment files efficiently.

- Run a command (eg. `samtools flagstat`) without any arguments to **return the syntax ('usage') and the options**.

#### Sorting and indexing

`samtools flagstat file_name`: returns **simple statistics about the alignments** in the BAM/SAM file, eg. number of reads mapped, properly paired, duplicates, and other alignment categories. 

`samtools sort input_file_name output_file_name`: **sorts alignments by position**/genomic coordinates (default) or by read name (`-n`). 
- This operation can be very slow so again we can add `nohup` before it and `&` after it to make it run in the background.
- After it's sorted, when we view the header of the file, in the first line after '@HD' there will be an eg. 'SQ:coordinate' to denote how it's sorted. 

`samtools index sorted_file_name`: **after the file is *coordinate-sorted*, creates an index** (sorted_file_name.bam.bai file). A BAM index (.bai file) is *a small auxiliary file* that allows rapid random access to specific genomic regions in the .bam file. **After that, an indexed BAM file allows tools like `samtools view` to quickly jump to a region (eg. chr1:1,000,000-2,000,000) without loading the entire file**. 
- Index a FASTA file via `samtools faidx genome.fa`. Also creates a FASTA index (.fai file) and the purpose for indexing is the same (just consider how massive a FASTA file can be...).
- We have an example at the end of this module (using the command `samtools tview`) featuring both an indexed BAM and an indexed FASTA, which i think further explains this indexing matter. 

---

`samtools merge merged_file_name file1 file2 ... `: merges together multiple alignment files.

`samtools view`: three different operations
- **View alignment**. The default `samtools view file_name` does not include the header portion of the file (recall!!). Add `-h` to include header in output, or `-H` to only view the headers (and exclude the alignments portion).
- **Conversion between BAM and SAM formats**. Without options, the default output (also when the input is a .bam file) is a .sam file; to save it, just `> output.bam` to redirect it. 
    - To convert into the BAM format: add `-b` to output in binary. 
    - The option `-T` specifies a reference genome, so when that option is present it is followed by a FASTA (.fa) file.  
- **Allow extraction of alignments only within a range of the genome**. The input file must be **indexed**. We can either use the optional argument in the end of the command to specify the search region, eg. `samtools view input.bam "chr1:1000-2000" | head` and view the first ten lines to see whether the genomic coordinate matches our specified range (reminder: use dash for a range rather than a colon!!!!!), or use the option `-L bed_file.bed` to only include reads overlapping this BED file, which will contain intervals. 

#### BEDtools

BEDtools is a versatile suite that contains tools for **genome arithmetic** (especially based on overlaps), eg. finding overlapping intervals (intersections), clustering intervals, subtracting intervals, etc. It also allows for **format conversions** (eg. between BAM and BED), and contains subcommands for **sequence extractions** such as `getfasta` which uses the given intervals (.gtf or .bed) to extract the FASTA sequences.

- The subcommands usually require explicitly typing the option (eg. argument file A preceded by `-A`) before entering the corresponding argument.

- Subcommand + `--help`: returns the help page for syntax (usage) and option choices. 

`bedtools intersect`: Finds overlaps between two sets of genomic feature from input files eg. BED, BAM, etc. Enables us to extract reads/regions overlapping a reference.
- `bedtools intersect [OPTIONS] -a file_A -b file_B`
- Different types of **output**: `-wo` is like INNER JOIN, returning overlap-only lines of both A and B (concatenated vertically with B on the right); so `-wao` is more or less like LEFT JOIN (with 'LEFT' being A). `-wa` and `-wb` just return the overlapped lines in A *or* B, which means no 'JOIN's. `-v` only reports entries in A with no overlaps in B.
- `-f`: minimum overlap required, specified by a fraction of the length of A. There's also the option `-r` which requires the fraction of overlap be applied equally to A and B. 
- `-split`: in case we have a BED12 entry, in which one row corresponds to a whole cluster that contains multiple blocks. Then each one of those blocks or exons should be 'splitted' to be considered separately when finding the intersection, so we're finding overlap by exons rather than the entire feature (the whole gene).
- For a simple intersection analysis, the `intersection` command is often accompanied (via pipeline `|`) by other basic UNIX commands like `cut`.

`bedtools bamtobed`: converts BAM alignments to BED format (extracts genomic intervals). 
- `bedtools bamtobed -i input.bam > output.bed` (`-i` is again just an explicit flag that specifies the input file for the argument)
    - `-split`: Split BAM CIGAR strings into separate exons (**as separate *BED12* entries**). Eg. used for splitting spliced alignments in RNA-seq. 
    - `-cigar`: Include CIGAR strings in the output BED entry as a 7th column.

`bedtools bedtobam`: converts BED/GTF files (feature records) to BAM format, requires a header from the reference genome.
- `bedtools bedtobam [OPTIONS] -i input.bed -g genome.fa.fai > output.bam` (`-g`: specifies genome file)
- The genome file requires a special format, it's rather a *header* file than an actual genome sequence file.
- `-bed12`: input is in BED12 format (i.e., 12 columns, has multiple blocks (eg. exons) in one entry)

`bedtools getfasta`: extract sequences from an input FASTA file into an output FASTA file based on feature coordinates in the BED/GFF/VCF file which specifies the range to extract. 
- `-split` indicates that our **BED12** file contains multiple exonic blocks for entry. So with this option, instead of returning lengthy sequences based on *the entire gene feature in one single entry due to the construction of the BED format*, `getfasta -split` will extract and concat exon sequences *separately* according to the blocks denoted in the additional columns of BED12. 


## Alignment & Sequence Variation
// Suffering from the poor quality of the subtitles,,, Q_Q

The bioinformatics workflow for handling and for detecting sequence variations (due to polymorphism) would start with alignment of reads to the reference genome, followed by **variant calling**: identifying variations (detecting SNPs, indels and structual variants) between reads and the genome. 

**(SAMtools) mpileup** format: a type of format to represent **variants**. Each line represents information about *a single position* (one 1-based coordinate) along the genome at which we have reads aligned. We then have a number of fields that are again separated by tabs `\t`, each containing one distinct piece of info eg. the chromosome number, the specific position/coordinate, the *depth* (the coverage/number of reads sequenced and aligned) at that position, a string of characters that *shows all observed bases in each of these aligned reads at this position*, etc.

**VCF** format: Variant Code Format (**BCF stands for Binary Code Format so it's machine-readable, smaller in size, and faster to process for pipelining. VCF is text-based so human-readable**). This format also gives information about sequence variation. Unlike the mpileup format (one line for *every base alike* in the genome at which we have reads), the VCF format gives detailed information about **only those positions at which we identify variation**.
- Involves **genotype likelihoods**, useful for variant calling.
- A simple example of how to represent sequence variations in VCF: (ID is a dot if the variation hasn't been identified in databases such as dbSNP; REF shows the base in the reference while ALT shows how the base varied in the read.) Note how indels are represented by showing bases preceding the variation, which is why the 'POS' (position) starts earlier.
  
![VCF](https://freeimghost.net/i/Screenshot-2025-04-26-at-09.37.19.xo7ODP)

### Genomic tools

#### Alignment tools: Bowtie and BWA

These are **very fast tools for mapping large numbers of sequences**. They do so by using a compressed representation of the genome as an **index** (recall from previous sections!!), so that while mapping, each read is searched against the genome index to identify matches, and then more complex alignment algorithms are used to create the entire map. Both Bowtie and BWA search for contiguous matches that only have a small number of indels.


**Bowtie2**:
- If we're not working with a model organism (for which we can download a copy of its genome index) or we only want to map a strict number of regions, we have to first create an index by `bowtie2-build <reference_fasta_file> <index_file_prefix>`, which will create the genome index file with a .bt2 (for Bowtie2) extension. 
- After we've got the genome index, the next step is to actually **map reads to the genome index**. `bowtie2 [options] -x <index_file_prefix> <reads.fq> -S <output.sam>`
  - `--help` to print the usage/options message.
  - `-x` specifies the genome index file (ensure that the argument only include the **prefix of the file name with no file extensions** like .bt2)
  - Input files (reads) are FASTQ .fq/.fastq by default (`-q`); `-f` if input files are (multi-)FASTA .fa/.mfa.
  - `-S` specifies that the output is in SAM format.
  - We can choose the mapping to be **global alignment** (`--end-to-end`, default) or just match a portion of the input read, i.e., **local alignment** (`--local`. Recall, this one's basically finding the most similar *substring* in the two given strings rather than iterating over the entire strings (constructing the dynamic programming matrix) to compute the distance between them, as in 'global')
  - For both Bowtie2 and BWA (shown below) there's an option that specifies the number of threads (`-p <int>` in Bowtie2 and `-t <int>` in BWA MEM). This allows you to set the number of CPU cores for alignment (parallelizes; speeds up mapping especially for large datasets).

**BWA**: 
- Also starts with creating an index file. `bwa index <reference_fasta_file.fa>` to create a number of index files. 
- Mapping the reads to the index file. Three ways for alignment, but the most applicable one is `bwa mem [options] <reference_index.fa> <input_fastq> [input_fastq2]` and we can redirect the output to `> <output.sam>`. 
- For paired-end reads, the input files should be two (eg. R1.fq and R2.fq, as explained earlier in the section on SAM).

#### Sequence variation tools: SAMtools mpileup and BCFtools
**samtools mpileup**: Based on *aligned reads in the BAM* input file, generates a summary of sequence variations between the reads and the genome. `samtools mpileup [options] <input.bam> > output.mpileup` the default output format is just the **.mpileup format we've covered earlier for variation representation**.
- `-f <reference.fa>` specifies the *indexed* fasta file for the reference genome. (Note that the BAM file must also be indexed using `samtools index <input.bam>`, covered earlier.)
- `-v` outputs VCF file (text-based, for analysis by humans. add an option `-u` (forms `-uv`) so that the VCF isn't compressed); `-g` outputs BCF file (binary, compressed, thus faster processing in pipelines). 
  - When either option is used, we can generate **genotype likelihoods** (since this field is contained in the two formats). There are also a number of options for calculating genotype likelihoods. 

**BCFtools**
A suite including many commands (like SAMtools) for manipulating VCF/BCF files (eg. indexing, viewing, and other familiar operations), especially useful for variant calling. Usually operates on the output from `samtools mpileup`.

Some common options for all options are: 
- `-O <z|u|z|v>` specifies output type (compressed/uncompressed VCF (`z`/`v`)/BCF (`b`/`u`)). For final analysis in variant calling, output in VCF since it's text-based.
- `-o <file>` specifies output file name.

`bcftools view [options] <input> [region1 [...]]`: for VCF/BCF conversion, view (remember that BCF was binary so viewing in text-based human-friendly VCF), filter (via the specified regions). Sounds so similar to `samtools view`.

`bcftools call [options] <input>`: for SNP/indel variant calling from VCF/BCF, used in conjunction with `samtools mpileup`. 
- `-v` outputs only variant sites.
- `-m` for multi-allelic calls.

#### Simple application: variant calling
**Alignment** (via alignment tools eg. BWA MEM or Bowtie2, index first tho) -> **variant calling** (via SAMtools mpileup and then BCFtools `call`) -> **analysis**.

> Q: Why use `bcftools call` after `samtools mpileup` when the latter can already output a VCF file?

> A: The latter outputs raw data representing **candidate sites** of variation. The former, on the other hand, allows **further evaluation** by applying statistical models to *assign* (rather than simply calculating the *likelihoods*) **genotypes** (eg. 0/1!! OvO) into the 'GT' field, handles multiallelic sites, filters low-quality calls (since `samtools mpileup` *outputs all mismatches* including potential sequencing errors) by calculating confidence metrics, etc. 

Here's an example. Assume that we already have the sequenced reads in a FASTQ file aligned with the reference FASTA file into the 'example.bam' file through alignment tools. 

```bash
samtools index example.bam
samtools mpileup -g -f reference.fasta example.bam > example.bcf

bcftools call -v -m -O z -o example.vcf.gz example.bcf

zcat example.vcf.gz
```

**After analyzing variations in the VCF file** (analysis is prob the main purpose for a variant calling task), we can also visualize one variant (record its positional info from the VCF file) of interest by `samtools tview -p <chromosome_num:position> -d T example.bam reference.fasta`. This is a command for BAM visualizations given a specific position on the chromosome. Both the BAM file for the read and the FASTA file for the reference must be **indexed to allow the tool to quickly jump to a region according to it's genomic positions!!**
- `-p` specifies its position, while `-d T` displays the visualization as text (default is by HTML, i.e. `-d H`).

## Tools for Transcriptomics
There are usually a wide variety of transcripts due to alternative splicing, and some variants of them may cause diseases. In transcriptomics, we aim to answer:  
- What genes are transcribed in the sample? What transcripts are expressed for each gene? (**Transcript assembly**, in which the aligned transcripts are clustered in gene-oriented groups to construct a basic structure, represented by graph. We identify the most likely splice variants from the reads we have (usually even more varied) by enumerating and then assigning reads back and calculating the expression level, using algorithms, and selecting the top variants as being the most likely.)
- How many copies of the gene are transcribed? How many copies of each type of the transcript are there? (**Quantification** from the basic structure)
- How do the expression levels and splicing patterns differ between two different conditions (eg. disease vs. healthy tissues)? (**Differential analysis**)

So the RNA-seq workflow is: **sequencing -> alignment/mapping of reads -> transcript assembly and quantification [-> differential analysis]**

![Screenshot](https://freeimghost.net/images/2025/04/26/Screenshot-2025-04-26-at-14.46.55-copy.png)

// The source is old... (and so is this course i know...) but i'll cling to it rn.

#### Running a shell script
An alternative to just typing the command in the terminal, we can create a simple shell script. We can just open an editor eg. by `vi` with a file name, and edit the script to write the command(s). The commands are the same as what we would type in the command line, eg. `mkdir` and `tophat -o ...`, but we can also add shortcuts (? i'm not sure whether they are called so) to make the commands shorter and the arguments more organized. 

Type the name of the script in the command to run it. We can also set `nohup` and redirect stdout and stderr to a .log so it won't output anything~~~.

#### Alignment: TopHat 
Maps FASTQ sequences from spliced transcripts to whole genomes, requires a Bowtie2 index and uses Bowtie2 as its alignment engine. Returns a BAM file ('accepted_hits.bam') for aligned reads (also outputs detected splice junctions in BED format, etc.) in a directory. 
```bash
tophat -o output_directory -G reference_annotation.gtf genome_index_prefix R1.fastq R2.fastq
```
Recall that the GTF format represents annotated genomic features, containing the informative *9th column* and a separate record for each exon (of each transcript).

**Note that even though TopHat is doing an alignment task, it uses the `-G <reference_annotation.gtf>` to *improve splice junction detection and guide alignment of RNA-seq reads***, since exon-exon junctions are not contiguous (i.e. spliced) in the genome. The GTF file provides *known transcript structures* (exons, introns, splice sites), allowing TopHat to reduce false-positive novel junction predictions.

#### Transcript assembly: Cufflinks
 Assembles transcripts from RNA-seq reads aligned by TopHat ('accepted_hits.bam', the output of `tophat`) and estimates their abundances.
```bash
cufflinks -o output_directory -g reference_annotation.gtf /tophat_output_dir/accepted_hits.bam
```
Outputs `transcripts.gtf` for assembled transcripts, **sample specifc**.
- One issue is overassembly, is Cufflinks may predict many novel transcripts that are **artifacts** (i.e. technical false positives caused by PCR/sequencing/analysis methods lolol). If we use `-G` rather than `-g` (the latter uses the input reference_annotation.gtf as **a guide/anchor** but *still allows novel isoform discovery*), this flag forces Cufflinks to *only use the annotated transcripts from the input GTF* to avoid overassembly (and no novel isoform discovery).
- Other output files include two .fpkm files that give us the expression estimates.

#### Differntial analysis for expression and splicing: Cuffdiff
Before this operation, the full workflow also involves - `cuffcompare` that **filters individual sample assmebly (resulted from Cufflinks) before merging sample assemblies** (returns a .combined.gtf), and
- `cuffmerge` that **merges transcript assemblies from different samples** into a single, unified transcript catalog (merged.gtf) for `cufflinks`. 
  - Cuffmerge requires us to create a .txt file (eg. `assemblies.txt`) listing paths to all transcripts.gtf files **each** resulted from an **individual sample** Cufflinks assembly.   

```bash
cufflinks -o sample1_out accepted_hits.bam
cuffcompare -r reference_annotation.gtf -o output_prefix sample1_out/transcripts.gtf sample2_out/transcripts.gtf
cuffmerge -g compare_out.combined.gtf assemblies.txt
```

**Cufflinks** outputs a set of .diff files (also tab-delimited with headers, contains info like P-values (after all, it's comparing between two samples)) with two types: expressions and splicing. 
```bash
cuffdiff -o diff_out \
         -L cond1,cond2 \
         merged.gtf \
         cond1_rep1.bam,cond1_rep2.bam \
         cond2_rep1.bam,cond2_rep2.bam
```
- `-L`: labels for conditions.
- Input merged GTF file and BAM files grouped by condition (rep1 means replicate1). **Both these types of files are required for expression quantification**: the merged GTF defines the transcripts to quantify, while the BAM files provide the **raw read data** for actual quantification and statistical testing.

#### Visualization: IGV
We can use IGV (Integrated Genomics Viewer) to visualize and manually annotate the results of these prior differential analyses.

There are two types of files that we need to use to visualize in IGV: the alignment files (from TopHat) and GTF files after transcript assembly (Cufflinks/Cuffmerge). Again, both types of files need be **indexed** to allow easy access to random locations by the tool. 
- Recall that the files must be **sorted by genomic position** before they are able to be **indexed**. 
- For GTF files, this can be done simply by `igvtools sort <input_file.gtf> <sorted_file.gtf>`. And then we can index them; still simple: `igvtools index <sorted_file.gtf>`. This operation will also output a corresponding `sorted_file.gtf.idx` auxiliary index file.
- For the alignment files (BAM), we can first extract only part of the entire sequence by the [region] argument in `samtools view -b` (eg. "chr1:1000-2000". `-b` added since we still want the output to be in BAM; recall the three functions of `samtools view`!), and then use `samtools index` to index it.