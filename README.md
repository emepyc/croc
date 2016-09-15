

### CROC - Finding chromosomic clusters in genomes

CROC is a tool that aims to help in the identification and analysis
of genomic gene clusters. This method has been successfully used before in the
identification of chromosomal clusters in different eukaryotic species and make
use of serialized hypergeometric distribution tests in a sliding window approach



### Requirements

CROC is written in Perl. So it should run in every OS with Perl. If you experience
problems using it, please, send an email to emepyc[AT]gmail.es

The code for hypergeometric distribution tests is written in C and loaded in the
Perl script via the 'Inline::C' Perl module. If you encounter problems executing this
software, you can use an "all Perl" version of this module with the "--allPerl"
command line option when invoking the script.


### Installation


Download the file croc-<VERSION>.tar.gz
Unpack it and "cd" into the created directory
Altenatively, you can make a soft link of the executable file (croc.pl) in a path directory.

To view all available options type:
$ ./croc.pl -h
or
$ perl croc.pl -h

This distribution comes with already available genomes and samples to test the script in the
"Species" and "Datasets" folders.


### Options

```
-w | --window [integer]    : Width of the length to calculate the clusters (Default: 30000)
-o | --offset [integer]    : Offset between windows (Default: 10000). Only has effect with -w
-g | --gwindow [integer]   : Defines the window by number of genes instead of by number of nucleotides (and sets offset to 1 gene)
-p | --pval [float]        : Max allowed pvalue per cluster before multiple testing correction (Default: 0.05)
-m | --min_genes [integer] : Min number of genes that can for a cluster
-a | --genomeref           : Tells the program to use the entire genome as reference in the statistical calculations (by default, the references are each chromosome)
-t | --fdr [string]        : Method to apply to correct multiple testing. Possibilities are:
                                    BO => Bonferroni
                                    HO => Holms
                                    BH => Benjamini & Hochberg (default)
                                    BL => Benjamini & Liu
                                    NO => No correction
--allPerl                  : Use pure Perl code to perform the hypergeometric distribution tests
                             Using this option, the script will run slower, but will run in systems where you can't use Inline::C
-h | --help                : Displays this document
```


### Formats


+ Sample file: Should be a file containing a list of genes (1 per line).
+ Reference file: Should contain a first block of lines having the chromosome names and its length separated by a tab character. After that, a blank line.
Finally, another block of lines describing the genomic features (f.e. genes). These lines must have the following fields separated by tabs:

        <ID> <Chr> <Start> <End> <Name>



### Bug report and feedback

If you find any issue using this software, please, contact me at:
emepyc[AT]gmail.com
