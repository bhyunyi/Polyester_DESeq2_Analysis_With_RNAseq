# Installs Polyester if it is not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = "/path/R/libraries")
BiocManager::install("polyester", lib = "/path/R/libraries")

# Installs Biostrings if it is not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = "/path/R/libraries")
BiocManager::install("Biostrings", lib = "/path/R/libraries")

# Imports necessary libraries
library(polyester, lib.loc = "/path/R/libraries")
library(Biostrings)

# sets the working directory to output files into
setwd("/path/Chr22/")

# creates a small fasta file with 50 transcripts 
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)
small_fasta = fasta[1:50]
writeXStringSet(small_fasta, 'chr22_small.fa')
fold_changes = c(rep(1, 50))
outdir = 'simulated_reads'

# ~20x coverage ----> reads per transcript = length/readlength * 20
# "width" is operating on a DNAStringSet (from Biostrings)
readspertx = round(20 * width(small_fasta) / 100)
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
                    num_reps=20, fold_changes=fold_changes, outdir=outdir) 
