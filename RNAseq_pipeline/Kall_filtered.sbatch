#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:10:00
#SBATCH --mem=8GB
#SBATCH --job-name=all_filt
#SBATCH --output=slurm_%j.out

#removing all fastq files from pseudoalign folder
for file in /path/Kallisto_pseudoalign/*RR*
do
rm ${file}/*.fasta
done

#purge all and then load python
module purge
module load python/intel/3.8.6

#set pwd to directory containing all SRA directories
cd /path/Kallisto_pseudoalign


python << EOF
import pandas as pd
import os
import re

all_data = []

for dirname in os.listdir('.'):
    if re.match(r'^(sample)', dirname):
        tsv_file = f"{dirname}/{dirname}_filtered_abundance.tsv"
        if os.path.exists(tsv_file):
            df = pd.read_csv(tsv_file,delimiter = '\t')
            df['sample'] = dirname
            all_data.append(df)

final_df = pd.concat(all_data, ignore_index=True)

final_df.to_csv('Kall_filtered.csv', index=False)
EOF
