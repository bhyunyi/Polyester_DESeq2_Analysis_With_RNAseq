#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=Kall4one
#SBATCH --output=slurm_%j.out


#initialize variables
sra=""
cell_type=""

#parse command-line arguments for -idx and -sra flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -sra) sra="$2"; shift ;;
        -cell_type) cell_type="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

echo "sra file: $sra"
echo "cell_type: $cell_type"

cd /path/Kallisto_pseudoalign
mkdir ${sra}
cd /path/Kallisto_pseudoalign/${sra}

cp /path/simulated_reads/${sra}_1.fasta .
cp /path/simulated_reads/${sra}_2.fasta .

#purge all modules and load kallisto
module purge
module load kallisto/0.50.1

#run the kallisto pseudoalignment
kallisto quant -i /path/Chr22/chr22_small.idx -o ./ "$sra"_1.fasta "$sra"_2.fasta


# Filter the abundance file
awk -F'\t' 'NR==1 || $5>0 || $6>0' abundance.tsv > "$sra"_filtered_abundance.tsv


#deleting the fastas
rm "$sra"_1.fasta
rm "$sra"_2.fasta

