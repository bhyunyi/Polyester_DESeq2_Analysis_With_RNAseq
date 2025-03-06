#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=Sall4one
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

cd /scratch/by2372/Deseq2_analysis/Polyester_Test/Salmon_pseudoalign
mkdir ${sra}
cd /scratch/by2372/Deseq2_analysis/Polyester_Test/Salmon_pseudoalign/${sra}

cp /scratch/by2372/Deseq2_analysis/Polyester_Test/simulated_reads/${sra}_1.fasta .
cp /scratch/by2372/Deseq2_analysis/Polyester_Test/simulated_reads/${sra}_2.fasta .

#purge all modules and load salmon
module purge
module load salmon/1.10.0

#run the salmon pseudoalignment
salmon quant -i /scratch/by2372/Deseq2_analysis/Polyester_Test/Salmon_index -l A -1 ${sra}_1.fasta -2 ${sra}_2.fasta -p 4 --validateMappings -o quants

#purging modules cuz otherwise awk fails due to salmon
module purge

# Filter the abundance file
awk -F'\t' 'NR==1 || $4>0 || $5>0' quants/quant.sf > "$sra"_filtered_abundance.tsv


#deleting the fastqs
rm "$sra"_1.fasta
rm "$sra"_2.fasta
