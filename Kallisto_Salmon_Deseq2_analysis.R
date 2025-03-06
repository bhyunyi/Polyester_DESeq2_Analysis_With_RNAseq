# initialize the type of immune cell to run analysis on
sample_name = "M0"

# install DESeq2 package if package has not already been installed
if (!requireNamespace("DESeq2", quietly = TRUE))
  install.packages("DESeq2", lib = "/scratch/by2372/R/libraries")

# load libaries
library(DESeq2)
library(ggplot2)

# Function for generation of count dataframes from raw data inputs
generate_count_data = function(raw_data){
  
  # Will change condition label based on if the data is from Kallisto or salmon.
  # Kallisto and Salmon have read count data in different columns
  raw_colnames = colnames(raw_data)
  if (raw_colnames[4] == "est_counts"){
    # fourth column in kallisto output has est_counts
    count_num = 4
    SRA_label = "_K"
  } else {
    # fifth column in salmon output has numreads
    count_num = 5
    SRA_label = "_S"
  }
  
  # creates a list of dataframes by splitting the raw data based on the SRA ID column
  data_split = split(raw_data,raw_data$SRA)
  
  # creates a lists of all SRR IDs in the lists of dataframes
  samples = c()
  for (x in 1:length(data_split)){
    
    SRA = data_split[[x]][1,8]
    SRA_sample = paste0(SRA,SRA_label)
    samples = append(samples,SRA_sample)
    
  }
  
  # creates a list of lists of gene names for each SRA ID
  gene_names_list = c()
  for(x in 1:length(data_split)){
    
    g_names = data_split[[x]][,6]
    g_names = list(g_names)
    gene_names_list = append(gene_names_list,g_names)
    
  }
  
  # creates of a list of genes that overlap with all SRA IDs
  final_gene_names = Reduce(intersect,gene_names_list)
  
  # creates a count data frame with genes as rows and SRA IDs as columns
  count_df = data.frame(matrix(NA,nrow = length(final_gene_names),ncol = length(samples)))
  rownames(count_df) = final_gene_names
  colnames(count_df) = samples
  
  # fills in count dataframe with raw read counts
  for (y in 1:nrow(count_df)){
    
    for (z in 1:ncol(count_df)){
      
      count_subset = subset(data_split[[z]], external_gene_name == rownames(count_df)[y], select = count_num)
      count_df[y,z] = count_subset[1,1]
    }
    
  }
  
  # function returns a count dataframe for DESeq2 analysis
  return(count_df)
  
}

# set working directory and import raw transcription factor file
setwd("/scratch/by2372/Kallisto_pseudoalign_BY/")
Kal_raw = read.csv("Kall_filtered_tfs.csv", header = TRUE, sep = ',')
# Only choose samples with common cell identity for DESeq2 analysis
Kal_cell_ID = Kal_raw[grepl(sample_name,Kal_raw$cell_identity),]


# set working directory and import raw transcription factor file
setwd("/scratch/by2372/Salmon_pseudoalign_BY/")
Sal_raw = read.csv("Sall_filtered_tfs.csv", header = TRUE, sep = ',')
# Only choose samples with common cell identity for DESeq2 analysis
Sal_cell_ID = Sal_raw[grepl(sample_name,Sal_raw$cell_identity),]

# Generate count data 
Kal_count_data = generate_count_data(Kal_cell_ID)
Sal_count_data = generate_count_data(Sal_cell_ID)

# find intersecting expressed genes between Kallisto and Salmon outputs 
Kal_count_rows = rownames(Kal_count_data)
Sal_count_rows = rownames(Sal_count_data)
common_names = Reduce(intersect,list(Kal_count_rows,Sal_count_rows))

# filter both out for common genes and then column bind the two count dataframes
Kal_filtered_count = Kal_count_data[common_names,]
Sal_filtered_count = Sal_count_data[common_names,]
Kal_Sal_count_data = cbind(Kal_filtered_count,Sal_filtered_count)

Kal_num_sample = ncol(Kal_filtered_count)
Sal_num_sample = ncol(Sal_filtered_count)

# create column data for DESeq2 analysis
Kal_Sal_col_data = data.frame(Sample = colnames(Kal_Sal_count_data), 
                              Condition = c(rep("Kallisto",Kal_num_sample),rep("Salmon",Sal_num_sample)))

### DESeq2 Analysis ###

# performs the DESeq2 analysis
dd_countdata = as.matrix(round(Kal_Sal_count_data))  
dds = DESeqDataSetFromMatrix(countData=dd_countdata, colData=Kal_Sal_col_data, design=~Condition)
dds = DESeq(dds)
results_dds = results(dds)
results_dds = merge(as.data.frame(results_dds), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
#results_dds = na.omit(results_dds)

# makes log2foldchange column numeric and creates negative log of padj column
results_dds$log2FoldChange = as.numeric(results_dds$log2FoldChange)
results_dds$neglog_padj = -log10(ifelse(is.na(results_dds$padj), 1e-10, results_dds$padj))
# add a column to the results to signify statistical significance
results_dds$significant = ifelse(results_dds$neglog_padj > -log10(0.05), "Significant", "Not Significant")

# Finds all genes with greater than 5 log2FoldChange and with a signficant adjusted p-value
plot_labels = results_dds[abs(results_dds$log2FoldChange) > 5 & results_dds$padj < 0.05,]

# Volcano plot with labels for significant genes with differential expression
volcano_plot = ggplot(results_dds, aes(x = log2FoldChange, y = neglog_padj, color = significant)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 5, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = -5, linetype = "dotted", color = "blue") +
  geom_text(data = plot_labels, aes(label = Row.names, vjust = -0.5), size = 2) +
  annotate("text", x = Inf, y = Inf, label = paste0(length(results_dds$Row.names)," genes present"), hjust = 1.1, vjust = 1.1, size = 3) +
  xlim(-10,10)+
  labs(title = paste0("Volcano Plot for ",sample_name," Samples"), x = "Log2 Fold change", y = "-log(adjusted p-value)", 
       color = "Significance of Differential Expression")

# Saves the volcano plot to a directory
setwd("/scratch/by2372/Deseq2_analysis/")
ggsave(paste0(sample_name,"_volcanoplot.png"), plot = volcano_plot, width = 10, height = 6, dpi = 300)
