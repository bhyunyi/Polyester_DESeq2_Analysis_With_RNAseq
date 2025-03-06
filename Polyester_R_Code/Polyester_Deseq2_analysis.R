# install DESeq2 package if package has not already been installed
if (!requireNamespace("DESeq2", quietly = TRUE))
  install.packages("DESeq2", lib = "/path/R/libraries")

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
  data_split = split(raw_data,raw_data$sample)
  
  # creates a lists of all SRR IDs in the lists of dataframes
  samples = c()
  for (x in 1:length(data_split)){
    
    SRA = data_split[[x]][1,6]
    SRA_sample = paste0(SRA,SRA_label)
    samples = append(samples,SRA_sample)
    
  }
  
  # creates a list of lists of gene names for each SRA ID
  gene_names_list = c()
  for(x in 1:length(data_split)){
    
    g_names = data_split[[x]][,7]
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

load("/path/simulated_reads/sim_counts_matrix.rda")

# set working directory and import raw transcription factor file
setwd("/path/Kallisto_pseudoalign/")
Kal_raw = read.csv("Kall_filtered.csv", header = TRUE, sep = ',')
Kal_raw$external_gene_name = Kal_raw$target_id

# set working directory and import raw transcription factor file
setwd("/path/Salmon_pseudoalign/")
Sal_raw = read.csv("Sall_filtered.csv", header = TRUE, sep = ',')
Sal_raw$external_gene_name = Sal_raw$Name

# Generate count data 
Kal_count_data = generate_count_data(Kal_raw)
Sal_count_data = generate_count_data(Sal_raw)

# Alter the rownames from the counts_matrix
counts_rownames = rownames(counts_matrix)
new_counts_rows = c()
for (x in counts_rownames){
  x_split = strsplit(x, "\\|")[[1]]
  new_string = paste0(x_split[1:4], collapse = "|")
  new_string = paste0(new_string,"|")
  new_counts_rows = append(new_counts_rows,new_string)
}

rownames(counts_matrix) = new_counts_rows

# find intersecting expressed genes 
Kal_count_rows = rownames(Kal_count_data)
common_kallisto = Reduce(intersect,list(Kal_count_rows,rownames(counts_matrix)))
Sal_count_rows = rownames(Sal_count_data)
common_Salmon = Reduce(intersect,list(Sal_count_rows,rownames(counts_matrix)))

# filter both out for common genes and then column bind the two count dataframes
Kal_filtered_count = Kal_count_data[common_kallisto,]
common_K_filtered = counts_matrix[common_kallisto,]
final_Kal_count = cbind(Kal_filtered_count,common_K_filtered)

Sal_filtered_count = Sal_count_data[common_Salmon,]
common_S_filtered = counts_matrix[common_Salmon,]
final_Sal_count = cbind(Sal_filtered_count,common_S_filtered)

# create column data for DESeq2 analysis Kallisto
Kal_col_data = data.frame(Sample = colnames(final_Kal_count), 
                          Condition = c(rep("Kallisto",20),rep("Actual",20)))

Sal_col_data = data.frame(Sample = colnames(final_Sal_count), 
                          Condition = c(rep("Salmon",20),rep("Actual",20)))

### DESeq2 Analysis Kallisto ###

dd_countdata_K = as.matrix(round(final_Kal_count))  

# perform DESeq2 analysis and get the results
dds_K = DESeqDataSetFromMatrix(countData=dd_countdata_K, colData=Kal_col_data, design=~Condition)
dds_K = DESeq(dds_K)
results_dds_K = results(dds_K)
results_dds_K = merge(as.data.frame(results_dds_K), as.data.frame(counts(dds_K, normalized=TRUE)), by="row.names", sort=FALSE)

results_dds_K$log2FoldChange = as.numeric(results_dds_K$log2FoldChange)
results_dds_K$neglog_padj = -log10(ifelse(is.na(results_dds_K$padj), 1e-10, results_dds_K$padj))
# add a column to the results to signify statistical significance
results_dds_K$significant = ifelse(results_dds_K$neglog_padj > -log10(0.05), "Significant", "Not Significant")

# Finds all genes with greater than 5 log2FoldChange and with a signficant adjusted p-value
plot_labels_K = results_dds_K[abs(results_dds_K$log2FoldChange) > 5 & results_dds_K$padj < 0.05,]

# Volcano plot with labels for significant genes with differential expression
volcano_plot_K = ggplot(results_dds_K, aes(x = log2FoldChange, y = neglog_padj, color = significant)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 5, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = -5, linetype = "dotted", color = "blue") +
  geom_text(data = plot_labels_K, aes(label = Row.names, vjust = -0.5), size = 2) +
  annotate("text", x = Inf, y = Inf, label = paste0(length(results_dds_K$Row.names)," genes present"), hjust = 1.1, vjust = 1.1, size = 3) +
  xlim(-10,10)+
  labs(title = paste0("Volcano Plot for Kallisto Estmates vs. Actual Abundances"), x = "Log2 Fold change", y = "-log(adjusted p-value)", 
       color = "Significance of Differential Expression")

# Save the volcano plot to a directory
setwd("/path/Volcano_Plots/")
ggsave("Polyester_Kallisto_volcanoplot.png", plot = volcano_plot_K, width = 10, height = 6, dpi = 300)


### DESeq2 Analysis Salmon ###

dd_countdata_S = as.matrix(round(final_Sal_count))  
# perform DESeq2 analysis and get the results
dds_S = DESeqDataSetFromMatrix(countData=dd_countdata_S, colData=Sal_col_data, design=~Condition)
dds_S = DESeq(dds_S)
results_dds_S = results(dds_S)
results_dds_S = merge(as.data.frame(results_dds_S), as.data.frame(counts(dds_S, normalized=TRUE)), by="row.names", sort=FALSE)

results_dds_S$log2FoldChange = as.numeric(results_dds_S$log2FoldChange)
results_dds_S$neglog_padj = -log10(ifelse(is.na(results_dds_S$padj), 1e-10, results_dds_S$padj))
# add a column to the results to signify statistical significance
results_dds_S$significant = ifelse(results_dds_S$neglog_padj > -log10(0.05), "Significant", "Not Significant")

# Finds all genes with greater than 5 log2FoldChange and with a signficant adjusted p-value
plot_labels_S = results_dds_S[abs(results_dds_S$log2FoldChange) > 5 & results_dds_S$padj < 0.05,]

# Volcano plot with labels for significant genes with differential expression
volcano_plot_S = ggplot(results_dds_S, aes(x = log2FoldChange, y = neglog_padj, color = significant)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 5, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = -5, linetype = "dotted", color = "blue") +
  geom_text(data = plot_labels_S, aes(label = Row.names, vjust = -0.5), size = 2) +
  annotate("text", x = Inf, y = Inf, label = paste0(length(results_dds_S$Row.names)," genes present"), hjust = 1.1, vjust = 1.1, size = 3) +
  xlim(-10,10)+
  labs(title = paste0("Volcano Plot for Salmon Estmates vs. Actual Abundances"), x = "Log2 Fold change", y = "-log(adjusted p-value)", 
       color = "Significance of Differential Expression")
# Save the volcano plot to a directory
setwd("/scratch/by2372/Deseq2_analysis/Plots/")
ggsave("Polyester_Salmon_volcanoplot.png", plot = volcano_plot_S, width = 10, height = 6, dpi = 300)

