#Assignment 2: Transcriptomics:
#Libraries ====
library("readr")
library("DESeq2")
library("tximport")
library(ggplot2)
library(org.Sc.sgd.db) #Bioconductor annotation data package for yeast
library(apeglm)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(tidyr)
library(pheatmap)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(AnnotationDbi)
library(scales)

#TODO: ====
# - Fix headings and labels of sections
# - Fix plot labels
# - Talk about the whole dataset
#Change all gene symbols to gene names in the figures
#Add figure for DEseQ2
#Use facetwrap from ggplot
#Lecture 8 slides -> Do plot of up v down regulated genes for ORA
#FIX COMMENTS

#Set directory as source file location

# Load files ====


#Point to quantification files from Salmon
dir <- "quants"
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
names(files) <- samples

#Checking all the files are there
files #Should output nine files

# Creating transcript-to-gene map ====
#With tx2gene and the Ensembl database
#Transcript ID usually is the systematic name (ORF)
#Pull the transcript ID from database to make sure they match
#Create a table with col1 (transcript ID) and col2 (gene ID) for Salmon to summarize counts to gene level

#Extract the map with the Bioconductor Genome wide annotation for Yeast (https://www.bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html

#Extract identifiers from an annotation database object (org.Sc.sgd.db)
#Keytype = ENSEMBL to keep it consistent with gene IDs(e.g. not Uniprot)
keys <- keys(org.Sc.sgd.db, keytype = "ENSEMBL")
#See first few gene symbols ("YIL064W"...)
head(keys)

#The Salmon quant names have _mRNA in the name -> need to add to IDs to match Salmon outputs
tx_names_mRNA <- paste0(keys, "_mRNA")

#Manually construct the 2-column table [TranscriptID, GeneID]
#For yeast: Transcript ID = Gene ID for most of them, so we can repeat the column 
tx2gene <- data.frame(TXNAME = tx_names_mRNA,
                      GENEID = keys)
#Check that it has two columns
dim(tx2gene) #5801 x 2 --> Aligns with Yeast genome size
head(tx2gene)
all(file.exists(files)) #Should be true

#Import files to aggregate transcript counts to gene level
#txIn = default true, whether the incoming files are transcript level
#txOut = default false, whether it should output transcript level
#tx2gene = needs column order to be transcript id, gene id
#ignoreTxVersion = split the tx id on the . to remove version info for matching with tx id

txi <- tximport(files, type= "salmon", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene, ignoreTxVersion = TRUE)

#Check that it worked
head(txi$counts) #Should see output for 9 samples


#Statistical comparison with DESeq2 for results
#Need to tell R which samples belong to Early, Med, Mature groups to get p-vals

#Create metadata table for SRA IDs
#Check the order R loaded the files in: 
print(colnames(txi$counts)) #Sorted numerically

#Defining the conditions in that by biofilm:
#SRR...57, 58, 59 = Mature
#SRR...60, 61, 62 = Thin
#SRR...63, 64, 65 = Early 

conditions <- c(rep("Mature", 3), 
                rep("Thin", 3),
                rep("Early", 3))

#Adding a replicates column for shapes or labels in plots
replicates <- factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3))
#Create Sample Table
sampleTable <- data.frame(condition = factor(conditions), 
                          replicate = replicates)

#Assign row names
rownames(sampleTable) <- colnames(txi$counts)

#Set reference level as Early (baseline)
sampleTable$condition <- relevel(sampleTable$condition, ref = "Early")

#Verify the match and make sure table is in tidy format (one sample per row)
print(sampleTable)

# Running DESeq2 =====

#Comparing all the groups pairwise:
levels(sampleTable$condition) #Check that it gives all three

#The ~condition tells DESeq2 to model the data based on Mature, thin, early
#Fitting a model log2(count) = Intercept(Early) + B1(mature vs early) + B2 (thin vs early)
dds <- DESeqDataSetFromTximport(txi, 
                                colData = sampleTable, 
                                design = ~condition)

#Running analysis: Using the default Wald test -> pairwise test between two groups
dds <- DESeq(dds)


#Check comparison names (Early is baseline so see: Mature vs Early and Thin vs Early)
resultsNames(dds) 

#Mature vs. Early
mature_v_early <- results(dds, contrast = c("condition", "Mature", "Early"))

#Mature vs. Thin --> DESeq2 handles with contrast (it's not a direct coefficient for the GLM because Early is reference)
mature_v_thin <- results(dds, contrast = c("condition", "Mature", "Thin"))

#Early vs. Thin
early_v_thin <- results(dds, contrast = c("condition", "Early", "Thin"))

#Re-leveling to get Mature vs. Thin as a coefficient for shrinkage
sampleTable$condition <- relevel(sampleTable$condition, ref = "Thin")
dds2 <- DESeqDataSetFromTximport(txi, 
                                colData = sampleTable, 
                                design = ~condition)

#Running analysis: Using the default Wald test -> pairwise test between two groups
dds2 <- DESeq(dds2)


#Check comparison names (Thin is baseline so see: Mature vs Thin and Early vs Thin)
resultsNames(dds2) 


# Shrinkage & Filtering ====
#to reduce noise from low-count or highly variable genes
#Apply shrinkage to each pair

MvE_resLFC <- lfcShrink(dds, coef="condition_Mature_vs_Early", type="apeglm")
TvE_resLFC <- lfcShrink(dds, coef="condition_Thin_vs_Early", type="apeglm")
MvT_resLFC <- lfcShrink(dds2, coef="condition_Mature_vs_Thin", type="apeglm")

#With shrinkage example plot
plotMA(MvE_resLFC, ylim=c(-2,2))


#Finding_top_X_genes
process_res <- function(resLFC_object, top_x = 10) {
  #Filter data for most significant 5 genes based on adjusted p-value for the three pairs
  #We need to mark genes as upregulated, downregulated, or not significant to colour them in ggplot
  #Exclude genes with under 2-fold change (log2FoldChange < 1)
  df <- as.data.frame(resLFC_object)
  df$gene <- rownames(df) #Set as the systematic ORF names
  
  #Remove NA rows
  df <- na.omit(df)
  
  #Add significance column based on adjusted p-values
  #Exclude genes with under 2-fold change (log2FoldChange < 1)
  df$significant <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1,
                           ifelse(df$log2FoldChange > 0, "Up", "Down"),
                           "Not Sig")
  
  #Map official gene names (e.g. FLO11, COX1)
  #Switch the key types to be standard gene name instead of systematic ORF name
  gene_name_map <- mapIds(org.Sc.sgd.db,
                          keys = df$gene,
                          column = "COMMON",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  df$gene_name <- gene_name_map #Map the symbols to their common name
  df$gene_label <- ifelse(is.na(df$gene_name), #If the gene name didn't map just keep the symbol
                         df$gene,
                         df$gene_name)
  
  #Identify top N genes (sorting by p value)
  top_genes <- df %>% 
    filter(significant != "Not Sig") %>% 
    slice_min(order_by = pvalue, n = top_x)
  return(list(full = df, top = top_genes))
}

MvE <- process_res(MvE_resLFC)
MvT <- process_res(MvT_resLFC)
TvE <- process_res(TvE_resLFC)

# Plots ====
## PCA Plot ----

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)
# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display conditions by colour, and shape
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") + 
  theme_minimal()
## Volcano Plots ----
# We plot log2foldchange against -log10(adjusted p value)

#Combining the three datasets into one master
combined_df <- bind_rows(
  MvE$full %>% mutate(comparison = "Mature vs. Early"),
  MvT$full %>% mutate(comparison = "Mature vs. Thin"),
  TvE$full %>% mutate(comparison = "Thin vs. Early")
)

plot_volcano <- function(df) {
  
  ggplot(df, aes(x = log2FoldChange,
                 y = -log10(padj), #Use padjusted 
                 color = significant)) +
    geom_point(alpha = 0.4) + #Labels easier to read if points are transparent
    scale_color_manual(values = c("Down" = "blue", 
                                  "Not Sig" = "gray", 
                                  "Up" = "red")) +
    
    #Facet wrap splits the plot into 3 based on the comparison
    facet_wrap(~comparison, nrow = 1)+
    
    #Add labels using the new top_x subset
    geom_text_repel(data = df %>% 
                      group_by(comparison) %>% 
                      filter( significant != "Not Sig")%>% 
                      slice_min(order_by = padj, n = 10), #Use padjusted over pvalue
                    aes(label = gene_label),
                    color = "black",
                    size = 3.5, 
                    max.overlaps = Inf,
                    box.padding = 0.5) +
    
    labs(x = "Log2 Fold Change", 
         y = "-Log10 p-value", 
         title = "Differential Gene Expression analysis of top 10 genes") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines") #Add space between plots
    )
}
plot_volcano(combined_df)

#GO and ORA ====
## Map IDs to get ENTREZID ====
gene_map <- bitr(combined_df$gene,
                 fromType = "ORF", 
                 toType = "ENTREZID",
                 OrgDb = org.Sc.sgd.db) %>% 
  distinct(ORF, .keep_all = TRUE)
combined_df <- combined_df %>% left_join(gene_map, by = c("gene" = "ORF"))

## Create named list for compareCluster
gene_list <- list(
  "MvE Up"   = combined_df %>% filter(comparison == "Mature vs. Early", padj < 0.05, log2FoldChange > 1) %>% pull(ENTREZID) %>% na.omit() %>% as.character(),
  "MvE Down" = combined_df %>% filter(comparison == "Mature vs. Early", padj < 0.05, log2FoldChange < -1) %>% pull(ENTREZID) %>% na.omit() %>% as.character(),
  "MvT Up"   = combined_df %>% filter(comparison == "Mature vs. Thin", padj < 0.05, log2FoldChange > 1) %>% pull(ENTREZID) %>% na.omit() %>% as.character(),
  "MvT Down" = combined_df %>% filter(comparison == "Mature vs. Thin", padj < 0.05, log2FoldChange < -1) %>% pull(ENTREZID) %>% na.omit() %>% as.character(),
  "TvE Up"   = combined_df %>% filter(comparison == "Thin vs. Early", padj < 0.05, log2FoldChange > 1) %>% pull(ENTREZID) %>% na.omit() %>% as.character(),
  "TvE Down" = combined_df %>% filter(comparison == "Thin vs. Early", padj < 0.05, log2FoldChange < -1) %>% pull(ENTREZID) %>% na.omit() %>% as.character()
)

lengths(gene_list) #Should be more than zero for each list

# Define our background list of genes to compare to
# ORA needs a universe set
all_genes <- as.character(na.omit(unique(combined_df$ENTREZID))) #Length = 5721

## GO enrichment ====
ck_go <- compareCluster(geneCluster = gene_list,
                        fun = "enrichGO",
                        universe = all_genes,
                        OrgDb = org.Sc.sgd.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH", #Default multiple testing algorithm
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,#FDR control
                        readable = FALSE
                        )

head(as.data.frame(ck_go)) #Check there are > 0 in counts

#Run for KEGG ----
ck_kegg <- compareCluster(geneCluster = gene_list,
                          fun = "enrichKEGG",
                          universe = all_genes,
                          organism = "sce", #Yeast organism
                          keyType = "ncbi-geneid",
                          pAdjustMethod = "BH", #Default
                          pvalueCutoff = 0.5,
                          qvalueCutoff = 0.2)

## Make the plot with three panes ----

#Split the column into comparison and direction
ck_df <- as.data.frame(ck_kegg) %>% 
  separate(Cluster, into = c("Comparison", "Direction"), sep = " ") %>% 
  separate(GeneRatio, into = c("num", "den"), sep = "/", remove = FALSE) %>%
  mutate(GeneRatioNum = as.numeric(num) / as.numeric(den)) %>% 
  group_by(Comparison, Direction) %>% 
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>% 
  ungroup()

ggplot(ck_df, aes(x = Direction, y = Description)) +
  
  #Use GeneRatioNum for the size
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  scale_size(range = c(3, 8), breaks = breaks_pretty(n = 5), name = "Gene Ratio") + #Make the gene ratio legend less
  # This creates the 3-facet layout
  facet_wrap(~Comparison, scales = "free_y") + 
  theme_bw() +
  labs(title = "KEGG Pathway Enrichment: Up vs Down",
       x = "Regulation Direction", 
       y = "Pathway",
       size = "Gene Ratio") +
  theme(axis.text.x = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))




#ROUGH DRAFT CODE ====
#Define Up and Down gene sets
# We would need to split those out ourselves as our gene set of interest
# Here's an example for upregulation:
#BP = biological process
plot_combined_go <- function(df) {
  
  #Mapped IDs
  gene_map <- bitr(df$gene,
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID",
                   OrgDb = org.Sc.sgd.db) %>% 
              distinct(ENSEMBL, .keep_all = TRUE)
  df <- df %>% left_join(gene_map, by = c("gene" = "ENSEMBL"))
  
  #Define background universe (all genes)
  all_genes <- as.character(na.omit(unique(df$ENTREZID)))
  
  #Helper to process each pairwise comparison
  run_ora <- function(sub_df) {
    comp_name <- unique(sub_df$comparison)
  
    #Split up and downreg genes
    up_ids <- sub_df %>% filter(padj < 0.05, log2FoldChange > 1) %>% 
      pull(ENTREZID) %>% na.omit()
    dn_ids <- sub_df %>% filter(padj < 0.05, log2FoldChange < -1) %>% 
      pull(ENTREZID) %>% na.omit()
    
    #Run GO enrichment
    ego_bp_up <- enrichGO(gene = up_ids,
                          universe = all_genes,
                          OrgDb = org.Sc.sgd.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = FALSE)
    ego_bp_dn <- enrichGO(gene = dn_ids,
                          universe = all_genes,
                          OrgDb = org.Sc.sgd.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = FALSE)
    
    #Combine into a tidy dataframe
    res_up <- as.data.frame(ego_bp_up) %>% mutate(Direction = "Upregulated", Comparison = comp_name)
    res_dn <- as.data.frame(ego_bp_db) %>% mutate(Direction = "Downregulated", Comparison = comp_name)
    
    return(bind_rows(res_up, res_dn))
  }
  
  #Apply helper to each comparison pair and combine
  plot_data <- df %>% 
    group_split(comparison) %>% 
    map_dfr(run_ora)
  
  #Clean up for plotting
  plot_data <- plot_data %>% 
    separate(GeneRatio, into = c("num", "den"), sep = "/") %>% 
    mutate(GeneRatioNum = as.numeric(num) / as.numeric(den)) %>% 
    group_by(Comparison, Direction) %>% 
    slice_min(p.adjust, n = 10, with_ties = FALSE) %>% #Keep top 10 terms
    ungroup()
  
  #Make the faceted plot
  ggplot(plot_data, aes(x = Direction, y = Description)) +
    g
} 
#Already done the DESeq2 analysis, and that our results data is comparing all pairwise combos:
# Make a dataframe containing only our results table post-shrinkage, and with NA values trimmed
# Convert Ensembl IDs to Entrez IDs
# Remove version numbers (e.g., .9 from ENSG00000189221.9)
# Map to Entrez IDs
make_gene_map <- function(df) {
  gene_map <- bitr(df$gene,
                   fromType = "ENSEMBL", 
                   toType = c("ENTREZID", "COMMON", "SGD"),
                   OrgDb = org.Sc.sgd.db) %>% 
                distinct(ENSEMBL, .keep_all = TRUE)
  return(gene_map)
}
#ROUGH DRAFT CODE ====
gene_map <- bitr(res_df$gene, 
                 fromType = "ENSEMBL", 
                 toType = c("ENTREZID", "COMMON", "SGD"),
                 OrgDb = org.Sc.sgd.db) %>% 
                distinct(ENSEMBL, .keep_all = TRUE)

# Add the mapping to results, join based on systematic names
res_df <- res_df %>% 
  left_join(gene_map, by = c("gene" = "ENSEMBL"))

head(res_df)

# Convert to character before pulling the list
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  mutate(ENTREZID= as.character(ENTREZID)) %>% # Force character type
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

# Define our background list of genes to compare to
# Rembember that ORA needs an "interesting gene" set and a background set.
all_genes <- res_df %>%
  pull(ENTREZID) %>% #CHECK THIS ====
  na.omit() %>%
  unique() %>% 
  as.character()

# Here we'll do a GO analysis for only Biological Process #Most useful for analysis
# (Try Molecular Function or Cellular Component instead and see what you get!)
ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH", #CHECK THIS ====
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE) #So it doesn't convert results into SYMBOL column (keep COMMON or GENENAME for Yeast database)
head(as.data.frame(ego_bp))

# Now we'll do a KEGG analysis
#Keep keytypes to default (Kegg which works for EntrezIDs)
sig_genes_kegg <- res_df %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  pull(gene) %>% 
  unique()
kegg_enrich <- enrichKEGG(gene = sig_genes_kegg,
                          organism = 'sce',
                          keyType = 'kegg',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

#Check for results: 
head(as.data.frame(kegg_enrich))

# Dot plots ====
dotplot(ego_bp, showCategory = 20, title = "GO Biological Process")

dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")


# Bar plot ====
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")

barplot(kegg_enrich, showCategory = 15, title = "KEGG Biological Process")

# Enrichment map (Displays linked GO terms)
emapplot(pairwise_termsim(ego_bp), showCategory = 30)

# Notice that none of these plots split our genes by upregulated/downregulated?
# We would need to split those out ourselves as our gene set of interest
# Here's an example for upregulation:
#BP = biological process

upregulated_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

downregulated_genes <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange < -1) %>% 
  pull(ENTREZID) %>% 
  na.omit() %>% 
  unique()

ego_bp_up <- enrichGO(gene = upregulated_genes,
                      universe = all_genes,
                      OrgDb = org.Sc.sgd.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = FALSE)

ego_bp_down <- enrichGO(gene = downregulated_genes,
                        universe = all_genes,
                        OrgDb = org.Sc.sgd.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = FALSE)

# Plots of upregulated genes:
# Important note! Upregulated genes aren't the same thing as upregulated biological effects.
# What would happen if an important repressor of angiogenesis was upregulated? 
# Would it show up in upregulation or downregulation? What would happen to angiogenesis?
# Most of the time a gene related to a process is an enhancement but for some biological systems the gene related can be a negative regulator so you have to look carefully as to whether that ene is actually upregulating the process or not
# Normally you say regulation when discussing the genes and use "activated or repressed" for the physiological process
dotplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

dotplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")



barplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

#Should this be kegg_enrich or kegg_up?
barplot(kegg_up, showCategory = 15, title = "KEGG - Upregulated Genes")

heatplot(ego_bp_down, foldChange = ego_bp_down, showCategory = 20)


p1 <- ggplot(ego_bp_up, aes(x, y, fill = zscore)) + geom_tile() + labs(title = "Upregulated")
p2 <- ggplot(ego_bp_down, aes(x, y, fill = value)) + geom_tile() + labs(title = "Downregulated")

plot(p1)


### Preparing data ----

MvE_res_df <- as.data.frame(MvE_resLFC)
MvE_res_df$gene <- rownames(MvE_res_df)
MvE_res_df$significant <- ifelse(MvE_res_df$padj < 0.05 & abs(MvE_res_df$log2FoldChange) > 1, 
                                 ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
MvE_res_df <- na.omit(res_df)

#Switch the key types to be standard gene name instead of systematic ORF name
gene_name_map <- mapIds(org.Sc.sgd.db,
                        keys = res_df$gene,
                        column = "COMMON",
                        keytype = "ENSEMBL",
                        multiVals = "first")
res_df$gene_name <- gene_name_map

### Identify top N genes (sorting by p value) ----
top_x <- 5
top_genes <- res_df %>% 
  filter(significant != "Not Sig") %>% 
  slice_min(order_by = pvalue, n = top_x) %>% 
  as.data.frame()

#Heatmap ====

# Remove NA values
resLFC <- na.omit(resLFC)

# Select top 20 genes by p values from non-NA genes (could also sort by logfoldchange)
top_x <- 20
top_genes <- head(order(abs(resLFC$padj), decreasing = FALSE), top_x) #EntrezIDs
gene_names <- rownames(resLFC)[top_genes] #ORF Gene names

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)
# Store counts in a matrix for the heatmap
mat <- assay(vsd)[gene_names, ]

# Add annotation for condition and replicates
annotation_df <- sampleTable[, c("condition", "replicate")]
colnames(annotation_df) <- c("Condition", "Replicate")

# Create heatmap ---> Write the gene names legitametaly ---> ==== use gene symbols
pheatmap(mat, 
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = FALSE,
)

#PCA Plot ----

#Use same VST as heatmap
# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display conditions by colour, and shape
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples")
coord_fixed()
