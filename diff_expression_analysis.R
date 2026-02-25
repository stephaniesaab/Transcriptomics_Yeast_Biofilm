#Assignment 2: Transcriptomics:
library("readr")
library("DESeq2")
library("tximport")
library(ggplot2)
library(org.Sc.sgd.db) #Bioconductor annotation data package for yeast
library(apeglm)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(enrichplot)
library(DOSE)
library(clusterProfiler)

#TODO: ====
# - Fix headings and labels of sections
# - Fix plot labels
# - Talk about the whole dataset


#Set directory as source file location

#Load files


#Point to quantification files from Salmon
dir <- "quants"
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
names(files) <- samples

#Checking all the files are there
files #Should output nine files

#Creating transcript-to-gene map with tx2gene and the Ensembl database
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
## CHECK THIS ====
tx2gene <- data.frame(TXNAME = tx_names_mRNA,
                      GENEID = keys)
#Check that it has two columns
dim(tx2gene)
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

#Running DESeq2 =====

#The ~condition tells DESeq2 to model the data based on Mature, thin, early
dds <- DESeqDataSetFromTximport(txi, 
                                colData = sampleTable, 
                                design = ~condition)

#Running analysis: performing a likelihood ratio test
dds <- DESeq(dds)


#Check comparison names ==== #WHy does it say intercept and why not early vs thin?
resultsNames(dds) 

#Comparing groups pairwise:

#Mature vs. Early
mature_v_early <- results(dds, contrast = c("condition", "Mature", "Early"))

#Mature vs. Thin
mature_v_thin <- results(dds, contrast = c("condition", "Mature", "Thin"))

#Early vs. Thin
early_v_thin <- results(dds, contrast = c("condition", "Early", "Thin"))

#Check comparison names
resultsNames(dds)


#Run the Likelihood ratio test (LRT) to find genes that just differ among the conditions
dds_lrt <- DESeq(dds, test = "LRT",
                 reduced = ~1)
#Get the results
res_any_diff <- results(dds_lrt)

resultsNames(dds)



#Bioinformatics trinity: Visualization (MA/Volcano), clustering (heatmap), dimensionality reduction (PCA)

#Visualization ====

#Shrinkageto reduce noise from low-count or highly variable genes
resLFC <- lfcShrink(dds, coef="condition_Mature_vs_Early", type="apeglm")

#Filter data for most significant 10 genes

#MA plot without shrinkage
plotMA(res_any_diff, ylim=c(-2,2))
plotMA(mature_v_early, ylim=c(-2, 2))

#With shrinkage
plotMA(resLFC, ylim=c(-2,2))

#Volcano plot ====

### Preparing data ----
#We need to mark genes as upregulated, downregulated, or not significant to colour them in ggplot
#We'll also exclude genes with under 2-fold change (log2FoldChange < 1)
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                             ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
res_df <- na.omit(res_df)

### Identify top N genes (sorting by p value) ----
top_x <- 10
top_genes <- res_df %>% 
  filter(significant != "Not Sig") %>% 
  slice_min(order_by = pvalue, n = top_x) %>% 
  as.data.frame()

### Plot ----
# Here's some ggplot code for a volcano plot
# We plot log2foldchange against -log10(adjusted p value)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.4) + #Labels easier to read if points are transparent
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  
  #Add labels using the new top_x subset
  geom_text_repel(data = top_genes, 
                  aes(label = gene),
                  color = "black",
                  size = 3.5, 
                  max.overlaps = Inf,
                  box.padding = 0.5) +
  
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
       title = "Volcano Plot: Mature vs Early") +
  theme_minimal() +
  theme(legend.position = "right")

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

# Add annotation for cell lines and treatment
annotation_df <- sampleTable[, c("condition", "replicate")]
colnames(annotation_df) <- c("Condition", "Replicate")

# Create heatmap
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


#Tutorial 8 ====

#Already done the DESeq2 analysis, and that our results data is comparing mature vs. early (and all the other combos):
# Make a dataframe containing only our results table post-shrinkage, and with NA values trimmed

# Convert Ensembl IDs to Entrez IDs
# Remove version numbers (e.g., .9 from ENSG00000189221.9)

# Map to Entrez IDs
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


#====== tutorial 7 code =====

data2 -> de noising algorithm
#TALK about the whole dataset ====

#BP = biological processses --> Good one for the three categories for gene ontology 