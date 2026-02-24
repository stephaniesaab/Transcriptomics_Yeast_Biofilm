#Assignment 2: Transcriptomics:
library("readr")
library("DESeq2")
library("tximport")
library(ggplot2)
library(org.Sc.sgd.db) #Bioconductor annotation data package

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

#Verify the match
print(sampleTable)

#Running DESeq2


