# RNA-seq data analysis

# 1. Getting the data ----
#* git clone https://github.com/guigolab/rnaseq-course

# 2. Gene annotation ----
#* cd rnaseq-course
#* zless refs/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz
#* zgrep CSDE1 refs/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz > CSDE1.gff
#* awk '$3=="transcript"' CSDE1.gff | wc -l
#*   18
#* awk '$3=="exon"' CSDE1.gff | wc -l
#*   189

#* What is the biotype of the gene?
#*   gene_type "protein_coding"
#* What is the biotype of each transcript of this gene?
#*   transcript_type "protein_coding"

#* Advanced: Extract gene identifiers
#*   zcat refs/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz | awk -F"\t" '$3=="gene" && $9~/gene_type "protein_coding"/' | grep -Eo 'ENS[A-Z]+[0-9]{11}[.]?[0-9]*' | sort -u > protein_coding_IDs.txt

# 3. RNA-seq data analysis ----
setwd('rnaseq-course/analysis/')

# 3.1 Analysis setup ----
.libPaths("/software/rg/rnaseq/rpackages_2025/")

library(scales)
library(dplyr)
library(magrittr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vroom)
library(edgeR)
library(pheatmap)

# 3.2 Data preparation ----

# 3.2.1 Data import ----
raw_counts <- vroom("../quantification/raw_counts.tsv")
raw_counts

meta <- vroom("../quantification/metadata.tsv")

gene_annotation <- vroom("../refs/gene_annotation.tsv",
                         col_names = c(
                           "chrom", "start", "end", "strand",
                           "gene_id", "gene_type", "gene_name"
                         )
)

#* What are the dimensions of the metadata and gene annotation tables?
dim(meta)
#*   [1] 12 12

dim(gene_annotation)
#*   [1] 58780     7

raw_counts <-
  raw_counts %>%
  column_to_rownames("gene_id") %>%
  `[`(, meta$SampleID)

identical(colnames(raw_counts), meta$SampleID)

# 3.2.2 Data filtering ----
cpms <-
  raw_counts %>%
  as.matrix() %>%
  cpm()

#* Can you filter the raw_counts table to keep the genes that have more than 1 cpm in all our samples?
cpm_filter <-
  cpms %>%
  # Test which values pass the filter of 1 cpm
  `>`(1) %>%
  # Sum up in how many samples each gene passes the filter
  rowSums() %>%
  # Test which genes pass the filter for all 12 samples
  `==`(12)

raw_counts_filtered <-
  raw_counts[cpm_filter, ]

# An alternative
## raw_counts_filtered <- cpms[apply(cpms, 1, function(a) sum(a > 1) == 12), ]

# 3.2.3 Data Normalization ----
dge <-
  DGEList(counts = raw_counts_filtered,
          group = meta$Treatment)

print(dge)

dge <- calcNormFactors(object = dge, method = "TMM")

print(dge)

# 3.2.4 FPKM matrix ----

#* Can you use the information contained in the gene_annotation table to calculate the length of all the genes included in the dge object?
y_genes <-
  # Obtain ENSID of genes in PAR in chrY
  gene_annotation %>%
  filter(gene_id %in% rownames(dge$counts)) %>%
  group_by(gene_id) %>%
  summarise(n = n()) %>%
  filter(n == 2) %>%
  pull(gene_id)

# Get gene lengths
gene_lengths <-
  gene_annotation %>%
  # filter for genes in our dataset
  filter(gene_id %in% rownames(dge$counts) &
           # filter out genes in chrY,
           # whose chrX counterpart is in our dataset
           !(gene_id %in% y_genes & chrom == "chrY")) %>%
  # Order genes the same way they are ordered in our expression matrix
  mutate(order = match(gene_id, rownames(dge$counts))) %>%
  arrange(order) %>%
  # Calculate gene length
  mutate(length = end - start + 1) %>%
  pull(length)

#* Check the example solution above. Then, try to answer the following questions:

##* Why do we have +1 in the length calculation formula?
##*   Because these coordinates are 1-based

##* The length we used for the normalization is based on the genomic coordinates of the genes. Is this the right length to use?
##*   We should be using exon length

##* What does the y_genes object contain, and how do we use it for filtering? Tip: check the output of the following command: grep("_", rownames(raw_counts), value = TRUE)
##*   Genes in the pseudoautosomal regions of chromosome Y

fpkms <- rpkm(dge,
              log = TRUE,
              gene.length = gene_lengths
)

#* What is the base of the logarithm used by the rpkm function?
#*   ?rpkm --> log logical, if TRUE then log2 values are returned.

#* Does the rpkm function use a pseudocount? If it does, what is the value it assigns to it?
#*   prior.count = 2

# 3.3 Clustering and PCA ----

# 3.3.1 Clustering Analysis ----
fpkms.cor <- cor(fpkms)

head(fpkms.cor)

#* Which are the dimensions of the correlation matrix? How is this size determined?
dim(fpkms.cor)
#*   [1] 12 12 --> number of samples

fpkms.dist <- dist(x = fpkms.cor)

fpkms.hclust <- hclust(d = fpkms.dist)

pdf('hclust.pdf')
plot(fpkms.hclust)
dev.off()

#* How many clustering methods are available in the hclust function? Which is the default method?
#* ?hclust --> 8 methods, default is 'complete'

forAnnotation <-
  meta %>%
  dplyr::select(SampleID, Treatment_Duration) %>%
  column_to_rownames("SampleID")

pdf('heatmap.pdf')
pheatmap(
  mat = fpkms.cor,
  annotation_col = forAnnotation,
  annotation_row = forAnnotation
)
dev.off()

#* The default method for calculating the correlation in the cor function is pearson correlation. How does the output change if we use spearman correlation instead?
fpkms.cor <- cor(fpkms, method = 'spearman')

fpkms.dist <- dist(x = fpkms.cor)

fpkms.hclust <- hclust(d = fpkms.dist)

pdf('hclust_spearman.pdf')
plot(fpkms.hclust)
dev.off()

pdf('heatmap_spearman.pdf')
pheatmap(
  mat = fpkms.cor,
  annotation_col = forAnnotation,
  annotation_row = forAnnotation
)
dev.off()

#* Does the heatmap show what you would expect?
#*   Sample ENCFF443WJB should be nearer to the rest of 4h samples

# 3.3.2 Principal Component Analysis ----
pca <-
  fpkms %>%
  t() %>%
  prcomp(center = TRUE, scale. = TRUE)

#* What do the center and scale options do?
#*   center: shifted to be zero centered; 
#*   scale: scaled to have unit variance

summary(pca)

#* How many principal components explain 90% of the variation in our data?
#*   8 components

#* Can you make a scatter plot with the values of the first two principal components?
pdf('pca.pdf')
pca$x %>%
  # Change matrix to a data.frame
  as.data.frame() %>%
  # Add column with annotation about treatment duration
  cbind(Treatment_Duration = meta$Treatment_Duration) %>%
  # Make base ggplot
  ggplot(aes(PC1, PC2,
             col = Treatment_Duration,
             fill = Treatment_Duration
  )) +
  # Add ellipse layer
  stat_ellipse(
    geom = "polygon", level = 0.80,
    col = "black", alpha = 0.5
  ) +
  # Add scatter plot layer
  geom_point(shape = 21, col = "black")
dev.off()

#* Does the output of the PCA agree with that of the clustering analysis?
#*   Yes, the outlier from 4h remains there

# 3.4 Differential gene expression ----

# 3.4.1 Simple Comparison between two groups ----
dge <- estimateDisp(dge)

pdf('BCV.pdf')
plotBCV(dge)
dev.off()

#* Does the observed trend of the dispersion make sense?
#*   Yes, as lower CPM display higher dispersion

dge.result <-
  exactTest(object = dge, pair = c("control", "dexamethasone"))

summary(decideTests(object = dge.result, p.value = 0.01))
#       dexamethasone-control
#Down                     408
#NotSig                 11260
#Up                       592


#* What is FDR correction, why is it needed? Repeat the DGE analysis and set FDR to 5%. Count the number of hits.
#*   It corrects type I error (rejecting a true null hypothesis), i.e, finding more ositive results than the actual ones
#*   an FDR value of 0.05 means that 5% of “declared” positive results are truly negative
summary(decideTests(object = dge.result, p.value = 0.05))
#       dexamethasone-control
#Down                     671
#NotSig                 10754
#Up                       835

pdf('MAplot.pdf')
plotMD(dge.result)
dev.off()

# 3.4.2 Multiple comparisons ----
# Change group
dge <-
  DGEList(counts = raw_counts_filtered,
          group = meta$Treatment_Duration)

# Calculate again normalization factors
dge <- calcNormFactors(object = dge, method = "TMM")

dge.design <-
  model.matrix(object = ~ 0 + Treatment_Duration, data = meta)
dge.design

dge <- estimateDisp(y = dge, design = dge.design)

# Create contrasts for possible comparisons
dge.contrasts <-
  makeContrasts(
    "2vs0" = Treatment_Duration2hr - Treatment_Duration0hr,
    "4vs0" = Treatment_Duration4hr - Treatment_Duration0hr,
    "4vs2" = Treatment_Duration4hr - Treatment_Duration2hr,
    "TreatmentVsControl" =
      (Treatment_Duration4hr + Treatment_Duration2hr) / 2 -
      Treatment_Duration0hr,
    levels = dge.design
  )
dge.contrasts

#* Do you understand which is the comparison indicated by the last column of the contrast?
#*   The mean of treatments vs the control

dge.fit <- glmQLFit(dge, dge.design)

dge.result <-
  glmQLFTest(dge.fit,
             contrast = dge.contrasts[, c("2vs0", "4vs0")]
  )

summary(decideTests(object = dge.result, p.value = 0.01))
#       LR test on 2 degrees of freedom
#NotSig                           11731
#Sig                                529

topTags(object = dge.result)

allDeg <-
  # Get table of differential expression results
  topTags(object = dge.result, n = Inf, p.value = 0.01)$table %>%
  # push rownames to a column
  rownames_to_column("gene_id") %>%
  # make the table a tibble
  as_tibble() %>%
  # filter for genes with an aboslute logFC of 0.5
  # or greater in at least one of the two comparisons
  filter(abs(logFC.2vs0) >= 0.5 | abs(logFC.4vs0) >= 0.5) %>%
  # Include columns for gene_name and
  # gene_type from the annotation table
  left_join(gene_annotation %>%
              filter(gene_id %in% rownames(dge$counts) &
                       !(gene_id %in% y_genes & chrom == "chrY")) %>%
              dplyr::select(gene_id, gene_name, gene_type))

#* What does a value of 1 for the logFC mean in terms of fold-change?
#*   It indicates that one value is the double of the other (2**1 = 2)

#* Which time point, 2 or 4 hours, has more DEGs?
#*   2 hours
summary(decideTests(object = glmQLFTest(dge.fit, contrast = dge.contrasts[, c("2vs0")]), p.value = 0.01))
#       -1*Treatment_Duration0hr 1*Treatment_Duration2hr
#Down                                                124
#NotSig                                            11860
#Up                                                  276

#*   4 hours
summary(decideTests(object = glmQLFTest(dge.fit, contrast = dge.contrasts[, c("4vs0")]), p.value = 0.01))
#       -1*Treatment_Duration0hr 1*Treatment_Duration4hr
#Down                                                211
#NotSig                                            11617
#Up                                                  432

#*   4 hours has more DEGs

#* Are there any DEGs that are over-expressed in 2vs0 but under-expressed in 4vs0?
allDeg[allDeg$logFC.2vs0 > 0 & allDeg$logFC.4vs0 < 0, ]$gene_id
# [1] "ENSG00000198618.5" "ENSG00000233276.4"

#* Another method for correcting p-values is the Bonferroni correction. Try to use this instead, and compare the distribution of the adjusted p-values generated in the two cases.
allDeg_bonf <-
  # Get table of differential expression results
  topTags(object = dge.result, n = Inf, p.value = 0.01, adjust.method = 'bonferroni')$table %>%
  # push rownames to a column
  rownames_to_column("gene_id") %>%
  # make the table a tibble
  as_tibble() %>%
  # filter for genes with an aboslute logFC of 0.5
  # or greater in at least one of the two comparisons
  filter(abs(logFC.2vs0) >= 0.5 | abs(logFC.4vs0) >= 0.5) %>%
  # Include columns for gene_name and
  # gene_type from the annotation table
  left_join(gene_annotation %>%
              filter(gene_id %in% rownames(dge$counts) &
                       !(gene_id %in% y_genes & chrom == "chrY")) %>%
              dplyr::select(gene_id, gene_name, gene_type))

pdf('pval_dist.pdf')
hist(allDeg$FDR, main = 'FDR correction')
hist(allDeg_bonf$FWER, main = 'Bonferroni correction')
dev.off()

#* Can you identify the genes that are differentially expressed only in one of the two comparisons?
two_hr_specific <-
  allDeg %>%
  filter(abs(logFC.2vs0) > 0.5) %>%
  pull(gene_id) %>%
  setdiff(allDeg %>%
            filter(abs(logFC.4vs0) > 0.5) %>%
            pull(gene_id))

four_hr_specific <-
  allDeg %>%
  filter(abs(logFC.4vs0) > 0.5) %>%
  pull(gene_id) %>%
  setdiff(allDeg %>%
            filter(abs(logFC.2vs0) > 0.5) %>%
            pull(gene_id))

# 3.4.3 DEG visualization ----
fpkms.degs <- fpkms[allDeg %>% pull(gene_id), ]

pdf('DEG_heatmap.pdf')
pheatmap(
  mat = fpkms.degs,
  show_rownames = FALSE,
  annotation_col = forAnnotation,
  clustering_method = "single",
  scale = "row"
)
dev.off()

#* Can you make a volcano plot?
pdf('volcanoplot.pdf')
# Get table with logFC and p-value for all genes tested
topTags(object = dge.result, n = Inf)$table %>%
  # Make a label column based on our significance threshold
  mutate(
    Significant =
      ifelse(FDR < 0.01, "FDR < 0.01", "Not Sig")
  ) %>%
  # Merge logFC in one column, add another one indicating the comparison
  pivot_longer(c("logFC.2vs0", "logFC.4vs0"),
               names_to = "comparison",
               values_to = "logFC"
  ) %>%
  separate(comparison, into = c("prefix", "comparison")) %>%
  dplyr::select(-prefix) %>%
  # Make volcano plot using ggplot
  ggplot(aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  # Split plot for the two comparisons
  facet_wrap(~comparison, nrow = 1, scales = "free_x") +
  # Add labels for the top30 degs
  # geom_text_repel(
  geom_label(
    data =
      allDeg %>%
      # Do the same processing as before to the allDeg table
      pivot_longer(c("logFC.2vs0", "logFC.4vs0"),
                   names_to = "comparison",
                   values_to = "logFC"
      ) %>%
      separate(comparison, into = c("prefix", "comparison")) %>%
      dplyr::select(-prefix) %>%
      arrange(desc(abs(logFC))) %>%
      head(10),
    aes(label = gene_name),
    size = 3,
    position = position_jitter(width = 1,
                               height = 1,
                               seed = 150)
  )
dev.off()

# 3.5 Functional enrichment ----

# Set of overexpressed genes
overexpressed <-
  allDeg %>%
  filter(logFC.2vs0 > 0 & logFC.4vs0 > 0) %>%
  separate(gene_id, into = c("gene_id", "num")) %>%
  pull(gene_id)

# Set of underexpressed genes
underexpressed <-
  allDeg %>%
  filter(logFC.2vs0 < 0 & logFC.4vs0 < 0) %>%
  separate(gene_id, into = c("gene_id", "num")) %>%
  pull(gene_id)

# Set of all DEGs
degs <-
  allDeg %>%
  separate(gene_id, into = c("gene_id", "num")) %>%
  pull(gene_id)

# Set of all genes used initially in our analysis
background <-
  raw_counts_filtered %>%
  rownames_to_column("gene_id") %>%
  separate(gene_id, into = c("gene_id", "num")) %>%
  pull(gene_id)

# Save the identifiers of the different sets of genes
write.table(x = overexpressed, 
            file = "../functional_analysis/overexpressed_genes.txt", 
            sep = "\n", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(x = underexpressed, 
            file = "../functional_analysis/underexpressed_genes.txt", 
            sep = "\n", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(x = degs, 
            file = "../functional_analysis/degs.txt", 
            sep = "\n", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(x = background, 
            file = "../functional_analysis/background_genes.txt", 
            sep = "\n", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# #* What is the background? How should you choose it?
# #*    the background should include any gene that could have been positive
#
# #* Can a gene have more than 1 GO term associated? Can you think of any category of genes that are not suitable for this kind of analysis?
# #*   Yes, genes can be defined by more than one category
# 
# #* Which functions do you find enriched when considering all DEG genes? Do they make sense with the biological process we are looking at?
# #*   Functions associated with the immune system (white cell, TGF, TNF), as well as glucocorticoid-specific pathways and adipogenesis
# 
# #* Can you make a barplot illustrating the number of DEGs in our dataset that belong to the top enriched KEGG and/or REAC pathways?
# For bonferroni correction
go_result <- read.table('../functional_analysis/output/degs_bonferroni.csv',
                        sep = ',',
                        header = TRUE)

pdf('GOenrichment_barplot_bonferroni.pdf',
    width = 10)
go_result[go_result$source %in% c("KEGG", "REAC", "WP"), ] %>%
  # Sort enriched pathways by p-value
  arrange(adjusted_p_value) %>%
  # Keep top 15
  head(15) %>%
  ggplot(aes(y = intersection_size, x = term_name)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_col(aes(fill = negative_log10_of_adjusted_p_value)) +
  scale_fill_gradient(
    high = "firebrick",
    low = "azure4"
  ) +
  # Flip coordinates
  coord_flip() +
  # Split plot for the three sources (i.e. KEGG, REAC, WP)
  facet_wrap(~source, nrow = 1, scales = "free_x") +
  # Change theme
  theme_bw()
dev.off()

# For custom g:SCS correction
go_result <- read.table('../functional_analysis/output/degs_gSCS.csv',
                        sep = ',',
                        header = TRUE)

pdf('GOenrichment_barplot_gSCS.pdf',
    width = 10)
go_result[go_result$source %in% c("KEGG", "REAC", "WP"), ] %>%
  # Sort enriched pathways by p-value
  arrange(adjusted_p_value) %>%
  # Keep top 15
  head(15) %>%
  ggplot(aes(y = intersection_size, x = term_name)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_col(aes(fill = negative_log10_of_adjusted_p_value)) +
  scale_fill_gradient(
    high = "firebrick",
    low = "azure4"
  ) +
  # Flip coordinates
  coord_flip() +
  # Split plot for the three sources (i.e. KEGG, REAC, WP)
  facet_wrap(~source, nrow = 1, scales = "free_x") +
  # Change theme
  theme_bw()
dev.off()


# Alternative GO enrichment with gprofiler2 ----

# pdf('GOenrichment.pdf')
# Functional overrepresentation analysis for overexpressed genes
gost(
  query = overexpressed,
  custom_bg = background,
  organism = "hsapiens"
) %>%
  gostplot(capped = FALSE, interactive = F)


# Functional overrepresentation analysis for underexpressed genes
gost(
  query = underexpressed,
  custom_bg = background,
  organism = "hsapiens"
) %>%
  gostplot(capped = FALSE, interactive = F)


# Functional overrepresentation analysis for all DEGs
gost(
  query = degs,
  custom_bg = background,
  organism = "hsapiens"
) %>%
  gostplot(capped = FALSE, interactive = TRUE)
# dev.off()


# #* How does the result change if you change the method for p-value correction in the gost function?
# # pdf('GOenrichment_bonferroni.pdf')
# # Functional overrepresentation analysis for overexpressed genes
# gost(
#   query = overexpressed,
#   custom_bg = background,
#   organism = "hsapiens",
#   correction_method = 'bonferroni'
# ) %>%
#   gostplot(capped = FALSE, interactive = F)
# 
# 
# # Functional overrepresentation analysis for underexpressed genes
# gost(
#   query = underexpressed,
#   custom_bg = background,
#   organism = "hsapiens",
#   correction_method = 'bonferroni'
# ) %>%
#   gostplot(capped = FALSE, interactive = F)
# 
# 
# # Functional overrepresentation analysis for all DEGs
# gost(
#   query = degs,
#   custom_bg = background,
#   organism = "hsapiens",
#   correction_method = 'bonferroni'
# ) %>%
#   gostplot(capped = FALSE, interactive = F)
# # dev.off()
# 

# Alternative GO enrichment with topGO ----
library(topGO)
library(biomaRt)

# Load complete database from ENSEMBL
db <- useMart("ensembl",
              dataset = "hsapiens_gene_ensembl")

# Subset the GO terms from our background genes
go_db <- getBM(attributes = c("go_id", 
                              "ensembl_gene_id"), 
               filters = "ensembl_gene_id", 
               values = background, 
               mart = db) %>%
  unstack() # Convert it to a list 

# Functional overrepresentation analysis for overexpressed genes
## Keep only annotated genes 
over4go <- overexpressed[overexpressed %in% names(go_db)]

## And generate named factor to signal the candidate genes to be tested (1 for candidate genes)
over4go <- factor(as.integer(background %in% over4go))
names(over4go) <- background

## Generate topGO data object
over_GOdata <- new("topGOdata", 
                   ontology = "BP", # choosing the biological process ontology
                   allGenes = over4go, 
                   annot = annFUN.gene2GO, # annotation function when annotations are given in a gene-to-GO manner
                   gene2GO = go_db)

## Test for significant enrichments
over_result <- runTest(over_GOdata, 
                       algorithm = "weight01", # take into account GO hierarchy
                       statistic = "fisher") 

## Generate a table for the results
over_allGO <- usedGO(over_GOdata)
over_res <- GenTable(over_GOdata, 
                     weightFisher = over_result, 
                     orderBy = "weightFisher", 
                     topNodes = length(over_allGO))

# Functional overrepresentation analysis for underexpressed genes
## Keep only annotated genes 
under4go <- underexpressed[underexpressed %in% names(go_db)]

## And generate named factor to signal the candidate genes to be tested (1 for candidate genes)
under4go <- factor(as.integer(background %in% under4go))
names(under4go) <- background

## Generate topGO data object
under_GOdata <- new("topGOdata", 
                    ontology = "BP", # choosing the biological process ontology
                    allGenes = under4go, 
                    annot = annFUN.gene2GO, # annotation function when annotations are given in a gene-to-GO manner
                    gene2GO = go_db)

## Test for significant enrichments
under_result <- runTest(under_GOdata, 
                        algorithm = "weight01", # take into account GO hierarchy
                        statistic = "fisher") 

## Generate a table for the results
under_allGO <- usedGO(under_GOdata)
under_res <- GenTable(under_GOdata, 
                      weightFisher = under_result, 
                      orderBy = "weightFisher", 
                      topNodes = length(over_allGO))

# Functional overrepresentation analysis for all DEGs
## Keep only annotated genes 
degs4go <- degs[degs %in% names(go_db)]

## And generate named factor to signal the candidate genes to be tested (1 for candidate genes)
degs4go <- factor(as.integer(background %in% degs4go))
names(degs4go) <- background

## Generate topGO data object
degs_GOdata <- new("topGOdata", 
                   ontology = "BP", # choosing the biological process ontology
                   allGenes = degs4go, 
                   annot = annFUN.gene2GO, # annotation function when annotations are given in a gene-to-GO manner
                   gene2GO = go_db)

## Test for significant enrichments
degs_result <- runTest(degs_GOdata, 
                       algorithm = "weight01", # take into account GO hierarchy
                       statistic = "fisher") 

## Generate a table for the results
degs_allGO <- usedGO(degs_GOdata)
degs_res <- GenTable(degs_GOdata, 
                     weightFisher = degs_result, 
                     orderBy = "weightFisher", 
                     topNodes = length(degs_allGO))


