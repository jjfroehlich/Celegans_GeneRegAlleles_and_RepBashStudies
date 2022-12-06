##################################################/
## R script for: "Published C. elegans Gene-Regulatory Alleles and Reporter Bashing Studies" ####
## Jonathan J. Froehlich, Nikolaus Rajewsky
##
## The script is separated into five sections:
## 0 Setup
## 1 Load libraries
## 2 Load data
## 3 Transcriptome cleanup
## 4 Annotate the datasets
## 5 Analyze gene-regulatory alleles and reporter bashing studies
##
## see document outline for table of contents (press Ctrl+Shift+O in RStudio Windows)
## The script performs the following main steps: 
## -Downloads “classical alleles” and transcriptome/genome annotation from WormBase; 
## -Cleans up transcriptome annotation to remove overlaps conservatively (e.g., 3′UTR sequence that overlaps with coding sequence); 
## -Categorizes each classical allele by its overlap with genomic feature (e.g., coding, synonymous, 5′-/3′UTR, intergenic); 
## -Exports a spreadsheet with this information together with the original WormBase annotation; 
## -This table then has to be populated manually for example in a common spreadsheet software (we provide the annotated table in our Extended Data, the script should also find and download these); 
## -Exports .bed files of alleles according to their overlap with genomic feature to help in curation; 
## -Loads the manually populated spreadsheet; 
## -Exports .pdf with the plots shown in this publication, plus additional plots; 
## -Exports browser shots of each allele together with the transcriptome annotation to help in curation.
##
##
##################################################/

## set working directory:
working_directory <- "D:/R_projects/Froehlich_Celegans_GeneRegAlleles_and_RepBashStudies"

####################################/
## 0 Setup            ##############
####################################/

##+++++++++++++++++++++++++++++++++++
## creating directories
##+++++++++++++++++++++++++++++++++++

## set working directory, path was defined above
setwd(working_directory)

## defines sub directories
data_dir <- "data/"
output_dir <- "output/"
supplemental_data_dir <- "supplemental_data/"
output_sub_dir1 <- "output/classical_alleles_bed_files/"
output_sub_dir2 <- "output/classical_alleles_browser_shots/"
output_sub_dir3 <- "output/cleaned_regions_bed_files/"
output_sub_dir4 <- "output/cleaned_regions_browser_shots/"

## creates sub directories
if (!dir.exists(data_dir)) {dir.create(data_dir)}
if (!dir.exists(output_dir)) {dir.create(output_dir)}
if (!dir.exists(supplemental_data_dir)) {dir.create(supplemental_data_dir)}
if (!dir.exists(output_sub_dir1)) {dir.create(output_sub_dir1)} # add "recursive = TRUE" to dir.create if necessary
if (!dir.exists(output_sub_dir2)) {dir.create(output_sub_dir2)}
if (!dir.exists(output_sub_dir3)) {dir.create(output_sub_dir3)}
if (!dir.exists(output_sub_dir4)) {dir.create(output_sub_dir4)}

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## checking if data sets exists, if not downloads them, ~1 GB
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## transcriptome annotation, 8 MB
if (!file.exists("data/c_elegans.PRJNA13758.WS284.canonical_geneset.gtf.gz"))
{download.file("https://downloads.wormbase.org/releases/WS284/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS284.canonical_geneset.gtf.gz",
               "data/c_elegans.PRJNA13758.WS284.canonical_geneset.gtf.gz")}

## transcriptome annotation, 750 MB
if (!file.exists("data/c_elegans.PRJNA13758.WS284.annotations.gff3.gz"))
{download.file("https://downloads.wormbase.org/releases/WS284/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS284.annotations.gff3.gz",
               "data/c_elegans.PRJNA13758.WS284.annotations.gff3.gz")}

## conservation scores, 350 MB
if (!file.exists("data/ce11.phyloP135way.bw"))
{download.file("http://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP135way/ce11.phyloP135way.bw",
               "data/ce11.phyloP135way.bw")}

## Classical Alleles from WormBase, 3 MB
# Version used for the publication, WS284
if (!file.exists("data/VARIATIONS_CLASSICAL_ALLELES.gff3"))
{download.file("https://github.com/jonathanfroehlich/Celegans_GeneRegAlleles_and_RepBashStudies/blob/main/VARIATIONS_CLASSICAL_ALLELES_WS284.gff3",
               "data/VARIATIONS_CLASSICAL_ALLELES.gff3")}

# Alternatively, latest version directly from WormBase
# if (!file.exists("data/VARIATIONS_CLASSICAL_ALLELES.gff3"))
# {download.file("https://wormbase.org/tools/genome/gbrowse/c_elegans_PRJNA13758/?gbgff=1;l=VARIATIONS_CLASSICAL_ALLELES;s=0;f=save+gff3;format=gff3",
#                "data/VARIATIONS_CLASSICAL_ALLELES.gff3")}
# 
# Alternatively one can download this data manually one-by-one from the WormBase Jbrowse
# In case you do this you need to uncomment the code that merges these files, see below in the section "2-Load data" - "cele classical alleles" 
# "Classical alleles" from WormBase release WS284, downloaded one-by-one from the WormBase Jbrowse at https://wormbase.org/tools/genome/jbrowse-simple/?data=data/c_elegans_PRJNA13758. Click on little triangle next to track title, "save track data", select "whole reference sequence", format "gff3". Do this on each chromosome.
# "WS284_release_Classical alleles_I.gff3"
# "WS284_release_Classical alleles_II.gff3"
# "WS284_release_Classical alleles_II.gff3"
# "WS284_release_Classical alleles_IV.gff3"
# "WS284_release_Classical alleles_V.gff3"
# "WS284_release_Classical alleles_X.gff3"

## Extended Data from Froehlich & Rajewsky which contains manual annotations (see last step of 4c in this script)
if (!file.exists("supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx"))
{download.file("https://github.com/jonathanfroehlich/Celegans_GeneRegAlleles_and_RepBashStudies/blob/main/Extended_Data_Tables_S1_S2_S3.xlsx",
               "supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx")}


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## checking if libraries exists, if not installs them
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if(!require(dplyr)){install.packages("dplyr")}
if(!require(data.table)){install.packages("data.table")}
if(!require(data.table)){install.packages("patchwork")}
if(!require(stringr)){install.packages("openxlsx")}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# if(!require(biomaRt)){BiocManager::install("biomaRt")}
if(!require(GenomicRanges)){BiocManager::install("GenomicRanges")}
if(!require(GenomicFeatures)){BiocManager::install("GenomicFeatures")}
if(!require(rtracklayer)){BiocManager::install("rtracklayer")}
if(!require(plyranges)){BiocManager::install("plyranges")}
if(!require(Gviz)){BiocManager::install("Gviz")}

if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(forcats)){install.packages("forcats")}
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
if(!require(viridis)){install.packages("viridis")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(ggrepel)){install.packages("ggrepel")}
if(!require(ggpp)){install.packages("ggpp")}
if(!require(scales)){install.packages("scales")}
if(!require(see)){install.packages("see")}

####################################/
## 1 Load libraries   ##############
####################################/

## set working directory, path was defined above
setwd(working_directory)

## loads libraries
## data transformation
library(dplyr)
library(data.table)
library(patchwork)
# library(svglite)
library(stringr)
library(openxlsx)
# library(readxl)

## genomic data
# library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)
# library(Biostrings)

## browser shots
library(Gviz)

## plotting
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(ggpp)
library(ggrepel)
library(scales)
library(see) # for geom_violinhalf

########################################/
## 2 Load data             #############
########################################/

##++++++++++++++++++++++++++++
### transcriptome ###
##++++++++++++++++++++++++++++

## transcriptome from gtf file
gtf.gr <- import("data/c_elegans.PRJNA13758.WS284.canonical_geneset.gtf.gz") # from: https://downloads.wormbase.org/releases/WS284/species/c_elegans/PRJNA13758/
seqlevelsStyle(gtf.gr) <- "UCSC"

## add proper seqinfo # found in header of c_elegans.PRJNA13758.WS284.annotations.gff3.gz
seqinfo(gtf.gr) <- Seqinfo(seqnames=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"),
                           seqlengths=c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794),
                           isCircular=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                           genome="ce11")

## filter transcriptome for "protein_coding" vs "ncRNAs" in data.table then re-export and load
## its hacky but I couldnt find a better way, several other approaches failed 
## (I had tried alternative approaches using 1. TxDb.Celegans.UCSC.ce11.refGene Annotation package, 2. makeTxDbFromEnsembl, and 3. makeTxDbFromBiomart, each with filtering for protein coding when loading)

## coerce to data.table for this
gtf.dt <- as.data.table(gtf.gr)
gtf_protein.dt <- gtf.dt[gene_biotype == "protein_coding"]
gtf_ncRNA.dt <- gtf.dt[gene_biotype != "protein_coding"]

## back to granges
gtf_protein.gr <- as_granges(gtf_protein.dt)
gtf_ncRNA.gr <- as_granges(gtf_ncRNA.dt)

## assign seqinfo # found in header of c_elegans.PRJNA13758.WS284.annotations.gff3.gz
seqinfo(gtf_protein.gr) <- Seqinfo(seqnames=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"),
                                   seqlengths=c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794),
                                   isCircular=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                                   genome="ce11")
seqinfo(gtf_ncRNA.gr) <- Seqinfo(seqnames=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"),
                                 seqlengths=c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794),
                                 isCircular=c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                                 genome="ce11")

## make txdb
txdb <- makeTxDbFromGRanges(gtf.gr, drop.stop.codons=FALSE)
txdb_protein <- makeTxDbFromGRanges(gtf_protein.gr, drop.stop.codons=FALSE)
txdb_ncRNA <- makeTxDbFromGRanges(gtf_ncRNA.gr, drop.stop.codons=FALSE)

## remove unnecessary files
rm(gtf.dt)
rm(gtf_protein.dt)
rm(gtf_ncRNA.dt)
rm(gtf_protein.gr)
rm(gtf_ncRNA.gr)

##++++++++++++++++++++++++++++++++++
### classical alleles from WormBase ###
##++++++++++++++++++++++++++++++++++

classical_alleles.gr <- import("data/VARIATIONS_CLASSICAL_ALLELES.gff3")
seqlevelsStyle(classical_alleles.gr) <- "UCSC"

# ## alternatively, if you have downloaded from the WormBase jbrowse one-by-one, "classical alleles" from WormBase release WS284
# ## then uncomment these steps, that merge the separate files
# a_1 <- import("data/WS284_release_Classical alleles_I.gff3")
# a_2 <- import("data/WS284_release_Classical alleles_II.gff3")
# a_3 <- import("data/WS284_release_Classical alleles_II.gff3")
# a_4 <- import("data/WS284_release_Classical alleles_IV.gff3")
# a_5 <- import("data/WS284_release_Classical alleles_V.gff3")
# a_6 <- import("data/WS284_release_Classical alleles_X.gff3")
# 
# ## merge
# classical_alleles.gr <- c(a_1, a_2, a_3, a_4, a_5, a_6)
# seqlevelsStyle(classical_alleles.gr) <- "UCSC"
# 
# ## remove unnecessary 
# rm(a_1, a_2, a_3, a_4, a_5, a_6)


##################################################/
## 3 Transcriptome cleanup ####
## removing any overlaps
## e.g. results in all 3'UTR nucleotides that do not overlap with any CDS.
## This is to remove any confounding CDS overlap when analyzing non-coding regions. 
##################################################/

##+++++++++++++++++++++++++++++++++++++++++++++++++
### 3.1 - extract UTRs/CDS/introns per gene ####
##+++++++++++++++++++++++++++++++++++++++++++++++++

##-------------------------------------------------.
## 3UTR granges per GENE
## source: Herve Pages post on Bioconductor, 2019 https://support.bioconductor.org/p/124164/
##-------------------------------------------------.

utr3_by_tx.gr <- threeUTRsByTranscript(txdb_protein)

## Get all 3' UTRs in a GRanges object:
all_utr3.gr <- unlist(utr3_by_tx.gr, use.names=FALSE)

## The names on 'utr3_by_tx' are the **internal** transcript ids (the
## same you see in the 'tx_id' metadata column of the object returned
## by 'transcripts(txdb_protein)', note that they're stored as integer values
## in the SQLite db). Let's propagate them to 'all_utr3.gr':
mcols(all_utr3.gr)$tx_id <- rep(as.integer(names(utr3_by_tx.gr)), lengths(utr3_by_tx.gr))

## Map transcripts to genes:
tx2gene <- mcols(transcripts(txdb_protein, columns=c("tx_id", "tx_name", "gene_id")))
tx2gene$gene_id <- as.character(tx2gene$gene_id)

## Add "tx_name" and "gene_id" metadata columns to 'all_utr3.gr':
m <- match(mcols(all_utr3.gr)$tx_id, tx2gene$tx_id)
mcols(all_utr3.gr) <- cbind(mcols(all_utr3.gr), tx2gene[m, -1L, drop=FALSE])

## Finally we can get the 3' UTRs grouped by gene by splitting all_utr3.gr by the gene_id metadata column:
utr3_by_gene.gr <- split(all_utr3.gr, mcols(all_utr3.gr)$gene_id)

## collapse ("reduce") overlapping isoforms
utr3_by_gene_collapsed.gr <- reduce(utr3_by_gene.gr)

## export .bed to check in browser
# seqlevelsStyle(utr3_by_gene_collapsed.gr) <- "UCSC" # adjust style to fit UCSC browser (e.g. chrI vs chr1 vs Chr1 etc.)
# export.bed(utr3_by_gene_collapsed.gr, "utr3_by_gene_collapsed.bed")

##-------------------------------------------------.
## 5UTR granges per GENE 
## (like above for 3'UTRs) 
##-------------------------------------------------.

utr5_by_tx.gr <- fiveUTRsByTranscript(txdb_protein)
all_utr5.gr <- unlist(utr5_by_tx.gr, use.names=FALSE)
mcols(all_utr5.gr)$tx_id <- rep(as.integer(names(utr5_by_tx.gr)), lengths(utr5_by_tx.gr))
tx2gene <- mcols(transcripts(txdb_protein, columns=c("tx_id", "tx_name", "gene_id")))
tx2gene$gene_id <- as.character(tx2gene$gene_id)
m <- match(mcols(all_utr5.gr)$tx_id, tx2gene$tx_id)
mcols(all_utr5.gr) <- cbind(mcols(all_utr5.gr), tx2gene[m, -1L, drop=FALSE])
utr5_by_gene.gr <- split(all_utr5.gr, mcols(all_utr5.gr)$gene_id)
utr5_by_gene_collapsed.gr <- reduce(utr5_by_gene.gr)

##-------------------------------------------------.
## CDS granges per GENE 
## (like above for 3'UTRs) 
##-------------------------------------------------.

cds_by_tx.gr <- cdsBy(txdb_protein, by="tx")
all_cds.gr <- unlist(cds_by_tx.gr, use.names=FALSE)
mcols(all_cds.gr)$tx_id <- rep(as.integer(names(cds_by_tx.gr)), lengths(cds_by_tx.gr))
tx2gene <- mcols(transcripts(txdb_protein, columns=c("tx_id", "tx_name", "gene_id")))
tx2gene$gene_id <- as.character(tx2gene$gene_id)
m <- match(mcols(all_cds.gr)$tx_id, tx2gene$tx_id)
mcols(all_cds.gr) <- cbind(mcols(all_cds.gr), tx2gene[m, -1L, drop=FALSE])
cds_by_gene.gr <- split(all_cds.gr, mcols(all_cds.gr)$gene_id)
cds_by_gene_collapsed.gr <- reduce(cds_by_gene.gr)

rm(m) #remove unnecessary 

##-------------------------------------------------.
## intronic granges per GENE 
##-------------------------------------------------.
introns_by_tx.gr <- intronicParts(txdb_protein, linked.to.single.gene.only=TRUE)
introns_by_gene_collapsed_u.gr <- unlist(reduce(split(introns_by_tx.gr, introns_by_tx.gr$gene_id)))


##+++++++++++++++++++++++++++++++++++++++++++++++++
### 3.2 - remove overlapping nucleotides to obtain unique 5UTR/3UTR/CDS ####
##+++++++++++++++++++++++++++++++++++++++++++++++++

## unlist
utr3_by_gene_collapsed_u.gr <- unlist(utr3_by_gene_collapsed.gr)
utr5_by_gene_collapsed_u.gr <- unlist(utr5_by_gene_collapsed.gr)
cds_by_gene_collapsed_u.gr <- unlist(cds_by_gene_collapsed.gr)

## Make "negative" of regions:

## 1.Get chromosome info from the db and make a GRanges out of it
chroms.gr <- as(seqinfo(txdb_protein),"GRanges")

# ## in case txdb_protein doesnt contain proper seqinfo
# ## "chroms.gr" can be constructed manually  
# chroms.gr <- as_granges(data.frame(seqnames = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"),
#                                  start = c(1, 1, 1, 1, 1, 1, 1),
#                                  end = c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794), 
#                                  width = c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794), 
#                                  strand = c("*", "*", "*", "*", "*", "*", "*")))
# ## if there is an issue also check seqlevels style (e.g. chrI,...,chrM, versus: I,...,MtDNA etc.)

## 2.Take the set difference between chromosomes and regions, using setdiff_ranges from library(plyranges)
utr3_gaps.gr <- setdiff_ranges(chroms.gr,utr3_by_gene_collapsed_u.gr)
utr5_gaps.gr <- setdiff_ranges(chroms.gr,utr5_by_gene_collapsed_u.gr)
cds_gaps.gr <- setdiff_ranges(chroms.gr,cds_by_gene_collapsed_u.gr)

## clean up each type of region by removing bp that overlap with other type of region (regions: 5UTR/3UTR/CDS)
## only keep what overlaps with "negative" ("*_gaps"), using join_overlap_intersect from library(plyranges)
utr3_by_gene_collapsed_clean.gr <- join_overlap_intersect(utr3_by_gene_collapsed_u.gr, utr5_gaps.gr)
utr3_by_gene_collapsed_clean.gr <- join_overlap_intersect(utr3_by_gene_collapsed_clean.gr, cds_gaps.gr)

## same for utr5 and cds
utr5_by_gene_collapsed_clean.gr <- join_overlap_intersect(utr5_by_gene_collapsed_u.gr, utr3_gaps.gr)
utr5_by_gene_collapsed_clean.gr <- join_overlap_intersect(utr5_by_gene_collapsed_clean.gr, cds_gaps.gr)
cds_by_gene_collapsed_clean.gr <- join_overlap_intersect(cds_by_gene_collapsed_u.gr, utr3_gaps.gr)
cds_by_gene_collapsed_clean.gr <- join_overlap_intersect(cds_by_gene_collapsed_clean.gr, utr5_gaps.gr)

## same for introns
introns_by_gene_collapsed_clean.gr <- join_overlap_intersect(introns_by_gene_collapsed_u.gr, utr3_gaps.gr)
introns_by_gene_collapsed_clean.gr <- join_overlap_intersect(introns_by_gene_collapsed_clean.gr, utr5_gaps.gr)
introns_by_gene_collapsed_clean.gr <- join_overlap_intersect(introns_by_gene_collapsed_clean.gr, cds_gaps.gr)

## copy-paste gene_id as new metadata column
utr3_by_gene_collapsed_clean.gr$gene_id  <- names(utr3_by_gene_collapsed_clean.gr)
utr5_by_gene_collapsed_clean.gr$gene_id  <- names(utr5_by_gene_collapsed_clean.gr)
cds_by_gene_collapsed_clean.gr$gene_id  <- names(cds_by_gene_collapsed_clean.gr)
introns_by_gene_collapsed_clean.gr$gene_id  <- names(introns_by_gene_collapsed_clean.gr)

# ## split into GrangesList by gene, to have separated ranges grouped by gene
# utr3_by_gene_collapsed_clean_list.gr <- split(utr3_by_gene_collapsed_clean.gr, names(utr3_by_gene_collapsed_clean.gr))

## check in browser, export as .bed
seqlevelsStyle(utr3_by_gene_collapsed.gr) <- "UCSC"
export.bed(utr3_by_gene_collapsed.gr, "output/cleaned_regions_bed_files/utr3_by_gene_collapsed.bed")

seqlevelsStyle(utr3_by_gene_collapsed_clean.gr) <- "UCSC"
export.bed(utr3_by_gene_collapsed_clean.gr, "output/cleaned_regions_bed_files/utr3_by_gene_collapsed_clean.bed")

seqlevelsStyle(introns_by_gene_collapsed_clean.gr) <- "UCSC"
export.bed(introns_by_gene_collapsed_clean.gr, "output/cleaned_regions_bed_files/introns_by_gene_collapsed_clean.bed")


##+++++++++++++++++++++++++++++++++++++++++++++++++
### 3.3 - extract ups/downs regions per gene ####
##+++++++++++++++++++++++++++++++++++++++++++++++++

##--------------------------------------.
## define genic and intergenic regions
##--------------------------------------.
## code snippets from Chris Seidel at https://research.stowers.org/cws/CompGenomics/Tutorial/peak_assignment.html

## txdb_protein: thus only protein coding genes
allTx.gr <- transcriptsBy(txdb_protein, by='gene')

## collapse all transcripts into a single range
allTx.gr <- reduce(allTx.gr)

## Then, let's get rid of the list format
genes.gr <- unlist(allTx.gr)
genes.gr$gene_id <- names(genes.gr)

# make a copy 
genic.gr <- genes.gr

# remove strand information
strand(genic.gr) <- '*'

## reduce gene overlaps to one feature
genic <- reduce(genic.gr)

## Now that we have a definition of all "genic" feature space
## we can take the inverse of this by subtracting it from the
## chromosome definitions.

## Get the chromosome info from the db and make a GRanges out of it
chroms.gr <- as(seqinfo(txdb_protein),'GRanges')

## Take the set difference between chromosomes and genes
intergenic.gr <- setdiff_ranges(chroms.gr,genic.gr)

## check in browser, export as .bed
seqlevelsStyle(intergenic.gr) <- "UCSC"
export.bed(intergenic.gr,"output/cleaned_regions_bed_files/intergenic.bed")

##----------------------------------------------.
## define intergenic upstream ("ups") ranges until next gene
##----------------------------------------------.

## define ups region coordinates
upstream_genes <- follow(genes.gr, unstrand(genes.gr))
upstream_genes[is.na(upstream_genes)] <- 1
genes.gr$ups_gene_id <- names(genes.gr[upstream_genes])
genes.gr$ups_dist <- distance(genes.gr, genes.gr[upstream_genes], ignore.strand=TRUE)
genes.gr$ups_dist[is.na(genes.gr$ups_dist)] <- 0 # 1 doesnt make a difference

## create list of genes that do not have an upstream region (because they overlap with another gene)
## Because: follow and precede ignore overlapping ranges, and will find the next range, thus resulting in a region extending across the overlapping gene. Therefore we will have to filter these out.
genes_start_nt.gr <- mutate(anchor_5p(genes.gr), width = 1) # reduce genes to 1 bp from the start, library(plyranges)
# genes_start_nt.gr <- resize(genes.gr, 1, fix="start") # GRanges alternative
genes_start_nt.gr <- shift_upstream(genes_start_nt.gr, 1) # shift by 1 bp upstream
overlaps <- findOverlaps(genes.gr, genes_start_nt.gr, ignore.strand = TRUE) 
genes.gr$ups_dist[subjectHits(overlaps)] <- 0
ups_regions_to_remove <- as.character(genes.gr$gene_id[subjectHits(overlaps)]) # find ups regions that overlap with another gene

## create upstream regions
ups_regions.gr <- trim(flank_upstream(genes.gr, width=genes.gr$ups_dist)) 

## remove upstream regions that shouldnt exist
ups_regions.dt <- as.data.table(ups_regions.gr)
ups_regions.dt <- ups_regions.dt[!(gene_id %in% ups_regions_to_remove)] # coerce to data.table for this :(
ups_regions.gr <- as_granges(ups_regions.dt)
names(ups_regions.gr) <- ups_regions.gr$gene_id

## subtract ncRNAs
ncRNAs.gr <- transcripts(txdb_ncRNA)
ncRNAs_gaps.gr <- setdiff_ranges(chroms.gr, ncRNAs.gr)
ups_regions_clean.gr <- join_overlap_intersect(ups_regions.gr, ncRNAs_gaps.gr)

## export as .bed to check in browser
seqlevelsStyle(ups_regions.gr) <- "UCSC"
export.bed(ups_regions.gr,"output/cleaned_regions_bed_files/ups_regions.bed")

seqlevelsStyle(ups_regions_clean.gr) <- "UCSC"
export.bed(ups_regions_clean.gr,"output/cleaned_regions_bed_files/ups_regions_clean.bed")

## remove unnecessary objects
rm(upstream_genes)
rm(genes_start_nt.gr)
rm(ups_regions_to_remove)
rm(ups_regions.dt)

##----------------------------------------------.
## define intergenic DOWNS ranges until next gene
##----------------------------------------------.

## define downs region coordinates
downstream_genes <- precede(genes.gr, unstrand(genes.gr))
downstream_genes[is.na(downstream_genes)] <- 1
genes.gr$downs_gene_id <- names(genes.gr[downstream_genes])
genes.gr$downs_dist <- distance(genes.gr, genes.gr[downstream_genes], ignore.strand=TRUE)
genes.gr$downs_dist[is.na(genes.gr$downs_dist)] <- 0 # 1 doesnt make a difference

## create list of genes that do not have an downs region (because they overlap with another gene)
## Because: follow and precede ignore overlapping ranges, and will find the next range, thus resulting in a region extending across the overlapping gene. Therefore we will have to filter these out.
genes_end_nt.gr <- mutate(anchor_3p(genes.gr), width = 1) # reduce genes to 1 bp from the start, library(plyranges)
# genes_end_nt.gr <- resize(genes.gr, 1, fix="end") # GRanges alternative
genes_end_nt.gr <- shift_downstream(genes_end_nt.gr, 1) # shift by 1 bp downs
overlaps <- findOverlaps(genes.gr, genes_end_nt.gr, ignore.strand = TRUE) 
genes.gr$downs_dist[subjectHits(overlaps)] <- 0
downs_regions_to_remove <- as.character(genes.gr$gene_id[subjectHits(overlaps)]) # find downs regions that overlap with another gene

## create downs regions
downs_regions.gr <- trim(flank_downstream(genes.gr, width=genes.gr$downs_dist)) 

## filter out downs regions that overlap with another gene
## coerce to data.table for this
downs_regions.dt <- as.data.table(downs_regions.gr)
downs_regions.dt <- downs_regions.dt[!(gene_id %in% downs_regions_to_remove)]
downs_regions.gr <- as_granges(downs_regions.dt)
names(downs_regions.gr) <- downs_regions.gr$gene_id

## subtract ncRNAs
ncRNAs.gr <- transcripts(txdb_ncRNA)
ncRNAs_gaps.gr <- setdiff_ranges(chroms.gr, ncRNAs.gr)
downs_regions_clean.gr <- join_overlap_intersect(downs_regions.gr, ncRNAs_gaps.gr)

## export as .bed to check in browser
seqlevelsStyle(downs_regions.gr) <- "UCSC"
export.bed(downs_regions.gr,"output/cleaned_regions_bed_files/downs_regions.bed")

seqlevelsStyle(downs_regions_clean.gr) <- "UCSC"
export.bed(downs_regions_clean.gr,"output/cleaned_regions_bed_files/downs_regions_clean.bed")

## remove unnecessary objects
rm(downstream_genes)
rm(genes_end_nt.gr)
rm(downs_regions_to_remove)
rm(downs_regions.dt)

##----------------------------------------------.
## only 500 bp ups / downs ranges 
##----------------------------------------------.
ups_regions.gr$ups_dist_500 <- ups_regions.gr$ups_dist
ups_regions.gr$ups_dist_500[ups_regions.gr$ups_dist > 500] <- 500

downs_regions.gr$downs_dist_500 <- downs_regions.gr$downs_dist
downs_regions.gr$downs_dist_500[downs_regions.gr$downs_dist > 500] <- 500

ups_500.gr <- mutate(anchor_3p(ups_regions.gr), width = (ups_regions.gr$ups_dist_500)) 
downs_500.gr <- mutate(anchor_5p(downs_regions.gr), width = (downs_regions.gr$downs_dist_500))

## subtract ncRNAs
ncRNAs.gr <- transcripts(txdb_ncRNA)
ncRNAs_gaps.gr <- setdiff_ranges(chroms.gr, ncRNAs.gr)

ups_500_clean.gr <- join_overlap_intersect(ups_500.gr, ncRNAs_gaps.gr)
downs_500_clean.gr <- join_overlap_intersect(downs_500.gr, ncRNAs_gaps.gr)

## export as .bed to check in browser
seqlevelsStyle(ups_500.gr) <- "UCSC"
export.bed(ups_500.gr,"output/cleaned_regions_bed_files/ups_500.bed")

seqlevelsStyle(downs_500.gr) <- "UCSC"
export.bed(downs_500.gr,"output/cleaned_regions_bed_files/downs_500.bed")

seqlevelsStyle(ups_500_clean.gr) <- "UCSC"
export.bed(ups_500_clean.gr,"output/cleaned_regions_bed_files/ups_500_clean.bed")

seqlevelsStyle(downs_500_clean.gr) <- "UCSC"
export.bed(downs_500_clean.gr,"output/cleaned_regions_bed_files/downs_500_clean.bed")

##----------------------------------------------.
## remove ups500 from ups regions 
##----------------------------------------------.
ups_500_gaps.gr <- setdiff_ranges(chroms.gr, ups_500_clean.gr)
ups_clean.gr <- join_overlap_intersect(ups_regions_clean.gr, ups_500_gaps.gr)

## export as .bed to check in browser
seqlevelsStyle(ups_clean.gr) <- "UCSC"
export.bed(ups_clean.gr,"output/cleaned_regions_bed_files/ups_clean.bed")

##+++++++++++++++++++++++++++++++++++++++++++++++++
### 3.4 - plot browser shots to check results ####
##+++++++++++++++++++++++++++++++++++++++++++++++++

## fetch all genes
genes.dt <- as.data.table(genes.gr)
## sample 5 random genes
genes_to_plot.dt <- genes.dt[sample(.N, 5)]
## add some chosen genes
genes_to_plot.dt <- rbind(genes_to_plot.dt, genes.dt[gene_id %in% c("WBGene00302984", "WBGene00011144")])
genes_to_plot.list <- genes_to_plot.dt$gene_id

## modify txdb_protein so that gene_name is shown as well (and not only gene_id e.g. WBGene00014472, ...)
gtf2.gr<-gtf.gr
gtf2.gr$gene_id[!is.na(gtf.gr$gene_name)] <- gtf.gr$gene_name[!is.na(gtf.gr$gene_name)] ## alternative1: replace gene_id with gene_name, except for fields were gene_name is NA
# gtf2.gr$gene_id <- paste(gtf.gr$gene_name, gtf.gr$gene_id) ## alternative2: concatenate both columns (because some NAs in gene_name) and replace gene_id
txdb2 <- makeTxDbFromGRanges(gtf2.gr, drop.stop.codons=FALSE)

## open one output file
## disable if desiring individual output files, and enable relevant options in the loop e.g. jpeg()/pdf(), and dev.off()
pdf(file = "output/cleaned_regions_browser_shots/browser_shots_all.pdf")

## for loop
for(g in genes_to_plot.list)
{
  ## extract coordinates (chr, start, end)
  g_coord <- genes_to_plot.dt[genes_to_plot.dt$gene_id==g,]

  gen <- "ce11"
  chr <- as.character(g_coord$seqnames)
  start <- (g_coord$start)-500
  end <- (g_coord$end)+500

  ## enable if desiring individual output files
  # jpeg(file = paste0("output/cleaned_regions_browser_shots/browser_shot_", g, ".jpg"), res = 150)
  # pdf(file = paste0("output/cleaned_regions_browser_shots/browser_shot_", g, ".pdf"))

  ## import only part of bigwig file to speed things up:
  GR <- makeGRangesFromDataFrame(data.frame(chr=chr, start=start, end=end, score=1))
  phyloP_135.bw <- import("data/ce11.phyloP135way.bw", which=GR) #grab only relevant region. Otherwise, if loading the whole .bw before the loop, process is 1.5-1.8x slower
  # phyloP_26.bw <- import("data/ce11.phyloP26way.bw", which=GR) #grab only relevant region

  ## define all the tracks to plot
  ideoTrack <- IdeogramTrack(genome = gen,
                             chromosome = chr,
                             lwd=0.7,
                             fontface=1,
                             fontsize=9,
                             fontcolor="dark grey",
                             col="red",
                             col.border.title="dark grey")

  axisTrack1 <- GenomeAxisTrack(lwd=1, fontface=1, fontsize=9, col = "dark grey", labelPos = "below", fontcolor="dark grey")

  axisTrack2 <- GenomeAxisTrack(lwd=1, fontface=1, fontsize=9, col = "dark grey", labelPos = "below", scale = 0.1)

  Track1<-GeneRegionTrack(range=ups_regions.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          # fontface.group=1,
                          # fontsize.group=12,
                          # fontcolor.group="light grey",
                          # fontcolor="light grey",
                          # fontcolor.title="light grey",
                          col="light grey",
                          col.line="light grey",
                          fill="light grey",
                          shape="box",
                          name="")

  Track2<-GeneRegionTrack(range=ups_regions_clean.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="black",
                          fontcolor="black",
                          fontcolor.title="black",
                          col="black",
                          col.line="black",
                          fill="dark grey",
                          shape="box",
                          name="ups")

  Track3<-GeneRegionTrack(range=downs_regions.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="light grey",
                          fontcolor="light grey",
                          col="light grey",
                          col.line="light grey",
                          fill="light grey",
                          shape="box",
                          name="")

  Track4<-GeneRegionTrack(range=downs_regions_clean.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="black",
                          col="black",
                          col.line="black",
                          fill="dark grey",
                          shape="box",
                          name="downs")

  Track5<-GeneRegionTrack(range=utr3_by_gene_collapsed_u.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="light grey",
                          fontcolor="light grey",
                          col="light grey",
                          col.line="light grey",
                          fill="light grey",
                          shape="box",
                          name="")

  Track6<-GeneRegionTrack(range=utr3_by_gene_collapsed_clean.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="black",
                          col="black",
                          col.line="black",
                          fill="dark grey",
                          shape="box",
                          name="utr3")

  Track7<-GeneRegionTrack(range=utr5_by_gene_collapsed_u.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="light grey",
                          fontcolor="light grey",
                          col="light grey",
                          col.line="light grey",
                          fill="light grey",
                          shape="box",
                          name="")

  Track8<-GeneRegionTrack(range=utr5_by_gene_collapsed_clean.gr,
                          genome=gen,
                          chromosome=chr,
                          start=start,
                          end=end,
                          showId=FALSE,
                          geneSymbols=FALSE,
                          transcriptAnnotation="gene",
                          lwd=0.7,
                          fontface.group=1,
                          fontsize.group=12,
                          fontcolor.group="black",
                          col="black",
                          col.line="black",
                          fill="dark grey",
                          shape="box",
                          name="utr5")

  txdbTrack<-GeneRegionTrack(range=txdb2,
                             genome=gen,
                             chromosome=chr,
                             start=start,
                             end=end,
                             showId=TRUE,
                             geneSymbols=TRUE,
                             transcriptAnnotation="gene",
                             # collapseTranscripts = "longest",
                             lwd=0.7,
                             fontface.group=1,
                             fontsize.group=12,
                             fontcolor.group="black",
                             col="black",
                             col.line="black",
                             lex=8,
                             fill="dark grey",
                             name="txdb")

  phyloP_135Track<-DataTrack(range=phyloP_135.bw,
                             genome=gen,
                             chromosome=chr,
                             start=start,
                             end=end,
                             type="histogram",
                             # window=-1,
                             ylim = c(-4,20),
                             col="grey",
                             col.line="grey",
                             baseline=0,
                             col.baseline="dark grey",
                             col.histogram="grey",
                             fill.histogram="grey",
                             name="135")
  ## plot
  plotTracks(list(ideoTrack, axisTrack1, axisTrack2,
                  Track1, Track2, Track3, Track4, Track5, Track6, Track7, Track8,
                  txdbTrack, phyloP_135Track),
             chromosome=chr,
             from=start,
             to=end,
             # extend.right=1500,
             # extend.left=1500,
             sizes=c(0.4, 0.4, 0.4,
                     0.3, 0.7, 0.3, 0.7, 0.3, 0.7, 0.3, 0.7,
                     1, 1),
             fontface.title=1,
             fontsize.title=12,
             fontcolor.title="black",
             col.border.title="transparent",
             background.title="transparent",
             lwd.title=1, #
             col.axis="black")

  ## enable if using individual output files
  # dev.off()
}
dev.off()



#######################################################/
## 4 Annotate the datasets ####
## 
#######################################################/

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 4.1 - subset alleles according to overlaps  ####
##++++++++++++++++++++++++++++++++++++++++++++++++++

##  Wormbase
##  classical
##   alleles 
##     │
##     ├── 1. CDS overlap ("_CDS") 
##     │        │
##     │        ├── synonymous
##     │        └── others
##     │
##     └── 2. others ("_notCDS")
##              │
##              ├── 3. ncRNA overlap ("_ncRNAs")
##              └── 4. others ("_notncRNAs")
##                      │
##                      ├── 5. within introns ("_introns")
##                      └── 6. 5'/3'-UTR & intergenic ("_in53")


## remove transposon insertion sites from classical alleles
## because they are mostly not forward genetics alleles with phenotype (many are MosSCI insertion sites)
classical_alleles.gr <- filter(classical_alleles.gr, type != "transposable_element_insertion_site")

## 1. alleles that overlap CDS
seqlevelsStyle(cds_by_gene_collapsed_u.gr) <- "UCSC"
classical_alleles_CDS.gr <- subsetByOverlaps(classical_alleles.gr, cds_by_gene_collapsed_u.gr, type="any", ignore.strand=TRUE)
classical_alleles_CDS.dt <- as.data.table(classical_alleles_CDS.gr)
classical_alleles_CDS.dt <- distinct(classical_alleles_CDS.dt) #remove duplicate rows

## 2. alleles that are not CDS and thus within [intron+intergenic+5UTR+3UTR](=cds_gaps.gr)
seqlevelsStyle(classical_alleles.gr) <- "UCSC"
seqlevelsStyle(cds_gaps.gr) <- "UCSC"
classical_alleles_notCDS.gr <- subsetByOverlaps(classical_alleles.gr, cds_gaps.gr, type="within", ignore.strand=TRUE)
classical_alleles_notCDS.dt <- as.data.table(classical_alleles_notCDS.gr)
classical_alleles_notCDS.dt <- distinct(classical_alleles_notCDS.dt) #remove duplicate rows

## 3. alleles that overlap ncRNAs (piRNAs, miRNAs, pseudogenes, tRNAs, etc.)
## first make ncRNA transcript granges
ncRNAs.gr <- transcripts(txdb_ncRNA)
ncRNAs.dt <- as.data.table(ncRNAs.gr)
## then
seqlevelsStyle(classical_alleles_notCDS.gr) <- "UCSC"
seqlevelsStyle(ncRNAs.gr) <- "UCSC"
classical_alleles_ncRNAs.gr <- subsetByOverlaps(classical_alleles_notCDS.gr, ncRNAs.gr, type="any", ignore.strand=TRUE)
classical_alleles_ncRNAs.dt <- as.data.table(classical_alleles_ncRNAs.gr)
classical_alleles_ncRNAs.dt <- distinct(classical_alleles_ncRNAs.dt) #remove duplicate rows

## 4. alleles that do not overlap ncRNAs
## first make ncRNA_gaps 
ncRNA_gaps.gr <- setdiff_ranges(chroms.gr,ncRNAs.gr)
## then
seqlevelsStyle(classical_alleles_notCDS.gr) <- "UCSC"
seqlevelsStyle(ncRNA_gaps.gr) <- "UCSC"
classical_alleles_notncRNAs.gr <- subsetByOverlaps(classical_alleles_notCDS.gr, ncRNA_gaps.gr, type="within", ignore.strand=TRUE)
classical_alleles_notncRNAs.dt <- as.data.table(classical_alleles_notncRNAs.gr)
classical_alleles_notncRNAs.dt <- distinct(classical_alleles_notncRNAs.dt) #remove duplicate rows

## 5. within introns
seqlevelsStyle(classical_alleles_notncRNAs.gr) <- "UCSC"
seqlevelsStyle(introns_by_gene_collapsed_u.gr) <- "UCSC"
classical_alleles_introns.gr <- subsetByOverlaps(classical_alleles_notncRNAs.gr, introns_by_gene_collapsed_u.gr, type="within", ignore.strand=TRUE)
classical_alleles_introns.dt <- as.data.table(classical_alleles_introns.gr)
classical_alleles_introns.dt <- distinct(classical_alleles_introns.dt) #remove duplicate rows

## 6. within intergenic & 5'/3'UTR 
## first make introns+CDS 
introns_CDS.gr <- append(introns_by_gene_collapsed_u.gr, cds_by_gene_collapsed_u.gr)
## then make intergenic+5UTR+3UTR (which is introns+CDS gaps)
intron_cds_gaps.gr <- setdiff_ranges(chroms.gr, introns_CDS.gr)
## then
seqlevelsStyle(classical_alleles_notncRNAs.gr) <- "UCSC"
seqlevelsStyle(intron_cds_gaps.gr) <- "UCSC"
classical_alleles_in53.gr <- subsetByOverlaps(classical_alleles_notncRNAs.gr, intron_cds_gaps.gr, type="within", ignore.strand=TRUE)
classical_alleles_in53.dt <- as.data.table(classical_alleles_in53.gr)
classical_alleles_in53.dt <- distinct(classical_alleles_in53.dt) #remove duplicate rows

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 4.2 - export as .bed to check in UCSC browser ####
##++++++++++++++++++++++++++++++++++++++++++++++++++
seqlevelsStyle(classical_alleles.gr) <- "UCSC"
classical_alleles.gr$score <- "1"
classical_alleles.gr$score <- as.numeric(classical_alleles.gr$score)
export.bed(classical_alleles.gr, "output/classical_alleles_bed_files/classical_alleles.bed")

seqlevelsStyle(classical_alleles_CDS.gr) <- "UCSC"
classical_alleles_CDS.gr$score <- "1"
classical_alleles_CDS.gr$score <- as.numeric(classical_alleles_CDS.gr$score)
export.bed(classical_alleles_CDS.gr, "output/classical_alleles_bed_files/classical_alleles_CDS.bed")

seqlevelsStyle(classical_alleles_notCDS.gr) <- "UCSC"
classical_alleles_notCDS.gr$score <- "1"
classical_alleles_notCDS.gr$score <- as.numeric(classical_alleles_notCDS.gr$score)
export.bed(classical_alleles_notCDS.gr, "output/classical_alleles_bed_files/classical_alleles_notCDS.bed")

seqlevelsStyle(ncRNAs.gr) <- "UCSC"
ncRNAs.gr$score <- "1"
ncRNAs.gr$score <- as.numeric(ncRNAs.gr$score)
export.bed(ncRNAs.gr, "output/classical_alleles_bed_files/ncRNAs.bed")

seqlevelsStyle(classical_alleles_ncRNAs.gr) <- "UCSC"
classical_alleles_ncRNAs.gr$score <- "1"
classical_alleles_ncRNAs.gr$score <- as.numeric(classical_alleles_ncRNAs.gr$score)
export.bed(classical_alleles_ncRNAs.gr, "output/classical_alleles_bed_files/classical_alleles_ncRNAs.bed")

seqlevelsStyle(classical_alleles_notncRNAs.gr) <- "UCSC"
classical_alleles_notncRNAs.gr$score <- "1"
classical_alleles_notncRNAs.gr$score <- as.numeric(classical_alleles_notncRNAs.gr$score)
export.bed(classical_alleles_notncRNAs.gr, "output/classical_alleles_bed_files/classical_alleles_notncRNAs.bed")

seqlevelsStyle(classical_alleles_introns.gr) <- "UCSC"
classical_alleles_introns.gr$score <- "1"
classical_alleles_introns.gr$score <- as.numeric(classical_alleles_introns.gr$score)
export.bed(classical_alleles_introns.gr, "output/classical_alleles_bed_files/classical_alleles_introns.bed")

seqlevelsStyle(classical_alleles_in53.gr) <- "UCSC"
classical_alleles_in53.gr$score <- "1"
classical_alleles_in53.gr$score <- as.numeric(classical_alleles_in53.gr$score)
export.bed(classical_alleles_in53.gr, "output/classical_alleles_bed_files/classical_alleles_in53.bed")

## remove unused
rm(classical_alleles_CDS.gr, classical_alleles_notCDS.gr, ncRNAs.gr, classical_alleles_ncRNAs.gr, classical_alleles_notncRNAs.gr, classical_alleles_introns.gr)

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 4.3 - export as .xlsx, add manual annotation, import ####
##++++++++++++++++++++++++++++++++++++++++++++++++++
classical_alleles_CDS.dt$overlap_feature_R <- "CDS"
classical_alleles_ncRNAs.dt$overlap_feature_R <- "ncRNAs"
classical_alleles_introns.dt$overlap_feature_R <- "introns"
classical_alleles_in53.dt$overlap_feature_R <- "in53"

## modify classical_alleles_CDS.dt: add "CDS_synonymous" where X_Consequence contains "synonymous_variant"
classical_alleles_CDS.dt[consequence == "synonymous_variant", overlap_feature_R := "CDS_synonymous"]

## combine
classical_alleles_export.dt <- rbind(classical_alleles_CDS.dt, classical_alleles_ncRNAs.dt, classical_alleles_introns.dt, classical_alleles_in53.dt)

## remove unused
rm(classical_alleles_CDS.dt, classical_alleles_notCDS.dt, classical_alleles_ncRNAs.dt, classical_alleles_notncRNAs.dt, classical_alleles_introns.dt, classical_alleles_in53.dt)

## add columns for manual annotation
classical_alleles_export.dt$name <- classical_alleles_export.dt$public_name
classical_alleles_export.dt$overlap_feature_R_postManAnnot <- classical_alleles_export.dt$overlap_feature_R
classical_alleles_export.dt$consequence_postManAnnot <- classical_alleles_export.dt$consequence
classical_alleles_export.dt$ManAnnotSyn_truly_synonymous_and_directly_affects_phenotype <- "-"
classical_alleles_export.dt$ManAnnotSyn_mechanism_suggested <- "-"
classical_alleles_export.dt$width_postManAnnot <- classical_alleles_export.dt$width
classical_alleles_export.dt$type_postManAnnot <- classical_alleles_export.dt$type
classical_alleles_export.dt$production_method_postManAnnot <- classical_alleles_export.dt$production_method
classical_alleles_export.dt$engineered_postManAnnot <- classical_alleles_export.dt$engineered
classical_alleles_export.dt$ManAnnotIn53_distance_to_coding <- "-"
classical_alleles_export.dt$ManAnnotIn53_location <- "-"
classical_alleles_export.dt$ManAnnotIn53_gene <- "-"
classical_alleles_export.dt$ManAnnotIn53_variants_tested <- "-"
classical_alleles_export.dt$ManAnnotIn53_mechanism_explained <- "-"
classical_alleles_export.dt$ManAnnotIn53_mechanism <- "-"
classical_alleles_export.dt$ManAnnotIn53_consequence_on_expression <- "-"
classical_alleles_export.dt$ManAnnotIn53_mutation_expression_consequence <- "-"
classical_alleles_export.dt$ManAnnotIn53_phenotype <- "-"
classical_alleles_export.dt$ManAnnot_comment <- "-"
classical_alleles_export.dt$ManAnnot_ref_regions <- "-"
classical_alleles_export.dt$ManAnnot_reference <- "-"
classical_alleles_export.dt$ManAnnot_ref_year <- "-"

## subset for export
synon_annot <- classical_alleles_export.dt[overlap_feature_R == "CDS_synonymous"]
in53_annot <- classical_alleles_export.dt[overlap_feature_R == "in53" | consequence == "3_prime_UTR_variant" | consequence == "5_prime_UTR_variant"]
classical_alleles_export.dt[is.na(classical_alleles_export.dt)] <- "-" # to prevent dropping rows in next line
rest_annot <- classical_alleles_export.dt[!(overlap_feature_R == "CDS_synonymous" | overlap_feature_R == "in53" | consequence == "3_prime_UTR_variant" | consequence == "5_prime_UTR_variant")]

## prep sheet to manually add reporter bashing studies
reporter_bashing_annot <- data.table()
reporter_bashing_annot$overlap_feature_R <- "-"
reporter_bashing_annot$gene <- "-"
reporter_bashing_annot$min_width <- "-"
reporter_bashing_annot$max_width <- "-"
reporter_bashing_annot$min_distance_to_coding <- "-"
reporter_bashing_annot$max_distance_to_coding <- "-"
reporter_bashing_annot$expression_technique <- "-"
reporter_bashing_annot$approach <- "-"
reporter_bashing_annot$mutation <- "-"
reporter_bashing_annot$ManAnnotIn53_variants_tested <- "-"
reporter_bashing_annot$mechanism_explained <- "-"
reporter_bashing_annot$evidence <- "-"
reporter_bashing_annot$mechanism <- "-"
reporter_bashing_annot$consequence_on_expression <- "-"
reporter_bashing_annot$ManAnnotIn53_mutation_expression_consequence <- "-"
reporter_bashing_annot$comment <- "-"
reporter_bashing_annot$ManAnnot_reference <- "-"
reporter_bashing_annot$ManAnnot_ref_year <- "-"

## export as excel file 
dataset_names <- list('all' = classical_alleles_export.dt, 'synon_annot' = synon_annot, 'in53_annot' = in53_annot, 'reporter_annot' = reporter_bashing_annot)
write.xlsx(dataset_names, file = "output/table_for_manual_annotation.xlsx")
rm(dataset_names)

## then add manual annotations in excel by using the file "output/table_for_manual_annotation.xlsx"
## sheet "synon_annot" for synonymous mutations 
## sheet "in53_annot" for promoter/enhancer/5'UTR/3'UTR
## populate empty fields in excel file, annotate for example using Wormbase and search for variation by its name, and perform general literature mining
## sheet "reporter_annot" for reporter bashing studies
## use literature searches to find relevant studies and populate empty fields
##
## important! after the manual annotation, save the annotated file as "supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx"


#######################################################/
## 5 Analyze gene-regulatory alleles and reporter bashing publications ####
## 
#######################################################/

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 5.1 - load annotated data (e.g. the Extended Data from Froehlich & Rajewsky)    ####
##++++++++++++++++++++++++++++++++++++++++++++++++++
## load annotated classical alleles 
synon_annot <- read.xlsx("supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx", sheet = "synon_annot")
in53_annot <-  read.xlsx("supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx", sheet = "in53_annot")

## recombine the now annotated variants with rest
classical_alleles_post_ManAnnot.dt <- rbind(rest_annot, synon_annot, in53_annot, fill=TRUE)
classical_alleles_post_ManAnnot.dt[classical_alleles_post_ManAnnot.dt == "-"] <- NA # put NAs back
rm(synon_annot, in53_annot)

## load reporter bashing studies 
reporter_annot <- setDT(read.xlsx("supplemental_data/Extended_Data_Tables_S1_S2_S3.xlsx", sheet = "reporter_annot"))

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 5.2 - plots    ####
##++++++++++++++++++++++++++++++++++++++++++++++++++
## code to add percentage labels to stacked bar charts from MikeM https://stackoverflow.com/questions/66631575/ggplot-stacked-bar-plot-with-percentage-labels

## dataset PRE manual annotation 
#---------------------------------/
## as data.frame for plots
classical_alleles.dt <- as.data.table(classical_alleles.gr)
classical_alleles.dt <- distinct(classical_alleles.dt) #remove duplicate rows
## add own manual annotation (2 new columns that group several consequences together, con_sum1, con_sum2)
summarize_consequence_types.dt <- data.table(consequence  = c("start_lost", "stop_gained", "stop_lost", "frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "coding_sequence_variant", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "prom_enh_variant", "core_prom_variant", "prom_variant", "enh_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", NA),
                                             con_sum1 = c("coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "splicing", "splicing", "splicing", "splicing", "ncRNA", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", NA),
                                             con_sum2 = c("coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", NA)
)
classical_alleles.dt <- merge.data.table(classical_alleles.dt, summarize_consequence_types.dt, by="consequence")
rm(summarize_consequence_types.dt)
## custom order (to group coding, splicing, non-coding consequences)
classical_alleles.dt$consequence <- factor(classical_alleles.dt$consequence, levels = c("start_lost", "stop_gained", "stop_lost", "frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "coding_sequence_variant", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "prom_enh_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", NA))
classical_alleles.dt$con_sum1 <- factor(classical_alleles.dt$con_sum1, levels = c("coding", "splicing", "ncRNA", "non-coding", NA))
classical_alleles.dt$con_sum2 <- factor(classical_alleles.dt$con_sum2, levels = c("coding/splicing/ncRNA", "non-coding", NA))
classical_alleles_export.dt$overlap_feature_R <- factor(classical_alleles_export.dt$overlap_feature_R, levels = c("CDS", "CDS_synonymous", "introns", "ncRNAs", "in53"))

## dataset POST manual annotation
#---------------------------------/
## add own manual annotation (2 new columns that group several consequences together, con_sum1, con_sum2)
summarize_consequence_types.dt <- data.table(consequence_postManAnnot  = c("start_lost", "stop_gained", "stop_lost", "frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "coding_sequence_variant", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "prom_enh_variant", "core_prom_variant", "prom_variant", "enh_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", NA),
                                             con_sum1 = c("coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "coding", "splicing", "splicing", "splicing", "splicing", "ncRNA", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", NA),
                                             con_sum2 = c("coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "coding/splicing/ncRNA", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", "non-coding", NA)
)
classical_alleles_post_ManAnnot.dt <- merge.data.table(classical_alleles_post_ManAnnot.dt, summarize_consequence_types.dt, by="consequence_postManAnnot")
rm(summarize_consequence_types.dt)
## custom order (to group coding, splicing, non-coding consequences)
classical_alleles_post_ManAnnot.dt$consequence_postManAnnot <- factor(classical_alleles_post_ManAnnot.dt$consequence_postManAnnot, levels = c("start_lost", "stop_gained", "stop_lost", "frameshift_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "coding_sequence_variant", "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant", "enh_variant", "prom_variant", "core_prom_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", NA))
classical_alleles_post_ManAnnot.dt$con_sum1 <- factor(classical_alleles_post_ManAnnot.dt$con_sum1, levels = c("coding", "splicing", "ncRNA", "non-coding", NA))
classical_alleles_post_ManAnnot.dt$con_sum2 <- factor(classical_alleles_post_ManAnnot.dt$con_sum2, levels = c("coding/splicing/ncRNA", "non-coding", NA))
classical_alleles_post_ManAnnot.dt$overlap_feature_R_postManAnnot <- factor(classical_alleles_post_ManAnnot.dt$overlap_feature_R_postManAnnot, levels = c("CDS", "CDS_synonymous", "introns", "ncRNAs", "enh_variant", "prom_variant", "core_prom_variant", "5_prime_UTR_variant", "3_prime_UTR_variant"))
classical_alleles_post_ManAnnot.dt$ManAnnot_ref_regions <- factor(classical_alleles_post_ManAnnot.dt$ManAnnot_ref_regions, levels = c("prom_enh", "prom_enh_3UTR", "5UTR", "3UTR"))

## dataset reporter bashing studies
#---------------------------------/
## custom order 
reporter_annot$overlap_feature_R <- factor(reporter_annot$overlap_feature_R, levels = c("prom_enh", "prom_enh/introns", "prom_enh/introns/5/3UTR", "3'UTR", "3'UTR/5'UTR"))

## plots
#-------/

## page 1
## to focus on forward genetics alleles we exclude CRISPR-Cas9 generated alleles
p1 <- ggplot(classical_alleles.dt[which(is.na(classical_alleles.dt$engineered)),] %>% dplyr::count(consequence) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=consequence)) +
  geom_bar(stat="identity", alpha=0.8, colour="black", linewidth=0.35, width=0.4) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.25,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.22), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(sprintf("%1.1f", pct*100),"%"))) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.45))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Alleles (count)", x="", title="molecular \n consequence \n - before curation", fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  scale_fill_viridis(discrete = TRUE) 
plots_p1 <- (p1)

## page 2
p2 <- ggplot(classical_alleles_post_ManAnnot.dt[which(is.na(classical_alleles_post_ManAnnot.dt$engineered_postManAnnot)),] %>% dplyr::count(consequence_postManAnnot) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=consequence_postManAnnot)) +
  geom_bar(stat="identity", alpha=0.8, colour="black", linewidth=0.35, width=0.4) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.25,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.22), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(sprintf("%1.2f", pct*100),"%"))) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.45))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Alleles (count)", x="", title="molecular \n consequence \n - after curation", fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  scale_fill_viridis(discrete = TRUE) 
plots_p2 <- (p2)

## page 3
p1 <- ggplot(classical_alleles_post_ManAnnot.dt[which(is.na(classical_alleles_post_ManAnnot.dt$engineered_postManAnnot)),] %>% dplyr::count(type) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=forcats::fct_reorder(type, pct, .desc=TRUE))) +
  geom_bar(stat="identity", alpha=0.5, colour="black", linewidth=0.35, width=0.2) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.20,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.13), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(sprintf("%1.1f", pct*100),"%"))) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.55))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Alleles (count)", x="", title="Mutation types", fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"), 
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  scale_fill_brewer(palette = "RdYlGn") # https://colorbrewer2.org
p2 <- ggplot(classical_alleles_post_ManAnnot.dt[which(is.na(classical_alleles_post_ManAnnot.dt$engineered_postManAnnot)),] %>% dplyr::count(con_sum1) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=con_sum1)) +
  geom_bar(stat="identity", alpha=0.5, colour="black", linewidth=0.35, width=0.2) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.20,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.13), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(sprintf("%1.1f", pct*100),"%"))) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.55))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Alleles (count)", x="", title="Alleles \n", fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  scale_fill_viridis(discrete = TRUE)
p3 <- ggplot(classical_alleles_post_ManAnnot.dt[which(classical_alleles_post_ManAnnot.dt$con_sum1 == "non-coding" & is.na(classical_alleles_post_ManAnnot.dt$engineered_postManAnnot)),] 
             %>% dplyr::count(consequence_postManAnnot) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=consequence_postManAnnot)) +
  geom_bar(stat="identity", alpha=0.8, colour="black", linewidth=0.35, width=0.2) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.13,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.13), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(pct*sum(n), "/", sum(n), "  (", sprintf("%1.1f", pct*100),"%)"))) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.55))) +
  scale_y_continuous(n.breaks=7, expand = expansion(mult = c(0, 0.05))) + # scale_y_continuous(breaks=c(0,25,50)) +
  labs(y="Alleles (count)", x="", title=paste0("Alleles, \n non-coding \n"), fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  # scale_fill_manual(values = c("#ffa600","#ffbf00","#ffe100","#ffff00")) # scale_fill_brewer(palette = "YlOrBr")
  scale_fill_manual(values = c("#c7eba1","#d9ee95","#d9ee95","#d9ee95","#ecf18c","#fef392")) # scale_fill_brewer(palette = "YlOrBr")

## now plot number of studies (+include CRISPR-Cas9 studies)
p4 <- ggplot(unique(classical_alleles_post_ManAnnot.dt[,c("ManAnnot_ref_regions","con_sum1","ManAnnot_reference", "ManAnnot_ref_year")][which(classical_alleles_post_ManAnnot.dt$con_sum1 == "non-coding" & classical_alleles_post_ManAnnot.dt$consequence_postManAnnot != c("synonymous_variant")),]) 
             %>% dplyr::count(ManAnnot_ref_regions) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=ManAnnot_ref_regions)) +
  geom_bar(stat="identity", alpha=0.8, colour="black", linewidth=0.35, width=0.2) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.13,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.13), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(pct*sum(n), "/", sum(n), "  (", sprintf("%1.1f", pct*100),"%)"))) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.55))) +
  # scale_y_continuous(breaks=c(0,25,50)) +
  scale_y_continuous(n.breaks=8, expand = expansion(mult = c(0, 0.05))) +
  labs(y="Studies (count)", x="", title=paste0("Allele \n studies (incl. CRISPR) \n"), fill="type:") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  # scale_fill_manual(values = c("#ffbf00","#ffe100","#ffff00"))
  scale_fill_manual(values = c("#d9ee95","#d9ee95","#ecf18c","#fef392"))
plots_p3 <- (p1|p2)/(p3|p4)

## page 4, alleles & reporter studies - variants/study, sequence function, distance etc.
## panel: reporter studies overview
p1 <- ggplot(reporter_annot %>% dplyr::count(overlap_feature_R) %>% mutate(pct=n/sum(n)), aes(x="", n, fill=overlap_feature_R)) +
  geom_bar(stat="identity", alpha=0.8, colour="black", linewidth=0.35, width=0.2) +
  geom_text_repel(box.padding=0.07, size=3, xlim=c(1.17,5), max.overlaps=Inf, 
                  position=position_stacknudge(vjust = 0.5, x = 0.11), # seems to not work after latest updates
                  aes(y=pct*sum(n), label=paste0(pct*sum(n), "/", sum(n), "  (", sprintf("%1.1f", pct*100),"%)"))) +
  scale_x_discrete(expand = expansion(mult = c(0.17, 0.55))) +
  scale_y_continuous(n.breaks=8, expand = expansion(mult = c(0, 0.05))) +
  labs(y="Studies (count)", x="", title=paste0("reporter bashing \n studies "), fill="non-coding sequence(s):") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5)) + 
  scale_fill_manual(values = c("#d9ee95","#d9ee95","#d9ee95","#fef392","#fef392"))

## panel: variants/study
## classical alleles: sum up variants/study
variants_per_study.dt <- classical_alleles_post_ManAnnot.dt[overlap_feature_R_postManAnnot %in% c("enh_variant","prom_variant","core_prom_variant","5_prime_UTR_variant","3_prime_UTR_variant")][
  ,c("ManAnnotIn53_variants_tested", "ManAnnot_reference", "production_method_postManAnnot")][
    ,(ManAnnotIn53_variants_tested = sum(as.numeric(ManAnnotIn53_variants_tested))), by = c("ManAnnot_reference", "production_method_postManAnnot") ]
variants_per_study.dt <- setnames(variants_per_study.dt, "V1", "variants_per_study")
variants_per_study.dt$type <- "endogenous"
variants_per_study.dt[production_method_postManAnnot == "CRISPR_Cas9"]$type <- "endogenous_CRISPR"

## subset reporter studies
DT <- reporter_annot[,c("ManAnnotIn53_variants_tested", "ManAnnot_reference")]
DT <- setnames(DT, "ManAnnotIn53_variants_tested", "variants_per_study")
DT$type <- "reporter"

## merge classical alleles & reporter studies
variants_per_study.dt <- rbind(variants_per_study.dt[,-c("production_method_postManAnnot")], DT)
rm(DT)

p2 <- ggplot(variants_per_study.dt, aes(y=variants_per_study, x=factor(type), group=type)) +
  geom_violin(trim=TRUE, color=NA, fill="#d1d1d1") + # bw=1
  geom_boxplot(width=0.07, color="black", size=0.4, outlier.shape=NA, alpha=0.5) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, stackratio=2, color=NA, fill="black", alpha=0.25) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(y="Tested variants per study (count)",
       x="",
       title="Variants per study", fill="category") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line.x = element_line(lineend = "square"),
        axis.line.y = element_line(lineend = "square"),
        plot.title=element_text(hjust = 0.5))

# install.packages("see")
# library(see)
# p2 <- ggplot(variants_per_study.dt, aes(y=variants_per_study, x=factor(type), group=type, fill=type)) +
#   geom_violinhalf(trim=FALSE, color=NA, fill="#d1d1d1", position = position_nudge(x=0.1,y=0)) + # bw=1
#   geom_boxplot(width=0.07, color="black", size=0.7, outlier.shape=NA) +
#   geom_dotplot(binaxis='y', stackdir='up', dotsize=0.3, color=NA, fill="black", alpha=0.6, position=position_nudge(x=0.1)) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   labs(y="variants per study",
#        x="type",
#        title="variants per study", fill="category") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.key.size = unit(0.3, 'cm'),
#         text = element_text(size=10),
#         axis.text = element_text(color = "black"),
#         axis.ticks = element_line(color = "black"),
#         axis.line.x = element_line(lineend = "square"),
#         axis.line.y = element_line(lineend = "square"),
#         plot.title=element_text(hjust = 0.5))

rm(variants_per_study.dt)

plots_p4 <- (p1|plot_spacer())/(p2)


## panel: sequence function
sequence_function.dt <- classical_alleles_post_ManAnnot.dt[overlap_feature_R_postManAnnot %in% c("enh_variant","prom_variant","core_prom_variant","5_prime_UTR_variant","3_prime_UTR_variant")][
  ManAnnotIn53_mutation_expression_consequence %in% c("down","up","both")][
    ,c("ManAnnotIn53_mutation_expression_consequence", "ManAnnot_reference", "ManAnnotIn53_gene", "overlap_feature_R_postManAnnot")]
sequence_function.dt[overlap_feature_R_postManAnnot %in% c("enh_variant","prom_variant","core_prom_variant")]$overlap_feature_R_postManAnnot <- "prom_enh_variant"
sequence_function.dt <- unique(sequence_function.dt)
sequence_function.dt$type <- "A" # for "Allele"

DT <- reporter_annot[
  ,c("ManAnnotIn53_mutation_expression_consequence", "ManAnnot_reference", "overlap_feature_R", "gene")]
DT <- setnames(DT, "gene", "ManAnnotIn53_gene")
DT <- DT[overlap_feature_R %in% c("prom_enh", "prom_enh/introns/5/3UTR", "prom_enh/introns"), overlap_feature_R_postManAnnot := "prom_enh_variant"][overlap_feature_R %in% c("3'UTR", "3'UTR/5'UTR"),overlap_feature_R_postManAnnot := "3_prime_UTR_variant"][, -c("overlap_feature_R")]
DT$type <- "R" # for "Reporter"

sequence_function.dt <- rbind(sequence_function.dt, DT)
rm(DT)
sequence_function.dt <- sequence_function.dt[overlap_feature_R_postManAnnot %in% c("prom_enh_variant"), overlap_feature_R_postManAnnot := "prom/enh"][overlap_feature_R_postManAnnot %in% c("3_prime_UTR_variant"), overlap_feature_R_postManAnnot := "3'UTR"][overlap_feature_R_postManAnnot %in% c("5_prime_UTR_variant"), overlap_feature_R_postManAnnot := "5'UTR"]

## count number of studies (n= ?)
number_of_studies.dt <- sequence_function.dt[, .N, by=.(overlap_feature_R_postManAnnot, type)]
number_of_studies.dt
## count rows
sequence_function.dt <- sequence_function.dt[, .N, by=.(ManAnnotIn53_mutation_expression_consequence, overlap_feature_R_postManAnnot, type)]
## in pct
sequence_function.dt <- sequence_function.dt[ , pct := N / sum( N ) * 100 , by = .(overlap_feature_R_postManAnnot, type) ]

## custom order 
sequence_function.dt$overlap_feature_R_postManAnnot <- factor(sequence_function.dt$overlap_feature_R_postManAnnot, levels = c("prom/enh", "5'UTR", "3'UTR"))
sequence_function.dt$ManAnnotIn53_mutation_expression_consequence <- factor(sequence_function.dt$ManAnnotIn53_mutation_expression_consequence, levels = c("up", "down", "both"))

p1 <- ggplot(sequence_function.dt, aes(x=type, y=pct, fill=ManAnnotIn53_mutation_expression_consequence, group=N)) +
  geom_bar(stat="identity", position = "stack", alpha=0.5, colour="black", size=0.4, width=1) +
  facet_grid( ~ overlap_feature_R_postManAnnot, scales = "free_x", space = "free") + 
  # facet_grid( ~ overlap_feature_R_postManAnnot) + # this keeps an empty space were no values are present (5'UTR, reporter)
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Studies concluding each (%)", x="", title="Mutation: effect \n on gene expression", fill="Gene expression") +
  theme_classic() +
  theme(legend.position="right", 
        legend.key.size = unit(0.3, 'cm'), 
        text = element_text(size=10), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"), 
        axis.line.x = element_line(lineend = "square"), 
        axis.line.y = element_line(lineend = "square"), 
        plot.title=element_text(hjust=0.5), 
        strip.background = element_rect(size=NA),
        strip.text.x = element_text(color = "black")) + 
  scale_fill_manual(values = c("#f48483","#8fc7d6","#cccccc"))
rm(sequence_function.dt)
rm(number_of_studies.dt)

## panel: distance
DT <- reporter_annot[
  ,c("min_distance_to_coding", "max_distance_to_coding", "overlap_feature_R")]
DT <- DT[overlap_feature_R %in% c("prom_enh", "prom_enh/introns/5/3UTR", "prom_enh/introns"), overlap_feature_R_postManAnnot := "prom_enh_variant"][overlap_feature_R %in% c("3'UTR", "3'UTR/5'UTR"),overlap_feature_R_postManAnnot := "3_prime_UTR_variant"][, -c("overlap_feature_R")]
DT <- DT[overlap_feature_R_postManAnnot == "prom_enh_variant"]
DT2 <- DT[, c("min_distance_to_coding")]
DT2 <- setnames(DT2, "min_distance_to_coding", "distance")
DT2$type <- "Reporter_closest"
DT3 <- DT[, c("max_distance_to_coding")]
DT3 <- setnames(DT3, "max_distance_to_coding", "distance")
DT3$type <- "Reporter_furthest"
DT4 <- classical_alleles_post_ManAnnot.dt[overlap_feature_R_postManAnnot %in% c("enh_variant","prom_variant","core_prom_variant")][,c("ManAnnotIn53_distance_to_coding", "ManAnnotIn53_gene", "ManAnnot_reference" )]
DT4 <- unique(DT4)
DT4 <- DT4[,-c("ManAnnotIn53_gene", "ManAnnot_reference")]
DT4 <- setnames(DT4, "ManAnnotIn53_distance_to_coding", "distance")
DT4$type <- "Alleles"
distance.dt <- rbind(DT2, DT3, DT4)

## set 0 values to 1 (because log10 for plot)
distance.dt[distance == 0, distance := 0.1] 

# ## order for plot
# distance.dt$type <- factor(distance.dt$type, levels = c("endogenous", "reporter_min", "reporter_max"))
# 
# p4 <- ggplot(distance.dt, aes(y=as.numeric(distance), x=factor(type), group=type, fill=type)) +
#   geom_violin(trim=TRUE, color=NA, fill="#d1d1d1") + # bw=1
#   geom_boxplot(width=0.07, color="black", size=0.4, outlier.shape=NA, alpha=0.5) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, stackratio=2, color=NA, fill="black", alpha=0.25) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=7),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   labs(y="distance to gene (bp)",
#        x="",
#        title="distance to gene", fill="category") +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.key.size = unit(0.3, 'cm'),
#         text = element_text(size=10),
#         axis.text = element_text(color = "black"),
#         axis.ticks = element_line(color = "black"),
#         axis.line.x = element_line(lineend = "square"),
#         axis.line.y = element_line(lineend = "square"),
#         plot.title=element_text(hjust = 0.5)) + 
#   coord_flip()

## order for plot
distance.dt$type <- factor(distance.dt$type, levels = c("Reporter_furthest", "Reporter_closest", "Alleles"))

p2 <- ggplot(distance.dt, aes(y=as.numeric(distance), x=factor(type), group=type)) +
  geom_violinhalf(trim=TRUE, color=NA, fill="#d1d1d1", position = position_nudge(x=0.1,y=0)) + # bw=1
  geom_boxplot(width=0.15, color="black", size=0.4, outlier.shape=NA, alpha = 0.5, position = position_nudge(x=-0.05,y=0), fill="grey") +
  geom_dotplot(binaxis='y', stackdir='up', dotsize=0.4, stackratio=1.7, color=NA, fill="black", alpha=0.2, position=position_nudge(x=0.1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=7),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(y="Distance to gene (bp)",
       x="",
       title="Distance to gene", fill="category") +
  theme_classic() +
  theme(legend.position = "none",
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line.x = element_line(lineend = "square"),
        axis.line.y = element_line(lineend = "square"),
        plot.title=element_text(hjust = 0.5)) + 
  coord_flip()

plots_p5 <- (p1|p2)/(plot_spacer())

rm(DT, DT2, DT3, DT4)
rm(distance.dt)



pdf("output/plot_classical_alleles.pdf")
print(plots_p1)
print(plots_p2)
print(plots_p3)
print(plots_p4)
print(plots_p5)

dev.off()
rm(plots_p1, plots_p2, plots_p3, plots_p4, plots_p5)
rm(p1,p2,p3,p4)

##++++++++++++++++++++++++++++++++++++++++++++++++++
### 5.3 - automatic browser shots of alleles    ####
##++++++++++++++++++++++++++++++++++++++++++++++++++

## to inspect allele location manually
## with transcriptome annotation, +/- 500 bp 
## using library(Gviz)
## code snippets from van Arensbergen et al. 2019, at https://github.com/vansteensellab/SuRE-ms-scripts

## modify granges so that name of allele can be displayed
classical_alleles_in53.gr$gene <- classical_alleles_in53.gr$Name
classical_alleles_in53.df <- as.data.frame(classical_alleles_in53.gr)

## modify txdb so that gene_name is shown as well (and not only gene_id e.g. WBGene00014472, ...)
gtf2.gr<-gtf.gr
gtf2.gr$gene_id[!is.na(gtf.gr$gene_name)] <- gtf.gr$gene_name[!is.na(gtf.gr$gene_name)] ## alternative1: replace gene_id with gene_name, except for fields were gene_name is NA
# gtf2.gr$gene_id <- paste(gtf.gr$gene_name, gtf.gr$gene_id) ## alternative2: concatenate both columns (because some NAs in gene_name) and replace gene_id
txdb2 <- makeTxDbFromGRanges(gtf2.gr, drop.stop.codons=FALSE)
rm(gtf2.gr)

## list of the alleles for which to plot browser shots 
alleles_to_plot <- classical_alleles_in53.df$Name

## for loop
for(g in alleles_to_plot)
{
  ## extract coordinates (chr, start, end) for allele
  g_coord <- classical_alleles_in53.df[classical_alleles_in53.df$Name==g,]
  
  gen <- "ce11"
  chr <- as.character(g_coord$seqnames)
  start <- (g_coord$start)-1500
  end <- (g_coord$end)+1500
  
  ## define output file
  jpeg(file=paste0("output/classical_alleles_browser_shots/classical_alleles_in53_", g, ".jpg"), res = 150) 
  # pdf(file=paste0("output/classical_alleles_browser_shots/classical_alleles_in53_", g, ".pdf")) 
  
  ## import only part of bigwig file to speed things up: 
  GR<-makeGRangesFromDataFrame(data.frame(chr=chr, start=start, end=end, score=1))
  phyloP_135.bw <- import("data/ce11.phyloP135way.bw", which=GR) #grab only relevant region #bigwig file 29-Jun-2022 09:38)from http://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP135way/
  # phyloP_26.bw <- import("data/ce11.phyloP26way.bw", which=GR) #grab only relevant region #bigwig file 29-Jun-2022 09:38)from http://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP26way/
  
  ## define all the tracks to plot
  ideoTrack <- IdeogramTrack(genome = gen, 
                             chromosome = chr, 
                             lwd=0.7,
                             fontface=1,
                             fontsize=9,
                             fontcolor="dark grey", 
                             col="red",
                             col.border.title="dark grey")
  
  axisTrack <- GenomeAxisTrack(lwd=1, fontface=1, fontsize=9, col = "dark grey", labelPos = "below", fontcolor="dark grey")
  
  axisTrack2 <- GenomeAxisTrack(lwd=1, fontface=1, fontsize=9, col = "dark grey", labelPos = "below", scale = 0.1)
  
  allelesTrack<-GeneRegionTrack(range=classical_alleles_in53.gr,
                                genome=gen,
                                chromosome=chr,
                                start=start,
                                end=end,
                                showId=TRUE,
                                geneSymbols=TRUE,
                                transcriptAnnotation="gene",
                                lwd=0.7,
                                fontface.group=1,
                                fontsize.group=12,
                                fontcolor.group="black",
                                col="black",
                                col.line="black",
                                fill="grey",
                                name="alleles")
  
  txdbTrack<-GeneRegionTrack(range=txdb2,
                             genome=gen,
                             chromosome=chr,
                             start=start,
                             end=end,
                             # showId=TRUE,
                             geneSymbols=TRUE,
                             transcriptAnnotation="gene",
                             lwd=0.7, 
                             fontface.group=1,
                             fontsize.group=12,
                             fontcolor.group="black",
                             col="black",
                             col.line="black",
                             lex=8, 
                             fill="grey",
                             name="txdb")
  
  phyloP_135Track<-DataTrack(range=phyloP_135.bw,
                             genome=gen,
                             chromosome=chr,
                             start=start,
                             end=end,
                             type="l",
                             # window=-1,
                             ylim = c(-4,20),
                             col="grey",
                             col.line="grey",
                             baseline=0,
                             col.baseline="dark grey",
                             fill="grey",
                             name="135")
  
  ## plot
  plotTracks(list(ideoTrack, axisTrack, axisTrack2, allelesTrack, txdbTrack, phyloP_135Track),
             chromosome=chr,
             from=start,
             to=end,
             sizes=c(0.3, 0.3, 0.3, 0.8, 0.8, 0.4),
             fontface.title=1,
             fontsize.title=12,
             fontcolor.title="black",
             col.border.title="transparent",
             background.title="transparent",
             col.axis="black")
  
  dev.off() 
}


## show warnings
warnings()

