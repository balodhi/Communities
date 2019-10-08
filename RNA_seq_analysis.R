### RNA seq analysis from N2 and ep2 worm strains

# In this script we will analyse the RNA seq data collected by Yifan Wu, with 
# the following conditions:
# 	- N2 with OP50
# 	- N2 with GCB
# 	- ep2 with OP50 
# 	- ep2 with GCB

# N2 - OP50 will be considered as a control for our study, and the comparisons will be:
#	 - N2OP50 vs N2 GCB
#	 - N2OP50 vs ep2 OP50
#	 - N2OP50 vs ep2 GCB
#	 - ep2 OP50 vs ep2 GCB

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# maybe use ulimit -s 16384 before start R console

# analysis directory:  "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/RNA_seq"

library(tximport)
library(tidyverse)
library(DESeq2)
# notice that DESeq2 library masks 'rename' function from dplyr 
# library(ensembldb) # use only if you are going to deal with db
library(here)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(openxlsx)

# session options
options(width = 220)



# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

# load kegg tables from wormenrichr
# tidy it a bit, it's a mess when you load it in R
kegg = read.delim("KEGG_2019.txt", header = FALSE) 
kegg = kegg[,-2] 
kegg = as_tibble(kegg)

kegg = kegg %>% 
	gather(PathwayID, gene, -V1) %>% 
	select(-PathwayID) %>%
	filter(gene != '')

names(kegg) = c('PathwayID', 'gene_name')

# prepare a list with file names
files = file.path(dir,"quants", samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist

# create an object with all the reference genes and transcripts
# THIS DOES NOT SEEM TO BE WORKING
# txdb = TxDb.Celegans.UCSC.ce11.refGene::TxDb.Celegans.UCSC.ce11.refGene
# k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
# tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")


# let's make our database from ensembldb 
ah = AnnotationHub::AnnotationHub()
ahDb = AnnotationHub::query(ah, pattern = c("Caenorhabditis elegans", "EnsDb", 97))
ahEdb = ahDb[[1]]
# generate the database 
tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>% 
	tbl_df() %>%
	dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
	unnest

# join transcription info with gene ids and entrezids
info.join = tx2gene.complete %>% 
	tbl_df() %>%
	dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
	left_join(info)

write.csv(info, here('summary','gene_ids_mapping.csv'))

# subset to have tx_id in first column, and gene_id in second
tx2gene = tx2gene.complete[,c(1,7)]

# import quantification data 
txi = tximport(files, type = "salmon", tx2gene = tx2gene)




### starting analysis with DESeq2
# create DESeq data type to be analysed
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, design = ~  Sample)

# # prefilter, but that might not be necessary
# keep = rowSums(counts(ddsTxi)) >= 10
# ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "N2_OP50")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)


### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% tbl_df()

gene_counts = gene_counts %>% 
	gather(Name, counts, N2_GCB_1:ep2_OP50_4) %>% 
	select(gene_id = gene_list, Name, counts) %>%
	mutate(Name = as.factor(Name)) %>%
	left_join(tbl_df(samples), by = 'Name') %>%
	mutate(counts = as.double(counts),
		   gene_id = as.factor(gene_id),
		   Replicate = as.factor(Replicate)) %>% 
	left_join(info) %>%
	mutate(Sample = factor(Sample, levels = c('N2_OP50', 'N2_GCB', 'ep2_OP50', 'ep2_GCB'))) # refactor levels to show them in this order in the plots


# get results and tidy it
res = results(dds) 


# results with different shape of contrasts
# comparisons against control
res.ct.N2GCB = results(dds,   contrast = c("Sample", "N2_OP50" , "N2_GCB"))  
res.ct.N2GCB = lfcShrink(dds, contrast = c("Sample", "N2_OP50" , "N2_GCB"), res = res.ct.N2GCB, type = 'ashr')

res.ct.ep2GCB = results(dds,  contrast = c("Sample",  "N2_OP50", "ep2_GCB")) 
res.ct.ep2GCB = lfcShrink(dds, contrast = c("Sample", "N2_OP50", "ep2_GCB"), res = res.ct.ep2GCB, type = 'ashr')

res.ct.ep2OP50 = results(dds,  contrast = c("Sample",  "N2_OP50", "ep2_OP50"))   
res.ct.ep2OP50 = lfcShrink(dds, contrast = c("Sample", "N2_OP50", "ep2_OP50"), res = res.ct.ep2OP50, type = 'ashr')

# comparisons against ep2 OP50
res.ep2 = results(dds, contrast = c("Sample",   "ep2_OP50", "ep2_GCB")) 
res.ep2 = lfcShrink(dds, contrast = c("Sample", "ep2_OP50", "ep2_GCB"), res = res.ep2, type = 'ashr')

res.ep2.N2GCB = results(dds, contrast = c("Sample",   "ep2_OP50", "N2_GCB")) 
res.ep2.N2GCB = lfcShrink(dds, contrast = c("Sample", "ep2_OP50", "N2_GCB"), res = res.ep2.N2GCB, type = 'ashr')







# tidying the results
res.ct.N2GCB.tidy = as_tibble(res.ct.N2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'N2_GCB',
	Contrast_description = 'Comparison of N2_OP50 vs N2_GCB') %>%
	left_join(info) %>%
	left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.ct.ep2GCB.tidy = as_tibble(res.ct.ep2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'ep2_GCB',
	Contrast_description = 'Comparison of N2_OP50 vs ep2_GCB') %>%
	left_join(info) %>%
	left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.ct.ep2OP50.tidy = as_tibble(res.ct.ep2OP50, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'ep2_OP50',
	Contrast_description = 'Comparison of N2_OP50 vs ep2_OP50') %>%
	left_join(info) %>%
	left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())




res.ep2.tidy = as_tibble(res.ep2, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'ep2',
	Contrast = 'ep2_GCB',
	Contrast_description = 'Comparison of ep2_OP50 vs ep2_GCB') %>%
	left_join(info) %>%
	left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())


res.ep2.N2GCB.tidy = as_tibble(res.ep2.N2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'ep2',
	Contrast = 'N2_GCB',
	Contrast_description = 'Comparison of ep2_OP50 vs N2_GCB') %>%
	left_join(info) %>%
	left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())


results.complete = res.ct.N2GCB.tidy %>% rbind(res.ct.ep2GCB.tidy, res.ct.ep2OP50.tidy, res.ep2.tidy, res.ep2.N2GCB.tidy)

# write results in excel files
list_of_datasets = list('N2_OP50 vs N2_GCB' = res.ct.N2GCB.tidy, 
						'N2_OP50 vs ep2_GCB' = res.ct.ep2GCB.tidy, 
						'N2_OP50 vs ep2_OP50' = res.ct.ep2OP50.tidy,
						'ep2_OP50 vs ep2_GCB' = res.ep2.tidy,
						'ep2_OP50 vs N2_GCB' = res.ep2.N2GCB.tidy)

write.xlsx(list_of_datasets, here('analysis', 'complete_stats.xlsx'), colNames = T, rowNames = F) 



### MA plots for every comparison
plotMA(res.N2,  ylim=c(-3,3),  alpha = 0.05)
# idx <- identify(res.N2$baseMean, res.N2$log2FoldChange)
# rownames(res.N2)[idx]
quartz.save(file = here('summary', 'MAplot_N2.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.ep2,  ylim=c(-3,3),  alpha = 0.05)
quartz.save(file = here('summary', 'MAplot_ep2.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.GCB,  ylim=c(-3,3),  alpha = 0.05)
quartz.save(file = here('summary', 'MAplot_GCB.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.OP50,  ylim=c(-3,3),  alpha = 0.05)
quartz.save(file = here('summary', 'MAplot_OP50.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)



# transofrm data

vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)


# plot differences between different transformation data methods
df = bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
colnames(df)[1:2] <- c("mean", "sd")  

ggplot(df, aes(x = mean, y = sd)) + 
	geom_hex(bins = 100) +
	# coord_fixed() + 
	facet_grid( . ~ transformation) +
	theme_light()


quartz.save(file = here('summary', 'transformation_comparison.pdf'),
    type = 'pdf', dpi = 300, height = 4.5, width = 12)
