### RNA seq analysis from N2 and ep2 worm strains

# In this script we will analyse the RNA seq data collected by Yifan Wu, with 
# the following conditions:
# 	- N2 with OP50
# 	- N2 with GCB
# 	- ep2 with OP50 
# 	- ep2 with GCB

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
library(biobroom)
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
tx2gene.complete = transcripts(ahEdb, return.type="DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>% 
	tbl_df() %>%
	dplyr::select(width, gene_id, gene_name, gene_biotype, description, symbol)

# join transcription info with gene ids
info.join = tx2gene.complete %>% 
	tbl_df() %>%
	dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
	left_join(info)


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



## ploting some data
###
# quick check on a couple genes of interest
# WBGene00045401 = eol-1
# WBGene00018997 = F57B9.3
gns = c('WBGene00045401', 'WBGene00018997')

# acdh-1 and acdh-2
# WBGene00016943 = acdh-1
# WBGene00015894 = acdh-2
gns = c('WBGene00016943', 'WBGene00015894')

# WBGene00009221 = acs-2
# WBGene00010759 = cysl-2 

gns = c('WBGene00009221', 'WBGene00010759')

gene_counts %>%
	filter(gene_id %in% gns) %>%
	ggplot(aes(y = counts, x = Sample)) +
	geom_boxplot(aes(fill = Worm)) +
	# scale_y_continuous(trans='log2') +
	# facet_wrap(~symbol) +
	facet_wrap(~symbol, scales = 'free_y') +
	labs(x = 'Sample',
		 y = 'Normalised counts') +
	theme_classic()

quartz.save(file = here('summary', 'PCA_main.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)



# results with different shape of contrasts, tidy
res.N2 = results(dds,   contrast = c("Sample", "N2_OP50" , "N2_GCB"))  
res.N2 = lfcShrink(dds, contrast = c("Sample", "N2_OP50" , "N2_GCB"), res = res.N2)

res.ep2 = results(dds,  contrast = c("Sample",  "ep2_OP50", "ep2_GCB")) 
res.ep2 = lfcShrink(dds, contrast = c("Sample", "ep2_OP50", "ep2_GCB"), res = res.ep2)

res.GCB = results(dds,  contrast = c("Sample",  "N2_GCB", "ep2_GCB"))   
res.GCB = lfcShrink(dds, contrast = c("Sample", "N2_GCB", "ep2_GCB"), res = res.GCB)

res.OP50 = results(dds, contrast = c("Sample",   "N2_OP50", "ep2_OP50")) 
res.OP50 = lfcShrink(dds, contrast = c("Sample", "N2_OP50", "ep2_OP50"), res = res.OP50)

# tidying the results
res.N2.tidy = as_tibble(res.N2, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast = 'N2',
	Contrast_description = 'Comparison of N2_OP50 vs N2_GCB') %>%
	left_join(info) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything(), -symbol)

res.ep2.tidy = as_tibble(res.ep2, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast = 'ep2',
	Contrast_description = 'Comparison of ep2_OP50 vs ep2_GCB') %>%
	left_join(info) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything(), -symbol)

res.GCB.tidy = as_tibble(res.GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast = 'GCB',
	Contrast_description = 'Comparison of N2_GCB vs ep2_GCB') %>%
	left_join(info) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything(), -symbol)

res.OP50.tidy = as_tibble(res.OP50, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast = 'OP50',
	Contrast_description = 'Comparison of N2_OP50 vs ep2_OP50') %>%
	left_join(info) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything(), -symbol)

results.complete = res.N2.tidy %>% rbind(res.ep2.tidy, res.GCB.tidy, res.OP50.tidy)

# write results in excel files
list_of_datasets = list('N2 OP50_vs_GCB' = res.N2.tidy, 
						'ep2 OP50_vs_GCB' = res.ep2.tidy, 
						'OP50 N2_vs_ep2' = res.OP50.tidy,
						'GCB N2_vs_ep2' = res.GCB.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'), colNames = T, rowNames = F) 



### MA plots for every comparison
plotMA(res.N2,  ylim=c(-2,2),  alpha = 0.05)
# idx <- identify(res.N2$baseMean, res.N2$log2FoldChange)
# rownames(res.N2)[idx]
quartz.save(file = here('summary', 'MAplot_N2.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.ep2,  ylim=c(-2,2),  alpha = 0.05)
quartz.save(file = here('summary', 'MAplot_ep2.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.GCB,  ylim=c(-2,2),  alpha = 0.05)
quartz.save(file = here('summary', 'MAplot_GCB.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.OP50,  ylim=c(-2,2),  alpha = 0.05)
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



###
# Sample distances

sampleDists = dist(t(assay(rld)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
names = colnames(sampleDists) %>%
	str_split('_', simplify = T) %>%
	data.frame %>% tbl_df() %>%
	unite(sample, X1, X2, sep = " - ") %>%
	select(sample) %>%
	t %>% as.vector

colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
		col = colors)

quartz.save(file = here('summary', 'Euclidean_distances_samples.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 9)



##########
# PCA data

# plotPCA(rld, intgroup = c("Bacteria", "Worm"))

pcaData = plotPCA(rld, intgroup = c("Bacteria", "Worm"), returnData = TRUE)
pcaData


# get info for the ellipses
ell = pcaData %>% group_by(Worm, Bacteria) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Bacteria, group = interaction(Worm, Bacteria))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Worm, Bacteria), linetype = Worm), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Worm, Bacteria), 
    linetype = Worm, fill = Bacteria), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme_classic()

quartz.save(file = here('summary', 'PCA_main_rld.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)



###
# PCA data manual way

# # pca_data = t(counts(rld, normalized = TRUE))
# pca_data = t(assay(rld))

# # lets compute the PCA
# res.pca = PCA(pca_data, scale.unit = FALSE, ncp = 5, graph = F)

# # metadata 
# meta_var = samples 

# # extract info about the individuals
# ind = get_pca_ind(res.pca)
# ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Worm, 
#   meta_var$Bacteria)

# colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Worm', 'Bacteria')

# # get ellipses based on the correlation
# getellipse = function(x, y, sc = 1) {
#   as.data.frame(ellipse::ellipse(cor(x, y),
#                                   scale = c(sd(x) * sc, sd(y) * sc),
#                                   centre = c(mean(x), mean(y))))
# }

# # make a data frame from ellipses
# ell = ind_df %>% group_by(Worm, Bacteria) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# # plot!
# ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Bacteria, group = interaction(Worm, Bacteria))) + 
#   geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
#   geom_path(data = ell, aes(x = x, y = y, group = interaction(Worm, Bacteria), linetype = Worm), size = 1) +
#   geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Worm, Bacteria), 
#     linetype = Worm, fill = Bacteria), size = 1, alpha = 0.3) +
#   xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
#   ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
#   ggtitle(paste("PCA of Bacteria:", paste(unique(ell$Bacteria), collapse = '_'))) + # paste inside a paste to collapse strain names
#   theme(plot.title = element_text(hjust = 0.5),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()) +
#   theme_classic()





########## Venn diagrams
### 

require(UpSetR)

# test
# listInput <- list(
# 	one = c('perro', 'casa', 'jardin', 'hogar', 'puto'), 
# 	two = c('hogar', 'puto', 'marte'), 
# 	three = c('casa', 'jardin', 'hogar'))

# upset(fromList(listInput), order.by = "freq")

N2.gns = res.N2.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
GCB.gns = res.GCB.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
OP50.gns = res.OP50.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector

listInput = list(
	N2 = N2.gns,
	ep2 = ep2.gns,
	GCB = GCB.gns,
	OP50 = OP50.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")

quartz.save(file = here('summary', 'UpSet_complete.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)




N2.gns = res.N2.tidy %>%   	 filter(padj < 0.05, Direction == 'Up') %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% 	 filter(padj < 0.05, Direction == 'Up') %>% select(gene_id) %>% t %>% as.vector
GCB.gns = res.GCB.tidy %>% 	 filter(padj < 0.05, Direction == 'Up') %>% select(gene_id) %>% t %>% as.vector
OP50.gns = res.OP50.tidy %>% filter(padj < 0.05, Direction == 'Up') %>% select(gene_id) %>% t %>% as.vector

listInput = list(
	N2 = N2.gns,
	ep2 = ep2.gns,
	GCB = GCB.gns,
	OP50 = OP50.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")

quartz.save(file = here('summary', 'UpSet_Upregulated.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)



N2.gns = res.N2.tidy %>%   	 filter(padj < 0.05, Direction == 'Down') %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% 	 filter(padj < 0.05, Direction == 'Down') %>% select(gene_id) %>% t %>% as.vector
GCB.gns = res.GCB.tidy %>% 	 filter(padj < 0.05, Direction == 'Down') %>% select(gene_id) %>% t %>% as.vector
OP50.gns = res.OP50.tidy %>% filter(padj < 0.05, Direction == 'Down') %>% select(gene_id) %>% t %>% as.vector

listInput = list(
	N2 = N2.gns,
	ep2 = ep2.gns,
	GCB = GCB.gns,
	OP50 = OP50.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")

quartz.save(file = here('summary', 'UpSet_Downregulated.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)






