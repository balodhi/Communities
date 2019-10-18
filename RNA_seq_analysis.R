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
c("Sample",  "N2_GCB", "ep2_OP50")

# results with different shape of contrasts
# comparisons against control
res.ct.N2GCB = results(dds,   contrast = c("Sample", "N2_GCB", "N2_OP50"))  
res.ct.N2GCB = lfcShrink(dds, contrast = c("Sample", "N2_GCB", "N2_OP50"), res = res.ct.N2GCB, type = 'ashr')

res.ct.ep2GCB = results(dds,  contrast = c("Sample", "ep2_GCB", "N2_OP50")) 
res.ct.ep2GCB = lfcShrink(dds, contrast = c("Sample", "ep2_GCB", "N2_OP50"), res = res.ct.ep2GCB, type = 'ashr')

res.ct.ep2OP50 = results(dds,  contrast = c("Sample", "ep2_OP50", "N2_OP50"))   
res.ct.ep2OP50 = lfcShrink(dds, contrast = c("Sample", "ep2_OP50", "N2_OP50"), res = res.ct.ep2OP50, type = 'ashr')

# comparisons against ep2 OP50
res.ep2 = results(dds, contrast = c("Sample", "ep2_GCB", "ep2_OP50")) 
res.ep2 = lfcShrink(dds, contrast = c("Sample", "ep2_GCB", "ep2_OP50"), res = res.ep2, type = 'ashr')

res.ep2.N2GCB = results(dds, contrast = c("Sample",  "N2_GCB", "ep2_OP50")) 
res.ep2.N2GCB = lfcShrink(dds, contrast = c("Sample",  "N2_GCB", "ep2_OP50"), res = res.ep2.N2GCB, type = 'ashr')







# tidying the results
res.ct.N2GCB.tidy = as_tibble(res.ct.N2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'N2_GCB',
	Contrast_description = 'Comparison of N2_OP50 vs N2_GCB') %>%
	left_join(info) %>%
	# left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.ct.ep2GCB.tidy = as_tibble(res.ct.ep2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'ep2_GCB',
	Contrast_description = 'Comparison of N2_OP50 vs ep2_GCB') %>%
	left_join(info) %>%
	# left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.ct.ep2OP50.tidy = as_tibble(res.ct.ep2OP50, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'Control',
	Contrast = 'ep2_OP50',
	Contrast_description = 'Comparison of N2_OP50 vs ep2_OP50') %>%
	left_join(info) %>%
	# left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())




res.ep2.tidy = as_tibble(res.ep2, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'ep2',
	Contrast = 'ep2_GCB',
	Contrast_description = 'Comparison of ep2_OP50 vs ep2_GCB') %>%
	left_join(info) %>%
	# left_join(kegg) %>%
	select(Contrast_description, Contrast, gene_id, gene_name, everything())


res.ep2.N2GCB.tidy = as_tibble(res.ep2.N2GCB, rownames = 'gene_id') %>% mutate(
	p_adj_stars = gtools::stars.pval(padj),
	Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
	Contrast_type = 'ep2',
	Contrast = 'N2_GCB',
	Contrast_description = 'Comparison of ep2_OP50 vs N2_GCB') %>%
	left_join(info) %>%
	# left_join(kegg) %>%
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
plotMA(res.ct.N2GCB,  ylim = c(-3,3),  alpha = 0.05)
# idx <- identify(res.N2$baseMean, res.N2$log2FoldChange)
# rownames(res.N2)[idx]
quartz.save(file = here('exploration', 'MAplot_Ct_vs_N2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.ct.ep2GCB,  ylim = c(-3,3),  alpha = 0.05)
quartz.save(file = here('exploration', 'MAplot_Ct_vs_ep2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.ct.ep2OP50,  ylim = c(-3,3),  alpha = 0.05)
quartz.save(file = here('exploration', 'MAplot_Ct_vs_ep2OP50.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)


plotMA(res.ep2,  ylim = c(-3,3),  alpha = 0.05)
quartz.save(file = here('exploration', 'MAplot_ep2OP50_vs_ep2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 11)

plotMA(res.ep2.N2GCB,  ylim = c(-3,3),  alpha = 0.05)
quartz.save(file = here('exploration', 'MAplot_ep2OP50_vs_N2GCB.pdf'),
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

quartz.save(file = here('exploration', 'Euclidean_distances_samples.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 8.5)



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
  scale_color_manual(values = c('#FFB531', '#1461B3')) +
  scale_fill_manual(values = c('#FFB531', '#1461B3')) +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme_classic()

quartz.save(file = here('exploration', 'PCA_main_rld.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)







##################################
### Exploratory plots/analysis ###
##################################


# barplot that shows the % of DE genes in the main comparisons
results.complete %>%
	filter(!is.na(padj)) %>%
	mutate(Sig = ifelse(padj < 0.05, 1, 0)) %>%
	group_by(Contrast_description, Contrast, Sig) %>%
	summarise(N = n()) %>%
	mutate(Total = sum(N),
		   Fraction = (N/Total)*100) %>%
	filter(Sig == 1) %>%
	ungroup %>%
	mutate(Contrast_description = str_replace(Contrast_description, 'Comparison of ', '')) %>%
	ggplot(aes(x = Contrast_description, y = Fraction)) +
	geom_bar(stat = 'identity', width = 0.9, aes(fill = Contrast_description)) +
	scale_fill_brewer(palette = "Dark2") + 
	# scale_x_discrete(limits = c('ep2-OP50 vs ep2-GCB', 'ep2-OP50 vs N2-GCB', 'N2-OP50 vs ep2-GCB', 'N2-OP50 vs ep2-OP50', 'N2-OP50 vs N2-GCB')) +	
	scale_y_continuous(limits = c(0, 70), breaks = c(0, 20, 40, 60, 70)) +
	labs(y = '% of DE genes',
         x = 'Condition') +
	theme_light() 


quartz.save(file = here('exploration', 'DEgenes_barplot.pdf'),
    type = 'pdf', dpi = 300, height = 10, width = 10)


# barplot that shows the % of DE genes in the 4 main comparisons


results.complete %>%
	filter(!is.na(padj)) %>%
	# filter(abs(log2FoldChange) >= 1) %>%
	# mutate(Sig = ifelse(padj < 0.05 && abs(log2FoldChange > 1), 1, 0)) %>%
	mutate(Sig = ifelse(padj < 0.05, 1, 0)) %>%
	group_by(Contrast_description, Contrast, Direction, Sig) %>%
	summarise(N = n()) %>%
	group_by(Contrast_description, Contrast) %>%
	mutate(Total = sum(N),
		   Fraction = round((N/Total)*100,1)) %>%
	group_by(Contrast_description, Contrast,  Sig) %>%
	arrange(desc(Direction)) %>%
	mutate(label_ypos = cumsum(Fraction)) %>%
	filter(Sig == 1) %>%
	ungroup %>%
	mutate(Contrast_description = str_replace(Contrast_description, 'Comparison of ', '')) %>%
	ggplot(aes(x = Contrast_description, y = Fraction, fill = Direction)) +
	geom_bar(stat = 'identity', width = 0.9) +
	geom_text(aes(y = label_ypos, label = Fraction), vjust = 1.6, size = 3.5) +
	scale_fill_brewer(palette = "Dark2") + 
	# scale_x_discrete(limits = c('N2', 'ep2', 'OP50', 'GCB')) +	
	scale_y_continuous(limits = c(0, 70), breaks = c(0, 20, 40, 60, 70)) +
	labs(y = '% of DE genes',
         x = 'Condition') +
	theme_light()

quartz.save(file = here('exploration', 'DEgenes_direction_barplot.pdf'),
    type = 'pdf', dpi = 300, height = 10, width = 10)






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


# WBGene00020343 = acs-3
# WBGene00003623 = nhr-25
gns = c('WBGene00016943','WBGene00009221', 'WBGene00020343', 'WBGene00000133')

# WBGene00001399 = fat-7
# WBGene00001477 = fmo-2
# WBGene00001500 = ftn-1
# WBGene00001851 = hif-1
# WBGene00011737 = sqst-1
# WBGene00001393 = fat-1
# WBGene00020930 = hlh-30
# WBGene00001178 = egl-9
# WBGene00007653 = cysl-1

gns = c('WBGene00001399', 'WBGene00001477', 'WBGene00001500', 'WBGene00001851', 'WBGene00011737', 'WBGene00001393', 'WBGene00020930', 'WBGene00001178', 'WBGene00007653')


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

quartz.save(file = here('exploration', 'boxplot_hlh30_path.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)




#####################
### Volcano plots ###
##################### 

require(EnhancedVolcano)

# datasets to plot
# res.ct.N2GCB 
# res.ct.ep2GCB
# res.ct.ep2OP50
# res.ep2 
# res.ep2.N2GCB

EnhancedVolcano(res.ct.N2GCB,
    lab = rownames(res.ct.N2GCB),
    x = 'log2FoldChange',
    y = 'padj',
    # xlim = c(-8, 8),
    FCcutoff = 2,
    colAlpha = 1,
    subtitle = '',
    title = 'Control vs N2-GCB')

quartz.save(file = here('exploration', 'Volcano_Ct_vs_N2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 9, width = 8)


EnhancedVolcano(res.ct.ep2GCB,
    lab = rownames(res.ct.ep2GCB),
    x = 'log2FoldChange',
    y = 'padj',
    # xlim = c(-8, 8),
    FCcutoff = 1,
    colAlpha = 1,
    subtitle = '',
    title = 'Control vs ep2-GCB')

quartz.save(file = here('exploration', 'Volcano_Ct_vs_ep2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 9, width = 8)

EnhancedVolcano(res.ct.ep2OP50,
    lab = rownames(res.ct.ep2OP50),
    x = 'log2FoldChange',
    y = 'padj',
    # xlim = c(-8, 8),
    FCcutoff = 1,
    colAlpha = 1,
    subtitle = '',
    title = 'Control vs ep2-OP50')

quartz.save(file = here('exploration', 'Volcano_Ct_vs_ep2OP50.pdf'),
    type = 'pdf', dpi = 300, height = 9, width = 8)

EnhancedVolcano(res.ep2,
    lab = rownames(res.ep2),
    x = 'log2FoldChange',
    y = 'padj',
    # xlim = c(-8, 8),
    FCcutoff = 1,
    colAlpha = 1,
    subtitle = '',
    title = 'ep2-OP50 vs ep2-GCB')

quartz.save(file = here('exploration', 'Volcano_ep2OP50_vs_ep2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 9, width = 8)


EnhancedVolcano(res.ep2.N2GCB,
    lab = rownames(res.ep2.N2GCB),
    x = 'log2FoldChange',
    y = 'padj',
    # xlim = c(-8, 8),
    FCcutoff = 1,
    colAlpha = 1,
    subtitle = '',
    title = 'Control vs N2-GCB')

quartz.save(file = here('exploration', 'Volcano_ep2OP50_vs_N2GCB.pdf'),
    type = 'pdf', dpi = 300, height = 9, width = 8)





#####################
### Venn diagrams ###
##################### 
 

require(UpSetR)

# original comparation
ct.N2GCB.gns = res.ct.N2GCB.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ct.ep2GCB.gns = res.ct.ep2GCB.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ct.ep2OP50.gns = res.ct.ep2OP50.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ep2.N2GCB.gns = res.ep2.N2GCB.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector



listInput = list(
						'Ctr vs N2_GCB' = ct.N2GCB.gns, 
						'Ctr vs ep2_GCB' = ct.ep2GCB.gns, 
						'Ctr vs ep2_OP50' = ct.ep2OP50.gns,
						'ep2_OP50 vs ep2_GCB' = ep2.gns,
						'ep2_OP50 vs N2_GCB' = ep2.N2GCB.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")
# upset(fromList(listInput), group.by = "sets")
quartz.save(file = here('exploration', 'UpSet_complete.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



# original comparation, 1 log2FC
ct.N2GCB.gns = res.ct.N2GCB.tidy %>%     filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ct.ep2GCB.gns = res.ct.ep2GCB.tidy %>%   filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ct.ep2OP50.gns = res.ct.ep2OP50.tidy %>% filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% 				 filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ep2.N2GCB.gns = res.ep2.N2GCB.tidy %>%   filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector



listInput = list(
						'Ctr vs N2_GCB' = ct.N2GCB.gns, 
						'Ctr vs ep2_GCB' = ct.ep2GCB.gns, 
						'Ctr vs ep2_OP50' = ct.ep2OP50.gns,
						'ep2_OP50 vs ep2_GCB' = ep2.gns,
						'ep2_OP50 vs N2_GCB' = ep2.N2GCB.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")
# upset(fromList(listInput), group.by = "sets")
quartz.save(file = here('exploration', 'UpSet_log2FC1_complete.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)




# original comparation, 1 log2FC
ct.N2GCB.gns = res.ct.N2GCB.tidy %>%     filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ct.ep2GCB.gns = res.ct.ep2GCB.tidy %>%   filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ct.ep2OP50.gns = res.ct.ep2OP50.tidy %>% filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% 				 filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector
ep2.N2GCB.gns = res.ep2.N2GCB.tidy %>%   filter(padj < 0.05) %>% select(gene_id) %>% t %>% as.vector



listInput = list(
						'Ctr vs N2_GCB' = ct.N2GCB.gns, 
						'Ctr vs ep2_GCB' = ct.ep2GCB.gns, 
						'Ctr vs ep2_OP50' = ct.ep2OP50.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")
# upset(fromList(listInput), group.by = "sets")
quartz.save(file = here('exploration', 'UpSet_controls.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)





# original comparation, 1 log2FC
ct.N2GCB.gns = res.ct.N2GCB.tidy %>%     filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ct.ep2GCB.gns = res.ct.ep2GCB.tidy %>%   filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ct.ep2OP50.gns = res.ct.ep2OP50.tidy %>% filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ep2.gns = res.ep2.tidy %>% 				 filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector
ep2.N2GCB.gns = res.ep2.N2GCB.tidy %>%   filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% select(gene_id) %>% t %>% as.vector



listInput = list(
						'Ctr vs N2_GCB' = ct.N2GCB.gns, 
						'Ctr vs ep2_GCB' = ct.ep2GCB.gns, 
						'Ctr vs ep2_OP50' = ct.ep2OP50.gns)

# generate overall differences plot
upset(fromList(listInput), order.by = "freq")
# upset(fromList(listInput), group.by = "sets")
quartz.save(file = here('exploration', 'UpSet_log2FC1_controls.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



















###########################
### Enrichment analysis ###
###########################


library(KEGG.db)
library(goseq)
# BiocManager::install("TxDb.Celegans.UCSC.ce11.refGene")
library(TxDb.Celegans.UCSC.ce11.refGene)
# BE CAREFUL: it masks select from dplyr
library(org.Ce.eg.db)

### this is a test to see the differentially enriched genes in https://amp.pharm.mssm.edu/WormEnrichr/
OP50.gns = res.OP50.tidy %>% dplyr::filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% dplyr::select(gene_name) %>% t %>% as.vector
write.table(OP50.gns, 'OP50_genes.txt', quote = FALSE, col.names = F, row.names = F)

N2.gns = res.N2.tidy %>% dplyr::filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% dplyr::select(gene_name) %>% t %>% as.vector
write.table(N2.gns, 'N2_genes.txt', quote = FALSE, col.names = F, row.names = F)




# KEGG databases
# get databases for genes and pathways from KEGG
kegg.links.entrez = limma::getGeneKEGGLinks('cel', convert = TRUE) 
kegg.links.ids = limma::getGeneKEGGLinks('cel')
path.ids = limma::getKEGGPathwayNames('cel', remove.qualifier = TRUE)
kegg.links = cbind(kegg.links.entrez, kegg.links.ids[,1])
colnames(kegg.links) = c('entrezid', 'PathwayID', 'KEGG_genes')

kegg.links = kegg.links %>% 
	as_tibble %>% 
	mutate(entrezid = as.integer(entrezid)) %>% 
	left_join(path.ids) %>%
	mutate(PathwayID = str_replace(PathwayID, 'path:cel', ''))



###
# FOR GOSEQ ANALYSIS
###

# remove the NAs
resdat = res.OP50[complete.cases(res.OP50$padj),]

degenes = as.integer(resdat$padj<0.05)
degenes = as.integer(resdat$padj<0.05 & abs(resdat$log2FoldChange) > 1)
names(degenes) = rownames(resdat)

# remove duplicate gene names
degenes = degenes[match(unique(names(degenes)),names(degenes))]

# table(degenes)

# the database we are using is this, let's extract the gene length
# ahEdb

# get the gene lengths
# txsByGene = transcriptsBy(ahEdb,"gene")
# lengthData = median(width(txsByGene))

# lengthData should be the same length as degenes
# test for gene list length
# length(intersect(names(lengthData), names(degenes)))

lengthDataBias = lengthData[names(degenes)]

# Fitting the probability weighting function (PWF) 
pwf = nullp(degenes, bias.data = lengthDataBias, plot.fit = FALSE)
# plot the PWF
plotPWF(pwf)


# Wallenius approximation for enrichment
GO.wall = goseq(pwf, "ce11", "ensGene")
GO.wall.kegg = goseq(pwf, "ce11", "ensGene", test.cats = "KEGG") # kegg categories


# random samplig approx., it takes a bit longer
GO.samp = goseq(pwf, "ce11", "ensGene", method = "Sampling", repcnt = 1000)


# compare the p-values from both methods
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),
 xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
 xlim=c(-3,0))
abline(0,1,col=3,lty=2)

# extract over-represented GOs after adjusting p-value with fdr
enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "fdr") < .05]

#
capture.output(for(go in enriched.GO[1:length(enriched.GO)]){
	print(GO.db::GOTERM[[go]])
	cat("--------------------------------------\n")
})


# this code needs this
ahEdb = ahDb[[1]]
# function for enrichment analysis
GO.enrich = function(data, FC = 0) {
	print('Starting analysis!')
	# remove the NAs
	resdat = data[complete.cases(data$padj),]
	degenes = as.integer(resdat$padj < 0.05 & abs(resdat$log2FoldChange) > FC)
	names(degenes) = rownames(resdat)

	# remove duplicate gene names
	degenes = degenes[match(unique(names(degenes)),names(degenes))]

	# get the gene lengths
	txsByGene = transcriptsBy(ahEdb,"gene")
	lengthData = median(width(txsByGene))
	# filter names
	lengthDataBias = lengthData[names(degenes)]

	print('Fitting the probability weighting function:')
	# Fitting the probability weighting function (PWF) 
	pwf = nullp(degenes, bias.data = lengthDataBias, plot.fit = FALSE)
	# plot the PWF
	plotPWF(pwf)
	print('Applying goseq analysis to both GOs and KEGG')
	# Wallenius approximation for enrichment
	GO.wall = goseq(pwf, "ce11", "ensGene"); GO.wall = as_tibble(GO.wall)
	GO.wall.kegg = goseq(pwf, "ce11", "ensGene", test.cats = "KEGG"); GO.wall.kegg = as_tibble(GO.wall.kegg)
	# extract over-represented GOs after adjusting p-value with fdr
	enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "fdr") < .05]
	enriched.KEGG = GO.wall.kegg$category[p.adjust(GO.wall.kegg$over_represented_pvalue, method = "fdr") < .05]
	print('Giving some human input to the GO categories, and saving everything!')
	# enriched GOs with info
	enriched.GO.info = capture.output(for(go in enriched.GO[1:10]){
		print(GO.db::GOTERM[[go]])
		cat("--------------------------------------\n")
	})
	out.list = list(
		Enriched.KEGG = GO.wall.kegg,
		Enriched.GO = GO.wall,
		Enriched.GO.info = enriched.GO.info)
	return(out.list)
}

##
# Calculate the enriched funcions in each result!!
# this takes a bit of time 

N2.enrich = GO.enrich(res.N2, FC = 1)
ep2.enrich = GO.enrich(res.ep2, FC = 1)
OP50.enrich = GO.enrich(res.OP50, FC = 1)
GCB.enrich = GO.enrich(res.GCB, FC = 1)


### PATHWAY ENRICHMENT ANALYSIS
# tidy results a bit
N2.kegg = N2.enrich$Enriched.KEGG %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'N2') %>%
	left_join(kegg.links %>% 
			dplyr::select(PathwayID, Description) %>% 
			unique, by = c('category' = 'PathwayID'))

ep2.kegg = ep2.enrich$Enriched.KEGG %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'ep2') %>%
	left_join(kegg.links %>% 
			dplyr::select(PathwayID, Description) %>% 
			unique, by = c('category' = 'PathwayID'))

OP50.kegg = OP50.enrich$Enriched.KEGG %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'OP50') %>%
	left_join(kegg.links %>% 
			dplyr::select(PathwayID, Description) %>% 
			unique, by = c('category' = 'PathwayID'))

GCB.kegg = GCB.enrich$Enriched.KEGG %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'GCB') %>%
	left_join(kegg.links %>% 
			dplyr::select(PathwayID, Description) %>% 
			unique, by = c('category' = 'PathwayID'))

# join everything in the same place
res.kegg = rbind(N2.kegg, ep2.kegg, OP50.kegg, GCB.kegg)

# how many enriched pathways we have
res.kegg %>%
	dplyr::filter(padj <= 0.05) %>%
	group_by(Contrast) %>%
	summarise(N = n())

res.kegg %>% filter(padj <= 0.05)

# filter the names of the enriched pathways
de.path = res.kegg %>%
	dplyr::filter(padj <= 0.05, !is.na(Description)) %>%
	select(Description) %>% unique %>% t %>% as.vector %>% sort


# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(), panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
            axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
            strip.text = element_text(colour = "black", face = "bold", 
                size = 10), axis.text.x = element_text(face = "bold", 
                colour = "black", size = 10, angle = 90, hjust = 1))


# plot enrichment p-values
res.kegg %>% 
	filter(Description %in% de.path) %>%
	mutate(logFDR = ifelse(-log10(padj)<0,0,-log10(padj)),
		   logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
		   Description = factor(Description, levels = de.path)) %>%
	ungroup %>%
	ggplot(aes(x = Contrast, y = Description)) +
		geom_tile(aes(fill = logFDRbin)) +
		scale_fill_manual(values = enrcols) + 
		labs(fill = 'FDR') + 
		p.theme


quartz.save(file = here('summary', 'KEGG_pathways_enrichment.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 7)


# # Get the mapping from ENSEMBL 2 Entrez
# en2eg = as.list(org.Ce.egENSEMBL2EG)
# # Get the mapping from Entrez 2 KEGG
# eg2kegg = as.list(org.Ce.egPATH)
# # Define a function which gets all unique KEGG IDs
# # associated with a set of Entrez IDs
# grepKEGG = function(id,mapkeys) {
# 	unique(unlist(mapkeys[id],use.names=FALSE))
# }
# # Apply this function to every entry in the mapping from
# # ENSEMBL 2 Entrez to combine the two maps
# kegg = lapply(en2eg,grepKEGG,eg2kegg)
# head(kegg)


# gns = c()
# for (i in 1:length(kegg)){
# 	if('00903' %in% kegg[[i]]){
# 		gns = c(gns, names(kegg[i]))
# 	}
# }

# info %>% filter(gene_id %in% gns)



###
### GO TERMS ANALYSIS



### PATHWAY ENRICHMENT ANALYSIS
# tidy results a bit
N2.go = N2.enrich$Enriched.GO %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'N2')

ep2.go = ep2.enrich$Enriched.GO %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'ep2') 

OP50.go = OP50.enrich$Enriched.GO %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'OP50') 

GCB.go = GCB.enrich$Enriched.GO %>%
	mutate(padj = p.adjust(over_represented_pvalue, method = 'fdr'),
		   Contrast = 'GCB') 




# join everything in the same place
res.go = rbind(N2.go, ep2.go, OP50.go, GCB.go)

# how many enriched pathways we have
res.go %>%
	dplyr::filter(padj <= 0.05) %>%
	group_by(Contrast) %>%
	summarise(N = n())

res.go %>% filter(padj <= 0.05)

# filter the names of the enriched pathways
de.go = res.go %>%
	dplyr::filter(padj <= 0.05, !is.na(term)) %>%
	select(term) %>% unique %>% t %>% as.vector %>% sort


# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)


# plot enrichment p-values
ont = 'BP'
res.go %>% 
	filter(term %in% de.go, ontology %in% ont) %>%
	mutate(logFDR = ifelse(-log10(padj)<0,0,-log10(padj)),
		   logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
		   term = factor(term, levels = de.go)) %>%
	ungroup %>%
	ggplot(aes(x = Contrast, y = term)) +
		geom_tile(aes(fill = logFDRbin)) +
		scale_fill_manual(values = enrcols) + 
		labs(fill = 'FDR') + 
		p.theme


quartz.save(file = here('summary', 'GO_BP_enrichment.pdf'),
    type = 'pdf', dpi = 300, height = 15, width = 8)



# plot enrichment p-values
ont = 'MF'
res.go %>% 
	filter(term %in% de.go, ontology %in% ont) %>%
	mutate(logFDR = ifelse(-log10(padj)<0,0,-log10(padj)),
		   logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
		   term = recode(term, 'oxidoreductase activity, acting on paired donors, with oxidation of a pair of donors resulting in the reduction of molecular oxygen to two molecules of water' = 'oxidoreductase activity, red. of molecular oxygen to 2 x H2O'),
		   term = recode(term, 'oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen' = 'oxidoreductase activity, with incorporation or red of molecular oxygen'),
		   term = factor(term, levels = de.go)) %>%
	ungroup %>%
	ggplot(aes(x = Contrast, y = term)) +
		geom_tile(aes(fill = logFDRbin)) +
		scale_fill_manual(values = enrcols) + 
		labs(fill = 'FDR') + 
		p.theme


quartz.save(file = here('summary', 'GO_MF_enrichment.pdf'),
    type = 'pdf', dpi = 300, height = 10, width = 7)



# plot enrichment p-values
ont = 'CC'
res.go %>% 
	filter(term %in% de.go, ontology %in% ont) %>%
	mutate(logFDR = ifelse(-log10(padj)<0,0,-log10(padj)),
		   logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
		   term = factor(term, levels = de.go)) %>%
	ungroup %>%
	ggplot(aes(x = Contrast, y = term)) +
		geom_tile(aes(fill = logFDRbin)) +
		scale_fill_manual(values = enrcols) + 
		labs(fill = 'FDR') + 
		p.theme


quartz.save(file = here('summary', 'GO_CC_enrichment.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 5)










