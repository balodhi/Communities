##
## script to analyse data from Biolog plates
##


library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(PFun)
library(forcats)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here)
library(Rtsne)

options(width = 220)

# setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/biolog")

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


dir.create('summary', showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create('exploration', showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create('analysis', showWarnings = TRUE, recursive = FALSE, mode = "0777")


info = read_csv('Biolog_metabolites.csv',quote = "\"")


# load data
pm1 = read_csv('pm1/Output/Summary.csv', quote = "\"") %>%
    rename(AUC_raw = `750nm_f_AUC`) %>% 
    filter(Strain != 'Marburg') %>%
    mutate(AUC_raw = replace_na(AUC_raw, 0)) %>%
    mutate(Strain = as.character(Strain),
           Strain = recode(Strain, '2783.0' = 'Ecloacae'),
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8]),
           Strain = as.factor(Strain)) %>% #Change Type column coding
    group_by(Strain, Type, Plate, Replicate, Group) %>%
    mutate(AUC = AUC_raw - AUC_raw[Metabolite == 'Negative Control']) %>% # Normalisation agains negative control in each plate
    select(Strain, Type, Replicate, Index, Plate, Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw, AUC) %>%
    ungroup 

    # load data
pm2 = read_csv('pm2/Output/Summary.csv', quote = "\"") %>%
    rename(AUC_raw = `750nm_f_AUC`) %>% 
    filter(Strain != 'MAR') %>%
    mutate(AUC_raw = replace_na(AUC_raw, 0)) %>%
    mutate(Strain = as.character(Strain),
           Strain = recode(Strain, '2783.0' = 'Ecloacae'),
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8]),
           Strain = as.factor(Strain)) %>% #Change Type column coding
    group_by(Strain, Type, Plate, Replicate, Group) %>%
    mutate(AUC = AUC_raw - AUC_raw[Metabolite == 'Negative Control']) %>% # Normalisation agains negative control in each plate
    select(Strain, Type, Replicate, Index, Plate, Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw, AUC) %>%
    ungroup 


data.b = rbind(pm1, pm2) 

data.sum = data.b %>% filter(Metabolite != 'Negative Control') %>%
    group_by(Strain, Class) %>%
    summarise(Score = sum(AUC)) 

# PCA

#Generate infor tables which describes samples
bioinfo = data.b %>% 
    group_by(Strain, Type) %>%
    summarise %>%
    data.frame

rownames(bioinfo) = bioinfo$SampleID


# filter by strain and Plate if you want
pca_b_data = data.b %>%
    filter(Metabolite != 'Negative Control')



# some data frame transformations
pca_b_data = pca_b_data %>%
    select(Strain, MetaboliteU, AUC) %>%
    spread(MetaboliteU, AUC) %>%
    data.frame(check.names = F)

rownames(pca_b_data) = pca_b_data[,1]
pca_b_data[,1] = NULL

# lets compute the PCA
res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 5, graph = F)

# metadata 
meta_var = bioinfo 

# extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], ind$coord[,4], ind$coord[,5], meta_var$Type, 
    meta_var$Strain)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Type', 'Strain')



# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim5, color = Strain)) + 
    geom_point(size = 6, show.legend = NA, alpha = 1) + 
    xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
    ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    theme_classic()

quartz.save(type = 'pdf', 
    file = here('analysis', 'PCA_MAIN_Dim14.pdf'), 
    width = 9, height = 9)



###
### secondary analyses
###

# variables information
var = get_pca_var(res.pca)

### secondary plots
# plot variables
fviz_pca_var(res.pca, col.var = "black")

quartz.save(type = 'pdf', 
    file = here('analysis', 'PCA_variables.pdf'), 
    width = 9, height = 9)

# quality of representation
corrplot::corrplot(var$cos2, is.corr = FALSE,  tl.col = "black", tl.cex = 0.3, tl.srt = 70)

quartz.save(type = 'pdf', 
    file = here('analysis', 'PCA_quality.pdf'), 
    width = 10, height = 49)

p1 = fviz_cos2(res.pca, choice = "var", axes = 1,  top = 50)
p2 = fviz_cos2(res.pca, choice = "var", axes = 2,  top = 50)
p3 = fviz_cos2(res.pca, choice = "var", axes = 3,  top = 50)
p6 = fviz_cos2(res.pca, choice = "var", axes = 4,  top = 50)
p7 = fviz_cos2(res.pca, choice = "var", axes = 5,  top = 50)
p4 = fviz_cos2(res.pca, choice = "var", axes = 1:2,top = 50)
p5 = fviz_cos2(res.pca, choice = "var", axes = 2:3,top = 50)


pdf(file = here('analysis', 'PCA_quality_dimensions.pdf'))
p1 
p2
p3
p6
p7
p4
p5
dev.off()


# Contributions of variables to PC1
p1 = fviz_contrib(res.pca, choice = "var", axes = 1, top = 50)

# Contributions of variables to PC2
p2 = fviz_contrib(res.pca, choice = "var", axes = 2, top = 50)

# Contributions of variables to PC3
p3 = fviz_contrib(res.pca, choice = "var", axes = 3, top = 50)

# Contributions of variables to PC3
p4 = fviz_contrib(res.pca, choice = "var", axes = 4, top = 50)

# Contributions of variables to PC3
p5 = fviz_contrib(res.pca, choice = "var", axes = 5, top = 50)


pdf(file = here( 'analysis' ,"/PCA_contrib.pdf"))
p1 
p2
p3
p4
p5
dev.off()




# contribution of individuals to PCA
p1 = fviz_contrib(res.pca, choice = "ind", axes = 1)
p2 = fviz_contrib(res.pca, choice = "ind", axes = 2)
p3 = fviz_contrib(res.pca, choice = "ind", axes = 3)
p4 = fviz_contrib(res.pca, choice = "ind", axes = 4)
p5 = fviz_contrib(res.pca, choice = "ind", axes = 5)

pdf(file = here('analysis', 'PCA_ind_contr.pdf'))
p1 
p2
p3
p4
p5
dev.off()


p1 = fviz_cos2(res.pca, choice = "ind", axes = 1,  top = 50)
p2 = fviz_cos2(res.pca, choice = "ind", axes = 2,  top = 50)
p3 = fviz_cos2(res.pca, choice = "ind", axes = 3,  top = 50)
p6 = fviz_cos2(res.pca, choice = "ind", axes = 4,  top = 50)
p7 = fviz_cos2(res.pca, choice = "ind", axes = 5,  top = 50)
p4 = fviz_cos2(res.pca, choice = "ind", axes = 1:2,top = 50)
p5 = fviz_cos2(res.pca, choice = "ind", axes = 2:3,top = 50)


pdf(file = here('analysis', 'PCA_ind_quality_dimensions.pdf'))
p1 
p2
p3
p6
p7
p4
p5
dev.off()









###


# filter by strain and Plate if you want
pca_b_data.sum = data.sum



# some data frame transformations
pca_b_data.sum = pca_b_data.sum %>%
    select(Strain, Class, Score) %>%
    spread(Class, Score) %>%
    data.frame(check.names = F)

rownames(pca_b_data.sum) = pca_b_data.sum[,1]
pca_b_data.sum[,1] = NULL

# lets compute the PCA
res.pca = PCA(pca_b_data.sum, scale.unit = TRUE, ncp = 5, graph = F)

# metadata 
meta_var = bioinfo 

# extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], ind$coord[,4], ind$coord[,5], meta_var$Type, 
    meta_var$Strain)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Type', 'Strain')



# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Strain)) + 
    geom_point(size = 6, show.legend = NA, alpha = 1) + 
    xlim(-4, 4) +
    ylim(-4, 4) +
    xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
    ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    theme_classic()

quartz.save(type = 'pdf', 
    file = here('analysis', 'PCA_MAIN_Dim12_byClass.pdf'), 
    width = 9, height = 9)





###############
### Heatmap ###
###############

heat.data.sum = t(scale((pca_b_data.sum), center = TRUE, scale = TRUE))

library(circlize)
col_fun = colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
col_fun(seq(-3, 3))

Heatmap(heat.data.sum, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Metabolite classes",
        column_title_side = "bottom",
        cluster_rows = FALSE)


quartz.save(type = 'pdf', 
    file = here('analysis', 'Heatmap_class.pdf'), 
    width = 5, height = 8)



# 

data.b

heat.data = t(scale((pca_b_data), center = TRUE, scale = TRUE))

library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun(seq(-3, 3))

Heatmap(heat.data, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Metabolite classes",
        column_title_side = "bottom",
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5))


quartz.save(type = 'pdf', 
    file = here('analysis', 'Heatmap_all.pdf'), 
    width = 7, height = 19)



heat.data = t(scale((pca_b_data), center = TRUE, scale = TRUE))

library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun(seq(-3, 3))

Heatmap(heat.data, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Metabolite classes",
        column_title_side = "bottom",
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5))


quartz.save(type = 'pdf', 
    file = here('analysis', 'Heatmap_all.pdf'), 
    width = 7, height = 19)




meta.list = data.b %>% filter(Strain == 'M131', Replicate == 1) %>%
    filter(Metabolite != 'Negative Control') %>%
    select(Metabolite, Class) %>%
    arrange(Metabolite)


# version with random colours

ha = rowAnnotation(Classes = meta.list$Class)

ht = Heatmap(heat.data, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Biolog metabolites",
        column_title_side = "bottom",
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5), 
        left_annotation = ha)





# specify colours for the classes

ha = HeatmapAnnotation(Classes = meta.list$Class,
    col = list(
               Classes = c("alcohol" = "#0C5BB0FF", "carbohydrate" = "#FA6B09FF", "carboxylic acid" = "#15983DFF",
                      "amide" = '#EC579AFF', "polymer" = '#16A08CFF', "amino acid" = '#149BEDFF', 
                      "ester" = '#A1C720FF', "amine" = '#FEC10BFF', "fatty acid" = '#EE0011FF')
    ), which = "row"
)

ht = Heatmap(heat.data, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Biolog metabolites",
        column_title_side = "bottom",
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5), 
        left_annotation = ha)
ht


# col_fun = colsel(9, palette = 'sat1')
# col_fun = yarrr::piratepal("basel")[1:9]

quartz.save(type = 'pdf', 
    file = here('analysis', 'Heatmap_all_classes.pdf'), 
    width = 7, height = 10)


##############################
### Heatmap, clean version ###
##############################

# filter by strain and Plate if you want
heat.data = data.b %>%
    filter(Metabolite != 'Negative Control') %>%
    arrange(Class) %>%
    select(Strain, MetaboliteU, AUC) %>%
    spread(MetaboliteU, AUC) %>%
    data.frame(check.names = F)


# remove first column, and replace matrix row names with it
rownames(heat.data) = heat.data[,1]
heat.data[,1] = NULL

# scale matrix with Z-scores
heat.data = t(scale((heat.data), center = TRUE, scale = TRUE))

# get Metabolite and Classes names, and arrange it by Class
meta.list = data.b %>% filter(Strain == 'M131', Replicate == 1) %>%
    filter(Metabolite != 'Negative Control') %>%
    select(Metabolite, Class) %>%
    arrange(Class)

# arrange matrix
heat.data = heat.data[meta.list$Metabolite,]

# matrix annotation
ha = HeatmapAnnotation(Classes = meta.list$Class,
    col = list(
               Classes = c("alcohol" = "#0C5BB0FF", "carbohydrate" = "#FA6B09FF", "carboxylic acid" = "#15983DFF",
                      "amide" = '#EC579AFF', "polymer" = '#16A08CFF', "amino acid" = '#149BEDFF', 
                      "ester" = '#A1C720FF', "amine" = '#FEC10BFF', "fatty acid" = '#EE0011FF')
    ), which = "row"
)
# heatmap
ht = Heatmap(heat.data, 
        col = col_fun,
        name = "Z-score",
        column_title = "Bacterial strains",
        row_title = "Biolog metabolites",
        column_title_side = "bottom",
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 5), 
        left_annotation = ha)
ht


quartz.save(type = 'pdf', 
    file = here('analysis', 'Heatmap_all_classes.pdf'), 
    width = 7, height = 10)





############################
### t-SNE transformation ###
############################



# transform the data, and then, plot
# change the BW_norm if you want to work with median instead of means

tsne_mat = pca_b_data

# t-SNE model
tsne_model_1 = Rtsne(as.matrix(tsne_mat), check_duplicates = FALSE, pca = TRUE, perplexity = 1, theta = 0.5, dims = 2, max_iter = 1000, num_threads = 6)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
d_tsne_1['Strain'] = rownames(tsne_mat)

## plotting the results without clustering
d_tsne_1 %>%
    ggplot(aes(x = V1, y = V2, colour = Strain)) +
    geom_point(size = 9, alpha = 1) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    ggtitle("t-SNE") +
    theme_light(base_size = 20) 

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/tSNE_worm_glucose.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)














