# community analysis

# libraries

library(tidyverse)
library(colorspace)
library(readxl)
library(SDMTools) # for weighted statistics
library(ggpubr)
library(dbscan)
library(ggnewscale)
library(openxlsx)
library(here)

# main path
# setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/mixes")

# session options
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)


# run only once
dir.create('summary', showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create('analysis', showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create('exploration', showWarnings = TRUE, recursive = FALSE, mode = "0777")



###
# user defined functions

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                  scale = c(sd(x) * sc, sd(y) * sc),
                                  centre = c(mean(x), mean(y))))
}

# function to extract clusters from DBSCAN
dbclust = function(data, eps = 200, minPts = 15) {
	db = dbscan(data, eps = eps, minPts = minPts)
	return(db$cluster)
}




### read files from raw_data

# N2

n2_1 = read_xlsx(here('raw_data' ,'/061219/copas/n2_filt.xlsx'), sheet = 'Summary')
n2_2 = read_xlsx(here('raw_data' ,'/061919/copas/n2_filt.xlsx'), sheet = 'Summary')
# n2_3 = read_xlsx(here('raw_data' ,'/062519/copas/n2_filt.xlsx'), sheet = 'Summary')
n2_4 = read_xlsx(here('raw_data' ,'/070219/copas/n2_filt.xlsx'), sheet = 'Summary')
n2_5 = read_xlsx(here('raw_data' ,'/070919/n2_filt.xlsx'), sheet = 'Summary')
n2_6 = read_xlsx(here('raw_data' ,'/071519/n2_filt.xlsx'), sheet = 'Summary')


# ep2

ep2_1 = read_xlsx(here('raw_data', '/061219/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_2 = read_xlsx(here('raw_data', '/061919/copas/ep2_filt.xlsx'), sheet = 'Summary')
# ep2_3 = read_xlsx(here('raw_data', '/062519/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_4 = read_xlsx(here('raw_data', '/070219/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_5 = read_xlsx(here('raw_data', '/070919/ep2_filt.xlsx'), sheet = 'Summary')
ep2_6 = read_xlsx(here('raw_data', '/071519/ep2_filt.xlsx'), sheet = 'Summary')



# n2_3 = n2_3 %>% filter(!Extinction > 750) %>%
#  	filter(!Extinction > 400 & TOF <1800) %>%
#  	filter(Bacteria %in% c('MG', 'Myb9', 'Myb131'))

# ep2_3 = ep2_3 %>%
# 	filter(Bacteria %in% c('MG', 'Myb9', 'Myb131'))


# sample n2_3 is a bit dirty

# n2_3 %>% filter(!Extinction > 750) %>%
# 	filter(!Extinction > 400 & TOF <1800) %>%
# 	ggplot(aes(x = TOF, y = Extinction, colour = Bacteria)) +
# 	geom_point() +
# 	facet_wrap(~Bacteria)

# n2_3 %>% filter(Bacteria == 'Mix') %>% 
# 	ggplot(aes(x = TOF, y = Extinction)) +
# 	geom_point()

#####################
### cleaning data ###
#####################


# filter wells without bacteria in dataset number 6
n2_6 = n2_6 %>% filter(!is.na(Bacteria))
ep2_6 = ep2_6 %>% filter(!is.na(Bacteria))

n2_6_M53 = n2_6 %>% filter(Bacteria == 'M53') %>% mutate(Bacteria = recode(Bacteria, 'M53' = 'Myb53'))
ep2_6_M53 = ep2_6 %>% filter(Bacteria == 'M53')%>% mutate(Bacteria = recode(Bacteria, 'M53' = 'Myb53'))


# it seems that I swaped M131 samples, so the N2 is in the ep2 dataset, and viceversa
m131_n2 = ep2_5 %>% filter(Bacteria == 'Myb131')
m131_ep2 = n2_5 %>% filter(Bacteria == 'Myb131')

n2_5 = n2_5 %>% filter(!Bacteria == 'Myb131') %>% rbind(m131_n2)
ep2_5 = ep2_5 %>% filter(!Bacteria == 'Myb131') %>% rbind(m131_ep2)



# join them together as a long table

n2 = rbind(n2_1, n2_2, n2_4, n2_5, n2_6_M53)
ep2 = rbind(ep2_1, ep2_2, ep2_4, ep2_5, ep2_6_M53)

# filter by thresholds, might be useful to be less stringent with ep2

n2 = n2 %>% 
	mutate(Group = ifelse(str_detect(Sample, '_g_'), 'GCB', 'Single')) %>%
	filter(TOF > 300, Extinction > 50)
ep2 = ep2 %>%
	mutate(Group = ifelse(str_detect(Sample, '_g_'), 'GCB', 'Single')) %>%
	filter(TOF > 300, Extinction > 50)

# n2[n2$Bacteria %in% c('Mix', 'Mix2783', 'MixGCB', 'MixMarb', 'MixGCBMarb'),]$Group = 'Mix'
# ep2[ep2$Bacteria %in% c('Mix', 'Mix2783', 'MixGCB', 'MixMarb', 'MixGCBMarb'),]$Group = 'Mix'

# establish Mix groups
n2[n2$Bacteria %in% c('Mix', 'MixGCB', 'MixMarb', 'MixGCBMarb'),]$Group = 'Mix'
ep2[ep2$Bacteria %in% c('Mix', 'MixGCB', 'MixMarb', 'MixGCBMarb'),]$Group = 'Mix'

# establish Triplet groups
n2[n2$Bacteria %in% c('GEOP50', 'GEMarb', 'GEM131', 'GEM71', 'GEM181', 'GEM9', 'GEMG'),]$Group = 'Triplet'
ep2[ep2$Bacteria %in% c('GEOP50', 'GEMarb', 'GEM131', 'GEM71', 'GEM181', 'GEM9', 'GEMG'),]$Group = 'Triplet'


# filter by dbscan method


## DBSCAN method

# df = (n2 %>% filter(Bacteria == 'OP50'))[,3:4]
# db = dbscan(df, eps = 200, minPts = 15)
# db.o = optics(df, minPts = 15)
# fviz_cluster(db, df, stand = FALSE, frame = FALSE, geom = "point")




# this loop code will calculate all the 
db.n2 = n2 %>% group_by(Bacteria, Group) %>%
	nest %>%
	mutate(group = map(data, ~dbclust(.x[,3:4]))) %>%
	unnest

db.ep2 = ep2 %>% group_by(Bacteria, Group) %>%
	nest %>%
	mutate(group = map(data, ~dbclust(.x[,3:4]))) %>%
	unnest



# plot how many points have been filtered by the method in monocultures
p1 = db.n2 %>% 
	filter(Group == 'Single') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = db.ep2 %>% 
	filter(Group == 'Single') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', 'DBSCAN_singles_n2_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 16)


# plot how many points have been filtered by the method in co_cultures
p1 = db.n2 %>% 
	filter(Group == 'GCB') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = db.ep2 %>% 
	filter(Group == 'GCB') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', 'DBSCAN_GCB_n2_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 16)



# plot how many points have been filtered by the method in co_cultures
p1 = db.n2 %>% 
	filter(Group == 'Mix') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = db.ep2 %>% 
	filter(Group == 'Mix') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', 'DBSCAN_Mix_n2_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 16)



# plot how many points have been filtered by the method in co_cultures
p1 = db.n2 %>% 
	filter(Group == 'Triplet') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = db.ep2 %>% 
	filter(Group == 'Triplet') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', 'DBSCAN_Triplet_n2_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 16)



#####################
### density plots ###
#####################

# density plots after filtering

p1 = db.n2 %>% 
	filter(group != 0, Group == 'Single') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	xlim(0, 3700) +
	ylim(0, 560) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Single') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	xlim(0, 3700) +
	ylim(0, 560) +
	theme_light()


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', '2d_density_single.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)



# density plots after filtering

p1 = db.n2 %>% 
	filter(group != 0, Group == 'GCB') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	xlim(0, 3700) +
	ylim(0, 560) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'GCB') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	xlim(0, 3700) +
	ylim(0, 560) +
	theme_light()


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', '2d_density_GCB.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)


# density plots after filtering

p1 = db.n2 %>% 
	filter(group != 0, Group == 'Mix') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	xlim(0, 3700) +
	ylim(0, 560) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Mix') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	xlim(0, 3700) +
	ylim(0, 560) +
	theme_light()


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', '2d_density_Mix.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)



# density plots after filtering

p1 = db.n2 %>% 
	filter(group != 0, Group == 'Triplet') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	xlim(0, 3700) +
	ylim(0, 560) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Triplet') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	xlim(0, 3700) +
	ylim(0, 560) +
	theme_light()


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)

quartz.save(file = here('exploration', '2d_density_Triplet.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)

####


p1 = db.n2 %>% 
	filter(group != 0 , Group == 'Single') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	# scale_fill_viridis_c(option = "D") +
	scale_fill_continuous_sequential(palette = "Purple-Blue") +
	facet_wrap(~Bacteria) +
	new_scale('fill') +
	stat_density_2d(data = db.ep2 %>% filter(group != 0, Group == 'Single'), 
		aes(x = TOF, y = Extinction, fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	scale_fill_continuous_sequential(palette = "Heat") +
	# scale_fill_viridis_c(option = "A") +
	# xlim(0, 3200) +
	# ylim(0, 530) +
	theme_light()

p2 = db.n2 %>% 
	filter(group != 0 , Group == 'GCB') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	# scale_fill_viridis_c(option = "D") +
	scale_fill_continuous_sequential(palette = "Purple-Blue") +
	facet_wrap(~Bacteria) +
	new_scale('fill') +
	stat_density_2d(data = db.ep2 %>% filter(group != 0, Group == 'GCB'), 
		aes(x = TOF, y = Extinction, fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	scale_fill_continuous_sequential(palette = "Heat") +
	# scale_fill_viridis_c(option = "A") +
	# xlim(0, 3200) +
	# ylim(0, 530) +
	theme_light()

p3 = db.n2 %>% 
	filter(group != 0 , Group == 'Mix') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	# scale_fill_viridis_c(option = "D") +
	scale_fill_continuous_sequential(palette = "Purple-Blue") +
	facet_wrap(~Bacteria) +
	new_scale('fill') +
	stat_density_2d(data = db.ep2 %>% filter(group != 0, Group == 'Mix'), 
		aes(x = TOF, y = Extinction, fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	scale_fill_continuous_sequential(palette = "Heat") +
	# scale_fill_viridis_c(option = "A") +
	# xlim(0, 3200) +
	# ylim(0, 530) +
	theme_light()


p4 = db.n2 %>% 
	filter(group != 0 , Group == 'Triplet') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	# scale_fill_viridis_c(option = "D") +
	scale_fill_continuous_sequential(palette = "Purple-Blue") +
	facet_wrap(~Bacteria) +
	new_scale('fill') +
	stat_density_2d(data = db.ep2 %>% filter(group != 0, Group == 'Triplet'), 
		aes(x = TOF, y = Extinction, fill = stat(nlevel)), geom = "polygon", alpha = 1) +
	scale_fill_continuous_sequential(palette = "Heat") +
	# scale_fill_viridis_c(option = "A") +
	# xlim(0, 3200) +
	# ylim(0, 530) +
	theme_light()

ggarrange(p1, p2, p3, p4,
          labels = c("Single", "Co-culture GCB", 'Mixes', 'Triplets'),
          ncol = 2, nrow = 2)


quartz.save(file = here('exploration', '2d_density_join.pdf'),
	type = 'pdf', dpi = 300, height = 16, width = 16)



##################################
### Scatter plots and barplots ###
##################################

# first of all, lets join both filtered datasets, and modify them to add more variables
# of interest: bio replicates, technical replicates, worm type...

db.n2['Worm'] = 'N2'
db.ep2['Worm'] = 'ep2'

copas = rbind(db.n2, db.ep2) %>% filter(group != 0) %>% select(-group)


# in single colonies, bacteria Myb27, Myb45, Myb53 and Myb71 lack data from ep2
# manually filtering them, not the best solution but it will have to work either way

m27 = ep2_1 %>% filter(Bacteria =='Myb27', Extinction > 10, TOF > 200)
m45 = ep2_1 %>% filter(Bacteria =='Myb45', Extinction > 10, TOF > 200) %>% filter(Extinction < 40, TOF < 700)
m53 = ep2_1 %>% filter(Bacteria =='Myb53', Extinction > 10, TOF > 200) %>% filter(Extinction < 50, TOF < 800)
m71 = ep2_2 %>% filter(Bacteria =='Myb71', Extinction > 10, TOF > 200, `Source well` %in% c('D9', 'E9', 'F9', 'D10', 'E10', 'F10')) %>%
	filter(Extinction < 50, TOF < 700)

m27['Group'] = 'Single'; m27['Worm'] = 'ep2'
m45['Group'] = 'Single'; m45['Worm'] = 'ep2'
m53['Group'] = 'Single'; m53['Worm'] = 'ep2'
m71['Group'] = 'Single'; m71['Worm'] = 'ep2'

copas = rbind(copas, m27, m45, m53, m71)

## and now, a silly movement: we are dropping m27 and m45 from the analysis

copas = copas %>% filter(!Bacteria %in% c('Myb27', 'Myb45')) %>% filter(!Sample %in% c('GCB_2_2'))


# lets make variables for the technical and biological replicates
# stupid loop, do it once

copas['biorep'] = 0
copas['techrep'] = 0
for (i in 1:dim(copas)[1]){
	copas$biorep[i] = as.integer(tail(strsplit(copas$Sample[i], '_')[[1]],2)[1])
	copas$techrep[i] = as.integer(tail(strsplit(copas$Sample[i], '_')[[1]],2)[2])
}

copas = copas %>% 
	select(Worm, Group, Bacteria, Sample, biorep, techrep, TOF, Extinction) %>%
	mutate(Worm = fct_relevel(Worm, 'N2', 'ep2'))
# save raw data into a file
list_of_datasets = list('copas' = copas)
write.xlsx(list_of_datasets, here('analysis', 'original_data.xlsx'), colNames = T, rowNames = F) 



# summarise data
copas.sum = copas %>%
	group_by(Worm, Group, Bacteria, Sample, biorep, techrep) %>%
	summarise(TOF_mean = mean(TOF),
			  TOF_sd = sd(TOF),
			  TOF_SEM = TOF_sd/sqrt(n()),
			  Ext_mean = mean(Extinction),
			  Ext_sd = sd(Extinction),
			  Ext_SEM = Ext_sd/sqrt(n())) %>%
	group_by(Worm, Group, Bacteria) %>%					# this is for weighted means and sd
	mutate(wTOF = 1/(TOF_sd**2),
		   wExt = 1/(Ext_sd**2)
		   )


# calculate means and sd (weighted and unweighted) of means

copas.sum2 = copas.sum %>%
	group_by(Worm, Group, Bacteria) %>%
	summarise(wTOF_mean = wt.mean(TOF_mean, wTOF),
			  wTOF_SD = wt.sd(TOF_mean, wTOF),
			  wExt_mean = wt.mean(Ext_mean, wExt),
			  wExt_SD = wt.sd(Ext_mean, wExt)) %>%
	ungroup


########################################################
### quick check that pooled sd are quite similar to normal SD

copas.sum3 = copas %>%
	group_by(Worm, Group, Bacteria) %>%
	summarise(TOF_mean = mean(TOF),
			  TOF_sd = sd(TOF),
			  Ext_mean = mean(Extinction),
			  Ext_sd = sd(Extinction)) %>%
	ungroup

dummy = copas %>% filter(Worm == 'ep2', Group == 'GCB', Bacteria == '2783') %>% select(TOF, Sample)
data = data.frame(y = ChickWeight$weight, g = ChickWeight$Diet)
ANOVAreplication::pooled.sd(dummy)

########################################################


# what I was representing before (wTOF_mean and wTOF_SD) is the MEAN and the SEM, but
# I think it will be more helpful to see the standard deviation of the samples rather
# than the standard error of the mean. But, take into account that the pooled sd could
# be a better indicator of the true variation of our samples (which will be kind of 
# similar, I believe)
#
# So, lets join the SD with copas.sum2

copas.sum2 = copas.sum2 %>% left_join(copas.sum3) 


# save summary statistics into a file
list_of_datasets = list('summary_by_sample' = copas.sum, 'summary_by_bacteria' = copas.sum2)
write.xlsx(list_of_datasets, here('analysis', 'sum_stats.xlsx'), colNames = T, rowNames = F) 



########################
# barplot with points and errorbars

# Single populations

gr = 'Single'
lvls = c('GCB', 'OP50', '2783', 'Marburg', 'MG', 'Myb131', 'Myb181', 'Myb9', 'Myb71', 'Myb53') # order in which the factors will appear in the plot
p1 = copas.sum2 %>% 
	filter(Group == gr) %>%
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wTOF_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = copas.sum %>% filter(Group == gr), 
		aes(x = Bacteria, y = TOF_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wTOF_mean - wTOF_SD, ymax = wTOF_mean + wTOF_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#8C5400', '#D68000')) +
	# scale_fill_manual(values = c('#8C5400', '#D68000')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
 	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) + # removes the spaces at the bottom of the barplot
	theme_classic()

p2 = copas.sum2 %>% 
	filter(Group == gr) %>%
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wExt_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = copas.sum %>% filter(Group == gr), 
		aes(x = Bacteria, y = Ext_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wExt_mean - wExt_SD, ymax = wExt_mean + wExt_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#FF260A', '#0656CC')) +
	# scale_fill_manual(values = c('#FF260A', '#0656CC')) +
	labs(x = 'Bacteria',
		y = 'Extinction') +
	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
	theme_classic()

ggarrange(p1, p2,
          labels = c("TOF", "Ext"),
          ncol = 1, nrow = 2)


quartz.save(file = here('exploration', paste0('barplot_means_',gr,'.pdf')),
	type = 'pdf', dpi = 300, height = 10, width = 9)



# GCB mixes

temp.sum2 = copas.sum2 %>% 
	filter(Group == gr) %>%
	rbind(copas.sum2 %>% filter(Group == 'Single', Bacteria == 'GCB'))
temp.sum = copas.sum %>%
	filter(Group == gr) %>%
	rbind(copas.sum %>% filter(Group == 'Single', Bacteria == 'GCB'))

gr = 'GCB'
lvls = c('GCB', 'OP50', '2783', 'Marburg', 'MG', 'Myb131', 'Myb181', 'Myb9', 'Myb71', 'Myb53') # order in which the factors will appear in the plot
p1 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wTOF_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum , 
		aes(x = Bacteria, y = TOF_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wTOF_mean - wTOF_SD, ymax = wTOF_mean + wTOF_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#8C5400', '#D68000')) +
	# scale_fill_manual(values = c('#8C5400', '#D68000')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 2200, label = 'Control') +
 	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) + # removes the spaces at the bottom of the barplot
	theme_classic()

p2 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wExt_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum, 
		aes(x = Bacteria, y = Ext_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wExt_mean - wExt_SD, ymax = wExt_mean + wExt_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#FF260A', '#0656CC')) +
	# scale_fill_manual(values = c('#FF260A', '#0656CC')) +
	labs(x = 'Bacteria',
		y = 'Extinction') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 350, label = 'Control') +
	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
	theme_classic()

ggarrange(p1, p2,
          labels = c("TOF", "Ext"),
          ncol = 1, nrow = 2)


quartz.save(file = here('exploration', paste0('barplot_means_',gr,'.pdf')),
	type = 'pdf', dpi = 300, height = 10, width = 9)


# All bacteria mixes

gr = 'Mix'
temp.sum2 = copas.sum2 %>% 
	filter(Group == gr) %>%
	rbind(copas.sum2 %>% filter(Group == 'Single', Bacteria == 'GCB'))
temp.sum = copas.sum %>%
	filter(Group == gr) %>%
	rbind(copas.sum %>% filter(Group == 'Single', Bacteria == 'GCB'))

lvls = c('GCB', 'Mix', 'MixGCB', 'MixMarb', 'MixGCBMarb') # order in which the factors will appear in the plot

p1 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wTOF_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum , 
		aes(x = Bacteria, y = TOF_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wTOF_mean - wTOF_SD, ymax = wTOF_mean + wTOF_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#8C5400', '#D68000')) +
	# scale_fill_manual(values = c('#8C5400', '#D68000')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 2200, label = 'Control') +
 	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) + # removes the spaces at the bottom of the barplot
	theme_classic()

p2 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wExt_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum, 
		aes(x = Bacteria, y = Ext_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wExt_mean - wExt_SD, ymax = wExt_mean + wExt_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#FF260A', '#0656CC')) +
	# scale_fill_manual(values = c('#FF260A', '#0656CC')) +
	labs(x = 'Bacteria',
		y = 'Extinction') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 350, label = 'Control') +
	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
	theme_classic()

ggarrange(p1, p2,
          labels = c("TOF", "Ext"),
          ncol = 1, nrow = 2)


quartz.save(file = here('exploration', paste0('barplot_means_',gr,'.pdf')),
	type = 'pdf', dpi = 300, height = 10, width = 9)

# Triplets

gr = 'Triplet'
temp.sum2 = copas.sum2 %>% 
	filter(Group == gr) %>%
	rbind(copas.sum2 %>% filter(Group == 'Single', Bacteria == 'GCB'))
temp.sum = copas.sum %>%
	filter(Group == gr) %>%
	rbind(copas.sum %>% filter(Group == 'Single', Bacteria == 'GCB'))

lvls = c('GCB', 'GEM131', 'GEM9', 'GEM71', 'GEMarb', 'GEOP50', 'GEMG', 'GEM181') # order in which the factors will appear in the plot

p1 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wTOF_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum , 
		aes(x = Bacteria, y = TOF_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wTOF_mean - wTOF_SD, ymax = wTOF_mean + wTOF_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#8C5400', '#D68000')) +
	# scale_fill_manual(values = c('#8C5400', '#D68000')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 2200, label = 'Control') +
 	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) + # removes the spaces at the bottom of the barplot
	theme_classic()

p2 = temp.sum2 %>% 
	mutate(Bacteria = factor(Bacteria, levels = lvls)) %>%
	ggplot(aes(x = Bacteria, y = wExt_mean, colour = Worm, group = Worm)) +
	geom_bar(aes(fill = Worm), stat = 'identity', color = 'black', position = position_dodge(), alpha = 0.7) +
	geom_point(data = temp.sum, 
		aes(x = Bacteria, y = Ext_mean, fill = Worm), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.99, size = 1, colour = 'gray10') +
	geom_errorbar(aes(ymin = wExt_mean - wExt_SD, ymax = wExt_mean + wExt_SD), position = position_dodge(0.9), width = 0.1, colour = 'gray20') +
	scale_fill_OkabeIto() +
	# scale_color_manual(values = c('#FF260A', '#0656CC')) +
	# scale_fill_manual(values = c('#FF260A', '#0656CC')) +
	labs(x = 'Bacteria',
		y = 'Extinction') +
	geom_vline(xintercept = 1.5, size = 0.8) + # adds a vertical line to separate control GCB from the other comparisons
	annotate('text', x = 1, y = 350, label = 'Control') +
	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
	theme_classic()

ggarrange(p1, p2,
          labels = c("TOF", "Ext"),
          ncol = 1, nrow = 2)


quartz.save(file = here('exploration', paste0('barplot_means_',gr,'.pdf')),
	type = 'pdf', dpi = 300, height = 10, width = 9)



###### 


# lets calculate what is the difference between ep2 and n2, that could be helpful
# for clarification:
# 	TOF_dif is the difference between N2 and ep2 values
# 	TOF_score is the N2 value divided by the TOF_dif
# we are interested in see which bacteria confer a better recovery from the
# wild phenotype, but also to correct against worms not well developed

dif.sum = copas.sum %>% 
	left_join(copas.sum2 %>% filter(Worm == 'N2') %>% ungroup %>% select(-Worm)) %>% 
	select(-biorep, -techrep, -wTOF, -wExt, -wTOF_SD, -wExt_SD) %>%
	rename(wTOF_N2 = wTOF_mean,
		   wExt_N2 = wExt_mean) %>%
	mutate(TOF_dif = wTOF_N2 - TOF_mean,
		   Ext_dif = wExt_N2 - Ext_mean,
		   TOF_score = wTOF_N2/TOF_dif,
		   Ext_score = wExt_N2/Ext_dif) %>%
	filter(Worm == 'ep2') %>%
	ungroup %>%
	select(-Worm)


# scatter plot of the scores

gr = 'Single'

ell = dif.sum %>%
	filter(Group == gr,!Bacteria %in% c('Myb27', 'Myb45'), !(Sample == 'GCB_2_2' & Group == 'Single')) %>%
	group_by(Bacteria) %>% 
	do(getellipse(.$TOF_score, .$Ext_score, 1)) %>% 
	data.frame

dif.sum %>%
	filter(Group == gr,!Bacteria %in% c('Myb27', 'Myb45'), !(Sample == 'GCB_2_2' & Group == 'Single')) %>%
	ggplot(aes(x = TOF_score, y = Ext_score)) +
	geom_point(aes(colour = Bacteria), size = 5) +
	geom_path(data = ell, aes(x = x, y = y, colour = Bacteria), size = 1) +
	geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Bacteria), fill = Bacteria), size = 1, alpha = 0.3) +
	scale_fill_manual(values = colsel(10, palette = 'sat1')) +
	scale_color_manual(values = colsel(10, palette = 'sat1')) +
	xlab('TOF score') + 
 	ylab('Extinction score') +
	theme_classic()

quartz.save(file = here('exploration', 'scatter_scores_Single.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 8)




####
# Is E cloacae different in ep2 when mixed with GCB?


ecloc = copas %>% filter(Bacteria == '2783', Worm == 'ep2') 
ecloc.sum = copas.sum %>% filter(Bacteria == '2783', Worm == 'ep2')
ecloc.sum2 = copas.sum2 %>% filter(Bacteria == '2783', Worm == 'ep2')


mg = copas %>% filter(Bacteria == 'MG', Worm == 'ep2') 


ggplot(ecloc, aes(y = TOF, x = Group, fill = Group)) +
	geom_boxplot()

model = lm(TOF ~ Group, data = ecloc)
summary(model)

wilcox.test(TOF ~ Group, data = ecloc)

# test over the means
ggplot(ecloc.sum, aes(y = TOF_mean, x = Group, fill = Group)) +
	geom_boxplot()


model = lm(TOF_mean ~ Group, data = ecloc.sum)
summary(model)

wilcox.test(TOF_mean ~ Group, data = ecloc.sum)


##
model = aov(TOF ~ Sample, data = ecloc %>% filter(Group == 'Single'))
summary(model)
TukeyHSD(model)



############
## general test for each condition
############



# function for statistical analyses
apply_lm = function(df, formula){
      lm(data = df, as.formula(formula))      
}

# REMEMBER to change the variable inside mcp function
apply_multilm = function(model){
	multcomp::glht(model, linfct = multcomp::mcp(Sample = 'Tukey'))
}

nested = copas %>%
	mutate(Sample = as.factor(Sample)) %>%
    group_by(Group, Worm, Bacteria) %>%
    nest()


# calculate linear models per sample
res = nested %>%
    mutate(models = map(.f = apply_lm, .x = data, formula = 'TOF ~ Sample')) %>%
	mutate(multcomp = map(.f = apply_multilm, .x = models)) %>%
	mutate(results = map(.f = summary, .x = multcomp)) %>%
	mutate(results = map(.f = tidy, .x = results)) %>%
	dplyr::select(Group, Bacteria, Worm, results) %>%
	unnest

# results summary for TOF, between samples
results.TOF.samples = res %>% rename(Contrast = lhs) %>% 
	mutate(FDR = p.adjust(p.value, method = 'fdr'),
		   p_stars = gtools::stars.pval(p.value),
		   FDR_stars = gtools::stars.pval(FDR))





# test between worm type, by bacteria
# REMEMBER to change the variable inside mcp function
apply_multilm = function(model){
	multcomp::glht(model, linfct = multcomp::mcp(Worm = 'Tukey'))
}
nested = copas %>%
	mutate(Worm = as.factor(Worm)) %>%
    group_by(Group, Bacteria) %>%
    nest()


# calculate linear models per sample
res = nested %>%
    mutate(models = map(.f = apply_lm, .x = data, formula = 'TOF ~ Worm')) %>%
	mutate(multcomp = map(.f = apply_multilm, .x = models)) %>%
	mutate(results = map(.f = summary, .x = multcomp)) %>%
	mutate(results = map(.f = tidy, .x = results)) %>%
	dplyr::select(Group, Bacteria, Worm, results) %>%
	unnest

# results summary for TOF, between samples
results.TOF.samples = res %>% rename(Contrast = lhs) %>% 
	mutate(FDR = p.adjust(p.value, method = 'fdr'),
		   p_stars = gtools::stars.pval(p.value),
		   FDR_stars = gtools::stars.pval(FDR))











