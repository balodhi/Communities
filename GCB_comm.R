# community analysis

# libraries

library(tidyverse)
library(colorspace)
library(readxl)
library(SDMTools) # for weighted statistics
library(ggpubr)
library(dbscan)
library(ggnewscale)
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



# function to extract clusters from DBSCAN
dbclust = function(data, eps = 200, minPts = 15) {
	db = dbscan(data, eps = eps, minPts = minPts)
	return(db$cluster)
}


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
# lets make variables for the technical and biological replicates
# stupid loop, do it once

copas['biorep'] = 0
copas['techrep'] = 0
for (i in 1:dim(copas)[1]){
	copas$biorep[i] = as.integer(tail(strsplit(copas$Sample[i], '_')[[1]],2)[1])
	copas$techrep[i] = as.integer(tail(strsplit(copas$Sample[i], '_')[[1]],2)[2])
}

# copas = copas %>% filter(Sample != 'GCB_2_2')

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
		   wExt = 1/(Ext_sd**2),
		   wTOF_norm = wTOF/(sum(1/(TOF_sd**2))),
		   wExt_norm = wExt/(sum(1/(Ext_sd**2))))


# calculate means and sd (weighted and unweighted) of means

copas.sum2 = copas.sum %>%
	group_by(Worm, Group, Bacteria) %>%
	summarise(wTOF_mean = wt.mean(TOF_mean, wTOF),
			  wTOF_SD = wt.sd(TOF_mean, wTOF),
			  wExt_mean = wt.mean(Ext_mean, wExt),
			  wExt_SD = wt.sd(Ext_mean, wExt))




# barplot with points and error

gr = 'Triplet'
p1 = copas.sum2 %>% 
	filter(Group == gr) %>%
	ggplot(aes(x = Bacteria, y = wTOF_mean, colour = Worm, group = Worm)) +
	geom_errorbar(aes(ymin = wTOF_mean - wTOF_SD, ymax = wTOF_mean + wTOF_SD), position = position_dodge(0.9), width = 0.1) +
	geom_bar(aes(fill = Worm), stat = 'identity', position = position_dodge(), alpha = 0.1) +
	geom_point(data = copas.sum %>% filter(Group == gr), 
		aes(x = Bacteria, y = TOF_mean), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.3), alpha = 0.8) +
	scale_color_manual(values = c('#0000FEFF', '#107F01FF')) +
	scale_fill_manual(values = c('#0000FEFF', '#107F01FF')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
 	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) + # removes the spaces at the bottom of the barplot
	theme_classic()

p2 = copas.sum2 %>% 
	filter(Group == gr) %>%
	ggplot(aes(x = Bacteria, y = wExt_mean, colour = Worm, group = Worm)) +
	geom_errorbar(aes(ymin = wExt_mean - wExt_SD, ymax = wExt_mean + wExt_SD), position = position_dodge(0.9), width = 0.1) +
	geom_bar(aes(fill = Worm), stat = 'identity', position = position_dodge(), alpha = 0.1) +
	geom_point(data = copas.sum %>% filter(Group == gr), 
		aes(x = Bacteria, y = Ext_mean), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.3), alpha = 0.8) +
	scale_color_manual(values = c('#0000FEFF', '#107F01FF')) +
	scale_fill_manual(values = c('#0000FEFF', '#107F01FF')) +
	labs(x = 'Bacteria',
		y = 'TOF') +
	scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
	theme_classic()

ggarrange(p1, p2,
          labels = c("TOF", "Extinction"),
          ncol = 1, nrow = 2)


quartz.save(file = here('exploration', 'barplot_means_Triplet.pdf'),
	type = 'pdf', dpi = 300, height = 10, width = 9)







































