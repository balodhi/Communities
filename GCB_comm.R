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
n2_3 = read_xlsx(here('raw_data' ,'/062519/copas/n2_filt.xlsx'), sheet = 'Summary')
n2_4 = read_xlsx(here('raw_data' ,'/070219/copas/n2_filt.xlsx'), sheet = 'Summary')

# ep2

ep2_1 = read_xlsx(here('raw_data', '/061219/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_2 = read_xlsx(here('raw_data', '/061919/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_3 = read_xlsx(here('raw_data', '/062519/copas/ep2_filt.xlsx'), sheet = 'Summary')
ep2_4 = read_xlsx(here('raw_data', '/070219/copas/ep2_filt.xlsx'), sheet = 'Summary')

n2_3 = n2_3 %>% filter(!Extinction > 750) %>%
 	filter(!Extinction > 400 & TOF <1800) %>%
 	filter(Bacteria %in% c('MG', 'Myb9', 'Myb131'))

ep2_3 = ep2_3 %>%
	filter(Bacteria %in% c('MG', 'Myb9', 'Myb131'))


# sample n2_3 is a bit dirty

# n2_3 %>% filter(!Extinction > 750) %>%
# 	filter(!Extinction > 400 & TOF <1800) %>%
# 	ggplot(aes(x = TOF, y = Extinction, colour = Bacteria)) +
# 	geom_point() +
# 	facet_wrap(~Bacteria)

# n2_3 %>% filter(Bacteria == 'Mix') %>% 
# 	ggplot(aes(x = TOF, y = Extinction)) +
# 	geom_point()




# join them together as a long table

n2 = rbind(n2_1, n2_2, n2_3, n2_4)
ep2 = rbind(ep2_1, ep2_2, ep2_3, ep2_4)

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
n2[n2$Bacteria %in% c('GEOP50', 'GEMarb', 'GEM131', 'GEM71'),]$Group = 'Triplet'
ep2[ep2$Bacteria %in% c('GEOP50', 'GEMarb', 'GEM131', 'GEM71'),]$Group = 'Triplet'


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
	# xlim(0, 3200) +
	# ylim(0, 530) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Single') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	# xlim(0, 3200) +
	# ylim(0, 530) +
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
	# xlim(0, 3200) +
	# ylim(0, 530) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'GCB') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	# xlim(0, 3200) +
	# ylim(0, 530) +
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
	# xlim(0, 3200) +
	# ylim(0, 530) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Mix') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	# xlim(0, 3200) +
	# ylim(0, 530) +
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
	# xlim(0, 3200) +
	# ylim(0, 530) +
	scale_fill_viridis_c() +
	theme_light()

p2 = db.ep2 %>% 
	filter(group != 0, Group == 'Triplet') %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
	facet_wrap(~Bacteria) +
	scale_fill_viridis_c() +
	# xlim(0, 3200) +
	# ylim(0, 530) +
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







































# ####
# #### For presentation
# ####


# p1 = db.n2 %>% 
# 	filter(group != 0, Bacteria %in% c('OP50', 'GCB') ) %>%
# 	ggplot(aes(x = TOF, y = Extinction)) +
# 	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
# 	facet_wrap(~Bacteria) +
# 	xlim(0, 3500) +
# 	ylim(0, 500) +
# 	scale_fill_viridis_c() +
# 	theme_light()

# p2 = db.ep2 %>% 
# 	filter(group != 0, Bacteria %in% c('OP50', 'GCB'),
# 		TOF > 300, Extinction > 50 ) %>%
# 	ggplot(aes(x = TOF, y = Extinction)) +
# 	stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
# 	facet_wrap(~Bacteria) +
# 	scale_fill_viridis_c() +
# 	xlim(0, 3500) +
# 	ylim(0, 400) +
# 	theme_light()


# ggarrange(p1, p2, 
#           labels = c("N2", "ep2"),
#           ncol = 1, nrow = 2)

# quartz.save(file = here('summary', 'presentation_plot.pdf'),
# 	type = 'pdf', dpi = 300, height = 8, width = 8)











