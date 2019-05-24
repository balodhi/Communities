# test for Yifan's data


library(tidyverse)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(readxl)
library(SDMTools) # for weighted statistics
library(ggpubr)
library(dbscan)
library(here)

source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)


odir = 'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


n2 = read_xlsx('data/n2_1day.xlsx', sheet = 'n2_1day_adult_9_bacteria_1-1')
ep2 = read_xlsx('data/ep21day.xlsx', sheet = 'ep21day')

n2 = n2 %>% filter(TOF > 350, Extinction > 50)
ep2 = ep2 %>% filter(TOF > 350, Extinction > 50)


# plot data to see if it is similar to Yifan's data
ep2 %>% mutate(Bacteria = as.factor(Bacteria)) %>%
	# filter(Bacteria == 'Myb71') %>%
	ggplot(aes(x = TOF, y = Extinction, colour = Bacteria)) +
	geom_point() +
	facet_wrap(~Bacteria)



gcb.n2 = n2 %>% filter(Bacteria == 'GCB') %>% select(-Replicator, -Bacteria)
mg.n2 = n2 %>% filter(Bacteria == 'MG1655') %>% select(-Replicator, -Bacteria)


# lets compute the PCA
res1 = PCA(gcb.n2, scale.unit = FALSE, ncp = 5, graph = F)
res2 = PCA(mg.n2, scale.unit = FALSE, ncp = 5, graph = F)

cosa = get_pca_ind(res1)
cosa = as.data.frame(cosa$coord)

cosa %>% ggplot(aes(x = Dim.1, y = Dim.2)) + geom_point()

quartz.save(file = here('exploration', 'GCB_N2_PC1.pdf'), 
	type = 'pdf', dpi = 300, height = 7, width = 9)

cosa2 = get_pca_ind(res2)
cosa2 = as.data.frame(cosa2$coord)

cosa2 %>% ggplot(aes(x = Dim.1, y = Dim.2)) + geom_point()

quartz.save(file = here('exploration', 'MG1665_N2_PC1.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 9)


d1 = gcb.n2; d1['set'] = 'original'; colnames(d1) = c('var1', 'var2', 'set')
d2 = cosa; d2['set'] = 'PCA'; colnames(d2) = c('var1', 'var2', 'set')

dset = rbind(d1, d2) 

ggplot(dset, aes(x = var1, y = var2, colour = set)) + 
	geom_point() + 
	facet_wrap(~set, scales = "free")

quartz.save(file = here('exploration', 'PC1vsOriginal.pdf'),
	type = 'pdf', dpi = 300, height = 5, width = 9)


# let's do some plots and statistics

# stupid function to extract variabililty explained in the first PC
pca_info = function(data){
	var1 = data$eig[1,2]
	return(var1)
} 


# nest and PCA calculation
n2.pca = n2 %>% group_by(Bacteria) %>% 
	nest %>%
	mutate(pca = map(data, ~ PCA(.x %>% select(-Replicator), 
		scale.unit = FALSE, ncp = 2, graph = F)))


n2.pca %>% 
	mutate(var_exp = map(pca, ~pca_info(.x))) %>%
	select(Bacteria, var_exp) %>%
	unnest %>%
	ggplot(aes(x = Bacteria, y = var_exp, fill = Bacteria)) + 
		geom_bar(stat = "identity", colour = 'black') +
		labs(title = '% of variability explained in PC1') +
		theme_light()


quartz.save(file = here('exploration', 'var_exp_N2.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 7)





# nest and PCA calculation
ep2.pca = ep2 %>% group_by(Bacteria) %>% 
	nest %>%
	mutate(pca = map(data, ~ PCA(.x %>% select(-Replicator), 
		scale.unit = FALSE, ncp = 2, graph = F)))


ep2.pca %>% 
	mutate(var_exp = map(pca, ~pca_info(.x))) %>%
	select(Bacteria, var_exp) %>%
	unnest %>%
	ggplot(aes(x = Bacteria, y = var_exp)) + 
		geom_bar(stat = "identity") +
		labs(title = '% of variability explained in PC1') +
		theme_light()


quartz.save(file = here('exploration', 'var_exp_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 7)


# create a table with everything

total = rbind(n2, ep2)
total['strain'] = 'ep2'
total$strain[1:dim(n2)[1]] = 'n2'



total.pca = total %>% group_by(Bacteria, strain) %>% 
	nest %>%
	mutate(pca = map(data, ~ PCA(.x %>% select(-Replicator), 
		scale.unit = FALSE, ncp = 2, graph = F)))



total.pca %>% 
	mutate(var_exp = map(pca, ~pca_info(.x))) %>%
	select(Bacteria, strain, var_exp) %>%
	unnest %>%
	ggplot(aes(x = Bacteria, y = var_exp, fill = strain, colour = strain)) + 
		geom_bar(stat = "identity", position = position_dodge(), color = "black") +
		labs(title = '% of variability explained in PC1') +
		scale_fill_manual(values=c('#999999','#E69F00')) +
		theme_light() 

quartz.save(file = here('exploration', 'var_exp_comparison.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 7)




# now total, but by replicator (replicate)

total.pca = total %>% group_by(Bacteria, strain, Replicator) %>% 
	nest %>%
	mutate(pca = map(data, ~ PCA(.x , 
		scale.unit = FALSE, ncp = 2, graph = F)))



total.pca %>% 
	mutate(var_exp = map(pca, ~pca_info(.x))) %>%
	select(Bacteria, strain, Replicator, var_exp) %>%
	unnest %>%
	# filter(Bacteria == 'Myb131', strain == 'n2') %>% 
	ggplot(aes(x = Replicator, y = var_exp, fill = Replicator, colour = Replicator)) + 
		geom_bar(stat = "identity", position = position_dodge(), color = "black") +
		labs(title = '% of variability explained in PC1') +
		facet_wrap(Bacteria ~ strain, scales = "free_x") +
		# scale_fill_manual(values=c('#999999','#E69F00')) +
		theme_light() +
		theme(legend.position = "none")

quartz.save(file = here('exploration', 'var_exp_comparison_replicate.pdf'),
	type = 'pdf', dpi = 300, height = 13, width = 15)




################
## statistics ##
################

# let's start with something simple: comparison of three bacteria in the two strains

sub.n2 = n2.pca %>% filter(Bacteria %in% c('GCB', 'op50', 'MG1655'))

rescale = function(data){
	vector = data$ind$coord[,1] + abs(min(data$ind$coord[,1]))
	return(vector)
}

test.n2 = sub.n2 %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest

ggplot(test.n2, aes(x = Bacteria, y = coord, colour = Bacteria)) + geom_boxplot()



model = aov(coord ~ Bacteria, data = test.n2)




# boxplot of sample comparison


n2.coord = n2.pca %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest
n2.coord['strain'] = 'n2'


ep2.coord = ep2.pca %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest
ep2.coord['strain'] = 'ep2'


total.coord = rbind(n2.coord, ep2.coord)


total.coord %>% ggplot(aes(x = strain, y = coord, colour = strain)) +
	geom_boxplot(position = position_dodge()) +
	facet_wrap(~Bacteria, scales = "free")

quartz.save(file = here('exploration', 'boxplot_PC1.pdf'),
	type = 'pdf', dpi = 300, height = 13, width = 15)





total.coord['replicate'] = total$Replicator

## boxplot by replicate

total.coord %>% ggplot(aes(x = replicate, y = coord, fill = strain)) +
	geom_boxplot(position = position_dodge()) +
	facet_wrap(~Bacteria, scales = "free") +
	theme_light() +
	labs(title = 'Replicates Boxplot') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))


quartz.save(file = here('exploration', 'boxplot_PC1_replicates.pdf'),
	type = 'pdf', dpi = 300, height = 10, width = 12)



total.sum = total.coord %>% 
	group_by(Bacteria, strain, replicate) %>%
	summarise(
			Mean = mean(coord), 
			SD = sd(coord),
			Median = median(coord)) %>%
	group_by(Bacteria, strain) %>%
	mutate(w = 1/(SD**2),
		   w_norm = w/(sum(1/(SD**2)))) %>%
	ungroup

# calculate means and sd (weighted and unweighted) of means
w.sum = total.sum %>%
	group_by(Bacteria, strain) %>%
	summarise(aMean = mean(Mean),
			  aSD = sd(Mean),
			  w_Mean = wt.mean(Mean, w),
			  w_SD = wt.sd(Mean, w))


# barplot with points and error

w.sum %>% 
	ggplot(aes(x = Bacteria, y = w_Mean, colour = strain, group = strain)) +
	geom_errorbar(aes(ymin = w_Mean - w_SD, ymax = w_Mean + w_SD), position = position_dodge(0.9), width = 0.1) +
	geom_bar(aes(fill = strain), stat = 'identity', position = position_dodge(), alpha = 0.1) +
	geom_point(data = total.sum, aes(x = Bacteria, y = Mean), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1), alpha = 0.8) +
	scale_color_manual(values = c('#0000FEFF', '#107F01FF')) +
	scale_fill_manual(values = c('#0000FEFF', '#107F01FF')) +
	labs(x = 'Bacteria',
		y = 'PC1 mean value') +
	theme_light()

quartz.save(file = here('exploration', 'barplot_PC1_replicates.pdf'),
	type = 'pdf', dpi = 300, height = 6, width = 9)





######################
# filtering outliers #



cosa = total.coord %>% filter(Bacteria == 'Myb131', strain == 'n2')


# inter quartile range
IQR(cosa$coord) * 1.5 

# testing quartiles and IQR
unname(quantile(cosa$coord, probs=c(.25, .75))[2])
IQR(cosa$coord) * 1.5 + quantile(cosa$coord, probs=c(.25, .75))[2] 

# create a data frame with the limit values with two methods:
# method 1: calculate outilers by Q3 + 1.5*IQR | Q1 - 1.5*IQR
# method 2: calculate mean ± 3 * sigma

limits = total.coord %>% 
	group_by(Bacteria, strain, replicate) %>%
	summarise(
	iqr = IQR(coord),
	uplimit = IQR(coord) * 1.5 + unname(quantile(coord, probs=c(.25, .75))[2]),
	lowlimit = unname(quantile(coord, probs=c(.25, .75))[1]) - IQR(coord) * 1.5,
	upsigma = mean(coord) + (sd(coord) * 3),
	lowsigma = mean(coord) - (sd(coord) * 3)
	)


# join to total data frame

total = total.coord %>% left_join(limits)

### method 1
total_iqr = total %>% 
	filter(coord < uplimit & coord > lowlimit) %>%
	select(Bacteria:replicate)


# let's see how many values have been filtered

temp1 = total %>% 
	group_by(Bacteria, replicate, strain) %>%
	summarise(N = n())

temp2 = total_iqr %>%
	group_by(Bacteria, replicate, strain) %>%
	summarise(N2 = n())

temp1 %>% left_join(temp2) %>%
	mutate(diff = N - N2,
		   prop = 100 - (N2/N * 100)) %>%
	ggplot(aes(x = replicate, y = prop, fill = Bacteria)) +
	geom_bar(stat = "identity", position = position_dodge(), color = 'black') +
	facet_wrap(~ strain, ncol = 1) +
	theme_light() +
	labs(title = '% of outliers for each replicate (IQR method)',
		x = 'Replicates',
		y = '% of outliers') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))

quartz.save(file = here('exploration', 'outliers_iqr.pdf'),
	type = 'pdf', dpi = 300, height = 9, width = 11)




### method 2
total_sig = total %>% 
	filter(coord < upsigma & coord > lowsigma) %>%
	select(Bacteria:replicate)


# let's see how many values have been filtered

temp1 = total %>% 
	group_by(Bacteria, replicate, strain) %>%
	summarise(N = n())

temp2 = total_sig %>%
	group_by(Bacteria, replicate, strain) %>%
	summarise(N2 = n())

temp1 %>% left_join(temp2) %>%
	mutate(diff = N - N2,
		   prop = 100 - (N2/N * 100)) %>%
	ggplot(aes(x = replicate, y = prop, fill = Bacteria)) +
	geom_bar(stat = "identity", position = position_dodge(), color = 'black') +
	facet_wrap(~ strain, ncol = 1) +
	theme_light() +
	labs(title = '% of outliers for each replicate (3 sigma method)',
		x = 'Replicates',
		y = '% of outliers') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))

quartz.save(file = here('exploration', 'outliers_sigma.pdf'),
	type = 'pdf', dpi = 300, height = 9, width = 11)


## after cleaning everything, let's create new boxplots


total_iqr %>% ggplot(aes(x = replicate, y = coord, fill = strain)) +
	geom_boxplot(position = position_dodge(1)) +
	facet_wrap(~Bacteria, scales = "free") +
	theme_light() +
	labs(title = 'Replicates Boxplot (IQR outliers)') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))


quartz.save(file = here('exploration', 'boxplot_PC1_iqr.pdf'),
	type = 'pdf', dpi = 300, height = 10, width = 12)



total_sig %>% ggplot(aes(x = replicate, y = coord, fill = strain)) +
	geom_boxplot(position = position_dodge(1)) +
	facet_wrap(~Bacteria, scales = "free") +
	theme_light() +
	labs(title = 'Replicates Boxplot (3_Sigma outliers)') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))


quartz.save(file = here('exploration', 'boxplot_PC1_sig.pdf'),
	type = 'pdf', dpi = 300, height = 10, width = 12)








##########################################
### Outlier dealing with original data ###
##########################################

# In this piece of code we will be testing several ways of dealing with the 
# outliers in our datasets. Two of them will be sample independent (quartiles
# and 3*sigma methods), whereas the other will depend on the sample itsef (isolation
# forests, local outiler, dbscan... among other possible methods).


# original datasets
# n2, ep2

### Method 1

# create limit values with IQR and 3sigma
n2_lim = n2 %>% 
	group_by(Replicator, Bacteria) %>%
	summarise(tof_iqr = IQR(TOF),
			  tof_uplimit = (IQR(TOF) * 1.5) + unname(quantile(TOF, probs = c(.25, .75))[2]),
			  tof_lowlimit = unname(quantile(TOF, probs = c(.25, .75))[1]) - (IQR(TOF) * 1.5),
			  tof_upsigma = mean(TOF) + (sd(TOF) * 3),
			  tof_lowsigma = mean(TOF) - (sd(TOF) * 3),
			  ext_iqr = IQR(Extinction),
			  ext_uplimit = IQR(Extinction) * 1.5 + unname(quantile(Extinction, probs = c(.25, .75))[2]),
			  ext_lowlimit = unname(quantile(Extinction, probs = c(.25, .75))[1]) - IQR(Extinction) * 1.5,
			  ext_upsigma = mean(Extinction) + (sd(Extinction) * 3),
			  ext_lowsigma = mean(Extinction) - (sd(Extinction) * 3)) %>%
	ungroup


# filter and plot
n2_filt = n2 %>% 
	left_join(n2_lim) %>%
	filter(TOF > tof_lowsigma) %>%
	filter(TOF < tof_upsigma) %>%
	filter(Extinction > ext_lowsigma) %>%
	filter(Extinction < ext_upsigma)

p1 = ggplot(n2_filt, aes(x = TOF, y = Extinction, colour = Bacteria)) +
	geom_point(alpha = 0.7, size = 2) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = ggplot(n2, aes(x = TOF, y = Extinction, colour = Bacteria)) +
	geom_point(alpha = 0.7, size = 2) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p2, p1, 
          labels = c("Original", "Filtered (3 * sigma)"),
          ncol = 2, nrow = 1)


quartz.save(file = here('exploration', 'Sigma_original_comparison.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)



#plot outliers selected in this method

n2_filt['Group'] = 1
test = n2 %>% left_join(n2_filt)

test %>% 
	mutate(Group = ifelse(is.na(Group), 0, 1)) %>%
	mutate(Group = as.factor(Group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = Group), alpha = 0.7, size = 1) + 
	theme_light() +
	facet_wrap(~Bacteria)

quartz.save(file = here('exploration', 'Sigma_outliers.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 10)


## DBSCAN method


df = (n2 %>% filter(Bacteria == 'op50'))[,3:4]
db = dbscan(df, eps = 200, minPts = 15)

fviz_cluster(db, df, stand = FALSE, frame = FALSE, geom = "point")



# function to extract clusters from DBSCAN
dbclust = function(data, eps = 230, minPts = 15) {
	db = dbscan(data, eps = eps, minPts = minPts)
	return(db$cluster)
}


# this loop code will calculate all the 
db.n2 = n2 %>% group_by(Bacteria) %>%
	nest %>%
	mutate(group = map(data, ~dbclust(.x[,2:3]))) %>%
	unnest

db.ep2 = ep2 %>% group_by(Bacteria) %>%
	nest %>%
	mutate(group = map(data, ~dbclust(.x[,2:3]))) %>%
	unnest



# plotting results
p1 = db.n2 %>% 
	# filter(Bacteria == 'op50') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)

p2 = db.ep2 %>% 
	# filter(Bacteria == 'op50') %>%
	mutate(group = as.factor(group)) %>%
	ggplot(aes(x = TOF, y = Extinction)) +
	geom_point(aes(colour = group)) +
	theme_light() +
	facet_wrap(~Bacteria)


ggarrange(p1, p2, 
          labels = c("N2", "ep2"),
          ncol = 2, nrow = 1)


quartz.save(file = here('exploration', 'DBSCAN_n2_ep2.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 17)


# plot number of outliers per sample
db.n2 %>%
	group_by(Bacteria, Replicator, group) %>%
	summarise(N = n()) %>%
	filter(group == 0) %>%
	ggplot(aes(x = Replicator, y = N, fill = Bacteria)) +
	geom_bar(stat = 'identity', color = 'black') +
	theme_light() +
	theme(
	axis.text.x = element_text(angle = 45))





# test again the PCR explanation with filtered data


# nest and PCA calculation
n2.pca.db = db.n2 %>% 
	filter(group != 0) %>%
	select(-group) %>%
	group_by(Bacteria) %>% 
	nest %>%
	mutate(pca = map(data, ~ PCA(.x %>% select(-Replicator), 
		scale.unit = FALSE, ncp = 2, graph = F)))


n2.pca.db %>% 
	mutate(var_exp = map(pca, ~pca_info(.x))) %>%
	select(Bacteria, var_exp) %>%
	unnest %>%
	ggplot(aes(x = Bacteria, y = var_exp, fill = Bacteria)) + 
		geom_bar(stat = "identity", colour = 'black') +
		labs(title = '% of variability explained in PC1') +
		theme_light()


quartz.save(file = here('exploration', 'var_exp_N2_db.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 7)





# let's start with something simple: comparison of three bacteria in the two strains

sub.n2 = n2.pca.db %>% filter(Bacteria %in% c('GCB', 'op50', 'MG1655'))

rescale = function(data){
	vector = data$ind$coord[,1] + abs(min(data$ind$coord[,1]))
	return(vector)
}

test.n2 = sub.n2 %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest

ggplot(test.n2, aes(x = Bacteria, y = coord, fill = Bacteria)) + 
	geom_boxplot(alpha = 0.9) +
	geom_jitter(position = position_jitter(0.3), alpha = 0.3) +
	theme_light()



model = aov(coord ~ Bacteria, data = test.n2)
TukeyHSD(model)

quartz.save(file = here('exploration', 'boxplot_test_db.pdf'),
	type = 'pdf', dpi = 300, height = 7, width = 7)




db.n2 %>% filter(Bacteria %in% c('GCB', 'MG1655')) %>%
	ggplot(aes(x = TOF, y = Extinction, colour = Bacteria)) +
	geom_point(size = 1, alpha = 0.5) +
	theme_light()




summary(res.man)



thing = db.n2 %>% filter(Bacteria %in% c('GCB', 'MG1655')) %>% mutate(Bacteria = as.factor(Bacteria))
res.man <- manova(cbind(TOF, Extinction) ~ Bacteria, data = thing)













# boxplot of sample comparison


n2.coord = n2.pca %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest
n2.coord['strain'] = 'n2'


ep2.coord = ep2.pca %>% mutate(coord = map(pca, ~rescale(.x))) %>% 
	select(Bacteria, coord) %>% unnest
ep2.coord['strain'] = 'ep2'


total.coord = rbind(n2.coord, ep2.coord)


total.coord %>% ggplot(aes(x = strain, y = coord, colour = strain)) +
	geom_boxplot(position = position_dodge()) +
	facet_wrap(~Bacteria, scales = "free")

quartz.save(file = here('exploration', 'boxplot_PC1.pdf'),
	type = 'pdf', dpi = 300, height = 13, width = 15)





total.coord['replicate'] = total$Replicator

## boxplot by replicate

total.coord %>% ggplot(aes(x = replicate, y = coord, fill = strain)) +
	geom_boxplot(position = position_dodge()) +
	facet_wrap(~Bacteria, scales = "free") +
	theme_light() +
	labs(title = 'Replicates Boxplot') +
	theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.text.x = element_text(angle = 45))


quartz.save(file = here('exploration', 'boxplot_PC1_replicates.pdf'),
	type = 'pdf', dpi = 300, height = 10, width = 12)



total.sum = total.coord %>% 
	group_by(Bacteria, strain, replicate) %>%
	summarise(
			Mean = mean(coord), 
			SD = sd(coord),
			Median = median(coord)) %>%
	group_by(Bacteria, strain) %>%
	mutate(w = 1/(SD**2),
		   w_norm = w/(sum(1/(SD**2)))) %>%
	ungroup

# calculate means and sd (weighted and unweighted) of means
w.sum = total.sum %>%
	group_by(Bacteria, strain) %>%
	summarise(aMean = mean(Mean),
			  aSD = sd(Mean),
			  w_Mean = wt.mean(Mean, w),
			  w_SD = wt.sd(Mean, w))


# barplot with points and error

w.sum %>% 
	ggplot(aes(x = Bacteria, y = w_Mean, colour = strain, group = strain)) +
	geom_errorbar(aes(ymin = w_Mean - w_SD, ymax = w_Mean + w_SD), position = position_dodge(0.9), width = 0.1) +
	geom_bar(aes(fill = strain), stat = 'identity', position = position_dodge(), alpha = 0.1) +
	geom_point(data = total.sum, aes(x = Bacteria, y = Mean), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1), alpha = 0.8) +
	scale_color_manual(values = c('#0000FEFF', '#107F01FF')) +
	scale_fill_manual(values = c('#0000FEFF', '#107F01FF')) +
	labs(x = 'Bacteria',
		y = 'PC1 mean value') +
	theme_light()

quartz.save(file = here('exploration', 'barplot_PC1_replicates.pdf'),
	type = 'pdf', dpi = 300, height = 6, width = 9)
















