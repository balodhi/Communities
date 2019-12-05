###
library(here)
library(tidyverse)
library(ggrepel)
library(openxlsx)
library(readxl)
library(RVenn)

options(width = 220)

# load reactions and compounds from KBbase to translate Sphin model to VMH names

KB_reacts =     read_xlsx('KBase2VMH_reactions.xlsx', sheet = 'RxnMetTranslation') 
KB_comp =     read_xlsx('KBase2VMH_reactions.xlsx', sheet = 'Compounds') 


Acro =     read_xlsx('Acro_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Acro')
GCB =       read_xlsx('GCB_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'GCB')
Arthro = read_xlsx('Arthro_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Arthro')
Marb =     read_xlsx('Marb_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Marb')
Com =       read_xlsx('Com_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Com')
EC =     	 read_xlsx('EC_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'EC')
MG =     	 read_xlsx('MG_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'MG')
EcoB =     read_xlsx('EcoB_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'EcoB')
Micro =   read_xlsx('Micro_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Micro')
Ochro =   read_xlsx('Ochro_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Ochro')
Rod =       read_xlsx('Rod_reacts.xlsx', sheet = 'React_abrev') %>% mutate(Strain = 'Rod')
Sphin =   read_xlsx('Sphin_reacts.xlsx', sheet = 'Reaction List') %>% mutate(Strain = 'Sphin')




#### fix Sphin reaction names

Sphin.rxn = Sphin$Abbreviation

# some string transformations to match values
# CHECK IF YOU CHANGE THE METABOLIC MODEL OF SPHINGOBACTERIUM
Sphin.rxn = str_replace(Sphin.rxn, '_e0', '(e)')
Sphin.rxn = str_replace(Sphin.rxn, 'R_', '')
Sphin.rxn = str_replace(Sphin.rxn, '_c0', '')
Sphin.rxn = str_replace(Sphin.rxn, 'bio1', 'biomass')
Sphin.rxn = str_replace(Sphin.rxn, 'EX_cpd15302', 'EX_cpd15302(c)')
Sphin.rxn = str_replace(Sphin.rxn, 'EX_cpd02701', 'EX_cpd02701(c)')
Sphin.rxn = str_replace(Sphin.rxn, 'EX_cpd11416', 'EX_cpd11416(c)')
Sphin.rxn = str_replace(Sphin.rxn, 'rxn03487', 'rxn00222') # manual curation 
Sphin.rxn = str_replace(Sphin.rxn, 'rxn02009', 'rxn02008') # manual curation 


Sphin.rxn = as.character(Sphin.rxn)


# putting back the vector with changed names, and adding names from the list

Sphin$Abbreviation = Sphin.rxn

Sphin = Sphin %>% 
	left_join(KB_reacts, by = c('Abbreviation' = 'DraftReactionID')) %>% 
	select(-Subsystem, -Notes) %>% 
	mutate(VMHreactionID = ifelse(is.na(VMHreactionID), Abbreviation, VMHreactionID)) # replaces empty entries by the original name

Sphin = Sphin %>% select(-Abbreviation) %>%
	rename(Abbreviation = VMHreactionID) %>%
	select(Abbreviation, everything())


# number of reactions by strain

n_react = c(dim(Acro)[1], 
			dim(GCB)[1], 
			dim(Arthro)[1], 
			dim(Marb)[1], 
			dim(Com)[1], 
			dim(EC)[1], 
			dim(MG)[1], 
			dim(EcoB)[1], 
			dim(Micro)[1], 
			dim(Ochro)[1], 
			dim(Rod)[1], 
			dim(Sphin)[1])

strains = c('Achromobacter \n xylosoxidans', 
			'Acinetobacter \n lwoffii',
			'Arthrobacter \n castelli',
			'Bacillus \n subtilis',
			'Comamonas \n testeroni',
			'Enterobacter \n cloacae',
			'Escherichia \n coli MG1655',
			'Escherichia \n coli SE11',
			'Microbacterium \n paraoxydans',
			'Ochrobactrum \n anthropi',
			'Rhodococcus \n erythropolis',
			'Sphingobacterium \n faecium')

strains_short = c(
	'Achromobacter',
	'Acinetobacter',
	'Arthrobacter',
	'Bacillus',
	'Comamonas',
	'Enterobacter',
	'Escherichia K12',
	'Escherichia B',
	'Microbacterium',
	'Ochrobactrum',
	'Rhodococcus',
	'Sphingobacterium')

df1 = data.frame(strains, n_react)
names(df1) = c('Species', 'N_reactions')

# reactions per strain
df1 %>% arrange(N_reactions) %>% 
	ggplot(aes(x = reorder(Species, - N_reactions), y = N_reactions)) +
	geom_bar(stat = "identity", fill = 'grey30', colour = 'black') +
	ylim(0, 2800) +
	labs(x = 'Species', y = 'Number of reactions') +
	theme_classic() 

quartz.save(file = here('summary', 'React_per_strain.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)





# store all reactions per strain in different vectors
Acro_react = unique(Acro$Abbreviation)
GCB_react = unique(GCB$Abbreviation)
Arthro_react = unique(Arthro$Abbreviation)
Marb_react = unique(Marb$Abbreviation)
Com_react = unique(Com$Abbreviation)
EC_react = unique(EC$Abbreviation)
MG_react = unique(MG$Abbreviation)
EcoB_react = unique(EcoB$Abbreviation)
Micro_react = unique(Micro$Abbreviation)
Ochro_react = unique(Ochro$Abbreviation)
Rod_react = unique(Rod$Abbreviation)
Sphin_react = unique(Sphin$Abbreviation)


list_react = list(
	Acro = Acro_react,
	GCB = GCB_react,
	Arthro = Arthro_react,
	Marb = Marb_react,
	Com = Com_react,
	EC = EC_react,
	MG = MG_react,
	EcoB = EcoB_react,
	Micro = Micro_react,
	Ochro = Ochro_react,
	Rod = Rod_react,
	Sphin = Sphin_react)


# union of all sets
mod = Venn(list_react)
# number of total and unique reactions of the 12 strains
total.rxns = sort(unite(mod))


# function that calculates the number of total and different reactions in a list of organisms
n.react = function(list){
	mod = Venn(list)
	total.rxns = unite(mod)

	## curate the final list
	# change _e for (e) and take unique list
	# remove every entry with '__'
	# remove every entry with '_copy1', '_copy2'...
	# remove every entry with 'biomass'
	
	total.rxns = str_replace(total.rxns, '_e', '(e)')
	total.rxns = total.rxns[!str_detect(total.rxns, 'biomass')]
	total.rxns = total.rxns[!str_detect(total.rxns, '__')]
	total.rxns = total.rxns[!str_detect(total.rxns, '_copy')]
	total.rxns = total.rxns[!str_detect(total.rxns, 'DM_atp_c_')]
	
	n = length(total.rxns)
	return(n)
}

# test
n.react(list_react)

# how to slice a list to get a list
cosa = c(list_react[1], list_react[2])



# x = 1:12
# comb = combn(x, 2)

# comb[,1]


##
# this function calculates all the possible combinations of n elements in a list, and computes 
# the average and std of the unique number of reactions for each combination
# it gives back a tibble 

comb.react = function(list, n){
	x = 1:12
	comb = combn(x, n)
	void = c()
	for (i in 1:dim(comb)[2]){
		j = n.react(list[comb[,i]])
		void = c(void, j)
	}
	Mean = mean(void)
	SD = sd(void)
	df = tibble(n, Mean, SD)
	return(df)
}
# test
comb.react(list_react, 3)


# create an empty tibble as a primer
reactxspc = tibble(n = 0, Mean = 0, SD = 0)

for (i in 2:12){
	tib = comb.react(list_react, i)
	reactxspc = bind_rows(reactxspc, tib)
}


# calculate the mean and sd of only one organism
a = c()
for (i in 1:12){
	n = length(list_react[[i]])
	a = c(a, n)
}
m1 = mean(a)
sd1 = sd(a)

uni.reacts = tibble(n = 1, Mean = m1, SD = sd1)

reactxspc = reactxspc %>%
	bind_rows(uni.reacts) %>%
	filter(n >= 1) %>%
	mutate(SD = ifelse(is.na(SD), 0, SD))


# plot results
reactxspc %>%
	ggplot(aes(x = n, y = Mean)) +
	geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2) +
	geom_line() +
	scale_x_continuous(breaks = seq(1, 12, by = 1), expand = c(0, 0)) +
	scale_y_continuous(limits = c(1000, 3400), expand = c(0,0)) +
	ylab("Number of \n unique reactions") +
	xlab("Community size") +
	# ylim(1000, 3400) +
	theme_bw()


quartz.save(file = here('summary', 'total_reactions.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 9)




### Mock community for Filipe

# Acinetobacter
# Achromobacter
# Comamonas
# Enterobacter cloacae
# Escherichia coli
# Rhodococcus
# Sphingobacterium





F_list_react = list(
	Acro = Acro_react,
	GCB = GCB_react,
	# Com = Com_react,
	EC = EC_react,
	MG = MG_react,
	Rod = Rod_react,
	Sphin = Sphin_react
	)

# union of all sets
mod = Venn(F_list_react)
# number of total and unique reactions of the 12 strains
total.rxns = sort(unite(mod))

# clean list
total.rxns = str_replace(total.rxns, '_e', '(e)')
total.rxns = total.rxns[!str_detect(total.rxns, 'biomass')]
total.rxns = total.rxns[!str_detect(total.rxns, '__')]
total.rxns = total.rxns[!str_detect(total.rxns, '_copy')]
total.rxns = total.rxns[!str_detect(total.rxns, 'DM_atp_c_')]



pos = length(total.rxns)

# plot results
reactxspc %>%
	ggplot(aes(x = n, y = Mean)) +
	geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2) +
	geom_line() +
	scale_x_continuous(breaks = seq(1, 12, by = 1), expand = c(0, 0)) +
	scale_y_continuous(limits = c(1000, 3400), expand = c(0,0)) +
	geom_segment(aes(x = 6, xend = 6, y = 1000, yend = pos), colour = 'grey50', alpha = .2) +
	# geom_segment(aes(x = 1, xend = 12, y = 3343 * 0.95, yend = 3343 * 0.95), colour = 'grey50', alpha = .2) +
	geom_hline(yintercept = 3343 * 0.95, linetype = 'longdash', colour = 'grey20', alpha = .5) +
	geom_hline(yintercept = 3343 * 0.90, linetype = 'longdash', colour = 'grey20', alpha = .5) +
	geom_point(aes(x = 6, y = pos), shape = 21, size = 2, fill = 'blue') +
	# geom_label_repel(x = 6,§ y = pos, label = pos) +
	annotate('label', x = 6.4, y = pos - 50, label = pos) +
	annotate('text', x = 1.5, y = (3343 * 0.95) + 25, label = '95%', colour = 'grey20') +
	annotate('text', x = 1.5, y = (3343 * 0.90) + 25, label = '90%', colour = 'grey20') +
	ylab("Number of \n unique reactions") +
	xlab("Community size") +
	# ylim(1000, 3400) +
	theme_bw()


quartz.save(file = here('summary', 'total_reactions_subcommunity.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 9)





# saving tables in an excel file
list_of_datasets = list('Reactions per strain' = df1, 'Summary stats' = reactxspc)

write.xlsx(list_of_datasets, here('summary', 'stats.xlsx'), colNames = T, rowNames = F) 






#################################################
#################################################































































