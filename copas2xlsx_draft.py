# -*- coding: utf-8 -*-
''' 
This python code will process COPAS output and transform it into an easier 
dataset to be analysed by other means 

Daniel Martinez, Jun - 2019
'''




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string
import math
import openpyxl
from optparse import OptionParser


# Release information
__version__ = '0.1'
_scriptname = 'copas2xlsx'
_verdata = 'Jun 2019'
_devflag = True



# Option parser

parser = OptionParser()
parser.add_option("-i", "--input", 
					dest = "inputfile", 
					metavar = "INPUT",
                  	help = "input file")
parser.add_option("-d", "--design", 
					dest = "designfile", 
					metavar = "INPUT",
                  	help = "input file")
parser.add_option("-o", "--output", 
					dest = "outputfile",  
					metavar = "OUTPUT", 
                  	help = "output file")

(options, args) = parser.parse_args()


# Functions

# creates Bacteria column in dataframe
def bact(df):
	v = []
	for i in range(df.shape[0]):
		v.append(df['Sample'][i].split('_')[0])
	return(v)

# outputs the number of rows for the scatter plot
def plot_nrow(df, var):
	elm = len(set(df[str(var)]))
	ncols = math.ceil(elm/3)
	return(ncols)






# Program Header
print('\n====================================================\n')
print(_scriptname + 'script, v' + __version__ , _verdata + 
	'\n =-= by Daniel Martinez =-=')
if(_devflag):
    print('\nWARNING! THIS IS JUST A DEVELOPMENT SUBRELEASE.' +
          '\nUSE IT AT YOUR OWN RISK!\n')
print('\n====================================================\n')


### Design

# reads excel file
ds = pd.read_excel('Design.xlsx', sheet_name = 0, index_col = 0)

# create a dictionary of values
ind = list(string.ascii_uppercase)[0:8]

# this loop associates each entry to an index, and cleans
# all the NaN values
wells = dict()
for index in ind:
	for n in range(1,13):
		if pd.isnull(ds.loc[str(index), n]) == False:
			wells[(str(index) + str(n))] = ds.loc[str(index), n]


### COPAS data

# read input file
copas = pd.read_csv('n2_1day_adult_9bacteria_1.txt', sep = "\t", header = 0)
copas2 = pd.read_csv('ep2_1day_adult_9bacteria_1_1.txt', sep = "\t", header = 0)
# remove rubish rows at the end of COPAS file
# I should do it in a more robust way
limit = copas.index[copas['Id'] == 'INSTRUMENT SETTINGS  HAVE BEEN MODIFIED DURING OR AFTER DATA ACQUISITION!'].tolist()
copas = copas.iloc[:limit[0], :]

# select columns to work
copas = copas[['Id', 'Source well', 'TOF', 'Extinction']]

#####
## copas 2
limit = copas2.index[copas2['Id'] == 'INSTRUMENT SETTINGS  HAVE BEEN MODIFIED DURING OR AFTER DATA ACQUISITION!'].tolist()
copas2 = copas2.iloc[:limit[0], :]

# select columns to work
copas2 = copas2[['Id', 'Source well', 'TOF', 'Extinction']]


# create a df from dict
wells_df = pd.DataFrame.from_dict(list(wells.items()))
wells_df.columns = ['Source well', 'Sample']
wells_df['Bacteria'] = bact(wells_df) # creates an extra column with bacterial names


# merge dfs
copas = pd.merge(copas, wells_df, on = 'Source well', how = 'left')
copas = copas[(copas['Extinction'] > 50) & (copas['TOF'] > 350)]

# merge dfs
copas2 = pd.merge(copas2, wells_df, on = 'Source well', how = 'left')
copas2 = copas2[(copas2['Extinction'] > 50) & (copas2['TOF'] > 350)]


copas.to_excel(str('copas') + '_filt.xlsx', sheet_name = 'Summary', index = False)



###############
## test zone ##
############### 




###
# scatter plots
# test


op50 = copas[copas['Bacteria'] == 'OP50']
Myb181 = copas[copas['Bacteria'] == 'Myb181']
Myb71 = copas[copas['Bacteria'] == 'Myb71']


elm = len(set(copas['Bacteria']))

fig, ax = plt.subplots(1, 3, sharex='col', sharey='row')

ax[0].scatter(op50['TOF'],	   op50['Extinction'], c = 'b', alpha = 0.7)
ax[1].scatter(Myb181['TOF'], Myb181['Extinction'], c = 'r', alpha = 0.7)
ax[2].scatter(Myb71['TOF'],   Myb71['Extinction'], c = 'k', alpha = 0.7)

plt.show()



# long version, but not relying on creating new variables
fig, ax = plt.subplots(2, 3, sharex = 'col', sharey = 'row')

ax[0,0].scatter(copas[copas['Bacteria'] == 'OP50']['TOF'], copas[copas['Bacteria'] == 'OP50']['Extinction'], c = 'b', alpha = 0.7)
ax[0,0].set_title('OP50')
ax[0,1].scatter(copas[copas['Bacteria'] == 'Myb181']['TOF'], copas[copas['Bacteria'] == 'Myb181']['Extinction'], c = 'r', alpha = 0.7)
ax[0,1].set_title('Myb181')
ax[0,2].scatter(copas[copas['Bacteria'] == 'Myb71']['TOF'],   copas[copas['Bacteria'] == 'Myb71']['Extinction'], c = 'k', alpha = 0.7)
ax[0,2].set_title('Myb71')
fig.suptitle('COPAS data', fontsize = 16)

plt.show()



nrows = plot_nrow(copas, 'Bacteria')

fig, ax = plt.subplots(nrows, 3, sharex = 'col', sharey = 'row')
	for name in set(copas[str('Bacteria')]):
		ax[0].scatter(copas[copas['Bacteria'] == name]['TOF'], copas[copas['Bacteria'] == name]['Extinction'], c = 'b', alpha = 0.7)






















