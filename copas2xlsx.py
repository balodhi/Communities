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
ds = pd.read_excel(options.designfile, sheet_name = 0, index_col = 0)

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
copas = pd.read_csv(options.inputfile, sep = "\t", header = 0)

# remove rubish rows at the end of COPAS file

limit = copas.index[copas['Id'] == 'INSTRUMENT SETTINGS  HAVE BEEN MODIFIED DURING OR AFTER DATA ACQUISITION!'].tolist()
copas = copas.iloc[:limit[0], :]

# select columns to work
copas = copas[['Id', 'Source well', 'TOF', 'Extinction']]


# create a df from dict
wells_df = pd.DataFrame.from_dict(list(wells.items()))
wells_df.columns = ['Source well', 'Sample']
wells_df['Bacteria'] = bact(wells_df) # creates an extra column with bacterial names


# merge dfs
copas = pd.merge(copas, wells_df, on = 'Source well', how = 'left')
# copas = copas[(copas['Extinction'] > 50) & (copas['TOF'] > 350)]


# save dataframe into an excel file
copas.to_excel(options.outputfile, sheet_name = 'Summary', index = False)



















