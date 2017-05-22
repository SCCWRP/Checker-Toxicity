import sys
import pandas as pd
import numpy as np
from scipy import stats
import xlrd
#from ordereddict import OrderedDict
import collections
import math
stdout = sys.stdout
reload(sys)
sys.setdefaultencoding('utf-8')
sys.stdout = stdout

## COMMON FUNCTIONS
def dcAddErrorToList(error_column, row, error_to_add,df):
	df.ix[int(row), 'row'] = str(row)
	if error_column in df.columns:
		# check if cell value is empty (nan) 
		if(pd.isnull(df.ix[int(row), error_column])):
			# no data exists in cell so add error
	      		df.ix[int(row), error_column] = error_to_add
			print("Row: %s, Error To Add: %s" % (int(row),error_to_add))
		else:
			# a previous error was recorded so append to it
			# even though there may be data lets check to make sure it is not empty
			if str(df.ix[int(row), error_column]):
				#print("There is already a previous error recorded: %s" % str(df.ix[int(row), error_column]))
				df.ix[int(row), error_column] = str(df.ix[int(row), error_column]) + "," + error_to_add
				print("Row: %s, Error To Add: %s" % (int(row),error_to_add))
			else:
				#print("No error is recorded: %s" % str(df.ix[int(row), error_column]))
	      			df.ix[int(row), error_column] = error_to_add
				print("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	else:
		df.ix[int(row), error_column] = error_to_add
		print("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	return df

### WORKSPACE START ###
# place the bight13 toxicity data in a location that the application can access
#df = pd.ExcelFile('/Users/pauls/Documents/Projects/Bight18/Training/clean.xlsx')
df = pd.ExcelFile('./summary.xlsx')

# get sheet names
df_tab_names = df.sheet_names

# create dictionary to hold dataframes
all_dataframes = collections.OrderedDict()

# loop through each sheet
count = 0
for tab in df_tab_names:
	tab_name = tab
    	### extract individual dataframes
    	tab = df.parse(tab)
        # if the sheet is blank skip to the next sheet
	if tab.empty:
		print('The application is skipping sheet "%s" because it is empty' % tab)
		continue
	# lowercase all column names
	tab.columns = [x.lower() for x in tab.columns]
    	### and put into dictionary object
    	all_dataframes[count] = tab
	### create tmp_row for tracking row numbers
	all_dataframes[count]['tmp_row'] = all_dataframes[count].index
	count = count + 1

### WORKSPACE END ###

### SUMMARY TABLE START ###

# summary must not be a groupby otherwise below functions wont work
# all_dataframes[1] is the toxicity results data
batch = all_dataframes[0]
result = all_dataframes[1]
wq = all_dataframes[2]
summary = all_dataframes[1]

def getCalculatedValues(grp):                                                                  
	grp['mean'] = grp['result'].mean()
	grp['n'] = grp['fieldreplicate'].sum()
	grp['stddev'] = grp['result'].std()
	grp['variance'] = grp['stddev'].apply(lambda x: x ** 2 )
	grp['coefficientvariance'] = ((grp['stddev']/grp['mean']) * 100)
	return grp
summary = summary.groupby(['stationid','toxbatch','fieldreplicate']).apply(getCalculatedValues)

# get all control records
cneg = summary[['stationid','toxbatch','sampletypecode','mean']].where(summary['sampletypecode'] == 'CNEG')
# get all non control records
nocneg = summary[['stationid','toxbatch','sampletypecode','mean']].where(summary['sampletypecode'] != 'CNEG')
# get all reference toxicant records just save them for now
reference_toxicants = summary.loc[summary['matrix'].isin(['reference toxicant','rt'])]
# drop all reference toxicants from the summary dataframe - not a part of summary results
summary = summary.loc[~summary['matrix'].isin(['reference toxicant','rt'])]

cneg = cneg.dropna()
nocneg = nocneg.dropna()

cneg['unique'] = np.nan
nocneg['unique'] = np.nan

control_mean = cneg.groupby(['stationid','toxbatch','mean', 'sampletypecode'])['unique'].nunique().reset_index()
result_mean = nocneg.groupby(['stationid','toxbatch','mean', 'sampletypecode'])['unique'].nunique().reset_index()

## prep code control_mean_stats_dict used in getPctControl function
#cneg_stats = summary[['stationid','toxbatch','sampletypecode','n','stddev','mean','variance']].where(summary['sampletypecode'] == 'CNEG')
#cneg_stats = cneg_stats.dropna()
#cneg_stats['unique'] = np.nan
#control_mean_stats = cneg_stats.groupby(['stationid','toxbatch','n','stddev','mean','variance'])['unique'].nunique().reset_index()
## create a dictionary lookup of toxbatch keys and corresponding important values
# drop unique column we used earlier
#control_mean_stats.drop('unique', axis=1, inplace=True)
# make toxbatch the index - we already group so it is unique
#control_mean_dict = control_mean.set_index('toxbatch')['mean'].to_dict()
#control_mean_stats.set_index("toxbatch", drop=True, inplace=True)
#control_mean_stats_dict = control_mean_stats.to_dict(orient="index")
## prep code control_mean_stats_dict end
# THE CODE ABOVE SEEMS TO BE UNNECESSARY - THE LINE BELOW CAN DO THE SAME THING - NOTE ADJUSTED LINE BELOW ALSO
## create a dictionary lookup of toxbatch keys and corresponding control mean values
control_mean_dict = control_mean.set_index('toxbatch')['mean'].to_dict()

def getPctControl(row):
    	## toxbatch control should always be 100
    	if(row['sampletypecode'] == 'CNEG'):
        	row['pctcontrol'] = 100
       	else:
            	if row['toxbatch'] in control_mean_dict:
                	# if the toxbatch is in the lookup dictionary then
                	# divide the result mean from the control mean and times by 100
                	# OLD LINE row['pctcontrol'] = ((row['mean']/control_mean_stats_dict[row['toxbatch']]['mean']) * 100)
			row['pctcontrol'] = ((row['mean']/control_mean_dict[row['toxbatch']]) * 100)
        return row
summary = summary.apply(getPctControl, axis=1)

## author - Tyler Vu
def getPValue(summary):
	for index, values in summary['toxbatch'].iteritems():
		station_code = summary.ix[index, 'stationid']
		cneg_result = summary[['result']].where((summary['sampletypecode'] == 'CNEG') & (summary['toxbatch'] == values))
		result_both = summary[['result']].where((summary['toxbatch'] == values) & (summary['stationid'] == station_code) )
		cneg_result = cneg_result.dropna()
		result_both = result_both.dropna()
		t, p = stats.ttest_ind(cneg_result, result_both, equal_var = False)
		summary.ix[index, 'tstat'] = t
		summary.ix[index, 'pvalue'] = p/2 #we divide by 2 to make it a 1 tailed
		if (t < 0):
			summary.ix[index, 'significance'] = 'NSC'
		else:
			if (p <= .05):
				summary.ix[index, 'significance'] = 'SC'
			else:
				if (summary.ix[index, 'sampletypecode'] == 'CNEG'):
					summary.ix[index, 'significance'] = 'X'
				else:
					summary.ix[index, 'significance'] = 'NSC'
getPValue(summary)

## author - Tyler Vu 
def getSQO(grp):
    if(grp['species'] == 'Eohaustorius estuarius'):
        if(grp['mean'] < 90):
            if (grp['pctcontrol'] < 82):
                if (grp['pctcontrol'] < 59):
                    grp['sqo'] = 'High Toxicity'
                else:
                    if (grp['significance'] == 'NSC'):
                        grp['sqo'] = 'Low Toxicity'
                    else:
                        grp['sqo'] = 'Moderate Toxicity'
            else:
                if (grp['significance'] == 'NSC'):
                    grp['sqo'] = 'Nontoxic'
                else:
                    grp['sqo'] = 'Low Toxicity'
        else:
            grp['sqo'] = 'Nontoxic'
    elif (grp['species'] == 'Mytilus galloprovincialis'):
        if (grp['mean'] < 80):
            if (grp['pctcontrol'] < 77):
                if (grp['pctcontrol'] < 42):
                    grp['sqo'] = 'High Toxicity'
                else:
                    if (grp['significance'] == 'NSC'):
                        grp['sqo'] = 'Low Toxicity'
                    else:
                        grp['sqo'] = 'Moderate Toxicity'
            else:
                if (grp['significance'] == 'NSC'):
                    grp['sqo'] = 'Nontoxic'
                else:
                    grp['sqo'] = 'Low Toxicity'
        else:
            grp['sqo'] = 'Nontoxic'
    return grp
summary = summary.apply(getSQO, axis=1)
summary.drop('result', axis=1, inplace=True)
summary.drop('labrep', axis=1, inplace=True)
# group on the following columns and reset as a dataframe rather than groupby object
#summary.groupby(['stationid','agency','sampletypecode','sampletypedesc','toxbatch','species','concentration','endpoint','resultunits','sqo','mean','n','stddev','pctcontrol','significance','qacode']).count().reset_index()
summary = summary.groupby(['stationid','agency','sampletypecode','sampletypedesc','toxbatch','species','concentration','endpoint','resultunits','sqo','mean','n','stddev','pctcontrol','significance','qacode']).size().to_frame(name = 'count').reset_index()
summary.to_csv('output.csv', sep='\t', encoding='utf-8')
print(summary['toxbatch'])

#summary.to_csv('output.csv', sep='\t', encoding='utf-8')
#print(list(summary))
#print(summary[['stationid','agency','sampletypecode','sampletypedesc','toxbatch','species','concentration','endpoint','resultunits','sqo','mean','n','stddev','pctcontrol','significance','qacode']])
#test = summary[['stationid','agency','sampletypecode','sampletypedesc','toxbatch','species','concentration','endpoint','resultunits','sqo','mean','n','stddev','pctcontrol','significance','qacode']]
#test.to_csv('output.csv', sep='\t', encoding='utf-8')
### SUMMARY TABLE END ###
