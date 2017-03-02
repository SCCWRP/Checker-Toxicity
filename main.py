import pandas as pd
import numpy as np
import xlrd
#from ordereddict import OrderedDict
import collections
import math
# place the bight13 toxicity data in a location that the application can access
#df = pd.ExcelFile('/Users/pauls/Documents/Projects/Bight18/Training/clean.xlsx')
df = pd.ExcelFile('./clean.xlsx')

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
    	### and put into dictionary object
    	all_dataframes[count] = tab
	count = count + 1

## SUMMARY RESULT CALCULATIONS

# summary must not be a groupby otherwise below functions wont work
# all_dataframes[1] is the toxicity results data
summary = all_dataframes[1]

# uses apply to get the mean and adds column to original dataframe
def get_mean(grp):
	grp['Mean'] = grp['Result'].mean()
        return grp
summary = summary.groupby(['StationID','ToxBatch','FieldReplicate']).apply(get_mean)

## uses apply to get the mean and adds column to original dataframe
def get_sum(grp):
    	grp['N'] = grp['FieldReplicate'].sum()
        return grp
summary = summary.groupby(['StationID','ToxBatch','FieldReplicate']).apply(get_sum)

## uses apply to get the mean and adds column to original dataframe
def get_std(grp):
	grp['StdDev'] = grp['Result'].std()
	return grp
summary = summary.groupby(['StationID','ToxBatch','FieldReplicate']).apply(get_std)

def get_variance(grp):
    	grp['Variance'] = grp['StdDev'].apply(lambda x: x ** 2 )
get_variance(summary)

# get all control records
cneg = summary[['StationID','ToxBatch','SampleTypeCode','Mean']].where(summary['SampleTypeCode'] == 'CNEG')
# get all non control records
nocneg = summary[['StationID','ToxBatch','SampleTypeCode','Mean']].where(summary['SampleTypeCode'] != 'CNEG')
# get all reference toxicant records just save them for now
reference_toxicants = summary.loc[summary['Matrix'].isin(['Reference Toxicant'])]
# drop all reference toxicants from the summary dataframe - not a part of summary results
summary = summary.loc[~summary['Matrix'].isin(['Reference Toxicant'])]

cneg = cneg.dropna()
nocneg = nocneg.dropna()

cneg['Unique'] = np.nan
nocneg['Unique'] = np.nan

control_mean = cneg.groupby(['StationID','ToxBatch','Mean', 'SampleTypeCode'])['Unique'].nunique().reset_index()
result_mean = nocneg.groupby(['StationID','ToxBatch','Mean', 'SampleTypeCode'])['Unique'].nunique().reset_index()

cneg_stats = summary[['StationID','ToxBatch','SampleTypeCode','N','StdDev','Mean','Variance']].where(summary['SampleTypeCode'] == 'CNEG')
cneg_stats = cneg_stats.dropna()
cneg_stats['Unique'] = np.nan
control_mean_stats = cneg_stats.groupby(['StationID','ToxBatch','N','StdDev','Mean','Variance'])['Unique'].nunique().reset_index()

## create a dictionary lookup of toxbatch keys and corresponding important values
# drop unique column we used earlier
control_mean_stats.drop('Unique', axis=1, inplace=True)
# make toxbatch the index - we already group so it is unique
control_mean_dict = control_mean.set_index('ToxBatch')['Mean'].to_dict()
control_mean_stats.set_index("ToxBatch", drop=True, inplace=True)
control_mean_stats_dict = control_mean_stats.to_dict(orient="index")

def getPctControl(row):
    	## toxbatch control should always be 100
    	if(row['SampleTypeCode'] == 'CNEG'):
        	row['PctControl'] = 100
       	else:
            	if row['ToxBatch'] in control_mean_dict:
                	# if the toxbatch is in the lookup dictionary then
                	# divide the result mean from the control mean and times by 100
                	#row['PctControl'] = ((row['Mean']/control_mean_dict[row['ToxBatch']]) * 100)
                	row['PctControl'] = ((row['Mean']/control_mean_stats_dict[row['ToxBatch']]['Mean']) * 100)
                else:
                    	# not sure what should happen with a reference toxicant
                    	# reference toxicants should be dropped - records are not include with summary table
                    	row['PctControl'] = 0
        return row

summary = summary.apply(getPctControl, axis=1)
print(summary)
