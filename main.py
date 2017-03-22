import pandas as pd
import numpy as np
from scipy import stats
import xlrd
#from ordereddict import OrderedDict
import collections
import math

## COMMON FUNCTIONS
def dcAddErrorToList(error_column, row, error_to_add,df):
	df.ix[int(row), 'row'] = str(row)
	if error_column in df.columns:
		# check if cell value is empty (nan) 
		if(pd.isnull(df.ix[int(row), error_column])):
			# no data exists in cell so add error
	      		df.ix[int(row), error_column] = error_to_add
		else:
			# a previous error was recorded so append to it
			# even though there may be data lets check to make sure it is not empty
			if str(df.ix[int(row), error_column]):
				print("There is already a previous error recorded: %s" % str(df.ix[int(row), error_column]))
				df.ix[int(row), error_column] = str(df.ix[int(row), error_column]) + "," + error_to_add
			else:
				print("No error is recorded: %s" % str(df.ix[int(row), error_column]))
	      			df.ix[int(row), error_column] = error_to_add
	else:
		df.ix[int(row), error_column] = error_to_add
	return df

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
	### create tmp_row for tracking row numbers
	all_dataframes[count]['tmp_row'] = all_dataframes[count].index
	count = count + 1

## OTHER TOXICITY CHECKS

# a reference toxicant in the matrix field must have data in the concentration field...cant be -88
reftox = (all_dataframes[1].loc[(all_dataframes[1]['Matrix'] == ('RT' or 'Reference Toxicant')) & (all_dataframes[1]['Concentration'] == -88)])
for item_number in (reftox.index):
	human_error = 'A "Reference Toxicant" record in the Matrix field can not have a -88 in the Concentration field'
	unique_error = '{ "column": "Concentration", "error_type": "Toxicity Error", "error": "%s" }' % (human_error)
	dcAddErrorToList("custom_errors_reftox",int(reftox['tmp_row'].loc[item_number]),unique_error,all_dataframes[1])

## SUMMARY RESULT TABLE

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

def get_p_value(summary):
	for index, values in summary['ToxBatch'].iteritems():
		station_code = summary.ix[index, 'StationID']
		cneg_result = summary[['Result']].where((summary['SampleTypeCode'] == 'CNEG') & (summary['ToxBatch'] == values))
		result_both = summary[['Result']].where((summary['ToxBatch'] == values) & (summary['StationID'] == station_code) )
		cneg_result = cneg_result.dropna()
		result_both = result_both.dropna()
		t, p = stats.ttest_ind(cneg_result, result_both, equal_var = False)
		summary.ix[index, 't - stat'] = t
		summary.ix[index, 'p - value'] = p/2 #we divide by 2 to make it a 1 tailed
	if (t < 0):
		summary.ix[index, 'Significance'] = 'NSC'
	else:
		if (p <= .05):
			summary.ix[index, 'Significance'] = 'SC'
		else:
			summary.ix[index, 'Significance'] = 'NSC'
get_p_value(summary)

#summary.to_csv('test.csv', sep='\t', encoding='utf-8')
print(summary)
