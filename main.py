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

### WORKSPACE START ###
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

### WORKSPACE END ###

### SUMMARY TABLE START ###

# summary must not be a groupby otherwise below functions wont work
# all_dataframes[1] is the toxicity results data
batch = all_dataframes[0]
result = all_dataframes[1]
wq = all_dataframes[2]
summary = all_dataframes[1]

def getCalculatedValues(grp):                                                                  
	grp['Mean'] = grp['Result'].mean()
	grp['N'] = grp['FieldReplicate'].sum()
	grp['StdDev'] = grp['Result'].std()
	grp['Variance'] = grp['StdDev'].apply(lambda x: x ** 2 )
	grp['CoefficientVariance'] = ((grp['StdDev']/grp['Mean']) * 100)
	return grp
summary = summary.groupby(['StationID','ToxBatch','FieldReplicate']).apply(getCalculatedValues)

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

## prep code control_mean_stats_dict used in getPctControl function
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
## prep code control_mean_stats_dict end

def getPctControl(row):
    	## toxbatch control should always be 100
    	if(row['SampleTypeCode'] == 'CNEG'):
        	row['PctControl'] = 100
       	else:
            	if row['ToxBatch'] in control_mean_dict:
                	# if the toxbatch is in the lookup dictionary then
                	# divide the result mean from the control mean and times by 100
                	row['PctControl'] = ((row['Mean']/control_mean_stats_dict[row['ToxBatch']]['Mean']) * 100)
                else:
                    	# not sure what should happen with a reference toxicant
                    	# reference toxicants should be dropped - records are not include with summary table
                    	row['PctControl'] = 0
        return row

summary = summary.apply(getPctControl, axis=1)

## author - Tyler Vu
def getPValue(summary):
	for index, values in summary['ToxBatch'].iteritems():
		station_code = summary.ix[index, 'StationID']
		cneg_result = summary[['Result']].where((summary['SampleTypeCode'] == 'CNEG') & (summary['ToxBatch'] == values))
		result_both = summary[['Result']].where((summary['ToxBatch'] == values) & (summary['StationID'] == station_code) )
		cneg_result = cneg_result.dropna()
		result_both = result_both.dropna()
		t, p = stats.ttest_ind(cneg_result, result_both, equal_var = False)
		summary.ix[index, 'TStat'] = t
		summary.ix[index, 'PValue'] = p/2 #we divide by 2 to make it a 1 tailed
		if (t < 0):
			summary.ix[index, 'Significance'] = 'NSC'
		else:
			if (p <= .05):
				summary.ix[index, 'Significance'] = 'SC'
			else:
				summary.ix[index, 'Significance'] = 'NSC'
getPValue(summary)

## author - Tyler Vu 
def getSQO(grp):
    if(grp['Species'] == 'Eohaustorius estuarius'):
        if(grp['Mean'] < 90):
            if (grp['PctControl'] < 82):
                if (grp['PctControl'] < 59):
                    grp['SQO'] = 'High Toxicity'
                else:
                    if (grp['Significance'] == 'NSC'):
                        grp['SQO'] = 'Low Toxicity'
                    else:
                        grp['SQO'] = 'Moderate Toxicity'
            else:
                if (grp['Significance'] == 'NSC'):
                    grp['SQO'] = 'Nontoxic'
                else:
                    grp['SQO'] = 'Low Toxicity'
        else:
            grp['SQO'] = 'Nontoxic'
    elif (grp['Species'] == 'Mytilus galloprovincialis'):
        if (grp['Mean'] < 80):
            if (grp['PctControl'] < 77):
                if (grp['PctControl'] < 42):
                    grp['SQO'] = 'High Toxicity'
                else:
                    if (grp['Significance'] == 'NSC'):
                        grp['SQO'] = 'Low Toxicity'
                    else:
                        grp['SQO'] = 'Moderate Toxicity'
            else:
                if (grp['Significance'] == 'NSC'):
                    grp['SQO'] = 'Nontoxic'
                else:
                    grp['SQO'] = 'Low Toxicity'
        else:
            grp['SQO'] = 'Nontoxic'
    return grp
summary = summary.apply(getSQO, axis=1)

### SUMMARY TABLE END ###

#print(summary)
#summary.to_csv('output.csv', sep='\t', encoding='utf-8')

## SUMMARY TABLE CHECKS ##
# the three blocks of code and corresponding for loops could be combined into one simpler function
def checkSummary(statement,column,warn_or_error,error_label,human_error):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,summary)

# 1 - Warning to check for data entry errors if the standard deviation for a sample exceeds 50 
checkSummary(summary.loc[(summary["StdDev"] > 50)].index.tolist(),'StdDev','Custom Warning','warning','Warning standard deviation exceeds 50.')
# 2 - Check that controls meet test acceptability requirements (can be done using summary table data).
### ***** missing needs to be done ******
# 3 - Mean should be greater than 90 where Species is equal to "Eohaustorius estuaries" or "EE" and SampleTypeCode is equal to "CNEG"
checkSummary(summary.loc[(summary['Species'] == 'EE') & (summary['SampleTypeCode'] == 'CNEG') & (summary['Mean'] < 90)].index.tolist(),'Mean','Custom Warning','warning','Does not meet control acceptability criterion; mean control value < 90')
checkSummary(summary.loc[(summary['Species'] == ('Eohaustorius estuarius' or 'EE')) & (summary['SampleTypeCode'] == 'CNEG') & (summary['Mean'] < 90)].index.tolist(),'Mean','Custom Warning','warning','Does not meet control acceptability criterion; mean control value < 90')
# 4 - Mean should be greater than 70 where Species is equal to "Mytilus galloprovinialis" or "MG" and SampleTypeCode is equal to "CNEG"
checkSummary(summary.loc[(summary['Species'] == ('Mytilus galloprovinialis' or 'MG')) & (summary['SampleTypeCode'] == 'CNEG') & (summary['Mean'] < 70)].index.tolist(),'Mean','Custom Warning','warning','Does not meet control acceptability criterion; mean control value < 70')
# 5 - Coefficient Variance should not be greater than 11.9 where Species is equal to "Eohaustorius estuaries" or "EE" and SampleTypeCode is equal to "CNEG" 
checkSummary(summary.loc[(summary['Species'] == ('Eohaustorius estuarius' or 'EE')) & (summary['SampleTypeCode'] == 'CNEG') & (summary['CoefficientVariance'] > 11.9)].index.tolist(),'CoefficientVariance','Custom Warning','warning','Does not meet control acceptability criterion; coefficient value > 11.9')
## END SUMMARY TABLE CHECKS ##

## LOGIC CHECKS ##
def checkLogic(statement,column,warn_or_error,error_label,human_error,dataframe):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,dataframe)
# Each toxbatch record must have a corresponding result record. Records are matched on QABatch and LabCode.
# 1 - All records for each table must have a corresponding record in the other tables due on submission. Join tables on Agency/LabCode and ToxBatch/QABatch
### first find matched rows based on toxbatch and result and put into a separate dataframe
brmatch = pd.merge(batch,result, on=['ToxBatch','QACode'], how='inner')
### check batch to see which combo toxbatch and labcode are not in the matched/merged dataframe above 
### check result to see which combo toxbatch and labcode are not in the matched/merged dataframe
checkLogic(batch[(~batch.ToxBatch.isin(brmatch.ToxBatch))&(batch.QACode.isin(brmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Batch Information record must have a corresponding Toxicity Result record. Records are matched on ToxBatch and LabCode.',batch)
checkLogic(result[(~result.ToxBatch.isin(brmatch.ToxBatch))&(result.QACode.isin(brmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Result record must have a corresponding Toxicity Batch record. Records are matched on ToxBatch and LabCode.',result)
### second find matched rows based on result and wq and put into a separate dataframe
rwmatch = pd.merge(result,wq, on=['ToxBatch','QACode'], how='inner')
checkLogic(result[(~result.ToxBatch.isin(rwmatch.ToxBatch))&(result.QACode.isin(rwmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Result Information record must have a corresponding Toxicity WQ record. Records are matched on ToxBatch and LabCode.',result)
checkLogic(wq[(~wq.ToxBatch.isin(rwmatch.ToxBatch))&(wq.QACode.isin(rwmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity WQ Information record must have a corresponding Toxicity Results record. Records are matched on ToxBatch and LabCode.',wq)
### third find matched rows based on batch and wq and put into a separate dataframe
bwmatch = pd.merge(batch,wq, on=['ToxBatch','QACode'], how='inner')
checkLogic(batch[(~batch.ToxBatch.isin(bwmatch.ToxBatch))&(batch.QACode.isin(bwmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Batch Information record must have a corresponding Toxicity WQ record. Records are matched on ToxBatch and LabCode.',batch)
checkLogic(wq[(~wq.ToxBatch.isin(bwmatch.ToxBatch))&(wq.QACode.isin(bwmatch.QACode))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity WQ Information record must have a corresponding Toxicity Batch record. Records are matched on ToxBatch and LabCode.',wq)
## END LOGIC CHECKS ##

## BATCH CHECKS ##
## END BATCH CHECKS ##

## RESULT CHECKS ##
def checkResults(statement,column,warn_or_error,error_label,human_error,dataframe):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,dataframe)
# 1. EACH BATCH WITH A MATRIX BS MUST HAVE ONE CNEG SAMPLE #
# first get unique cneg records from result dataframe
bsresult = result[['ToxBatch','SampleTypeCode']].where(result['SampleTypeCode'] == 'CNEG')
bsresult = bsresult.dropna() 
bsresult['Unique'] = np.nan
bsresult = bsresult.groupby(['ToxBatch','SampleTypeCode'])['Unique'].nunique().reset_index()
# second get unique batch records with a matrix of bs
bsbatch = batch[['ToxBatch','Matrix','tmp_row']].where(batch['Matrix'] == ("Bulk Sediment (whole sediment)" or "BS"))
bsbatch = bsbatch.dropna()
bsbatch['Unique'] = np.nan
bsbatch = bsbatch.groupby(['ToxBatch','Matrix'])['Unique'].nunique().reset_index()
# merge unique cneg and batch records on where they match
bsmerge = bsbatch.merge(bsresult, on='ToxBatch', how='inner')
# locate the rows in results that aren't in the merge of the two
checkResults(bsresult[(~bsresult.ToxBatch.isin(bsmerge.ToxBatch))].index.tolist(),'SampleTypeCode','Toxicity Error','error','Each batch record with a Matrix = BS must include at least one record with a SampleTypeCode = CNEG',result)


## END RESULT CHECKS ##
# old
'''
unequal_bs = br[(br['Matrix'].isin(['BS'])) & (~br['SampleTypeCode'].isin(["Grab","CNEG"]))]
checkResults(unequal_bs.index,'SampleTypeCode','Custom Error','error','You have a value in the Matrix column. This requires a corresponding value in batch record and the column SampleTypeCode must have or  you put ')

checkResults(all_dataframes[1].loc[(all_dataframes[1]['matrix_y'] == ('RT' or 'Reference Toxicant')) & (all_dataframes[1]['concentration'] == -88)].index.tolist(),'Concentration','Toxicity Error','error','A Reference Toxicant record in the Matrix field can not have a -88 in the Concentration field')
'''

## start of code to find out reference tox lab replicates that are out of range - warning only not error
'''
def get_labrep_max(grp):
	summary['Max'] = grp['LabRep'].max()
	return grp
summary = summary.groupby(['StationID','ToxBatch','Agency']).apply(get_labrep_max)

def get_reftox_date(row):
	if row['Species'] == ('Mytilus galloprovincialis' or 'MG' or 'Eohaustorius estuarius' or 'EE') and row['Max'] < 5:
		print(row)   

summary.apply(get_reftox_date,axis=1)

#print(summary)
summary.to_csv('output.csv', sep='\t', encoding='utf-8')

# "Reference Toxicant" in the Matrix field must have data in the Concentration field...cant be -88
#reftox = (all_dataframes[1].loc[(all_dataframes[1]['Matrix'] == ('RT' or 'Reference Toxicant')) & (all_dataframes[1]['Concentration'] == -88)])
#for item_number in (reftox.index):
#	human_error = 'A "Reference Toxicant" record in the Matrix field can not have a -88 in the Concentration field'
#	unique_error = '{ "column": "Concentration", "error_type": "Toxicity Error", "error": "%s" }' % (human_error)
#	dcAddErrorToList("custom_errors_reftox",int(reftox['tmp_row'].loc[item_number]),unique_error,all_dataframes[1])
'''
