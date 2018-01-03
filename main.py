import sys
import pandas as pd
import numpy as np
from scipy import stats
import xlrd
#from ordereddict import OrderedDict
import collections
import math
import urllib, json
stdout = sys.stdout
reload(sys)
sys.setdefaultencoding('utf-8')
sys.stdout = stdout

## COMMON FUNCTIONS
# PRINT all errors TO errorLog function placeholder used by web checker to print to multiple places
def errorLog(message):
        print(message)

def dcAddErrorToList(error_column, row, error_to_add,df):
	df.ix[int(row), 'row'] = str(row)
	if error_column in df.columns:
		# check if cell value is empty (nan) 
		if(pd.isnull(df.ix[int(row), error_column])):
			# no data exists in cell so add error
	      		df.ix[int(row), error_column] = error_to_add
			errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
		else:
			# a previous error was recorded so append to it
			# even though there may be data lets check to make sure it is not empty
			if str(df.ix[int(row), error_column]):
				#errorLog("There is already a previous error recorded: %s" % str(df.ix[int(row), error_column]))
				df.ix[int(row), error_column] = str(df.ix[int(row), error_column]) + "," + error_to_add
				errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
			else:
				#errorLog("No error is recorded: %s" % str(df.ix[int(row), error_column]))
	      			df.ix[int(row), error_column] = error_to_add
				errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	else:
		df.ix[int(row), error_column] = error_to_add
		errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	return df

### WORKSPACE START ###
# place the bight13 toxicity data in a location that the application can access
#df = pd.ExcelFile('/Users/pauls/Documents/Projects/Bight18/Training/clean.xlsx')
df = pd.ExcelFile('./errors.xlsx')

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
		errorLog('The application is skipping sheet "%s" because it is empty' % tab)
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
# result is the toxicity results data
batch = all_dataframes[0]
result = all_dataframes[1]
summary = all_dataframes[1]
wq = all_dataframes[2]
''' need working
for dataframe in all_dataframes.keys():
	df_sheet_and_table_name = dataframe.strip().split(" - ")
	table_name = str(df_sheet_and_table_name[2])
	if table_name == "tbltoxicitybatchinformation":
		batch = all_dataframes[dataframe]
	if table_name == "tbltoxicityresults":
		result = all_dataframes[dataframe]
		summary = all_dataframes[dataframe]
	if table_name == "tbltoxicitywq":
		wq = all_dataframes[dataframe]
'''

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
reference_toxicants = summary.loc[summary['matrix'].isin(['reference toxicant'])]
# drop all reference toxicants from the summary dataframe - not a part of summary results
summary = summary.loc[~summary['matrix'].isin(['reference toxicant'])]

cneg = cneg.dropna()
nocneg = nocneg.dropna()

cneg['unique'] = np.nan
nocneg['unique'] = np.nan

control_mean = cneg.groupby(['stationid','toxbatch','mean', 'sampletypecode'])['unique'].nunique().reset_index()
result_mean = nocneg.groupby(['stationid','toxbatch','mean', 'sampletypecode'])['unique'].nunique().reset_index()

# create controlvalue column
# copy control_mean dataframe column mean to controlvalue
control_mean['controlvalue'] = control_mean['mean']
summary = summary.merge(control_mean[['toxbatch','controlvalue']], how = 'left', on = ['toxbatch'])

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
		single_tail = p/2
		summary.ix[index, 'pvalue'] = single_tail #we divide by 2 to make it a 1 tailed
		if (t < 0):
			summary.ix[index, 'sigeffect'] = 'NSC'
		else:
			if (single_tail <= .05):
				summary.ix[index, 'sigeffect'] = 'SC'
			else:
				summary.ix[index, 'sigeffect'] = 'NSC'
getPValue(summary)

## author - Tyler Vu 
def getSQO(grp):
	#if(grp[grp.index.map(lambda x: x[0] in species)]):
    	#if(grp['species'].isin(['EE','Eohaustorius estuarius'])):
    	if(grp['species'] == 'EE'):
        	if(grp['mean'] < 90):
            		if (grp['pctcontrol'] < 82):
                		if (grp['pctcontrol'] < 59):
                    			grp['sqocategory'] = 'High Toxicity'
                		else:
                    			if (grp['sigeffect'] == 'NSC'):
                        			grp['sqocategory'] = 'Low Toxicity'
                    			else:
                        			grp['sqocategory'] = 'Moderate Toxicity'
            		else:
                		if (grp['sigeffect'] == 'NSC'):
                    			grp['sqocategory'] = 'Nontoxic'
              			else:
					grp['sqocategory'] = 'Low Toxicity'
                else:
            		grp['sqocategory'] = 'Nontoxic'
    	#elif (grp['species'].isin(['MG','Mytilus galloprovincialis'])):
    	elif (grp['species'] == 'MG'):
        	if (grp['mean'] < 80):
            		if (grp['pctcontrol'] < 77):
                		if (grp['pctcontrol'] < 42):
                    			grp['sqocategory'] = 'High Toxicity'
                		else:
                    			if (grp['sigeffect'] == 'NSC'):
                        			grp['sqocategory'] = 'Low Toxicity'
                    			else:
                        			grp['sqocategory'] = 'Moderate Toxicity'
            		else:
                		if (grp['sigeffect'] == 'NSC'):
                    			grp['sqocategory'] = 'Nontoxic'
                		else:
                    			grp['sqocategory'] = 'Low Toxicity'
        	else:
            		grp['sqocategory'] = 'Nontoxic'
    	return grp
summary = summary.apply(getSQO, axis=1)
summary.drop('result', axis=1, inplace=True)
summary.drop('labrep', axis=1, inplace=True)
# group on the following columns and reset as a dataframe rather than groupby object
#summary = summary.groupby(['stationid','lab','sampletypecode','toxbatch','species','concentration','endpoint','resultunits','sqocategory','mean','n','stddev','pctcontrol','sigeffect','qacode']).size().to_frame(name = 'count').reset_index()
### SUMMARY TABLE END ###

## SUMMARY TABLE CHECKS ##
# the three blocks of code and corresponding for loops could be combined into one simpler function
def checkSummary(statement,column,warn_or_error,error_label,human_error):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,summary)

# 1 - WARNING TO CHECK FOR DATA ENTRY ERRORS IF THE STANDARD DEVIATION FOR A SAMPLE EXCEEDS 50 
errorLog("## WARNING TO CHECK FOR DATA ENTRY ERRORS IF THE STANDARD DEVIATION FOR A SAMPLE EXCEEDS 50 ##")
errorLog(summary.loc[(summary["stddev"] > 50)])
checkSummary(summary.loc[(summary["stddev"] > 50)].index.tolist(),'StdDev','Custom Toxicity','error','Warning standard deviation exceeds 50.')
# 2 - MEAN SHOULD BE GREATER THAN 90 WHERE SPECIES IS EQUAL TO "EOHAUSTORIUS ESTUARIES" OR "EE" AND SAMPLETYPECODE IS EQUAL TO "CNEG"
errorLog("## MEAN SHOULD BE GREATER THAN 90 WHERE SPECIES IS EQUAL TO EOHAUSTORIUS ESTUARIES OR EE AND SAMPLETYPECODE IS EQUAL TO CNEG##")
errorLog(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 90)])
checkSummary(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 90)].index.tolist(),'Mean','Custom Toxicity','error','Does not meet control acceptability criterion; mean control value < 90')
# 3 - MEAN SHOULD BE GREATER THAN 70 WHERE SPECIES IS EQUAL TO "MYTILUS GALLOPROVINIALIS" OR "MG" AND SAMPLETYPECODE IS EQUAL TO "CNEG"
errorLog("## MEAN SHOULD BE GREATER THAN 70 WHERE SPECIES IS EQUAL TO MYTILUS GALLOPROVINIALIS OR MG AND SAMPLETYPECODE IS EQUAL TO CNEG ##")
errorLog(summary.loc[(summary['species'].isin(['Mytilus galloprovinialis','MG'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 70)])
checkSummary(summary.loc[(summary['species'].isin(['Mytilus galloprovinialis','MG'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 70)].index.tolist(),'Mean','Custom Toxicity','error','Does not meet control acceptability criterion; mean control value < 70')
# 4 - COEFFICIENT VARIANCE SHOULD NOT BE GREATER THAN 11.9 WHERE SPECIES IS EQUAL TO "EOHAUSTORIUS ESTUARIES" OR "EE" AND SAMPLETYPECODE IS EQUAL TO "CNEG" 
errorLog("## COEFFICIENT VARIANCE SHOULD NOT BE GREATER THAN 11.9 WHERE SPECIES IS EQUAL TO EOHAUSTORIUS ESTUARIES OR EE AND SAMPLETYPECODE IS EQUAL TO CNEG ##")
errorLog(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['coefficientvariance'] > 11.9)])
checkSummary(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['coefficientvariance'] > 11.9)].index.tolist(),'CoefficientVariance','Custom Toxicity','error','Does not meet control acceptability criterion; coefficient value > 11.9')
## END SUMMARY TABLE CHECKS ##
# rename a few columns to match with existing b13 column names
summary.rename(columns={"resultunits": "units"}, inplace=True)
# group on the following columns and reset as a dataframe rather than groupby object
summary = summary.groupby(['stationid','lab','sampletypecode','toxbatch','species','concentration','endpoint','units','sqocategory','mean','n','stddev','pctcontrol','pvalue','tstat','sigeffect','qacode','controlvalue']).size().to_frame(name = 'count').reset_index()
summary.to_csv('output.csv', sep='\t', encoding='utf-8')

## END SUMMARY TABLE CHECKS ##

## CHECKS ##
def checkData(statement,column,warn_or_error,error_label,human_error,dataframe):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,dataframe)
def checkLogic(statement,column,warn_or_error,error_label,human_error,dataframe):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,dataframe)

def dcValueAgainstMultipleValues (field,ucfield,listname,listfield,df):
    # author - Jordan Golemo
    # example: field = qacode, ucfield = QACode, listname = lu_toxtestacceptability, listfield = testacceptability, df = results
    # get published lookup list values
    url_search = "https://gis.sccwrp.org/arcgis/rest/services/bight2018%s/FeatureServer/0/query?where=1=1&returnGeometry=false&outFields=*&f=json" % listname
    jsonurl = urllib.urlopen(url_search).read()
    jsonlist = json.loads(jsonurl)
    list_of_codes = []
    # add lookup list values to array
    for i in range(len(jsonlist['features'])):
        list_of_codes.append(jsonlist['features'][i]['attributes'][listfield])
    # submitted field to check against lookup list
    df_search = df[field]
    df_search = df_search.fillna('empty')
    for i in range(len(df_search)):
        # submitted individual field/cell values are separated by commas like: A,B,C
        df_split = df_search[i].replace(',',' ').split()
        for j in range(len(df_split)):
            # find out if individual element is not in url lookup list
            if df_split[j] not in list_of_codes:
                invalid_code = df_search[i]
                # required so it cant be empty
                if df_split[j] == 'empty':
                    checkData(df[df[field].isnull()].tmp_row.tolist(),ucfield,'Toxicity Error','error','A code is required: <a href=http://checker.sccwrp.org/checker/scraper?action=help&layer=%s target=_blank>%s</a>' % (listname,listname),df)
                else:
                    checkData(df.loc[df[field]==invalid_code].tmp_row.tolist(),ucfield,'Toxicity Error','error','You have submitted an invalid code: %s. Please see lookup list: <a href=http://checker.sccwrp.org/checker/scraper?action=help&layer=%s target=_blank>%s</a>' % (df_split[j],listname,listname),df)

## LOGIC ##
# 1 - All records for each table must have a corresponding record in the other tables due on submission. Join tables on Agency/LabCode and ToxBatch/QABatch
### first find matched rows based on toxbatch and result and put into a separate dataframe
brmatch = pd.merge(batch,result, on=['toxbatch','lab'], how='inner')
### check batch to see which combo toxbatch and lab are not in the matched/merged dataframe above 
### check result to see which combo toxbatch and lab are not in the matched/merged dataframe
### make sure there are records that match between batch and result - otherwise big problem
if len(brmatch.index) != 0:
	# EACH TOXICITY BATCH INFORMATION RECORD MUST HAVE A CORRESPONDING TOXICITY RESULT RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE.
	errorLog("## EACH TOXICITY BATCH INFORMATION RECORD MUST HAVE A CORRESPONDING TOXICITY RESULT RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE. ##")
	errorLog(batch[(~batch.toxbatch.isin(brmatch.toxbatch))&(batch.lab.isin(brmatch.lab))])
	checkLogic(batch[(~batch.toxbatch.isin(brmatch.toxbatch))&(batch.lab.isin(brmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Batch Information record must have a corresponding Toxicity Result record. Records are matched on ToxBatch and LabCode.',batch)
else:
	# YOU HAVE ZERO MATCHING RECORDS BETWEEN TOXICITY BATCH AND RESULTS
	errorLog("## YOU HAVE ZERO MATCHING RECORDS BETWEEN TOXICITY BATCH AND RESULTS ##")
	unique_error = '{"column": "ToxBatch", "error_type": "Logic Error", "error": "Each Toxicity Batch Information record must have a corresponding Toxicity Result record. You have zero matching records between Toxicity Batch and Results"}'
	dcAddErrorToList('error',0,unique_error,batch)

errorLog("## EACH TOXICITY RESULT RECORD MUST HAVE A CORRESPONDING TOXICITY BATCH RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE. ##")
errorLog(result[(~result.toxbatch.isin(brmatch.toxbatch))&(result.lab.isin(brmatch.lab))])
checkLogic(result[(~result.toxbatch.isin(brmatch.toxbatch))&(result.lab.isin(brmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Result record must have a corresponding Toxicity Batch record. Records are matched on ToxBatch and LabCode.',result)

### second find matched rows based on result and wq and put into a separate dataframe
rwmatch = pd.merge(result,wq, on=['toxbatch','lab'], how='inner')
errorLog("## EACH TOXICITY RESULT INFORMATION RECORD MUST HAVE A CORRESPONDING TOXICITY WQ RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE. ##")
errorLog(result[(~result.toxbatch.isin(rwmatch.toxbatch))&(result.lab.isin(rwmatch.lab))])
checkLogic(result[(~result.toxbatch.isin(rwmatch.toxbatch))&(result.lab.isin(rwmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Result Information record must have a corresponding Toxicity WQ record. Records are matched on ToxBatch and LabCode.',result)
errorLog("## EACH TOXICITY WQ INFORMATION RECORD MUST HAVE A CORRESPONDING TOXICITY RESULTS RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE. ##")
errorLog(wq[(~wq.toxbatch.isin(rwmatch.toxbatch))&(wq.lab.isin(rwmatch.lab))])
checkLogic(wq[(~wq.toxbatch.isin(rwmatch.toxbatch))&(wq.lab.isin(rwmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity WQ Information record must have a corresponding Toxicity Results record. Records are matched on ToxBatch and LabCode.',wq)

### third find matched rows based on batch and wq and put into a separate dataframe
bwmatch = pd.merge(batch,wq, on=['toxbatch','lab'], how='inner')
errorLog("## EACH TOXICITY BATCH INFORMATION RECORD MUST HAVE A CORRESPONDING TOXICITY WQ RECORD. RECORDS ARE MATCHED ON TOXBATCH AND LABCODE. ##")
errorLog(batch[(~batch.toxbatch.isin(bwmatch.toxbatch))&(batch.lab.isin(bwmatch.lab))])
checkLogic(batch[(~batch.toxbatch.isin(bwmatch.toxbatch))&(batch.lab.isin(bwmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity Batch Information record must have a corresponding Toxicity WQ record. Records are matched on ToxBatch and LabCode.',batch)
errorLog("## Each Toxicity WQ Information record must have a corresponding Toxicity Batch record. Records are matched on ToxBatch and LabCode. ##")
errorLog(wq[(~wq.toxbatch.isin(bwmatch.toxbatch))&(wq.lab.isin(bwmatch.lab))])
checkLogic(wq[(~wq.toxbatch.isin(bwmatch.toxbatch))&(wq.lab.isin(bwmatch.lab))].index.tolist(),'ToxBatch','Logic Error','error','Each Toxicity WQ Information record must have a corresponding Toxicity Batch record. Records are matched on ToxBatch and LabCode.',wq)

# 2 - Check for the minimum number of replicates - ee and mg = 5 and na = 10
## first get a lab replicate count grouped on stationid, toxbatch, and species
dfrep = pd.DataFrame(result.groupby(['stationid','toxbatch','species']).size().reset_index(name='replicatecount'))
## merge the lab replicant group with results so that you can get the tmp_row - the lab rep count will be matched with each lab rep
## we will want to highlight them as a group rather than by row
dfrep = pd.merge(dfrep,result, on=['stationid','toxbatch','species'], how='inner')
errorLog("## A MINIMUM NUMBER OF 5 REPLICATES ARE REQUIRED FOR SPECIES EOHAUSTORIUS ESTUARIUS AND MYTILUS GALLOPROVINCIALIS ##")
errorLog(dfrep.loc[(dfrep['species'].isin(['Eohaustorius estuarius','EE','Mytilus galloprovincialis','MG'])) & (dfrep['replicatecount'] < 5)])
checkLogic(dfrep.loc[(dfrep['species'].isin(['Eohaustorius estuarius','EE','Mytilus galloprovincialis','MG'])) & (dfrep['replicatecount'] < 5)].tmp_row.tolist(),'LabRep','Logic Error','error','A minimum number of 5 replicates are required for species Eohaustorius estuarius and Mytilus galloprovincialis',result)
errorLog("## A MINIMUM NUMBER OF 10 REPLICATES ARE REQUIRED FOR SPECIES NEANTHES ARENACEODENTATA ##")
errorLog(dfrep.loc[(dfrep['species'] == 'NA') & (dfrep['replicatecount'] < 10)])
checkLogic(dfrep.loc[(dfrep['species'] == 'NA') & (dfrep['replicatecount'] < 10)].tmp_row.tolist(),'LabRep','Logic Error','error','A minimum number of 10 replicates are required for species Neanthes arenaceodentata',result)

# 3. EACH BS or SWI BATCH MUST HAVE A "REFERENCE TOXICANT" BATCH WITHIN A SPECIFIED DATE RANGE.
errorLog("## EACH BS or SWI BATCH MUST HAVE A REFERENCE TOXICANT BATCH WITHIN A SPECIFIED DATE RANGE. ##")
# get reference toxicant dataframe
batchrt = batch[['toxbatch','teststartdate', 'actualtestduration', 'actualtestdurationunits', 'referencebatch']].where(batch['matrix'].isin(['RT','Reference Toxicant']))
# drop emptys
batchrt = batchrt.dropna()
if len(batchrt.index) != 0:
	# get bs dataframe added swi on 21june17
	batchbs = batch[['toxbatch', 'matrix', 'species', 'teststartdate', 'actualtestduration', 'actualtestdurationunits', 'referencebatch','tmp_row']].where(batch['matrix'].isin(['BS','SWI','Bulk Sediment (whole sediment)','Sediment Water Interface']))
	# drop emptys
	batchbs = batchbs.dropna()
	# get bs dataframe
	if len(batchbs.index) != 0:
		# find any bs batch records with a missing rt 
		errorLog(batchbs[(~batchbs.referencebatch.isin(batchrt.toxbatch))])
		checkData(batchbs[(~batchbs.referencebatch.isin(batchrt.toxbatch))].tmp_row.tolist(),'Matrix','Toxicity Error','error','BS or SWI batch record is missing reference toxicant batch record',batch)
		# merge bs and rt
		bsmerge = pd.merge(batchbs, batchrt, how = 'inner', on = ['referencebatch'])
		# create date range column
		def checkRTDate(grp):
			grp['teststartdate_x'] = pd.to_datetime(grp['teststartdate_x'])
       			grp['teststartdate_y'] = pd.to_datetime(grp['teststartdate_y'])
			d = grp['teststartdate_x'] - grp['teststartdate_y']
       			grp['daterange'] = abs(d.days)
			return grp
		bsmerge = bsmerge.apply(checkRTDate, axis = 1)
		# checks by species and datarange
		errorLog(bsmerge.loc[(bsmerge['species'] == 'EE') & (bsmerge['daterange'] > 10)])
		checkLogic(bsmerge.loc[(bsmerge['species'] == 'EE') & (bsmerge['daterange'] > 10)].tmp_row.tolist(),'Matrix','Logic Error','toxicity_errors','Each BS or SWI batch must have a Reference Toxicant batch within a specified date range: EE less than 10 days',batch)
		errorLog(bsmerge.loc[(bsmerge['species'] == 'MG') & (bsmerge['daterange'] > 2)])
		checkLogic(bsmerge.loc[(bsmerge['species'] == 'MG') & (bsmerge['daterange'] > 2)].tmp_row.tolist(),'Matrix','Logic Error','toxicity_errors','Each BS or SWI batch must have a Reference Toxicant batch within a specified date range: MG less than 2 days',batch)
		errorLog(bsmerge.loc[(bsmerge['species'] == 'NA') & (bsmerge['daterange'] > 28)])
	else:
		unique_error = '{"column": "Matrix", "error_type": "Logic Error", "error": "A submission requires a Bulk Sediment record in batch submission"}'
		addErrorToList('toxicity_errors',0,unique_error,all_dataframes[store_keys[0]])
		errorsCount('custom')
else:
	unique_error = '{"column": "Matrix", "error_type": "Logic Error", "error": "A submission requires a Reference Toxicant record in batch submission"}'
	addErrorToList('toxicity_errors',0,unique_error,all_dataframes[store_keys[0]])
	errorsCount('custom')
## END LOGIC CHECKS ##

## BATCH CHECKS ##
# 1. EACH BATCH WITH A MATRIX OF BS MUST INCLUDE A CORRESPONDING RESULT CNEG SAMPLE
errorLog("## EACH BATCH WITH A MATRIX OF BS MUST INCLUDE A CORRESPONDING RESULT CNEG SAMPLE ##")
# first get unique cneg records from result dataframe
bsresult = result[['toxbatch','sampletypecode']].where(result['sampletypecode'] == 'CNEG')
bsresult = bsresult.dropna() 
bsresult['unique'] = np.nan
bsresult = bsresult.groupby(['toxbatch','sampletypecode'])['unique'].nunique().reset_index()
# second get unique batch records with a matrix of bs
bsbatch = batch[['toxbatch','matrix','tmp_row']].where(batch['matrix'].isin(['BS','SWI','Bulk Sediment (whole sediment)','Sediment Water Interface']))
bsbatch = bsbatch.dropna()
bsbatch['unique'] = np.nan
bsbatch = bsbatch.groupby(['toxbatch','matrix','tmp_row'])['unique'].nunique().reset_index()
# merge unique cneg and batch records on where they match
bsmerge = bsbatch.merge(bsresult, on='toxbatch', how='inner')
bslocate = bsbatch[(~bsbatch.toxbatch.isin(bsmerge.toxbatch))].toxbatch.tolist()
# label batch records
errorLog(bsbatch[(~bsbatch.toxbatch.isin(bsmerge.toxbatch))])
checkData(bsbatch[(~bsbatch.toxbatch.isin(bsmerge.toxbatch))].tmp_row.tolist(),'Result/SampleTypeCode','Toxicity Error','error','Each batch with a matrix of BS or SWI must include a corresponding result CNEG sample',batch)
# 2. EACH BATCH WITH A MATRIX OF RT MUST INCLUDE A CORRESPONDING RESULT WITH SAMPLETYPECODE = RFNH3.
errorLog("## EACH BATCH WITH A MATRIX OF RT MUST INCLUDE A CORRESPONDING RESULT WITH SAMPLETYPECODE = RFNH3. ##")
# first get unique rfnh3 records from result dataframe
rtresult = result[['toxbatch','sampletypecode']].where(result['sampletypecode'] == 'RFNH3')
rtresult = rtresult.dropna() 
rtresult['unique'] = np.nan
rtresult = rtresult.groupby(['toxbatch','sampletypecode'])['unique'].nunique().reset_index()
# second get unique batch records with a matrix of rt
rtbatch = batch[['toxbatch','matrix','tmp_row']].where(batch['matrix'].isin(["Reference Toxicant","RT"]))
rtbatch = rtbatch.dropna()
rtbatch['unique'] = np.nan
rtbatch = rtbatch.groupby(['toxbatch','matrix','tmp_row'])['unique'].nunique().reset_index()
# merge unique rt and batch records on where they match
rtmerge = rtbatch.merge(rtresult, on='toxbatch', how='inner')
errorLog(rtbatch[(~rtbatch.toxbatch.isin(rtmerge.toxbatch))])
checkData(rtbatch[(~rtbatch.toxbatch.isin(rtmerge.toxbatch))].tmp_row.tolist(),'Result/SampleTypeCode','Toxicity Error','error','Each batch with a matrix of RT must include a corresponding result SampleTypeCode = RFNH3',batch)

# 3. TESTACCEPTABILITY CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY) author - Jordan Golemo
errorLog("TESTACCEPTABILITY CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY)")
dcValueAgainstMultipleValues ('testacceptability','TestAcceptability','lu_toxtestacceptability','testacceptability',batch)

## END BATCH CHECKS ##

## RESULT CHECKS ##
# 1. CHECK IF SAMPLES WERE TESTED WITHIN 28 DAY HOLDING TIME
errorLog("## CHECK IF SAMPLES WERE TESTED WITHIN 28 DAY HOLDING TIME ##")
# merge result and batch on toxbatch but include teststartdate
df28 = pd.merge(result, batch[['toxbatch', 'teststartdate']], how = 'left', on = 'toxbatch')
# change the following field types to pandas datetime so they can be calculated (we arent changing submitted data)
df28['teststartdate'] = pd.to_datetime(df28['teststartdate'])
df28['samplecollectdate'] = pd.to_datetime(df28['samplecollectdate'])
# put day differences into own column
df28['checkdate'] = df28['teststartdate'] - df28['samplecollectdate']
# locate any records with a greater than 28 period
errorLog(df28.loc[df28['checkdate'].dt.days > 28])
checkData(df28.loc[df28['checkdate'].dt.days > 28].tmp_row.tolist(),'SampleTypeCode','Toxicity Error','error','Samples must be tested within a 28 day holding time.',result)

# 2. REFERENCE TOXICANT IN THE MATRIX FIELD MUST HAVE DATA IN CONCENTRATION FIELD. CAN'T BE -88.
errorLog("## REFERENCE TOXICANT IN THE MATRIX FIELD MUST HAVE DATA IN CONCENTRATION FIELD. CANT BE -88 ##")
errorLog(result.loc[result['matrix'].isin(['Reference Toxicant','RT']) & (result['concentration'] == -88)])
checkData(result.loc[result['matrix'].isin(['Reference Toxicant','RT']) & (result['concentration'] == -88)].tmp_row.tolist(),'Concentration','Toxicity Error','error','A "Reference Toxicant" record in the Matrix field can not have a -88 in the Concentration field',result)

# 3. STATION CHECK - A LAB IS ASSIGNED BOTH STATIONS AND TEST SPECIES. CHECK TO SEE IF THE SUBMISSION MATCHES BOTH.
errorLog("## STATION CHECK - A LAB IS ASSIGNED BOTH STATIONS AND TEST SPECIES. CHECK TO SEE IF THE SUBMISSION MATCHES BOTH. ##")
# concatenate station and species together - used below to match against whats returned from database
result['stationidspecies'] = result['stationid'] + "+" + result['species']
# lab list to search by
lab = result.lab.unique()
for l in lab:
	search_url = "https://gis.sccwrp.org/arcgis/rest/services/Bight18ToxicityAssignedSpecies/FeatureServer/0/query?where=lab=%27{0}%27&1=1&returnGeometry=false&outFields=stationid,lab,species&f=json".format(l)
	errorLog(search_url)
	response = urllib.urlopen(search_url)
	data = json.loads(response.read())
	# loop through json records and build station and species into a single string then add to list 
	search_list = []
	for i in data['features']:
		errorLog(i['attributes']['stationid']+ "+" + i['attributes']['species'])
		search_list.append(i['attributes']['stationid']+ "+" + i['attributes']['species'])
	errorLog(search_list)
	# find stations/species that dont match between submission and whats in database based on lab
	errorLog(result.loc[~result['stationidspecies'].isin(search_list)].stationid.tolist())
	checkData(result.loc[~result['stationidspecies'].isin(search_list)].tmp_row.tolist(),'StationID/Species','Toxicity Error','error','The station and species you submitted fails to match the lab assignment list',result)

# 4. QACODE CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY). author - Jordan Golemo
errorLog("QACODE CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY)")
dcValueAgainstMultipleValues ('qacode','QACode','lu_toxtestacceptability','testacceptability',result)

# drop temporary column
result.drop('stationidspecies', axis=1, inplace=True)
## END RESULT CHECKS ##

## START WQ CHECKS ##
# 1. CHECK THAT WATER QUALITY PARAMETERS ARE WITHIN ACCEPTABLE RANGES.
errorLog("## CHECK THAT WATER QUALITY PARAMETERS ARE WITHIN ACCEPTABLE RANGES. ##")
# merge wq and batch on toxbatch to get species from batch
dfwq = pd.merge(wq[['toxbatch','parameter','result']], batch[['toxbatch', 'species']], how = 'left', on = 'toxbatch')
errorLog(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','Mytilus galloprovincialis','EE','MG'])) & (dfwq['parameter'] == 'TEMP') & ((dfwq['result'] < 13) | (dfwq['result'] > 17))])
checkData(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','Mytilus galloprovincialis','EE','MG'])) & (dfwq['parameter'] == 'TEMP') & ((dfwq['result'] < 13) | (dfwq['result'] > 17))].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for TEMP not in acceptable range: must be between 13-17',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','Mytilus galloprovincialis','EE','MG'])) & (dfwq['parameter'] == 'SAL') & ((dfwq['result'] <= 30) | (dfwq['result'] >= 34))])
checkData(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','Mytilus galloprovincialis','EE','MG'])) & (dfwq['parameter'] == 'SAL') & ((dfwq['result'] <= 30) | (dfwq['result'] >= 34))].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for SAL not in acceptable range: must be between 30-34',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'DO') & (dfwq['result'] < 7.5)])
checkData(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'DO') & (dfwq['result'] < 7.5)].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for DO not in acceptable range: must be greater than 7.5',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'PH') & ((dfwq['result'] <= 7.7) | (dfwq['result'] >= 8.3))])
checkData(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'PH') & ((dfwq['result'] <= 7.7) | (dfwq['result'] >= 8.3))].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for PH not in acceptable range: must be between 7.7-8.3',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'NH3T') & (dfwq['result'] > 20)])
checkData(dfwq.loc[(dfwq['species'].isin(['Eohaustorius estuarius','EE'])) & (dfwq['parameter'] == 'NH3T') & (dfwq['result'] > 20)].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for NH3T not in acceptable range: must be less than 20',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Mytilus galloprovincialis','MG'])) & (dfwq['parameter'] == 'DO') & (dfwq['result'] < 4.0)])
checkData(dfwq.loc[(dfwq['species'].isin(['Mytilus galloprovincialis','MG'])) & (dfwq['parameter'] == 'DO') & (dfwq['result'] < 4.0)].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for DO not in acceptable range: must be greater than 4.0',wq)
errorLog(dfwq.loc[(dfwq['species'].isin(['Mytilus galloprovincialis','MG'])) & (dfwq['parameter'] == 'PH') & ((dfwq['result'] <= 7.6) | (dfwq['result'] >= 8.3))])
checkData(dfwq.loc[(dfwq['species'].isin(['Mytilus galloprovincialis','MG'])) & (dfwq['parameter'] == 'PH') & ((dfwq['result'] <= 7.6) | (dfwq['result'] >= 8.3))].index.tolist(),'Result','Toxicity WQ Error','error','Water quality parameter for paramter PH not in acceptable range: must be between 7.6-8.3',wq)
## END WQ CHECKS ##
