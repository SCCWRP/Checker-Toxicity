import sys, os, collections
import pandas as pd
import numpy as np

# custom functions
from functions import *

valid_flags = ('-f','--file','-h')

args = sys.argv

#VERBOSE = ('-v' in args) or ('--verbose' in args)


if '-h' in args:
	print(
	'''
	At the command line, run:
	
	python3 main.py [OPTIONAL -f <file path to tox data>]

	If the --file or -f flag is not provided, you will be prompted for the path to the excel file that you will be checking
	'''
	)
	sys.exit()

# if they provided the file path then get the excel file path, if not, prompt the user's input
if not ('-f' in args) or ('--file' in args):
	excel_path = input("Please enter the path to the excel file to check:")
else:
	excel_path = args[args.index('-f' if '-f' in args else '--file') + 1]
	

assert os.path.exists(excel_path), "File path {} not found".format(excel_path)
df = pd.ExcelFile(excel_path)

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

## create a dictionary lookup of toxbatch keys and corresponding control mean values
control_mean_dict = control_mean.set_index('toxbatch')['mean'].to_dict()

summary = summary.apply( lambda row: getPctControl(row, control_mean_dict), axis=1)

# initialize tstat column
summary['tstat'] = np.NaN

getPValue(summary)


summary = summary.apply(getSQO, axis=1)
summary.drop('result', axis=1, inplace=True)
summary.drop('labrep', axis=1, inplace=True)
# group on the following columns and reset as a dataframe rather than groupby object
#summary = summary.groupby(['stationid','lab','sampletypecode','toxbatch','species','concentration','endpoint','resultunits','sqocategory','mean','n','stddev','pctcontrol','sigeffect','qacode']).size().to_frame(name = 'count').reset_index()
### SUMMARY TABLE END ###

## SUMMARY TABLE CHECKS ##
# the three blocks of code and corresponding for loops could be combined into one simpler function


# 1 - WARNING TO CHECK FOR DATA ENTRY ERRORS IF THE STANDARD DEVIATION FOR A SAMPLE EXCEEDS 50 
errorLog("## WARNING TO CHECK FOR DATA ENTRY ERRORS IF THE STANDARD DEVIATION FOR A SAMPLE EXCEEDS 50 ##")
errorLog(summary.loc[(summary["stddev"] > 50)])
checkSummary(summary.loc[(summary["stddev"] > 50)].index.tolist(),'StdDev','Custom Toxicity','error','Warning standard deviation exceeds 50.',summary)
# 2 - MEAN SHOULD BE GREATER THAN 90 WHERE SPECIES IS EQUAL TO "EOHAUSTORIUS ESTUARIES" OR "EE" AND SAMPLETYPECODE IS EQUAL TO "CNEG"
errorLog("## MEAN SHOULD BE GREATER THAN 90 WHERE SPECIES IS EQUAL TO EOHAUSTORIUS ESTUARIES OR EE AND SAMPLETYPECODE IS EQUAL TO CNEG##")
errorLog(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 90)])
checkSummary(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 90)].index.tolist(),'Mean','Custom Toxicity','error','Does not meet control acceptability criterion; mean control value < 90',summary)
# 3 - MEAN SHOULD BE GREATER THAN 70 WHERE SPECIES IS EQUAL TO "MYTILUS GALLOPROVINIALIS" OR "MG" AND SAMPLETYPECODE IS EQUAL TO "CNEG"
errorLog("## MEAN SHOULD BE GREATER THAN 70 WHERE SPECIES IS EQUAL TO MYTILUS GALLOPROVINIALIS OR MG AND SAMPLETYPECODE IS EQUAL TO CNEG ##")
errorLog(summary.loc[(summary['species'].isin(['Mytilus galloprovinialis','MG'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 70)])
checkSummary(summary.loc[(summary['species'].isin(['Mytilus galloprovinialis','MG'])) & (summary['sampletypecode'] == 'CNEG') & (summary['mean'] < 70)].index.tolist(),'Mean','Custom Toxicity','error','Does not meet control acceptability criterion; mean control value < 70',summary)
# 4 - COEFFICIENT VARIANCE SHOULD NOT BE GREATER THAN 11.9 WHERE SPECIES IS EQUAL TO "EOHAUSTORIUS ESTUARIES" OR "EE" AND SAMPLETYPECODE IS EQUAL TO "CNEG" 
errorLog("## COEFFICIENT VARIANCE SHOULD NOT BE GREATER THAN 11.9 WHERE SPECIES IS EQUAL TO EOHAUSTORIUS ESTUARIES OR EE AND SAMPLETYPECODE IS EQUAL TO CNEG ##")
errorLog(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['coefficientvariance'] > 11.9)])
checkSummary(summary.loc[(summary['species'].isin(['Eohaustorius estuarius','EE'])) & (summary['sampletypecode'] == 'CNEG') & (summary['coefficientvariance'] > 11.9)].index.tolist(),'CoefficientVariance','Custom Toxicity','error','Does not meet control acceptability criterion; coefficient value > 11.9',summary)
## END SUMMARY TABLE CHECKS ##
# rename a few columns to match with existing b13 column names
summary.rename(columns={"resultunits": "units"}, inplace=True)
# group on the following columns and reset as a dataframe rather than groupby object
summary = summary.groupby(['stationid','lab','sampletypecode','toxbatch','species','concentration','endpoint','units','sqocategory','mean','n','stddev','pctcontrol','pvalue','tstat','sigeffect','qacode','controlvalue']).size().to_frame(name = 'count').reset_index()


## END SUMMARY TABLE CHECKS ##


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
# disable in standalone mode
# errorLog("TESTACCEPTABILITY CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY)")
# dcValueAgainstMultipleValues ('testacceptability','TestAcceptability','lu_toxtestacceptability','testacceptability',batch)

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
# temp disable in standalone mode
# errorLog("## STATION CHECK - A LAB IS ASSIGNED BOTH STATIONS AND TEST SPECIES. CHECK TO SEE IF THE SUBMISSION MATCHES BOTH. ##")
# # concatenate station and species together - used below to match against whats returned from database
# result['stationidspecies'] = "{}+{}".format(result['stationid'], result['species'])
# # lab list to search by
# lab = result.lab.unique()
# for l in lab:
# 	search_url = "https://gis.sccwrp.org/arcgis/rest/services/Bight18ToxicityAssignedSpecies/FeatureServer/0/query?where=lab=%27{0}%27&1=1&returnGeometry=false&outFields=stationid,lab,species&f=json".format(l)
# 	errorLog(search_url)
# 	response = urllib.urlopen(search_url)
# 	data = json.loads(response.read())
# 	# loop through json records and build station and species into a single string then add to list 
# 	search_list = []
# 	for i in data['features']:
# 		errorLog(i['attributes']['stationid']+ "+" + i['attributes']['species'])
# 		search_list.append(i['attributes']['stationid']+ "+" + i['attributes']['species'])
# 	errorLog(search_list)
# 	# find stations/species that dont match between submission and whats in database based on lab
# 	errorLog(result.loc[~result['stationidspecies'].isin(search_list)].stationid.tolist())
# 	checkData(result.loc[~result['stationidspecies'].isin(search_list)].tmp_row.tolist(),'StationID/Species','Toxicity Error','error','The station and species you submitted fails to match the lab assignment list',result)

# # 4. QACODE CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY). author - Jordan Golemo
# errorLog("QACODE CHECK - A SINGLE QACODE IS REQUIRED BUT MULTIPLE QACODES ARE POSSIBLE (MANY TO MANY)")
# dcValueAgainstMultipleValues ('qacode','QACode','lu_toxtestacceptability','testacceptability',result)

# drop temporary column
# result.drop('stationidspecies', axis=1, inplace=True)
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

if not os.path.exists('output'):
	os.makedirs('output')

writer = pd.ExcelWriter('output/report.xlsx', engine = 'xlsxwriter')

output_dfs = {
	'result_original': result.drop(columns = ['tmp_row','row','error']) if all([c in result.columns for c in ('tmp_row','row','error')]) else result,
	'result_errors': result[~result.error.isna()] if 'error' in result.columns else pd.DataFrame(),
	
	'batch_original': batch.drop(columns = ['tmp_row','row','error']) if all([c in batch.columns for c in ('tmp_row','row','error')]) else batch,
	'batch_errors': batch[~batch.error.isna()] if 'error' in batch.columns else pd.DataFrame(),
	
	'wq_original': wq.drop(columns = ['tmp_row','row','error']) if all([c in wq.columns for c in ('tmp_row','row','error')]) else wq,
	'wq_errors': wq[~wq.error.isna()] if 'error' in wq.columns else pd.DataFrame(),

	'summary': summary
}

for sheetname, df in output_dfs.items():
	df.to_excel(writer, sheet_name=sheetname, index = False)

writer.save()


