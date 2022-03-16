
import urllib, json
import pandas as pd
import numpy as np
from scipy import stats

## COMMON FUNCTIONS

# PRINT all errors TO errorLog function placeholder used by web checker to print to multiple places
def errorLog(message):
	print(message)

def dcAddErrorToList(error_column, row, error_to_add,df):
	df.at[int(row), 'row'] = str(row)
	if error_column in df.columns:
		# check if cell value is empty (nan) 
		if(pd.isnull(df.iat[int(row), df.columns.get_loc(error_column)])):
			# no data exists in cell so add error
			df.at[int(row), error_column] = error_to_add
			errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
		else:
			# a previous error was recorded so append to it
			# even though there may be data lets check to make sure it is not empty
			if str(df.at[int(row), error_column]):
				#errorLog("There is already a previous error recorded: %s" % str(df.iloc[int(row), df.columns.get_loc(error_column)]))
				df.at[int(row), error_column] = str(df.iloc[int(row), df.columns.get_loc(error_column)]) + "," + error_to_add
				errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
			else:
				#errorLog("No error is recorded: %s" % str(df.iloc[int(row), df.columns.get_loc(error_column)]))
				df.at[int(row), error_column] = error_to_add
				errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	else:
		df.at[int(row), error_column] = error_to_add
		errorLog("Row: %s, Error To Add: %s" % (int(row),error_to_add))
	return df

def checkSummary(statement,column,warn_or_error,error_label,human_error,df):
	for item_number in statement:
		unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
		dcAddErrorToList(error_label,item_number,unique_error,df)

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


# -- Stats functions -- #
def getCalculatedValues(grp):                                                                  
	grp['mean'] = grp['result'].mean()
	grp['n'] = grp['fieldreplicate'].sum()
	grp['stddev'] = grp['result'].std()
	grp['variance'] = grp['stddev'].apply(lambda x: x ** 2 )
	grp['coefficientvariance'] = ((grp['stddev']/grp['mean']) * 100)
	return grp


def getPctControl(row, control_mean_dict):
	## toxbatch control should always be 100
    if(row['sampletypecode'] == 'CNEG'):
        row['pctcontrol'] = 100
    else:
        if row['toxbatch'] in control_mean_dict:
			# if the toxbatch is in the lookup dictionary then
			# divide the result mean from the control mean and times by 100
			# OLD LINE row['pctcontrol'] = ((row['mean']/control_mean_stats_dict[row['toxbatch']]['mean']) * 100)
            row['pctcontrol'] = ((row['mean']/control_mean_dict[row['toxbatch']]) * 100)
        else:
            row['pctcontrol'] = np.NaN
    return row

## author - Tyler Vu
def getPValue(summary):
	for index, values in summary['toxbatch'].iteritems():
		station_code = summary.iloc[index, summary.columns.get_loc('stationid')]
		cneg_result = summary[['result']].where((summary['sampletypecode'] == 'CNEG') & (summary['toxbatch'] == values))
		result_both = summary[['result']].where((summary['toxbatch'] == values) & (summary['stationid'] == station_code) )
		cneg_result = cneg_result.dropna()
		result_both = result_both.dropna()
		t, p = stats.ttest_ind(cneg_result, result_both, equal_var = False)
		errorLog(summary.iloc[index])
		summary.at[index, 'tstat'] = t
		single_tail = p/2
		summary.at[index, 'pvalue'] = single_tail #we divide by 2 to make it a 1 tailed
		if (t < 0):
			summary.at[index, 'sigeffect'] = 'NSC'
		else:
			if (single_tail <= .05):
				summary.at[index,'sigeffect'] = 'SC'
			else:
				summary.at[index,'sigeffect'] = 'NSC'


## author - Tyler Vu 
def getSQO(grp):
	#if(grp[grp.index.map(lambda x: x[0] in species)]):
    #if(grp['species'].isin(['EE','Eohaustorius estuarius'])):
		
	if(grp['species'] == 'Eohaustorius estuarius'):
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
	elif (grp['species'] == 'Mytilus galloprovincialis'):
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
	else:
		grp['sqocategory'] = None
	return grp