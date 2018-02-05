import sys
import pandas as pd
import numpy as np
from scipy import stats
import xlrd
#from ordereddict import OrderedDict
import collections
import math
from math import sin, cos, sqrt, atan2, radians
import urllib, json
stdout = sys.stdout
reload(sys)
sys.setdefaultencoding('utf-8')
sys.stdout = stdout
import os
from sqlalchemy import create_engine
from pandas import DataFrame


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
batch = all_dataframes[0]
result = all_dataframes[1]
summary = all_dataframes[1]
wq = all_dataframes[2]
## CHECKS ##
def checkData(statement,column,warn_or_error,error_label,human_error,dataframe):
    for item_number in statement:
        unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
        dcAddErrorToList(error_label,item_number,unique_error,dataframe)
def checkLogic(statement,column,warn_or_error,error_label,human_error,dataframe):
    for item_number in statement:
        unique_error = '{"column": "%s", "error_type": "%s", "error": "%s"}' % (column,warn_or_error,human_error)
        dcAddErrorToList(error_label,item_number,unique_error,dataframe)

## LOGIC CHECKS ##

## END LOGIC CHECKS ##
# Eric - Station Check - A lab is assigned both stations and test species. Check to see if the submission matches both. If not, report back that the station and/or species doesn't match what was assigned. Possible error: "Unassigned Station" and "Unassigned Species". The html link for both station and species should point to a labs assigned values.
# Creates dataframe from database
eng = create_engine('postgresql://b18read:1969$Harbor@192.168.1.16:5432/bight2018')
spec_results = eng.execute("select stationid,lab,parameter from sample_assignment_table where datatype='Toxicity';")
db = DataFrame(spec_results.fetchall())
db.columns = spec_results.keys()

# Creates new field of unique stationID/species assignments 
db['unique_assignment'] = db.apply(lambda row: row.stationid + '|' + row.parameter, axis=1)
dblabs = db.groupby('lab')['unique_assignment'].unique().reset_index()

# Creates new dataframe from submitted data containing only pertinent fields and produces unique stationID/species submission
sublabs = result[['lab','stationid','species']]
sublabs = sublabs.drop(sublabs[sublabs['stationid']=='0000'].index)
sublabs['unique_submission'] = sublabs.apply(lambda row: row.stationid + '|' + row.species, axis=1)

# Merges assignment dataframe and submitted data dataframe by lab.
mer = pd.merge(sublabs[['lab','unique_submission']], dblabs, how = 'left', on = 'lab').set_index(sublabs.index)

# checks submitted stationID against assigned stationIDs
errorLog(mer.loc[mer.apply(lambda row: row.unique_submission  in row.unique_assignment, axis=1) == False])
checkData(mer.loc[mer.apply(lambda row: row.unique_submission  in row.unique_assignment, axis=1) == False].index.tolist(),'species','Custom Field','error','Unassigned Station or Species',result)

# Eric - A toxicity submission (batch, results,wq) requires that the field data be submitted first. To check all unique Result/StationID records should have a corresponding record in Field/Grab/StationID (make sure it wasn't abandoned also). This should be an error.
#### === THERE WAS NO STATIONID INCLUDED IN THE BATCH TAB === ####
eng = create_engine('postgresql://b18read:1969$Harbor@192.168.1.16:5432/bight2018')
tbl_stationoccupation_df = eng.execute("select tbl_stationoccupation.stationid, tbl_stationoccupation.abandoned from tbl_stationoccupation ")
db2 = DataFrame(tbl_stationoccupation_df.fetchall())
db2.columns = tbl_stationoccupation_df.keys()
# merged the databse(tbl_stationoccupation.stationid) with result (stationid)
merge_db2_result = result[['stationid']].merge(db2, how ='left', on = 'stationid')

#- Message should read "A comment is required when 'Yes' is listed in the Abandoned field
checkData(merge_db2_result.where(merge_db2_result['abandoned'].isin(['Yes'])).dropna(axis = 0, how = 'all').index.tolist(),'abandoned','occupation','Error','StationId is an Abandoned field',merge_db2_result)

merge_db2_wq = wq[['stationid']].merge(db2, how ='left', on = 'stationid')

#- Message should read "A comment is required when 'Yes' is listed in the Abandoned field
checkData(merge_db2_wq.where(merge_db2_wq['abandoned'].isin(['Yes'])).dropna(axis = 0, how = 'all').index.tolist(),'abandoned','occupation','Error','StationId is an Abandoned field',merge_db2_wq)


# db.to_csv('data_2', header=None, index=None, sep=' ', mode='a')
