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

eng = create_engine('postgresql://b18read:1969$Harbor@192.168.1.16:5432/bight2018')
spec_results = eng.execute("select * from sample_assignment_table where datatype='Toxicity';")
db = DataFrame(spec_results.fetchall())
db.columns = spec_results.keys()
# From database, determines lists of valid stationIDs and species
group1 = db.groupby('lab')['stationid'].unique()
group2 = db.groupby('lab')['parameter'].unique()
labs = pd.concat([group1,group2],axis=1).reset_index()

# Submitted stationID and species for specified lab
result_lab_data = result[['lab','stationid','species']]

# merge dataframe labs and dataframe result_lab_data
merge_labs_resultsData = pd.merge(result_lab_data[['lab','stationid','species']], labs[['lab','stationid','parameter']], how = 'left', on = 'lab')

# creates a new column of boolean series to check if every item in column species is contained in column paramter
merge_labs_resultsData['Check_species'] = merge_labs_resultsData.apply(lambda row: row['species'] in row['parameter'], axis=1)
checkData(merge_labs_resultsData.loc[merge_labs_resultsData['Check_species'] != True].index.tolist(),'species','Custom Field','warning','"Unassigned species"',merge_labs_resultsData)

# creates a new column of boolean series to check if every item in column stationid_x is contained in column stationsid_y
merge_labs_resultsData['Check_station'] = merge_labs_resultsData.apply(lambda row: row['stationid_x'] in row['stationid_y'], axis=1)
checkData(merge_labs_resultsData.loc[merge_labs_resultsData['Check_station'] != True].index.tolist(),'species','Custom Field','warning','"Unassigned station"',merge_labs_resultsData)
