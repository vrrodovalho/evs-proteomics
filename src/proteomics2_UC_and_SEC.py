# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho

These are the scripts for aggregation of data from qualitative proteomics for 
EVS purified by ultracentrifugation and chromatography.
It also verifies if there are any duplicated rows in ID column.
"""

import os
import sys
import pathlib
import pandas as pd
import shared_modules as md

def aggregate_proteomics_data(data_set={'df1':'','df2':''}, common_columns=[]):
    '''
    This function aggregates proteomics data from 2 subsets (e.g., two methods
    for EV purification).

    Parameters
    ----------
    data_set : dict
        Dictionary in the form {'df1_name' : df1, 'df2_name' : df2}
        containing 2 DataFrames that will be merged.
    common_columns : list
        List of strings, representing common columns whose content is equal 
        in both dataframes.
    Returns
    -------
    df : DataFrame
        A DataFrame containing data from the 2 subsets

    '''
    dataset_groups = list(data_set.keys())
    (df1, df2) = (data_set[i] for i in dataset_groups)
    df = df1.merge(df2, on=common_columns, how='outer', suffixes=dataset_groups)

    return df

def check_duplicates(df, column="query_name"):
    '''
    This fucntion checks if there are duplicated rows in a column of a 
    DataFrame.

    Parameters
    ----------
    df : DataFrame
        The Dataframe containing input data.
    column : str
        The name of the column that will be verified. 
        The default is "query_name".

    Returns
    -------
    dup : DataFrame
        A sub-DataFrame with duplicated rows.

    '''
    dup = df.loc[df[column].duplicated(),:]
    if not dup.empty:
        print("In column {} got {} duplicated rows: {}".format(
                column, dup.shape[0], dup.index))
    else:
        print("In column {}, got no duplicated rows".format(column))
    return dup

##############################################################################
    
# DIRECTORY SYSTEM
src_dir = 'C:\\Users\\venic\\Dropbox\\Scripts\\doutorado-proteomics-compare\\src'
src_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
main_dir = os.path.dirname(src_dir)
root_dir = os.path.dirname(main_dir)
data_dir = pathlib.Path(main_dir) / 'data'
input_dir = pathlib.Path(data_dir) / 'input'
output_dir = pathlib.Path(data_dir) / 'output'
sys.path.insert(0, root_dir)

# FILE PATHS
proteomics_SEC_file = input_dir / 'proteomics_SEC.xlsx'
proteomics_UC_file = input_dir / 'proteomics_UC.xlsx'

# READ FILES
proteomics_SEC_df = pd.read_excel(proteomics_SEC_file)
proteomics_UC_df = pd.read_excel(proteomics_UC_file)

# CHECK FILES
print("\nData for proteomics_SEC:")
md.print_basic_info(proteomics_SEC_df, columns=False, size=True, empties=True,
                      duplicates=False, check_columns=["query_name"])
print("\nData for proteomics_UC:")
md.print_basic_info(proteomics_UC_df, columns=False, size=True, empties=True,
                      duplicates=False, check_columns=["query_name"])

# PERFORM MERGE
data_set = {'_SEC' : proteomics_SEC_df,
            '_UC'  : proteomics_UC_df}
common_columns = ['query_name', 'gi', 'emb', 'description', 'Sequence', 
                  'Gene names', 'Protein names', 'locus_tag', 'Entry', 
                  'UniProtKB',  'Predicted protein name', 'eggNOG_description',
                  'GO_terms', 'EC_number', 'KEGG_KO', 'KEGG_Pathway', 
                  'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'KEGG_TC', 
                  'BRITE', 'BiGG_Reaction', 'COG', 'cello2go', 'predlipo', 
                  'CAZy']

proteomics_df = aggregate_proteomics_data(data_set=data_set, 
                                          common_columns=common_columns)

print("\nData for proteomics_SEC_and_UC:")
md.print_basic_info(proteomics_df, columns=False, size=True, empties=False,
                      duplicates=False, check_columns=["query_name"])
duplicates_df = check_duplicates(proteomics_df, column="query_name")

proteomics_df.to_excel(output_dir / 'proteomics_SEC_and_UC.xlsx', index=False)
duplicates_df.to_excel(output_dir / 'proteomics_SEC_and_UC_duplicates.xlsx', index=False)
