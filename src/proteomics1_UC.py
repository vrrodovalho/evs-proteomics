# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho

These are the scripts for processing and plotting data from
qualitative proteomics for EVS purified by ultracentrifugation.

"""

import os
import sys
import pathlib
import additional_info
import shared_modules as md
import numpy as np
import pandas as pd

def quantitative_proteomics(raw_df, info_df, replicates, fractions, media,
                            limits_factor=1.1, use=['XIC', 'PC', 'both'][2]):
    '''
    '''
    df = raw_df.copy()
    df['gi'] = df['accession'].str.extract(pat='gi\|(\d*)\|emb\S*')
    
    info_df = info_df.copy()
    info_df['gi'] = info_df['gi'].apply(str)
    
    mapping1 = dict(zip(info_df['gi'], info_df['locus_tag']))
    mapping2 = dict(zip(info_df['gi'], info_df['query_name']))
    mapping3 = dict(zip(info_df['gi'], info_df['description']))

    df['locus_tag'] = df['gi'].map(mapping1)
    df['query_name'] = df['gi'].map(mapping2)
    df['description'] = df['gi'].map(mapping3)
    
    new_mapping = {v : k for k, v in replicates.items()}
    fractions_x = {k + '.x' : [j + '.x' for j in v] 
                              for k, v in fractions.items()} 
    fractions_y = {k + '.y' : [j + '.y' for j in v] 
                              for k, v in fractions.items()}
    fractions = {**fractions_x, **fractions_y}
    
    media_x = {k + '.x' : [j + '.x' for j in v] 
                              for k, v in media.items()} 
    media_y = {k + '.y' : [j + '.y' for j in v] 
                              for k, v in media.items()}
    media = {**media_x, **media_y}
    
    df.columns = df.columns.to_series().replace(new_mapping, regex=True)
    
    # SELECT ACCORDING TO METHOD (XIC, PC, BOTH)
    if use == 'both':
        df['fold-change'] = df['Ratio_XIC$UF.YEL']
        df['fold-change'].fillna(df['Ratio_PC$UF.YEL'], inplace=True)
        df['p-value'] = df['pval_milieu.x']
        df['p-value'].fillna(df['pval_milieu.y'], inplace=True)
        df['p-value-adj'] = df['padj_milieu.x']
        df['p-value-adj'].fillna(df['padj_milieu.y'], inplace=True)
    elif use == 'XIC':
        df['fold-change'] = df['Ratio_XIC$UF.YEL']
        df['p-value'] = df['pval_milieu.x']
        df['p-value-adj'] = df['padj_milieu.x']
    elif use == 'PC':
        df['fold-change'] = df['Ratio_PC$UF.YEL']
        df['p-value'] = df['pval_milieu.y']
        df['p-value-adj'] = df['padj_milieu.y']
        
    # DEALING WITH ZERO VALUES FOR FC
    min_fc = df.loc[ ((df['fold-change'] != -np.inf) & 
                      (df['fold-change'] != 0)), 'fold-change'].min()
    min_fc_zeros = min_fc / limits_factor
    df.loc[:,'log-fold-change'] = np.log2( (df.loc[:,'fold-change']
                                              .replace(0, min_fc_zeros)))
    
    # DEALING WITH INFINITE VALUES FOR FC
    max_fc_log = df.loc[ df['log-fold-change'] !=  np.inf, 
                        'log-fold-change'].max()
    min_fc_log = df.loc[ df['log-fold-change'] != -np.inf, 
                        'log-fold-change'].min()
    df['log-fold-change'].replace( 
            { np.inf : max_fc_log * limits_factor , 
             -np.inf : min_fc_log * limits_factor }, inplace = True)
    
    # LOG OF P-VALUES
    df['log-p-value'] = np.log2(df['p-value'])
    df['log-p-value-adj'] = np.log2(df['p-value-adj'])
    return df

def quanti_subgroups(df, columns = {'FC':'log-fold-change', 'PV':'p-value-adj'},
                         cutoffs = {'p-value':0.05, 'logFC':1}):
    '''
    '''
    df = df.copy()
    FC = columns['FC']
    PV = columns['PV']
    up =   df.loc[ (df[PV] <= cutoffs['p-value']) & 
                            (df[FC] >= +(cutoffs['logFC'])) , : ]
    down = df.loc[ (df[PV] <= cutoffs['p-value']) & 
                            (df[FC] <= -(cutoffs['logFC'])) , : ]
    ns =  df.loc[ (df[PV] > cutoffs['p-value']) | 
                           ((df[FC] < +(cutoffs['logFC'])) &
                            (df[FC] > -(cutoffs['logFC']))), : ]
    up_proteins = up.loc[:,'query_name'].to_list()
    down_proteins = down.loc[:,'query_name'].to_list()
    ns_proteins = ns.loc[:,'query_name'].to_list()
    
    df['subgroup'] = np.where(df['query_name'].isin(up_proteins), 'up', 
                     np.where(df['query_name'].isin(down_proteins), 'down',
                     np.where(df['query_name'].isin(ns_proteins), 'not', '0')))
    
    total_i = len(df[FC].dropna())
    total_f = up.shape[0] + \
              down.shape[0] + \
              ns.shape[0]
    if total_i == total_f:
        print("\nQuantitative groups:\n")
        print("Got {} proteins divided into 3 subgroups".format(total_f))
        print("Up-regulated: {}".format(len(up_proteins)))
        print("Down-regulated: {}".format(len(down_proteins)))
        print("Not significative: {}".format(len(ns_proteins)))
        print("Cutoffs -> p-value: {}, logFC: {}".format(cutoffs['p-value'],
                                                         cutoffs['logFC']))
    else:
        print('Got a final sum of {}, expecting {}'.format(total_f, total_i))
    return df


def helper_creater_dict(media):
    '''
    '''
    media_dict = {}
    old_media_dict  = media.to_dict()
    for tuple_key in old_media_dict:
        values = old_media_dict[tuple_key]
        new_key = '_'.join(tuple_key)
        new_values = [new_key + '_' + x for x in values]
        media_dict[new_key] = new_values
    return media_dict

def extract_lists_from_dictionary(dictionary, l=2):
    '''
    2 levels of keys
    '''
    final_list = []
    if l == 1:
        for key1 in dictionary:
            value = dictionary[key1]
            final_list.extend(value)
    elif l == 2:
        for key1 in dictionary:
            for key2 in dictionary[key1]:
                value = dictionary[key1][key2]
                final_list.extend(value)
    return final_list

def quantitative2qualitative(df, replicats='', n_pep=5, pct=0.5, n_frac=1, 
                             output_dir='', print_info=False):
    '''
    '''
    quali = df.copy()
    cols = ['meth','med','frac','rep']
    fractions = pd.read_excel(replicats, header=None, names=cols)
    
    methods = fractions.groupby(['meth'])['rep'].apply(list)
    media = fractions.groupby(['meth','med'])['frac'].apply(set)  
    fractions = fractions.groupby(['meth','med','frac'])['rep'].apply(list)
    
    # for PC, if n <= n_pep, assign 0, otherwise 1
    PC_cols = methods['PC']
    quali.loc[:,PC_cols] = quali.loc[:,PC_cols].apply(
                                       lambda x: np.where(x < n_pep, 0, 1))
    
    # for XIC, if n is not null, assign 1
    XIC_cols = methods['XIC']
    quali.loc[:,XIC_cols] = quali.loc[:,XIC_cols].apply(
                                       lambda x: np.where(x.notnull(), 1, 0))
    
    # Agregate replicates per fraction in percentage
    fractions_dict = fractions.to_dict()
    for fraction in fractions_dict:
        reps = fractions_dict[fraction]
        fraction = '_'.join(fraction)
        quali[fraction] = quali.loc[:,reps].sum(axis=1) / len(reps)

    # COnsidere only pct above cutoff
    media_dict = helper_creater_dict(media)
    media_lists = extract_lists_from_dictionary(media_dict, l=1)
    quali.loc[:,media_lists] = quali.loc[:,media_lists].apply(
                                       lambda x: np.where(x >= pct, 1, 0))
    
    # Agregate fractions per medium if present in at least n_frac fractions
    for medium in media_dict:
        fracs = media_dict[medium]
        quali[medium]  = np.where( quali[fracs].sum(axis=1) >= n_frac, 1, 0)
        
    # If XIC is 1, then assign 1 for final value. Otherwise, look at PC. 
    # If PC is 1, then assign 1 for final value. Otherwise, assign 0.
    quali['UF']  = np.where(quali['XIC_UF'] == 1, 1,
                          np.where(quali['PC_UF']  == 1, 1, 0 ) )
    quali['YEL'] = np.where(quali['XIC_YEL'] == 1, 1,
                          np.where(quali['PC_YEL']  == 1, 1, 0 ) )
    quali['UF_unique']   = np.where((quali.loc[:,'UF'] == 1) & 
                                     (quali.loc[:,'YEL'] == 0), 1, 0)
    quali['YEL_unique']  = np.where((quali.loc[:,'UF'] == 0) & 
                                     (quali.loc[:,'YEL'] == 1), 1, 0)
    quali['both_UF_YEL'] = np.where((quali.loc[:,'UF'] == 1) & 
                                     (quali.loc[:,'YEL'] == 1), 1, 0)
    quali['medium'] = np.where(quali['UF_unique'] == 1, 'UF', 
                       np.where(quali['YEL_unique'] == 1, 'YEL', 
                       np.where(quali['both_UF_YEL'] == 1, 'UF/YEL', 0)))
    
    # Print info
    info_string = "Not median - n_pep {} - pct_rep {} - "
    info_string  = info_string.format(n_pep, pct)
    if print_info:
        print('\n\n')

        print(info_string)
        print("YEL: {}, UF: {}".format(quali[quali['YEL'] == 1]['YEL'].count(),
                                       quali[quali['UF'] == 1]['UF'].count()))
        
        uf_only= sorted(quali.loc[quali['medium'] == 'UF', 'query_name'].to_list())
        yel_only = sorted(quali.loc[quali['medium'] == 'YEL', 'query_name'].to_list())
        both = sorted(quali.loc[quali['medium'] == 'UF/YEL', 'query_name'].to_list())
        print("UF exclusive  ({}): {}".format(len(uf_only), ', '.join(uf_only)))
        print("YEL exclusive ({}): {}".format(len(yel_only), ', '.join(yel_only)))
        print("Bothe media ({}): {}".format(len(both), ', '.join(both)))

    # clean dataframe and export
    file_name = info_string + '.xlsx'
    file_path = output_dir / file_name
    quali = quali[(quali[['UF','YEL']] != 0).any(axis=1)]
    cols2keep = ['gi', 'locus_tag', 'query_name', 'description', 'UF', 'YEL', 
                 'medium']
    quali = quali.loc[:,cols2keep]
    quali.to_excel(file_path)
    
    return quali

def update_kegg_annotation(df, column, replace_dict):
    '''
    
    '''
    df = df.copy()
    df[column] = df[column].replace(replace_dict, regex=True)
    return df

def map_column(df1_target, df2_source, col1_ref, col2_new):
    '''
    '''
    df1_target = df1_target.copy()
    mapping_dict = dict(zip(df2_source[col1_ref], df2_source[col2_new]))
    df1_target[col2_new] = df1_target[col1_ref].map(mapping_dict)
    return df1_target

def aggregate(df, quanti_evs_CR, quali_evs_CR, 
              cutoffs = {'p-value':0.05, 'logFC':0}, print_info=True):
    '''
    Aggregates all proteomics data in a single dataframe
    Also filters out quantitative data that isn't validated qualitatively
    '''
    df = df.copy()
    cols_from_quali = ['gi', 'medium']
    cols_from_quanti = ['gi', 'fold-change', 'p-value', 'p-value-adj',
                        'log-fold-change', 'log-p-value', 'log-p-value-adj']
    
    sub_quali  = quali_evs_CR.loc[:,cols_from_quali]
    sub_quanti = quanti_evs_CR.loc[:,cols_from_quanti]
    
    df['gi'] = df['gi'].apply(str)
    sub_quali['gi'] = sub_quali['gi'].apply(str)
    sub_quanti['gi'] = sub_quanti['gi'].apply(str)
    
    for column in cols_from_quali:
        if column != "gi":
            df = map_column(df1_target=df, df2_source=sub_quali, 
                            col1_ref='gi', col2_new=column)
    for column in cols_from_quanti:
        if column != "gi":
            df = map_column(df1_target=df, df2_source=sub_quanti, 
                            col1_ref='gi', col2_new=column)
    
    df['medium'] = df['medium'].fillna('notEVs')
    cols_from_quanti_reduced = [i for i in cols_from_quanti if i != 'gi']
    df.loc[df.medium == 'notEVs', cols_from_quanti_reduced] = np.nan
    
    df = df.sort_values(by='medium')
    
    df = quanti_subgroups(df, cutoffs=cutoffs, columns={'FC':'log-fold-change', 
                                                        'PV':'p-value-adj'})

    if print_info:
        print("\nQualitative:")
        print(df.medium.value_counts())
        print("\nQuantitative:")
        print("Total of Fold Change values: {} ".format(
                len(df["log-fold-change"].dropna())))
        #df3["log-fold-change"].dropna().hist().plot()
            
    return df



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
proteins_info = input_dir / 'PF129_proteins.xlsx'
mass_dir = input_dir  / 'mass'
qualitative = mass_dir / 'qualitative'
quantitative = mass_dir / 'quantitative/merged'
mass_map_file_UC = qualitative / '20181221_all_combine_ultracentrifugation.xlsx'

# READ FILES
df = pd.read_excel(proteins_info)
replicats = input_dir / 'replicats.xlsx'
mass_map_raw_UC = pd.read_excel(open(mass_map_file_UC, 'rb'), 
                                sheet_name='compar unique sequence')
mass_uc_t =  pd.read_csv(quantitative / '18-053_Resultats_sans_seuil.csv')

#############################################################################

# PROTEOMICS - QUANTITATIVE AND QUALITATIVE5

method = 'both' # use XIC, PC or both data
n_pep = 5       # min number of peptides  (necessary to identify a protein)
pct_rep = 0.60  # min pct of replicates   (out of 100% of replicates)
n_frac = 1      # min number of fractions (out of 3 fractions)
quanti_cutoffs = {'p-value': 0.05, 'logFC': 0}

replicates_UC = additional_info.replicates_UC
fractions_UC = additional_info.fractions_UC
media_UC = additional_info.media_UC
map_quanti_UC = additional_info.map_quanti_UC

quanti_evs_UC = quantitative_proteomics(mass_uc_t, df, map_quanti_UC, 
                                        fractions_UC, media_UC, use=method)
quali_evs_UC = quantitative2qualitative(quanti_evs_UC, replicats=replicats, 
                                          n_pep=5, pct=0.60, n_frac=1,
                                          output_dir = output_dir)

df2_all = aggregate(df, quanti_evs_UC, quali_evs_UC, cutoffs=quanti_cutoffs)
df2_all = update_kegg_annotation(df2_all, column='KEGG_Pathway', 
                                  replace_dict={'ko01130':'ko01110',
                                                r",map\d+":""})


df3_EVs  = df2_all.loc[df2_all.medium != "notEVs",:]
df3_EVs_up   = df3_EVs[df3_EVs.subgroup == 'up'  ]
df3_EVs_down = df3_EVs[df3_EVs.subgroup == 'down']
df3_EVs_not  = df3_EVs[df3_EVs.subgroup == 'not' ]

md.print_basic_info(df2_all, columns=False, size=True, empties=False,
                    duplicates=True, check_columns=["query_name","gi"])
md.print_basic_info(df3_EVs, columns=False, size=True, empties=True,
                      duplicates=True, check_columns=["query_name","gi"])

df2_all.to_excel(output_dir / 'proteomics_UC.xlsx', index=False)


