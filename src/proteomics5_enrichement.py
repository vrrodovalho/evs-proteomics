# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho
"""

import os
import sys
import re
import pathlib
import pandas as pd
import numpy as np
import KEGG as kegg
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from tabulate import tabulate

'''
    
'''

def update_annotation(df, column, replace_dict):
    '''
    Updates a dataframe column based on a dictionary.

    Parameters
    ----------
    df : DataFrame
        DataFrame that will be modified.
    column : str
        Name of the column that will be modified.
    replace_dict : dict
        Dictionary whose keys will be replaced by its values in the selected
        column in df.

    Returns
    -------
    df : DataFrame
        The updated DataFrame.

    '''
    df = df.copy()
    df[column] = df[column].replace(replace_dict, regex=True)
    return df

def export_gmt(df, cat_col='KEGG_Pathway', cat_sep=',', genes_col='gi',
                description_map={}, replace_dict={}, cat_fill_na='?', 
                ref_org='', drop_unknowns=True, filter_by_size=False, 
                size_limit=(2,150), forbidden_set_prefixes=['map'], 
                output_dir='', file_name=''):
    '''
    Converts a df mapping gene -> categories to a df mapping category -> genes
    And creates a GMT file for using with gProfiler

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    cat_col : str, optional
        Name of the column with the categories annotation (KEGG, COG...). 
        The default is 'KEGG_Pathway'.
    cat_sep : str, optional
        The delimiter that separates multiple categories in a row. 
        The default is ','.
    genes_col : str, optional
        Name of the column with the genes. The default is 'gi'.
    description_map : dict, optional
        A dictionary that gives a description to each category. 
        That could be COG letters as keys and their meaning as values.
        The default is {}.
    replace_dict : dict, optional
        A dictionary to replace row values, useful to take out obsolete
        annotation. The default is {}.
    cat_fill_na : str, optional
        A value to fill missing values. The default is '?'.
    ref_org : str, optional
        A kegg organism code to be used as a reference for orthologues
        filtering. The default is ''.  
    drop_unknowns : bool, optional
        Whether to drop the functional category defined as unknown, that 
        previously had na values. The default is False.
    filter_by_size : bool, optional
        Whether to filter functional categories by min/max size.
        The default is False.
    size_limit : tuple, optional
        A tuple containing 2 integers, a min and a max size for the sets of
        functional categories. The default is (2,150).     
    forbidden_set_prefixes : list, optional
        If some gene sets are forbidden, they can be identified in a prefix
        list to be removed from the dataset. The default is ['map'].
    output_dir : str, optional
        Output directory. The default is ''.
    file_name : str, optional
        Output file name. The default is ''.

    Returns
    -------
    sub_df : DataFrame
        A DataFrame close to the GMT file produced.

    '''
    # simplify df
    sub_df = df.loc[:, [genes_col, cat_col]].copy()
    
    # make needed replacements
    if replace_dict:
        sub_df = update_annotation(sub_df, column=cat_col, 
                                   replace_dict=replace_dict)
    
    # fill na as specified
    sub_df[cat_col].fillna(cat_fill_na, inplace=True)
    
    # devide rows with multiple annotation based on delimiter
    if cat_sep == '':
        sub_df = (sub_df.set_index([genes_col])
                        .stack()
                        .apply(lambda x: pd.Series(list(x)))
                        .stack()
                        .unstack(-2)
                        .reset_index(-1, drop=True)
                        .reset_index()
                        )
    else:
        sub_df = (sub_df.set_index([genes_col])
                        .stack()
                        .str.split(cat_sep, expand=True)
                        .stack()
                        .unstack(-2)
                        .reset_index(-1, drop=True)
                        .reset_index()
                        )
    sub_df = ( sub_df.groupby(by=cat_col)[genes_col]
                      .apply(set)
                      .reset_index()
                      )
    
    # # filter by set size, to eliminate sets too short or too long
    if filter_by_size:
        if size_limit:
            min_size = min(size_limit)
            max_size = max(size_limit)
            sub_df['size'] = sub_df[genes_col].apply(len)
            sub_df = sub_df.sort_values(by=['size'], ascending=False)
            sub_df = sub_df.loc[ ( ( sub_df['size'] > min_size ) & \
                                   ( sub_df['size'] < max_size ) ),
                                [cat_col, genes_col]]
        else:
            str1 = "If filter_by_size is True, size_limit should be defined "
            str2 = "as a tuple containing 2 int: a min and a max limit."
            print(str1 + str2)
            
    sub_df = sub_df.set_index(cat_col)
    
    # take out unknown category (privously na)
    if drop_unknowns:
        sub_df = sub_df.drop([cat_fill_na])
            
    # take out forbidden gene sets 
    if forbidden_set_prefixes:
        for i in forbidden_set_prefixes:
            sub_df =  sub_df[~sub_df.index.str.startswith(i)]
            

    # Use a KEGG reference organism to drop unrelated pathways
    if ref_org:
        allowed_ids = search_allowed_pathways_ids(ref_org)
        sub_df =  sub_df[sub_df.index.isin(allowed_ids)]

    # change 1-column set to several columns and name them accordingly
    f = lambda x: 'element_{}'.format(x + 1)
    sub_df = pd.DataFrame(sub_df[genes_col].values.tolist(),
                            sub_df.index, dtype=object
                          ).rename(columns=f)
    sub_df = sub_df.reset_index()
    
    # Add description to gene sets, if available
    if description_map:
        description_map[cat_fill_na] = 'Unknown' 
        sub_df['description'] = sub_df[cat_col].map(description_map)
    else:
        sub_df['description'] = np.nan
    # reorder description column, to be in GMT style
    cols = list(sub_df.columns)
    cols.remove('description')
    cols.insert(1, 'description')
    sub_df = sub_df.loc[:,cols]
    
    # # save and return
    output_file = output_dir / file_name
    sub_df.to_csv(output_file, header=False, index=False, sep='\t')
    return sub_df

def generate_gprofiler_list(df, id_column='', category_filter={}, 
                            ordered_by='',  output_dir='', file_name=''):
    '''
    Returns a list of genes to use in GProfiler.

    Parameters
    ----------
    df : DataFrame
        The initial DataFrame.
    id_column : str
        The name of the column that contains gene IDs.
    category_filter : dict, optional
        A dictionary in which keys are column names (str) and values are
        allowed rows in that column (str). The default is {}.
    ordered_by : str, optional
        The name of the column that will be used to sort gene list.
        It could be a expression measure. The default is ''.
    output_dir : str
        Output directory.
    file_name : str
        Output file name.

    Returns
    -------
    string : str
         A list of genes to be used in GProfiler.

    '''
    df = df.copy()
    
    # Filter gene list by column values (such as a category)
    if category_filter:
        for col in category_filter:
            value = category_filter[col]
            df = df.loc[df[col] ==  value, :]
    
    # Sort gene list by column values (such as expression)
    if ordered_by:
        df[ordered_by] = df[ordered_by].abs()
        df = df.sort_values(by=ordered_by, ascending=False)
        
        min_value = df.iloc[0,  df.columns.get_loc(ordered_by)]
        max_value = df.iloc[-1, df.columns.get_loc(ordered_by)]
        string = "Ordered in {}, from {} to {}. ".format(ordered_by, 
                                                          min_value, max_value)
        print(string)
    
    # Make final string and files
    proteins = df.loc[:, id_column].astype(str).to_list()
    string = '\n'.join(proteins)
    output_file = output_dir / file_name
    with open(output_file, 'w') as output:
        output.write(string)
    
    return string

def merge_enrichment_sources(source_files={'name': 'dataframe'}, 
                             max_p_val=0.05, v_limit=6, replace_values={},
                             output_dir='', file_name=''):
    '''
    Merge enrichment results from different sources (KEGG, COG) into the
    same dataframe, corresponding to the same set of proteins.

    Parameters
    ----------
    source_files : dict
        A dictionary where the keys are string identifiers (KEGG, COG) and the
        values are the dataframes corresponding to the enrichment results
        corresponding to those strings.
    max_p_val : float, optional
        The p-value threshold of significance. The default is 0.05.
    v_limit : float, optional
        Vertical superior limit of log(p-value). Values exceeding that 
        threshold are capped to it. The default is 6.
    replace_values : dict, optional
        A dictionary where the keys are column names and the values are
        replacement dictionaries, containing key-value pairs for replacing
        values in that column. The default is {}.
    output_dir : str, optional
        Output directory. The default is ''.
    file_name : str, optional
        Output file name. The default is ''.
        
    Returns
    -------
    df : DataFrame
        A merged DataFrame.

    '''
    
    df = pd.DataFrame()
    for item_name in source_files:
        item = source_files[item_name]
        item['source'] = item_name
        df = pd.concat([df, item])
    
    df['log_p_value'] =  np.log10(df['adjusted_p_value']) * -1
    df['sig'] = np.where(df['adjusted_p_value'] <= max_p_val, 'sig.', 'not sig.')
    
    df = df.sort_values(by=['source', 'log_p_value'], ascending=False)
    df['log_p_value_capped'] = np.where(df['log_p_value'] >= v_limit, 
                            v_limit, df['log_p_value'])
    if replace_values:
        for col in replace_values:
            replace_dict = replace_values[col]
            df[col] = df[col].replace(replace_dict, regex=True)
    
    # save file
    df.to_excel(output_dir/file_name, index=False)
    
    return df

def plot_enrichment(df, data = {'x':'source', 'y':'log_p_value_capped', 
                   'label_col':'term_id', 'label_desc_col':'term_name'}, 
                    v_limit=6, max_p_val=0.05, 
                    significancy={'column':'sig','true':'sig.','false':'not sig.'},
                    jitter_val=0.3, s=4, reg_categories= {'column': 'sig', 
                    'true':'up', 'false':'down', 'true_color':'blue', 
                    'false_color':'red'}, title='Functional enrichment', 
                    save_fig=True,output_dir='',file_name='',file_format='tif', 
                    dpi=300):
    '''
    Plot enrichment 

    Parameters
    ----------
    df : DataFrame
        A dataframe containing the data to be plotted. Ideally generated by
        merge_enrichment_sources function.
    data : dict, optional
        A dictionary specifying column names in df for x, y and label values. 
        The default is {'x':'source', 'y':'log_p_value_capped', 
                        'label_col':'term_id', 'label_desc_col':'term_name'}.
    max_p_val : float, optional
        The p-value threshold of significance. The default is 0.05.
    v_limit : float, optional
        Vertical superior limit of log(p-value). Values exceeding that 
        threshold are capped to it. The default is 6.
    significancy : dict, optional
        A dictionary specifying which is the significancy column and what 
        values should be considered True and False. 
        The default is {'column':'sig','true':'sig.','false':'not sig.'}.
    jitter_val : float, optional
        Parameter for stripplot. Affects the points distribution.
        The default is 0.3.
    s : float, optional
        The size of the points in the graph. The default is 4.
    reg_categories : dict, optional
        A dictionary specifying regulation categories (up-regulated, 
        down-regulated), the column, their values and colors. 
        The default is {'column':'sig', 'true':'up', 'false':'down',
                        'true_color':'blue', 'false_color':'red'}.
    title : str, optional
        A title string to be plotted in the graph. 
        The default is 'Functional enrichment'.
    save_fig : bool, optional
        Wether to save figure or not. The default is True.
    output_dir : str, optional
        Output directory. The default is ''.
    file_name : str, optional
        Output file name. The default is ''.
    file_format : str, optional
        File format. The default is 'tif'.
    dpi : int, optional
        Resolution. The default is 300.

    Returns
    -------
    dict
        A dictionary containing the final DataFrame and a legend string.

    '''
    
    df = df.copy()
    fig = plt.figure()
    ax = plt.axes() 
    
    sub_df_sig = df.loc[ df[significancy['column']] == significancy['true'] ] 
    sub_df_not = df.loc[ df[significancy['column']] == significancy['false'] ]
    x = data['x']
    y = data['y']
    
    commons = {'ax':ax,'x':x,'y':y,'size':s,'marker':'s','jitter':jitter_val}
   
    # plot not significtives
    sns.stripplot(data=sub_df_not, linewidth=0.1, alpha=0.5,  color='grey', 
                  **commons)
    
    # plot significatives
    if reg_categories:
        palette = {reg_categories['true']:reg_categories['true_color'],
                   reg_categories['false']:reg_categories['false_color']}
        sns.stripplot(data=sub_df_sig,linewidth=0.5,alpha=1.0,palette=palette,
                      hue=reg_categories['column'],dodge=True, **commons)
    else:
        sns.stripplot(data=sub_df_sig,linewidth=0.5,alpha=1.0,color='blue',
                       **commons)
    
    # title?
    if title != '':
        plt.title(title, loc='center')
    
    # plot lines
    ax.set(ylim=(-0.2, v_limit+1))
    log_max_p_val = np.log10(max_p_val) * -1
    plt.axhline(y=log_max_p_val , color='grey', linewidth=0.5, linestyle='--')
    plt.axhline(y=v_limit , color='grey', linewidth=0.5, linestyle='--')
    
    # plot labels
    plt.xlabel('', fontsize=12, fontname="sans-serif")
    plt.ylabel('Statistical significance [-log10(P-value)]', fontsize=12, 
                                       fontname="sans-serif")
    
    # create a df with x-y coordinates only for significatives
    df_graph = pd.DataFrame({'x' : [], y : []})
    for i in range(len(ax.collections)):
        coll = ax.collections[i]
        x_values, y_values = np.array(coll.get_offsets()).T
        
        # look for significative y
        annotate = False
        for i in y_values:
            if i >= log_max_p_val:
                annotate = True
                break
        # if found significative y, add to df that will be used to annotate
        if annotate:
            sub_df = pd.DataFrame({'x':x_values, y:y_values})
            df_graph = pd.concat([df_graph, sub_df])
    
    # transfer id col to df_graph in order to have unique identifiers
    # and avoid label confusion
    unique_id = data['label_col']
    unique_id_desc = data['label_desc_col']
    df_graph[unique_id] = sub_df_sig[unique_id]

    # anottate significative y
    merged = sub_df_sig.merge(df_graph, on=[y, unique_id], how='left')  
    sig_x = merged['x']
    sig_y = merged[y]
    labels = merged[unique_id]
    coordinates = []
    for xi, yi, label in zip(sig_x, sig_y, labels):
        element = ax.annotate(label, xy=(xi,yi), xytext=(3,3), size=8,
                ha="center", va="top", textcoords="offset points")
        coordinates.append(element)
    
    # ajust labels to avoid superposition
    adjust_text(coordinates, autoalign='xy',
        arrowprops=dict(arrowstyle='<-, head_length=0.05, head_width=0.05',
                        color='black', alpha=0.6,  linewidth=0.5))
    plt.show()
    
    # return a legend string and file
    legend_df = sub_df_sig.loc[:,[unique_id, unique_id_desc]]
    legend = tabulate(legend_df, showindex=False)
    legend_file_name = '.'.join(file_name.split('.')[:-1]) + '.txt'
    output_legend = output_dir / legend_file_name
    with open(output_legend, 'w') as output:
        output.write(legend)
    
    # save
    if save_fig:
        fig.savefig(output_dir/file_name, format=file_format, dpi=dpi, 
                    bbox_inches="tight")

    return {'sub_df_sig':sub_df_sig, 'df_graph':df_graph,
            'df':merged, 'legend':legend}


def search_allowed_pathways_ids(ref_org, unknown='?'):
    '''
    Search in KEGG all the pathways ids for an organism

    Parameters
    ----------
    ref_org : str
        KEGG organism code.

    Returns
    -------
    allowed_ids : list
        List of allowed ids (with ko prefix).

    '''
    
    kegg_data = kegg.get_KEGG_data(org=ref_org, get_pathway_list=True, 
                get_genes_and_pathways=False,  format_conversion=False, 
                genes_names=False)
    org_pathways = kegg.parse_KEGG_pathways_description(kegg_data['pathways'])
    allowed_ref_ids = list(org_pathways.keys())
    allowed_ids = []
    p = '[a-z]+([0-9]+)'
    for ref_id in allowed_ref_ids:
        general_id = re.match(p,ref_id).groups()[0]
        general_id = 'ko' + general_id
        allowed_ids.append(general_id)
    allowed_ids.append(unknown)
    return allowed_ids

def export_tables(proteomics_df=None, proteomics_id_col='', enrichment_df=None, 
                  enrichment_id_col='', enrichment_src_col='', merge_all=False,
                  enrichment_desc_col='', split_ch=',', enrichment_filter={}, 
                  map_src2annot={}, output_dir='', file_name_prefix=''):
    '''
    Function to export merge proteomics annotation and functional enrichment
    table and filter based on specific rules.

    Parameters
    ----------
    proteomics_df : DataFrame
        A DataFrame containing proteomics annotation.
    proteomics_id_col : str
        The name of the column in proteomics_df where the protein ids are.
    enrichment_df : DataFrame
        A DataFrame containing enrichment results for proteins in proteomics_df.
    enrichment_id_col : str
        The name of the column where the functional category ids are specified.
    enrichment_src_col : str
        The name of the column where the source database is specified.
    enrichment_desc_col : str
        The name of the column where the description of id is specified.
    split_ch : str
        A character to split a string into a list of items in 
        enrichment_id_set_col. The default is ','.
    merge_all : bool
        Whether to merge all enriched categories elements in one single
        dataframe. Otherwise, they will be returned separated by category
        in a dictionary. The default is 'False'.
    enrichment_filter : dict, optional
        A dictionary describing a filter for enrichment_df the format 
        { col_name : [allowed_values] }. Only rows fulfilling these rules are
        accepted.
    map_src2annot : dict
        A dictionary with the relationship between 
        { col_name : [allowed_values] }. Only rows fulfilling these rules are
        accepted.
    output_dir : str, 
        The output directory.
    file_name_prefix : str
        A prefix for every output file name. The default is ''.

    Returns
    -------
    None.

    '''
    prot = proteomics_df.copy()
    enri = enrichment_df.copy()
    
    # get descritions
    desc = dict(zip( enri[enrichment_id_col], enri[enrichment_desc_col]))
    
    # filter enrichment data (significative)
    if enrichment_filter:
        for col in enrichment_filter:
            col_values = enrichment_filter[col]
            enri = enri.loc[enri[col].isin(col_values) ,:]
    
    # get dictionary of enriched categories by enrichment source
    enri_items = enri.loc[:,[enrichment_src_col, enrichment_id_col]]
    enri_items = ( enri_items.groupby(enrichment_src_col)[enrichment_id_col]
                  .apply(set).to_dict() )
    
    # search items in proteomics_df that correspond to enriched categories
    enriched_elements = {}
    appended_data = []
    prot = prot.fillna('?')
    for src in enri_items:
        where2look = map_src2annot[src]
        cats = enri_items[src]
        for cat in cats:
            description = desc[cat]
            sub_prot = prot.loc[prot[where2look].str.contains(cat) ,:]
            n_prot = sub_prot.shape[0]
            appended_data.append(sub_prot)
            print("{} \t{} \t(n={}) \t{}".format(src, cat, n_prot, 
                                                       description))
            
            enriched_elements[cat + ' : ' + description] = sub_prot
            file_name = '{}_{}_{}_{}.xlsx'.format(file_name_prefix, src, cat, 
                                                  description)
            sub_prot.astype(str).to_excel(output_dir / file_name, index=False)
    

    single_df = pd.concat(appended_data)
    single_df = single_df.drop_duplicates()
    file_name = '{}_merged.xlsx'.format(file_name_prefix)
    single_df.astype(str).to_excel(output_dir / file_name, index=False)

    # merge all enriched categories elements
    if merge_all:
        enriched_elements = single_df
    
    return enriched_elements
    

##############################################################################
    
# DIRECTORY SYSTEM
src_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
main_dir = os.path.dirname(src_dir)
root_dir = os.path.dirname(main_dir)
data_dir = pathlib.Path(main_dir) / 'data'
input_dir = pathlib.Path(data_dir) / 'input'
output_dir = pathlib.Path(data_dir) / 'output'
sys.path.insert(0, root_dir)

# FILE PATHS
proteomics_SEC_and_UC_file = input_dir / 'proteomics_SEC_and_UC_curated.xlsx'
proteomics_UC_file = input_dir / 'proteomics_UC.xlsx'
proteomics_core_file = input_dir / 'proteome_core.xlsx' 
proteomics_accessory_file = input_dir / 'proteome_accessory.xlsx'
proteomics_single_file = input_dir / 'proteome_single.xlsx'
proteomics_not_EVs_file = input_dir / 'proteome_not_EVs.xlsx'
cogs_file = input_dir / 'COGs.xlsx'
kegg_ko_storage_file = input_dir / 'kegg_ko.data'
gprofiler_core_kegg_file = input_dir / 'gProfiler_core_kegg.csv'
gprofiler_core_cog_file = input_dir / 'gProfiler_core_cog.csv'
gprofiler_accessory_kegg_file = input_dir / 'gProfiler_accessory_kegg.csv'
gprofiler_accessory_cog_file = input_dir / 'gProfiler_accessory_cog.csv'
gprofiler_single_kegg_file = input_dir / 'gProfiler_single_kegg.csv'
gprofiler_single_cog_file = input_dir / 'gProfiler_single_cog.csv'


# READ FILES
proteomics_SEC_and_UC = pd.read_excel(proteomics_SEC_and_UC_file)
proteomics_UC = pd.read_excel(proteomics_UC_file)
proteomics_core = pd.read_excel(proteomics_core_file)
proteomics_accessory = pd.read_excel(proteomics_accessory_file)
proteomics_single = pd.read_excel(proteomics_single_file)
proteomics_not_EVs = pd.read_excel(proteomics_not_EVs_file)
cogs_map = pd.read_excel(cogs_file)
kegg_data_ko = kegg.get_all_KEGG_data(kegg_ko_storage_file, org='ko')
gprofiler_core_kegg = pd.read_csv(gprofiler_core_kegg_file)
gprofiler_core_cog = pd.read_csv(gprofiler_core_cog_file)
gprofiler_accessory_kegg = pd.read_csv(gprofiler_accessory_kegg_file)
gprofiler_accessory_cog = pd.read_csv(gprofiler_accessory_cog_file)
gprofiler_single_kegg = pd.read_csv(gprofiler_single_kegg_file)
gprofiler_single_cog = pd.read_csv(gprofiler_single_cog_file)


# Common parameters
cog_names = dict(zip(cogs_map.l2, cogs_map.l2_description))
kegg_ko_pathways = kegg_data_ko['pathways_descriptions']
gmt_args = {'output_dir': output_dir, 'cat_fill_na':'?', 'genes_col':'gi',
            'replace_dict':{}, 'forbidden_set_prefixes':[],
            'drop_unknowns':False}


# EXPORT GMT FILES

# gmt_KEGG = export_gmt(df=proteomics_SEC_and_UC, 
#                       cat_col='KEGG_Pathway', ref_org='pfr', cat_sep=',', 
#                       description_map=kegg_ko_pathways, 
#                       file_name='gmt_KEGG.gmt', **gmt_args)
# gmt_COG = export_gmt(df=proteomics_SEC_and_UC, cat_col='COG', cat_sep='', 
#                      description_map=cog_names, 
#                      file_name='gmt_COG.gmt', **gmt_args)


# GENERATE LIST OF GIs TO USE IN GPROFILER
kargs = {'id_column':'gi', 'category_filter':{}, 'ordered_by':'',
          'output_dir': output_dir}
gprofiler_list_CORE = generate_gprofiler_list(df=proteomics_core, 
                            file_name='gprofiler_list_CORE.txt', **kargs)
# gprofiler_list_ACCESSORY = generate_gprofiler_list(df=proteomics_accessory, 
#                             file_name='gprofiler_list_ACCESSORY.txt', **kargs)
# gprofiler_list_SINGLE = generate_gprofiler_list(df=proteomics_single, 
#                             file_name='gprofiler_list_SINGLE.txt', **kargs)


# Common merge/plot parameters
merge_kargs = {'max_p_val':0.05, 'v_limit':6, 'output_dir':output_dir,
               'replace_values':{'term_id': {'ko':''}}}
plot_kargs = { 'v_limit':6, 'max_p_val':0.05, 'jitter_val':0.3, 's':4,
    'data': {'x':'source', 'y':'log_p_value_capped',
              'label_col':'term_id', 'label_desc_col':'term_name'},
    'significancy':{'column':'sig','true':'sig.','false':'not sig.'},
      'save_fig':True,'output_dir':output_dir, 'file_format':'tif', 'dpi':600,
      'title':'', 'reg_categories':{}
    }


# PLOT ENRICHMENT FOR CORE PROTEOME
source_files={'KEGG': gprofiler_core_kegg, 'COG' : gprofiler_core_cog}
enrichment_core = merge_enrichment_sources(source_files, 
                                           file_name='enrichment_core.xlsx', 
                                           **merge_kargs)
# plot_core = plot_enrichment(enrichment_core, 
#                             file_name='enrichment_core.tif', 
#                             **plot_kargs) 
                    

# PLOT ENRICHMENT FOR ACCESSORY PROTEOME
# source_files={'KEGG': gprofiler_accessory_kegg, 'COG' : gprofiler_accessory_cog}
# enrichment_accessory = merge_enrichment_sources(source_files, 
#                                          file_name='accessory_enrichment.xlsx',
#                                          **merge_kargs)
# plot_enrichment_accessory = plot_enrichment(enrichment_accessory, 
#                                         file_name='accessory_enrichment.tif',
#                                         **plot_kargs)

# PLOT ENRICHMENT FOR SINGLE PROTEOME
# source_files={'KEGG': gprofiler_single_kegg, 'COG' : gprofiler_single_cog}
# enrichment_single = merge_enrichment_sources(source_files, 
#                                          file_name='single_enrichment.xlsx',
#                                          **merge_kargs)
# plot_enrichment_single = plot_enrichment(enrichment_single, 
#                                          file_name='single_enrichment.tif',
#                                          **plot_kargs)

sig_tables_cor = export_tables(proteomics_df=proteomics_core, 
                         proteomics_id_col='gi', 
                         enrichment_df=enrichment_core, 
                         enrichment_id_col='term_id',
                         enrichment_desc_col='term_name',
                         enrichment_src_col='source', 
                         enrichment_filter={'sig':['sig.']},
                         map_src2annot = {'KEGG':'KEGG_Pathway', 'COG':'COG'},
                         output_dir=output_dir, file_name_prefix='core')