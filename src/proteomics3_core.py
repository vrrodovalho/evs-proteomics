# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho
"""

import os
import sys
import pathlib
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from venn import venn
import shared_modules as md

def filter_qualitative_data(df, id_col, qualitative_cols):
    '''
    FIlters only an ID column and the qualitative columns of a dataframe.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame containing all data to be filtered.
    id_col : STR
        Name of the column containing an unique ID.
    qualitative_cols : LIST
        List of names of columns containing qualitative data.

    Returns
    -------
    filtered : DataFrame
        Filtered Dataframe containing specified columns
    '''
    
    filtered = df.loc[:, [id_col] + qualitative_cols]
    return filtered

def filter_annotation_data(df, id_col, annotation_cols):
    '''
    FIlters only an ID column and the annotation columns of a dataframe.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame containing all data to be filtered.
    id_col : STR
        Name of the column containing an unique ID.
    annotation_cols : LIST
        List of names of columns containing annotation data.

    Returns
    -------
    filtered : DataFrame
        Filtered Dataframe containing specified columns
    '''
    
    filtered = df.loc[:, [id_col] + annotation_cols]
    return filtered

def create_binary_matrix(df, col_groups={}, med_groups={}, output_dir='',
                         file_name=''):
    '''
    Create a binary matrix of presence (0) / absence (1)
    based on groups of conditions. Fucntion designed to 
    2 conditions and 2 subconditions (the limiting line 
    is the where statement).

    Parameters
    ----------
    df : DataFrame
        A DataFrame containing conditions in string format ("UF", "UF/YEL").
    col_groups : DICT
        Dictionary specifying the relationship between columns in df and 
        contitions name that will appear in df2. 
        Ex: {'medium_SEC':['SEC_UF','SEC_YEL'],'medium_UC':['UC_UF','UC_YEL']}
    med_groups : DICT
        Dictionary specifying the relationship between groups of conditions
        that appear in df. The size of the list must be 2.
        Ex: {'UF' : ['UF', 'UF/YEL'], 'YEL' : ['YEL', 'UF/YEL']}

    Returns
    -------
    df : DataFrame
        Dataframe containing conditions columns and ) the presence (1) or
    aubsence (0) of each protein.

    '''
    
    df = df.copy()
    all_conditions = []
    
    # Iterate through groups and fill conditions columns with 0s and 1s
    print('\nCreating binary matrix')
    for group in col_groups:
        conditions = col_groups[group]
        all_conditions.extend(conditions)
        for condition in conditions:
            print('Processing group {}, condition {}'.format(
                                                          group, condition))
            for medium in med_groups:
                if medium in condition:
                    names = med_groups[medium]
                    df[condition] = np.where((df.loc[:,group] == names[0]) | 
                                             (df.loc[:,group] == names[1]), 
                                             1, 0)
    # Sort dataframe to have the core proteome in the beggining
    df2 = df.sort_values(by=all_conditions, ascending=[0 for i in range(
                                                        len(all_conditions))])
    # export to file
    df2.to_excel(output_dir / file_name, index=False)
    
    return df2

def bin_matrix2dict(df, true_condition=1, id_col='gi', include_cols=[]):
    '''
    Construct a dictionay of sets from a binary matrix

    Parameters
    ----------
    df : DataFrame
        Dataframe containing ID columns and SET columns.
    true_condition : STR / INT, optional
        The condition which will be considered true, in a way that the element
        belongs to the set. The true values will be replaced by ID values.
        The default is 1.
    id_col : STR, optional
        Column where the IDs are placed. They will be transferred to rows
        with the true condition and then the original column will be
        ignored. The default is 'gi'.
    include_cols : LIST, optional
        SET columns. Columns from which the sets will be created. 
        One set for each column.
        The default is [].

    Returns
    -------
    all_sets : Dictionary of lists
        A dictionary in the format {name : [set] } containing all the sets of
        elements constructed from the binnary matrix.
    '''
    
    df = df.copy()
    
    # replace true condition (1) by the ID string (gi)
    for condition in include_cols:
        df[condition] = np.where(df[condition] == true_condition, 
                                 df[id_col].astype(str), np.nan)
    
    # select cols that will make up the dictionaries
    sub_df = df.loc[:,include_cols]
    
    # create dict of sets
    all_sets = {i : set(sub_df.loc[:,i].dropna().to_list()) \
                for i in include_cols}
    
    return all_sets

def plot_venn(dictionary_of_elements, output_dir, file_name, file_format, dpi):
    '''
    Plot a venn diagram from a dictionary of sets.
    Uses: https://github.com/LankyCyril/pyvenn

    Parameters
    ----------
    dictionary_of_elements : DICT
        Dictionary of sets.
    output_dir : STR
        Output directory.
    file_name : STR
        Output file.
    file_format : TYPE
        Output file format.
    dpi : INT
        Resolution.

    Returns
    -------
    None.
    Plots a Venn Diagram and saves the figure in the specified directory.

    '''
    
    fig, ax = plt.subplots()
    cmaps = ["cool", "plasma", "viridis", "Set1"]
    # letters = iter(ascii_uppercase)
    # plt.title(title)
    
    #fmt = "{percentage:.1f}%"
    palette = sns.color_palette("colorblind", n_colors=4)
    fmt = "{size}"
    venn(dictionary_of_elements, fmt=fmt, cmap=cmaps[2], 
         fontsize=7, legend_loc="best", ax=ax)
    
    plt.draw()
    file_name += '.' + file_format
    plt.savefig(output_dir/file_name, format=file_format, dpi=dpi, 
                    bbox_inches="tight")
    return None

def plot_heatmap(bin_matrix, use_cols=[], id_col='', par_lw=0.5, 
                  palette_name="Blues", output_dir='', 
                  file_name='', file_format='tiff', dpi=600):
    '''
    Plot heatmap of presence / absence of the proteins in studied conditions

    Parameters
    ----------
    bin_matrix : DataFrame
        DESCRIPTION.
    id_col : STR, optional
        Column where the IDs are placed. The default is 'gi'.
    use_cols : LIST, optional
        List of columns to be included in the heatmap. 
        The default is [].
    palette_name : str, optional
        A Seaborn color palette name. The default is 'Blues'
    output_dir : STR
        Output directory.
    file_name : STR
        Output file.
    file_format : TYPE
        Output file format.
    dpi : INT
        Resolution.

    Returns
    -------
    None.
    Plots a heatmap and saves the figure in the specified directory.

    '''
    df = bin_matrix.copy()
    df = df.rename(columns = {id_col:'proteins'})
        
    # reorganize df to be in a estetically pleasant order
    df = df.loc[df[use_cols].sum(axis=1) >= 1, ['proteins'] + use_cols]
    df = df.reindex(df.sum().sort_values(ascending=False).index, axis=1)
    df = df.set_index('proteins')
    df = df.sort_values(by=list(df.columns),ascending=False)
    df['sum'] = df[use_cols].sum(axis=1)
    
    # change 1 values by the sum of the row
    for col in use_cols:
        df[col] = np.where(df[col] == 0, 0, df['sum'])  
    df = df.sort_values(by=['sum'] + list(df.columns),ascending=False)
    
    # keep sum col as series and use it to calculate the positions of 
    # vertical lines, then remove this col of df
    s = df.reset_index()['sum'].copy()
    sum_changes = s.diff()[s.diff() != 0].index.values
    sum_changes = sum_changes[1:]
    df = df.drop(['sum'], axis=1)
    df = df.transpose()
        
    # find min non zero
    min_value = int(np.nanmin(df.replace(0, np.nan).values))    
    # find max
    max_value = int(df.values.max())
    extreme_values = [0, min_value, max_value]
    
    # find intermediary values
    all_values = list(np.unique(df))
    int_values = [i for i in all_values if i not in extreme_values ]
    
    # create replace dict
    replace_dict = {        0 : np.nan,
                    min_value : 1,
                    max_value : 3   }
    for i in int_values:
        replace_dict[i] = 2
    print(replace_dict)
    
    # replace values in df
    df = df.replace(replace_dict)
    
    # use s to determine color scheme
    s = [0] + list(set(s))
    n_colors = 3
    # using start=2 for green, rot=0 to only have greens, dark=0.4 to have 
    # the darkest color not too dark and light=1 to have white as the 
    # lightest color
    palette = sns.palplot(sns.cubehelix_palette(start=0, rot=90, dark=0.5, 
                                                light=1, n_colors=n_colors))
    
    palette = sns.color_palette("colorblind", n_colors=n_colors)
   
    # close any existing plots
    plt.close("all")
    # fig, ax = plt.subplots(figsize=(35, 5))
    fig, ax = plt.subplots(figsize=(20, 5))
    sns.heatmap(df, ax=ax, vmin=min(s), vmax=max(s), robust=True, 
                annot=False, square=False, cmap=palette, linewidths=par_lw,
                 xticklabels=True, yticklabels=True)
    # ax.figure.tight_layout()
    
    
    # Manually specify colorbar labelling after it's been generated
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([i for i in range(1,n_colors+1)])
    cat = {1:'exclusive',2:'accessory',3:'core'}
    colorbar.set_ticklabels([cat[i] for i in range(1,n_colors+1)])
    colorbar.ax.tick_params(labelsize=20, axis='both', which='both', length=0)
    
    # make conditions labels bigger
    ax.tick_params(axis='y', which='both', labelsize=20, labelrotation=0)
    ax.set_ylabel('Conditions',size=30)
    
    # adjust proteins labels
    ticks = [i for i in range(0, df.shape[1], 25)]
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)
    ax.tick_params(axis='x', which='both', labelsize=15, labelrotation=0, 
                   reset=False)
    ax.set_xlabel('Proteins',size=30)
    
    # add separation horizontal lines
    ax.hlines([1, 2, 3], color='black', linewidth=2, linestyle='-',
              *ax.get_xlim())
    
    # add separation vertical lines
    # lines_indeces = [min(sum_changes), max(sum_changes)]
    # ax.vlines(lines_indeces, color='black', linewidth=2, linestyle='--',
    #           *ax.get_xlim())
    
    # save
    file_name += '.' + file_format
    plt.savefig(output_dir/file_name, format=file_format, dpi=dpi, 
                    bbox_inches="tight")
    return df

def extract_sub_proteomes(bin_matrix, true_value=1, groups_cols=[], annot='',
                          id_col='gi', output_dir='', filename='', func='core'):
    '''
    Extracts the proteins that are common to all conditions (core), not
    present in any condition (not), present in only 1 condition (single),
    or present in more than 1, but not all condition (accessory).
    The type of computation should be chosen with parameter func: 'core',
    'not', 'single' or 'accessory'.

    Parameters
    ----------
    bin_matrix : DataFrame
        Binnary matrix containing proteins present or absent in each 
        condition.
    true_value : INT/STR, optional
        The value that should be considered true in the binnary matrix.
        The default is 1.
    groups_cols : LIST, optional
        The list of columns that represent the different conditions. 
        The default is [].
    annot : DataFrame, optional
        The daframe supplied with the annotation to make sense of the proteins. 
        The default is ''.
    id_col : STR, optional
        The name of the column in annot whre the index is. 
        It should be also present in bin_matrix, in order to preform merge.
        The default is 'gi'.
    output_dir : STR
        Output directory.
    file_name : STR
        Output file.
    func: STR
        Kind of computation that will be executed to find a subset of proteins.
        Options: 'core', 'not', 'single', 'accessory' or 'all'.
        The default is 'core'.
        
    Returns
    -------
    merged : DataFrame
        DataFrame containing only the specified subdivision of the proteome 
        and its annotation.

    '''
    
    df = bin_matrix.copy()
    
    # Choose rules based on these rules and parameters supplied by the user
    n_groups = len(groups_cols)
    rules = {'core' :      {'lim_inf': n_groups, 'lim_sup': n_groups}, 
             'not' :       {'lim_inf': 0, 'lim_sup': 0}, 
             'yes' :       {'lim_inf': 1, 'lim_sup': n_groups},
             'single' :    {'lim_inf': 1, 'lim_sup': 1}, 
             'accessory' : {'lim_inf': 2, 'lim_sup': n_groups-1},
             'all'       : {'lim_inf': 0, 'lim_sup': n_groups} }
    lims = rules[func]
    
    # apply rules
    filtered = df.loc[ ((df[groups_cols].sum(axis=1) >= lims['lim_inf'] ) &
                        (df[groups_cols].sum(axis=1) <= lims['lim_sup'] )),:]
    
    if len(annot) > 0:
        merged = filtered.merge(annot, on=[id_col], how='inner')
    
    print("\nThe '{}' proteome contains {} proteins".format(
        func, filtered.shape[0]))
    
    excel_file = filename + '.xlsx'
    fasta_file = filename + '.fasta'
    merged.to_excel(output_dir / excel_file, index=False)
    
    md.export_fasta(merged, output_dir, output_file_name=fasta_file,
                filters = {}, id_column='gi', seq_column='Sequence')
    
    return merged

def extract_grouped_proteomes(bin_matrix, annot='', id_col='gi', 
                              output_dir='', filename='',
                              groups={('SEC_UF', 'SEC_YEL') :0, 
                                      ('UC_UF', 'UC_YEL')   :1 } ):
    
    df = bin_matrix.copy()
    
    # Choose rules based on these rules and parameters supplied by the user
    # apply rules
    for group in groups:
        group_val = groups[group]
        mask = (df[list(group)] == group_val).all(axis=1) 
        df = df[mask]
        print('filtered', df)
    
    if len(annot) > 0:
        merged = df.merge(annot, on=[id_col], how='inner')
    else:
        merged = df

    
    print("\nThe '{}' proteome contains {} proteins".format(
        str(groups), merged.shape[0]))
    
    excel_file = filename + '.xlsx'
    fasta_file = filename + '.fasta'
    merged.to_excel(output_dir / excel_file, index=False)
    
    md.export_fasta(merged, output_dir, output_file_name=fasta_file,
                filters = {}, id_column='gi', seq_column='Sequence')
    
    return merged

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
proteomics_file = input_dir / 'proteomics_SEC_and_UC_curated.xlsx'

# READ FILES
proteomics_df = pd.read_excel(proteomics_file)

# FILTER
id_col = 'gi'
quali_cols = ['medium_SEC','medium_UC']
annot_cols = ['query_name', 'emb', 'description', 'Sequence', 'Gene names',
       'Protein names', 'locus_tag', 'Entry', 'UniProtKB',
       'Predicted protein name', 'eggNOG_description', 'GO_terms', 
       'EC_number', 'KEGG_KO', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
       'KEGG_rclass', 'KEGG_TC', 'BRITE', 'BiGG_Reaction', 'COG', 'cello2go',
       'predlipo', 'CAZy']
proteomics_quali = filter_qualitative_data(proteomics_df, id_col, quali_cols)
proteomics_annot = filter_annotation_data(proteomics_df, id_col, annot_cols)

# BINNARY MATRIX
col_groups = {'medium_SEC' : ['SEC_UF', 'SEC_YEL'], 
              'medium_UC' : ['UC_UF',  'UC_YEL']}
med_groups = {'UF'  : ['UF', 'UF/YEL'],
              'YEL' : ['YEL', 'UF/YEL']}
bin_matrix = create_binary_matrix(proteomics_quali, col_groups=col_groups,
                                  med_groups=med_groups, output_dir=output_dir,
                                  file_name='proteome_bin_matrix.xlsx')

img_kargs = {'output_dir':output_dir, 'file_format':'tif', 'dpi':600}

# # PLOT VENN
include_cols = ['SEC_UF', 'SEC_YEL', 'UC_UF', 'UC_YEL']
elements = bin_matrix2dict(df=bin_matrix, true_condition=1, id_col=id_col, 
                            include_cols=include_cols)
plotted_venn = plot_venn(elements, file_name='venn_core', **img_kargs)


# PLOT HEATMAP
plotted_heatmap = plot_heatmap(bin_matrix, use_cols=include_cols, par_lw=0,
                          id_col=id_col, file_name='heatmap_all', **img_kargs)

kargs = {'true_value':1, 'groups_cols':include_cols, 
          'annot':proteomics_annot, 'output_dir':output_dir}

# # EXTRACT CORE PROTEOME

proteome_all = extract_sub_proteomes(bin_matrix, func='all', 
                                      filename='proteome_all', **kargs)

proteome_core = extract_sub_proteomes(bin_matrix, func='core', 
                                      filename='proteome_core', **kargs)

# EXTRACT NOT EV PROTEOME
proteome_not = extract_sub_proteomes(bin_matrix, func='not',
                                    filename='proteome_not_EVs', **kargs)

# EXTRACT YES EV PROTEOME
proteome_yes = extract_sub_proteomes(bin_matrix, func='yes', 
                                      filename='proteome_yes', **kargs)

# EXTRACT SINGLETONS
proteome_single = extract_sub_proteomes(bin_matrix, func='single',
                                      filename='proteome_single', **kargs)

# EXTRACT ACCESSORY PROTEOME
proteome_accessory = extract_sub_proteomes(bin_matrix, func='accessory',
                                  filename='proteome_accessory', **kargs)

# EXTRACT SEC-EXCLUSIVE PROTEOME
proteome_SEC_only = extract_grouped_proteomes(bin_matrix, id_col='gi', 
                                              annot=proteomics_annot, 
                                              output_dir=output_dir,
                              filename='proteome_SEC_only',
                              groups={('SEC_UF', 'SEC_YEL') :1, 
                                      ('UC_UF', 'UC_YEL')   :0 } )

# EXTRACT UC-EXCLUSIVE PROTEOME
proteome_UC_only = extract_grouped_proteomes(bin_matrix, id_col='gi', 
                                              annot=proteomics_annot, 
                                              output_dir=output_dir,
                              filename='proteome_UC_only',
                              groups={('SEC_UF', 'SEC_YEL') :0, 
                                      ('UC_UF', 'UC_YEL')   :1 } )

# GET IMMUNOMODULATORY PROTEINS
immunomodulatory_pf129 = ['CDP49125.1','CDP48267.1','CDP47860.1','CDP48574.1',
                          'CDP48241.1','CDP48736.1','CDP49252.1','CDP48242.1',
                          'CDP48273.1','CDP48858.1','CDP49687.1']
immunomodulatory_pf121_blasted_pf129 = ['659918109', '659917074',
                                        '659917660', '659917537',
                                        '659917391']

sub_core1 = proteome_core.loc[proteome_core['gi'].isin(
                             immunomodulatory_pf121_blasted_pf129),:]
                              
sub_core2 = proteome_core.loc[proteome_core['emb'].isin(
                             immunomodulatory_pf129, ),:]

sub_core = pd.concat([sub_core1, sub_core2])
sub_core_print = sub_core.loc[:,['query_name', 'gi', 'emb', 'description', 
                                 'KEGG_Pathway', 'COG', 'cello2go', 'predlipo']]
sub_core_print.to_excel(output_dir / 'selected_immuno.xlsx', index=False) 


# Export table for annexes
sel_cols = ['query_name', 'gi','emb','locus_tag','description', 'COG', 
            'cello2go','predlipo', 'SEC_UF', 'SEC_YEL', 'UC_UF', 'UC_YEL']
proteome_export = proteome_yes.loc[:,sel_cols]
replace_dict = {'query_name' : 'ID',
                'gi' : 'GI',
                'emb' : 'ACCESSION',
                'locus_tag' : 'LOCUS TAG',
                'description' : 'DESCRIPTION', 
                'COG' : 'COG',
                'cello2go' : 'LOCALIZATION',
                'predlipo' : 'LIPOPROTEIN'}
proteome_export = proteome_export.rename(columns=replace_dict)
proteome_export.to_excel(output_dir / 'proteome_export.xlsx', index=False)
 