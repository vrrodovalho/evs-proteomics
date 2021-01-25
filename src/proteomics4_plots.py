# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho

This module contains scripts for plotting graphs of Qualitative Proteomics
- Donuts for COG categories, with 2 hierachical levels
- Donuts for subcellular localizations with just 1 level
- Simple Venn diagramms with 2 categories and 1 intersection
- Heatmaps that shows frequencies distribution in 2 categories
 
"""

import os
import sys
import pathlib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib_venn

def calculate_cogs(cogs_annotation, cogs_hierarchy, annot_col='COG', 
                   hierarchy_col='l2'):
    '''
    Returns the COG hierarchy with the respective COG category frequencies.

    Parameters
    ----------
    cogs_annotation : DataFrame
        A dataframe containing a column with COG categories.
    cogs_hierarchy : DataFrame
        A dataframe containing a column with COG categories distribuition.
    annot_col : str, optional
        The name of the column in which there is the annotation for COG 
        categories. The default is 'COG'.
    hierarchy_col : str, optional
        The name of the column in which there is the desired level for COG 
        hierarchy representation. The default is 'l2'.

    Returns
    -------
    cogs : DataFrame
        A DataFrame containing COG categories frequencies in a given dataset.

    '''

    annot = cogs_annotation.copy()
    cogs = cogs_hierarchy.copy()
    
    # calculate freqs
    annot[annot_col].fillna('?', inplace=True)
    COGs = list(''.join(annot[annot_col].tolist()))
    COGs_freqs = {i:COGs.count(i) for i in set(COGs)}
    
    # complete cogs with freqs and create lists for graphs
    cogs['freqs'] = cogs[hierarchy_col].map(COGs_freqs)
    cogs.dropna(subset=['freqs'], inplace=True)
    
    # relative freqs
    cogs['rel'] = round(cogs['freqs'] / cogs['freqs'].sum() * 100,1)
    
    return cogs

def construct_list_of_colors(cogs_freqs, higher_hierarchy_col='l1'):
    '''
    Produces the external and internal colors for a COG nested donut plot.

    Parameters
    ----------
    cogs_freqs : DataFrame
        A DataFrame containing COG categories frequencies and the COG 
        2-level hierarchy system.
    higher_hierarchy_col : str
        The name of the column in cogs_freqs where there is the higer COG 
        category. Example: metabolism, information processing...
        The default is 'l1'.

    Returns
    -------
    A tuple of 2 lists of RGB colors. Each RGB color is a tuple of 3 values.
    These are the colors to be used internally and externally in the nested 
    donut plot.
    external_colors : list
        Lists of darker RGB colors that will be used externally.
    internal_colors : list
        Lists of gradually lighter RGB colors that will be used internally.

    '''
    palette = [plt.cm.Blues, plt.cm.Oranges, plt.cm.Purples, plt.cm.Greys]
    external_colors = [palette[0](0.7), palette[1](0.7), 
                       palette[2](0.7), palette[3](0.7)]
    internal_colors = []
    subgroups_sizes = cogs_freqs.groupby(higher_hierarchy_col).size().to_dict()
    group = 0
    higher_hierarchy_group = list(
        cogs_freqs.loc[:,higher_hierarchy_col].unique())
    for subgroup in higher_hierarchy_group:
        size = subgroups_sizes[subgroup]
        color_group = palette[group]
        group += 1
        decimals = np.arange(0.60, 0.05, -0.05)
        colors = [color_group(i) for i in decimals[:size]]
        internal_colors.extend(colors)
    return external_colors, internal_colors

def plot_nested_cog(cogs_annotation, annot_col, cogs_hierarchy_map, 
                    hierarchy_levels=['l1','l2'],  output_dir='', filename='', 
                    save_fig=True, file_format='tif', dpi=600):
    '''
    Plot a nested donut graph, with a 2-level hierarchy system for COG 
    categories.

    Parameters
    ----------
    cogs_annotation : DataFrame
        DataFrame containing one column with COG categories annotation.
    annot_col: str
        The name of the column in cogs_annotation where the annotation is.
    cogs_hierarchy_map : DataFrame
        A DataFrame contaning at least 2 columns with the 2-level hierarchy
        system.
    hierarchy_levels : list
        The list of 2 strings representing the names of the columns in 
        cogs_hierarchy_map that contains hierarchy level data. 
        The default is ['l1','l2'].
    output_dir : str
        The path where the file will be saved. 
    filename : str
        The name of the file that will be saved. 
    save_fig : bool
        Wether to save or not the image in a file. The default is True.

    Returns
    -------
    cogs_freqs : DataFrame
        A DataFrame containing COG categories frequencies in a given dataset.

    '''
    
    lvl = hierarchy_levels
    cogs_freqs = calculate_cogs(cogs_annotation=cogs_annotation, 
                                annot_col=annot_col,
                                cogs_hierarchy=cogs_hierarchy_map, 
                                hierarchy_col=lvl[1])
    

    l1 = cogs_freqs.groupby(lvl[0]).sum()['freqs'].reset_index()
    l2 = cogs_freqs[[lvl[1],'freqs']].reset_index()
    l1_names, l1_values = l1[lvl[0]].to_list(), l1['freqs'].to_list()
    l2_names, l2_values = l2[lvl[1]].to_list(), l2['freqs'].to_list()
    l1_names_lines = ['\nAND '.join(i.split('AND')) for i in l1_names]
    
    # Create colors
    external_colors, internal_colors = construct_list_of_colors(cogs_freqs, 
                                                   higher_hierarchy_col='l1')
    
    # initialize figure
    fig, ax = plt.subplots()
    ax.axis('equal')
     
    # First Ring (outside)
    mypie, _ = ax.pie(l1_values, radius=1.3, labels=l1_names_lines, 
                      labeldistance=1.1,
                      colors=external_colors )
    plt.setp( mypie, width=0.3, edgecolor='white')
    
    # Second Ring (Inside)
    mypie2, _ = ax.pie(l2_values, radius=1.3-0.3, 
                       labels=l2_names, labeldistance=0.8, 
                       colors=internal_colors)
    plt.setp( mypie2, width=0.4, edgecolor='white')
    plt.margins(0,0)
    
    if save_fig:
        output_file = output_dir / filename
        plt.savefig(output_file,format=file_format, dpi=dpi,
                    bbox_inches="tight")
    
    # show it
    plt.show()
    
    return cogs_freqs

def plot_donut(df, data_column='', title='', relative=True, absolute=True,
               palette=sns.color_palette("colorblind", 30), save_fig=True,
               output_dir='', file_name='',  file_format='tiff', dpi=600):
    '''
    Plots a DONUT graph, representing proportions of categories.

    Parameters
    ----------
    df : DataFrame
        A dataframe containing a column in which is the data whose frequencies
        will be plotted.
    data_column : str
        The name of the column in df where the data is placed. 
    title : str, optional
        A title that will be plotted above the gaph. The default is ''.
    relative : bool, optional
        Wether to include relative frequencies of the categories (%). 
        The default is True.
    absolute : bool, optional
        Wether to include absolute frequencies of the categories. 
        The default is True.
    palette : Seaborn palette object, optional
        Seaborn color palette object. The default is 
        sns.color_palette("colorblind", 30).
    save_fig : bool, optional
        Wether to save the image to a file. The default is True.
    output_dir : str, optional
        Directory where the image will be saved. 
    file_name : str, optional
        File name to save image. 
    file_format : str, optional
        File format to save image. The default is 'tiff'.
    dpi : int, optional
        Resolution to save image. The default is 600.

    Returns
    -------
    data : DataFrame
        DataFrame containing the calculated frequencies.

    '''
    df = df.copy()
    
    # renames
    rename_dict = {"Cytoplasmic Membrane":"Membrane"}
    df[data_column].replace(rename_dict, inplace=True)
    
    # get data, count them and put the results in X and Y lists
    df[data_column].fillna('?', inplace=True)
    df = df.reset_index()
    items = df[data_column].tolist()
    freqs = {i:items.count(i) for i in set(items)}
    data = pd.DataFrame.from_dict(freqs, orient='index', columns=['abs'])
    freqs_rel = {i: round(items.count(i)/len(items)*100,1) \
                          for i in set(items)}
    data['rel'] = data.index.to_series().map(freqs_rel)
    data = data.sort_values(by=['abs'], ascending=False)
    y = data['abs']
    
    # choose representation in absolute values and/or percentage
    if relative and absolute:
        data['formatted'] = data.index.astype(str) + ' ' + \
                            data['abs'].astype(str) + ' (' + \
                            data['rel'].astype(str) + '%)'
    elif relative and not absolute:
        data['formatted'] = data.index.astype(str) + ' ' + \
                            ' (' + \
                            data['rel'].astype(str) + '%)'
    elif not relative and absolute:
        data['formatted'] = data.index.astype(str) + ' (' + \
                            data['abs'].astype(str) + ')' 
    x = data['formatted']
    

    # start plotting
    fig, ax = plt.subplots()
    ax.axis('equal')  
    my_circle=plt.Circle( (0,0), 0.6, color='white')
    plt.pie(y, labels=x, colors=palette, startangle=0,counterclock=False)
    p = plt.gcf()
    p.gca().add_artist(my_circle)
        
    # draw!
    plt.draw()
    plt.title(title)    
    
    # save figure
    if save_fig:
        plt.savefig(output_dir/file_name, format=file_format, dpi=dpi, 
                    bbox_inches="tight")
    return data

def plot_venn(df, data_column="", weighted=False, intersection_cat='UF/YEL', 
              title='', palette=sns.color_palette("colorblind", 30), 
              save_fig=True, output_dir='', file_name='', file_format='tiff', 
              dpi=600):
    '''
    This is a function for plotting a simple venn diagramm with 2 categories.

    Parameters
    ----------
    df : DataFrame
        A dataframe containing a column in which there is the data to be
        plotted in a Venn diagramm.
    data_column : str
        The name of the column in df where the data is placed. 
    weighted : bool
        If False, venn circles have the same sizes. If True, the sizes are
        proportional to the numeric value.
    intersection_cat : str
        One of the 3 categories in data_column, that will be considered the
        intersection category in the venn diagramm.
    title : str, optional
        A title that will be plotted above the gaph. The default is ''.
    palette : Seaborn palette object, optional
        Seaborn color palette object. The default is 
        sns.color_palette("colorblind", 30).
    save_fig : bool, optional
        Wether to save the image to a file. The default is True.
    output_dir : str, optional
        Directory where the image will be saved. 
    file_name : str, optional
        File name to save image. 
    file_format : str, optional
        File format to save image. The default is 'tiff'.
    dpi : int, optional
        Resolution to save image. The default is 600.

    Returns
    -------
    data : DataFrame
        A DataFrame containing the values used to plot the venn.

    '''
    df = df.copy()

    # get data, count them and put the results in X and Y lists
    df[data_column].fillna('?', inplace=True)
    df = df.reset_index()
    items = df[data_column].tolist()
    freqs = {i:items.count(i) for i in set(items)}
    data = pd.DataFrame.from_dict(freqs, orient='index', columns=['freq'])
    data = data.sort_values(by=['freq'], ascending=False)

    # prepare subsets and set_labels for venn
    intersection_value = data.loc[ intersection_cat, 'freq']
    data = data.drop([intersection_cat])
    x = data.index.to_list()
    y = data['freq']    
    subsets = (y[0], y[1], intersection_value)
    set_labels = (x[0], x[1])
    
    # start plotting
    fig, ax = plt.subplots()
    if weighted:
        v = matplotlib_venn.venn2(subsets=subsets, set_labels=set_labels)
    else:
        v = matplotlib_venn.venn2_unweighted(subsets=subsets, 
                                             set_labels=set_labels)
        
    # set colors and alphas
    v.get_patch_by_id('10').set_color(palette[0])
    v.get_patch_by_id('10').set_alpha(0.8)
    v.get_patch_by_id('01').set_color(palette[1])
    v.get_patch_by_id('01').set_alpha(0.8)
    v.get_patch_by_id('11').set_color(palette[2])
    v.get_patch_by_id('11').set_edgecolor('none')
    v.get_patch_by_id('11').set_alpha(0.6)
    
    # set font sizes
    for text in v.set_labels:
        text.set_fontsize(14)
    for text in v.subset_labels:
        text.set_fontsize(16)
        
    # set label positions
    lbl_a = v.get_label_by_id("A")
    xa, ya = lbl_a.get_position()
    lbl_a.set_position((xa-0.2, ya+0.05))
    lbl_b = v.get_label_by_id("B")
    xb, yb = lbl_b.get_position()
    lbl_b.set_position((xb+0.25, yb+0.1))

    # draw!
    plt.draw()
    plt.title(title)    
    
    # save figure
    if save_fig:
        plt.savefig(output_dir/file_name, format=file_format, dpi=dpi, 
                    bbox_inches="tight")
    return data

def helper_frequencies(df, freq_column, split_char=False, 
                       forbidden_prefix='map'):
    '''
    This is a helper function that extracts qualitative data from a column
    of a DataFrame, counts them and return a dictionary of frequencies.
    If there is more than 1 value in a row, it is possible to split this 
    row and account for each value separately. It is also possible to
    exclude values based on a prefix.

    Parameters
    ----------
    df : DataFrame
        A dataframe containing the data to be analyzed.
    freq_column : str
        The column in df which contains the data whose frequencies will be
        calculated. 
    split_char : str, bool or NoneType, optional
        A character that will be used to split each row if there are multiple
        values in each row. If this is set to False, each row will be 
        considered a single value. If this is set to None, each row will be
        considered to contain multiple data, but no delimiter character.
        If this is set to a string, this string will be considered the split
        character that separates the values in each row. The default is False.
    forbidden_prefix : str, optional
        Values starting with the string set in this parameter will be excluded
        from analysis. The default is 'map'.

    Returns
    -------
    freqs : dict
        A dictionary of frequencies.

    '''
    values = df[freq_column].tolist()
    if split_char == False:
        new_values = values
    elif split_char == None:
        string = ''.join(values)
        new_values = list(string)
    else:
        string = split_char.join(values)
        new_values = string.split(split_char)
    
    if forbidden_prefix:    
        filtered_list = [x for x in new_values \
                         if not x.startswith(forbidden_prefix)]
    else:
        filtered_list = new_values
    freqs = {i:filtered_list.count(i) for i in set(filtered_list)}
    return freqs

def compare_n_frequencies(df, freq_column='COG', category_column='medium',
                          category_map={}, split_char=',', drop_empties=False):
    '''
    This functions compare the frequencies of values in 2 or more categories.
    Example: the frequencies of COG categories (the values) are compared in
             2 conditions (defined as categories), such as 2 culture media.
    It returns a dataframe with multiple frequencies columns, one column
    for each category that ws specified.

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    freq_column : str, optional
        The name of the column where are the values whose frequency will be 
        calculated. The default is 'COG'.
    category_column : str, optional
        The name of the column where the categories are specified. Each one
        of this categories will be represented as a column of frequencies in
        the final dataframe. The default is 'medium'.
    category_map : dict, optional
        A dictionary that maps multiple values for each category, grouping
        values in categories or changing their name. If each individual value 
        should be accounted as its own category, without change, an empty
        dictionary should be passed. Otherwise, categories should be organized 
        in the format {'':[]}. The default is {}.
   split_char : str, bool or NoneType, optional
        A character that will be used to split each row if there are multiple
        values in each row. If this is set to False, each row will be 
        considered a single value. If this is set to None, each row will be
        considered to contain multiple data, but no delimiter character.
        If this is set to a string, this string will be considered the split
        character that separates the values in each row. The default is False.
   drop_empties : bool, optional
        This allows to choose if empty annotations (nan) will not be 
        considered, being dropped (True) or considered (False) and atributted 
        to unknown function (S/?). The default is False.

    Returns
    -------
    new_df : TYPE
        DESCRIPTION.

    '''
    df = df.copy()
    
    # choose if empty annotations will not be considered or considered as 
    # unknown function (S/?)
    if drop_empties:
        df.dropna(subset=[freq_column], inplace=True)
    else:
        if freq_column == "COG":
            df[freq_column].fillna('S', inplace=True) 
        else:
            df[freq_column].fillna('?', inplace=True) 
    
    # if category mapping is supplied, get them. 
    # otherwise, consider each unique value in category column
    if category_map:
        categories = list(category_map.keys())
    else:
        categories = list(df[category_column].unique())
        category_map = { i:[i] for i in categories }
    # and a list of all COG categories found
    all_COGs = []
    
    # generate a dict  of frequencies for each category
    frequencies = {}
    for cat in categories:
        sub_df = df.loc[ df[category_column].isin( category_map[cat] ),:]
        freqs = helper_frequencies(sub_df, freq_column=freq_column, 
                                    split_char=split_char)
        all_COGs.extend(list(freqs.keys()))
        frequencies[cat] = freqs
    
    # generate a dataframe with all COG categories found and a column
    # with the frequencies for each category
    new_df = pd.DataFrame(list(set(all_COGs)), columns=[freq_column])
    
    for cat in categories:
        freqs = frequencies[cat]
        new_df[cat] = new_df[freq_column].map(freqs)
    # if nan was assigned to a category frequency, replace it by 0
    new_df.fillna(0, inplace=True)
    
    # add a column with the total sum to the dataframe
    new_df['Total'] = new_df.loc[:,categories].sum(axis=1)
    
    return new_df

def plot_heatmap(df, cat1_col='COG', cat1_split_char=False, cat2_col='', 
                 cat2_map='', sort_heatmap_by='Total', extra_map={},
                 replace_names={}, save_fig=True, output_dir='', filename='', 
                 file_format='tif', dpi=600, colors="Blues"):
    '''
    Plots a heatmap to show frequencies distribution towards 2 categories:
        a main category, such as COG category, and a secondary category,
        such as culture medium.

    Parameters
    ----------
    df : DataFrame
        A Dataframe containing at least 2 columns, representing the 2 
        variables that should be considered in the representation.
    cat1_col : str
        The name of the column of the main category. The default is 'COG'.
    cat1_split_char : str, bool or NoneType, optional
        A character that will be used to split each row if there are multiple
        values in each row. If this is set to False, each row will be 
        considered a single value. If this is set to None, each row will be
        considered to contain multiple data, but no delimiter character.
        If this is set to a string, this string will be considered the split
        character that separates the values in each row. The default is False.
    cat2_col : str
        The name of the column of the secondary category. 
        The default is 'medium'.
    cat2_map : dict, optional
        A dictionary that maps multiple values for each category, grouping
        values in categories or changing their name. If each individual value 
        should be accounted as its own category, without change, an empty
        dictionary should be passed. Otherwise, categories should be organized 
        in the format {'':[]}. The default is {}.
    sort_heatmap_by : str, optional
        The name of the column in final heatmap that will be used to sort it. 
        The default is 'Total', the column produced with the sum.
    extra_map : dict, optional
        A dictionary to replace some texts in heatmap, for a more complete
        annotation, if needed. The keys should be the the categories in 
        cat1_col. The default is {}.
    replace_names : dict, optional
        A dictionary to replace column names in the heatmap for a better 
        representation, id needed. The default is {}.
    save_fig : bool, optional
        Wether to save the image to a file. The default is True.
    output_dir : str, optional
        Directory where the image will be saved. 
    file_name : str, optional
        File name to save image. 
    file_format : str, optional
        File format to save image. The default is 'tiff'.
    dpi : int, optional
        Resolution to save image. The default is 600.
    colors : str, optional
        Colors for the heatmap. The default is "Blues".

    Returns
    -------
    df : DataFrame.
        A dataframe representing the final heatmap

    '''
    
    # calculate frequencies for main category
    column = cat1_col
    freqs = compare_n_frequencies(df, freq_column=column, drop_empties=False, 
                                  category_column=cat2_col,
                                  category_map=cat2_map, 
                                  split_char=cat1_split_char)
    
    # replace columns name in order to show in heatmap
    df = freqs.copy()
    if column in replace_names:
        column = replace_names[column]
    df.columns = df.columns.to_series().replace(replace_names)
    df = df.set_index(column)
    df = df.astype(int)
    df = df.reset_index()
    
    # if we have to complete the text of the main category that will be showed 
    # in the heatmap
    if extra_map:
        df[column] = df[column].map(extra_map)
    
    # if we need to sort the heatmap, for better presentation of the colors
    if len(sort_heatmap_by) > 0:
        df = df.sort_values(by=[sort_heatmap_by],ascending=False)
    
    # start plotting
    # choose sizes for heatmap depending if its a COG heatmap or smaller ones
    sns.set_style("dark")
    len_y, len_x = df.shape[0], df.shape[1]
    if column == 'COG category':
        sns.set(font_scale=0.8)
        cog_size=(19,5)
        len_x = cog_size[1] 
        len_y = cog_size[0]
        fig2 = plt.figure(figsize=[4, 4], constrained_layout=False)
    else:
        sns.set(font_scale=1.3)
        fig2 = plt.figure(constrained_layout=False)
    
    # add subplots
    gs = fig2.add_gridspec(len_y, len_x)
    ax1 = fig2.add_subplot(gs[:,:-1])
    ax2 = fig2.add_subplot(gs[:,-1:])
    df = df.set_index(column)
    df1 = df.iloc[:,:-1]
    df2 = df.iloc[:,-1:]
    
    # plot main heatmap
    sns.heatmap(df1, cmap=colors, center=1, annot=True, # annot_kws=squares_font, 
                fmt="d", ax=ax1, cbar = False, robust=True)
    df = df.reset_index()
    
    # y and x labels
    ax1.set_yticklabels(labels=df[column], va='center', rotation=0, 
                        position=(0,0.28)) 
    ax1.set_xticklabels(labels=df.columns[1:-1])
    
    # plot totals heatmap
    df2.index.name = None
    sns.heatmap(df2, cmap=colors, center=1, annot=True, #annot_kws=squares_font, 
                fmt="d", ax=ax2, cbar=True, xticklabels=True, yticklabels=False,
                robust=True)
    
    # save figure
    if save_fig:
        fig2.savefig(output_dir/filename, format=file_format, dpi=dpi, 
                    bbox_inches="tight")
    return df 


def hierarchical_grouped_barplots(df, cat1_col='COG', cat1_split_char=False, 
                                  cat2_col='', cat2_map='', 
                                  sort_heatmap_by='Total', extra_map={},
                                  replace_names={}, save_fig=True, 
                                  output_dir='', filename='',  
                                  file_format='tif', dpi=600, colors="Blues"):
    
    
    # calculate frequencies for main category
    column = cat1_col
    freqs = compare_n_frequencies(df, freq_column=column, drop_empties=False, 
                                  category_column=cat2_col,
                                  category_map=cat2_map, 
                                  split_char=cat1_split_char)
    
    # replace columns name in order to show in heatmap
    df = freqs.copy()
    if column in replace_names:
        column = replace_names[column]
    df.columns = df.columns.to_series().replace(replace_names)
    df = df.set_index(column)
    df = df.astype(int)
    df = df.reset_index()
   
    # if we have to complete the text of the main category that will be showed 
    # in the heatmap
    if extra_map:
        df[column] = df[column].map(extra_map)
    
    # if we need to sort the heatmap, for better presentation of the colors
    if len(sort_heatmap_by) > 0:
        df = df.sort_values(by=[sort_heatmap_by],ascending=False)
    
    df = pd.melt(df, id_vars=['COG category'], var_name='Condition', 
                 value_name='Absolute frequency',
                 value_vars=['UF-only','UF-YEL', 'YEL-only'])

    # penguins = sns.load_dataset("penguins")
    
    # fig = plt.figure(figsize=[4, 4], constrained_layout=False)
    
    # Draw a nested barplot by species and sex
    g = sns.catplot(data=df, kind="bar", x="COG category", y='Absolute frequency', 
                    hue="Condition", ci="sd", palette="colorblind", alpha=.9, 
                    height=6, aspect=2.0)
    # g.despine(left=True)

    # g.legend.set_title("")
    

    return df


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
proteomics_SEC_only_file = input_dir / 'proteome_SEC_only.xlsx'
proteomics_UC_only_file = input_dir / 'proteome_UC_only.xlsx'
cogs_file = input_dir / 'COGs.xlsx'

# READ FILES
proteomics_SEC_and_UC = pd.read_excel(proteomics_SEC_and_UC_file)
proteomics_UC = pd.read_excel(proteomics_UC_file)
proteomics_core = pd.read_excel(proteomics_core_file)
proteomics_accessory = pd.read_excel(proteomics_accessory_file)
proteomics_single = pd.read_excel(proteomics_single_file)
proteomics_not_EVs = pd.read_excel(proteomics_not_EVs_file)
proteomics_SEC_only = pd.read_excel(proteomics_SEC_only_file)
proteomics_UC_only = pd.read_excel(proteomics_UC_only_file)
cogs_map = pd.read_excel(cogs_file)


# common args
figs_kargs = {'output_dir':output_dir, 'save_fig':True, 'file_format':'tiff',
              'dpi':600}
cog_kargs = {'annot_col':'COG', 'cogs_hierarchy_map':cogs_map }
donut_kargs = {'palette':sns.color_palette("colorblind", 30), 
               'relative':True, 'absolute':True, }

# # PLOT CORE
# core_COG = plot_nested_cog(cogs_annotation=proteomics_core, 
#                             filename='core_COG.tif', **figs_kargs, **cog_kargs)
# core_PREDLIPO = plot_donut(df=proteomics_core, data_column='predlipo', 
#           title='', file_name='core_PREDLIPO.tif', **figs_kargs, **donut_kargs)
# core_CELLO = plot_donut(df=proteomics_core, data_column='cello2go', 
#             title='', file_name='core_CELLO.tif', **figs_kargs, **donut_kargs)

# # PLOT ACCESSORY
# accessory_COG = plot_nested_cog(cogs_annotation=proteomics_accessory, 
#                       filename='accessory_COG.tif', **figs_kargs, **cog_kargs)
# accessory_PREDLIPO = plot_donut(df=proteomics_accessory, data_column='predlipo', 
#                                 title='', file_name='accessory_PREDLIPO.tif', 
#                                 **figs_kargs, **donut_kargs)
# accessory_CELLO = plot_donut(df=proteomics_accessory, data_column='cello2go', 
#                               title='', file_name='accessory_CELLO.tif', 
#                               **figs_kargs, **donut_kargs)

# # # PLOT SINGLE
# single_COG = plot_nested_cog(cogs_annotation=proteomics_single, 
#                           filename='single_COG.tif', **figs_kargs, **cog_kargs)
# single_PREDLIPO = plot_donut(df=proteomics_single, data_column='predlipo', 
#                               title='', file_name='single_PREDLIPO.tif',  
#                               **figs_kargs, **donut_kargs)
# single_CELLO = plot_donut(df=proteomics_single, data_column='cello2go', 
#                           title='', file_name='single_CELLO.tif', 
#                           **figs_kargs, **donut_kargs)

# # PLOT SEC_only
SEC_only_COG = plot_nested_cog(cogs_annotation=proteomics_SEC_only, 
                          filename='SEC_only_COG.tif', **figs_kargs, **cog_kargs)
SEC_only_PREDLIPO = plot_donut(df=proteomics_SEC_only, data_column='predlipo', 
                              title='', file_name='SEC_only_PREDLIPO.tif',  
                              **figs_kargs, **donut_kargs)
SEC_only_CELLO = plot_donut(df=proteomics_SEC_only, data_column='cello2go', 
                          title='', file_name='SEC_only_CELLO.tif', 
                          **figs_kargs, **donut_kargs)

# # PLOT UC_only
UC_only_COG = plot_nested_cog(cogs_annotation=proteomics_UC_only, 
                          filename='UC_only_COG.tif', **figs_kargs, **cog_kargs)
UC_only_PREDLIPO = plot_donut(df=proteomics_UC_only, data_column='predlipo', 
                              title='', file_name='UC_only_PREDLIPO.tif',  
                              **figs_kargs, **donut_kargs)
UC_only_CELLO = plot_donut(df=proteomics_UC_only, data_column='cello2go', 
                          title='', file_name='UC_only_CELLO.tif', 
                          **figs_kargs, **donut_kargs)


# # QUALITATIVE PROTEOMICS 
df_UC_EVs = proteomics_UC.loc[proteomics_UC['medium'] != 'notEVs',:]

# # QUALITATIVE PROTEOMICS - PLOT UC VENN
# venn = plot_venn(df=df_UC_EVs, data_column="medium", weighted=False, 
#       intersection_cat='UF/YEL', title='', file_name='venn.tif', **figs_kargs)

# QUALITATIVE PROTEOMICS - PLOT UC HEATMAPS

# common args
heatmaps_kargs = { 'colors':"Blues", 'sort_heatmap_by':'Total', 
                  'replace_names': {'COG':'COG category', 
                                    'predlipo':'PRED-LIPO category',
                                    'cello2go':'Subcellular Localization', 
                                    'UF':'UF-only', 'YEL':'YEL-only'},
                  'cat2_col': 'medium', 
                  'cat2_map': {'UF':['UF'], 'UF-YEL':['UF/YEL'], 
                               'YEL':['YEL']}
                      }

# an extra map for the COGs
cog_names = dict(zip(cogs_map.l2, cogs_map.l2_description))
cog_names2 = { letter : letter + ' : ' + cog_names[letter]  \
              for letter in cog_names}
  
# # plot all heatmaps
# plot_hmp_COGs = plot_heatmap(df=df_UC_EVs, cat1_col='COG', 
#                              cat1_split_char=None, extra_map=cog_names2, 
#                              filename='quali_heatmap_COG.tif', 
#                              **figs_kargs, **heatmaps_kargs)    
# plot_hmp_lipo = plot_heatmap(df=df_UC_EVs, cat1_col='predlipo', 
#                              cat1_split_char=False, 
#                              filename='quali_heatmap_predlipo.tif', 
#                              **figs_kargs, **heatmaps_kargs)    
# plot_hmp_cello = plot_heatmap(df=df_UC_EVs, cat1_col='cello2go', 
#                               cat1_split_char=False, 
#                               filename='quali_heatmap_cello2go.tif', 
#                               **figs_kargs, **heatmaps_kargs)    
        
    
plot_grouped_bars = hierarchical_grouped_barplots(df=df_UC_EVs, cat1_col='COG',
                                                  cat1_split_char=None, 
                                                  **heatmaps_kargs)