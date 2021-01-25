#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 16:48:28 2019

@author: rodovalhovr
"""

import os, json, re, time, platform
import pandas as pd
import urllib.request


##############################################################################

def get_KEGG_data(org='hsa', get_pathway_list=True, get_genes_and_pathways=True,
                  format_conversion=True, format1='hsa', format2='uniprot',
                  genes_names=True):
    """
    Downlod updated data from KEGG.
    INPUT:  org: organism, which info should be downloaded: True/False,
            format1 and format2 for format conversion.
    OUTPUT: Dictionary with all info, not parsed.
    """
    
    # IN 1 ORGANISM: MAP PATHWAYS ID-DESCRIPTION
    if get_pathway_list:
        print("\nFrom KEGG, retrieving map PATHWAYS -> ID-DESCRIPTIONS...")
        url = 'http://rest.kegg.jp/list/pathway/{}'.format(org)
        with urllib.request.urlopen(url) as f:
            pathway_list_data = f.read().decode('utf-8').splitlines()
    else:
        pathway_list_data = None
    # IN 1 ORGANISM: MAP GENES-PATHWAYS
    if get_genes_and_pathways:
        print("From KEGG, retrieving map GENES -> PATHWAYS...")
        url = 'http://rest.kegg.jp/link/pathway/{}'.format(org)
        with urllib.request.urlopen(url) as f:
            genes_and_pathways_data = f.read().decode('utf-8').splitlines()
    else:
        genes_and_pathways_data = None
    # MAP FOR FORMAT CONVERSION (ex: HSA -> UNIPROT)
    if format_conversion and org != 'ko':
        print("From KEGG, retrieving map for FORMAT CONVERSION...")
        url = 'http://rest.kegg.jp/conv/{}/{}'.format(format1, format2)
        with urllib.request.urlopen(url) as f:
            format_map = f.read().decode('utf-8').splitlines()
    else:
        format_map = None
    if genes_names:
        print("From KEGG, retrieving map GENE IDS -> GENE NAMES...")
        url = 'http://rest.kegg.jp/list/{}'.format(org)
        with urllib.request.urlopen(url) as f:
            genes_data = f.read().decode('utf-8').splitlines()
    else:
        genes_data = None
    
    print("Done. Retrieved all information from KEGG.")
    return {'pathways':pathway_list_data, 'format_map':format_map,
            'genes_and_pathways':genes_and_pathways_data, 'genes':genes_data}

def parse_KEGG_pathways_description(raw_data):
    '''
    INPUT:  raw_data is a list of strings obtained from 
            get_KEGG_data['pathways']
    OUTPUT: dictionary mapping pathway id to description
    '''
    print("Parsing pathways descriptions...")
    df = pd.DataFrame(columns=['id','description','organism'])
    pattern = 'path:([\w\d]+)\t([\S\s]+)(?:- )*([\s\S]+)*'
    for line in raw_data:
        match_obj = re.match(pattern, line)
        matched = match_obj.groups()
        data = {'id':matched[0],'description':matched[1],'organism':matched[2]}
        df = df.append(data, ignore_index=True)
    return dict(zip(df.id, df.description))

def parse_KEGG_format_conversion(raw_data, org='hsa'):
    """
    INPUT:  raw_data is a list of strings obtained from 
            get_KEGG_data['format_map']
    OUTPUT: .   
    """
    if not raw_data:
        return None
    else:
        print("Parsing format conversion information...")
        pattern = '[\S]+:([\S]+)\t{}:([\d]+)'.format(org)
        conv = {}
        for string in raw_data:
            match_obj = re.match(pattern, string)
            (key, value) = (match_obj.group(1), match_obj.group(2)) 
            conv[key] = value
        return conv

def parse_KEGG_genes_and_pathways(raw_data, org='hsa'):
    """
    INPUT:  raw_data is a list of strings obtained from 
            get_KEGG_data['genes_and_pathways']
    OUTPUT: .   
    """
    print("Parsing genes-pathways mapping...")
    pattern = '{}:([\S]+)[\s\t]+path:({}[\S]+)'.format(org, org)
    gene2pathway_ids = {}
    for string in raw_data:
        if 'map' not in string:
            match_obj = re.match(pattern, string)
            (gene, pathway_id) = (match_obj.group(1), match_obj.group(2)) 
            if gene in gene2pathway_ids:
                gene2pathway_ids[gene].append(pathway_id)
            else:
                gene2pathway_ids[gene] = [pathway_id]
    return gene2pathway_ids

def parse_genes_names(raw_data, org='hsa'):
    '''
    INPUT:  raw_data is a list of strings obtained from 
            get_KEGG_data['pathways']
    OUTPUT: dictionary mapping pathway id to description
    '''
    print("Parsing genes IDs-genes names mapping...")
    df = pd.DataFrame(columns=['id','name'])
    pattern = '{}:([\w\d]+)\t([\S\s]+)'.format(org)
    for line in raw_data:
        match_obj = re.match(pattern, line)
        matched = match_obj.groups()
        gene_id = matched[0]
        gene_name = matched[1].split(';')[0]
        #.split(',')[0].replace('uncharacterized ','')
        data = {'id' : gene_id, 'name' : gene_name}
        df = df.append(data, ignore_index=True)
    return dict(zip(df.id, df.name))


def aggregate_all():
    '''
    '''
    
    return


def pathways_ids2pathways_descriptions(pathways_descriptions, genes_and_pathways):
    """
    INPUT:  pathways_descriptions: dict, ex: 
                          {'hsa00010':'Glycolysis / Gluconeogenesis'} 
            genes_and_pathways: dict, ex:
                          {10:['hsa00232', 'hsa00983', 'hsa01100', 'hsa05204']}
    OUTPUT: dict mapping genes to pathways description, ex:
                          {10: ['Glycolysis / Gluconeogenesis', ...]}
    """
    genes_and_pathways_description = {}
    for gene in genes_and_pathways:
        pathways_ids = genes_and_pathways[gene]
        for pathways_id in pathways_ids:
            pathway_description = pathways_descriptions[pathways_id]
            if gene in genes_and_pathways_description:
                genes_and_pathways_description[gene].append(pathway_description)
            else:
                genes_and_pathways_description[gene] = [pathway_description]
    return genes_and_pathways_description

def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file)
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return stat.st_mtime
        
def get_all_KEGG_data(storage_file='kegg.data', org='hsa'):
    storage_file_exists = os.path.isfile(storage_file)
    if storage_file_exists:
        get_from = 'file'
    else:
        get_from = 'web'
    print("\nProcessing KEGG data...")
    if get_from == 'web':
        print("Getting data from web...")
        parameters1= {'org':org, 'get_pathway_list':True, 
                  'get_genes_and_pathways':True, 'format_conversion':True, 
                  'genes_names':True,'format1':'hsa', 'format2':'uniprot'}
        raw_data = get_KEGG_data(**parameters1)
        pathways_descriptions = parse_KEGG_pathways_description(raw_data['pathways'])
        format_map = parse_KEGG_format_conversion(raw_data['format_map'], org='hsa')
        genes_and_pathways_ids = parse_KEGG_genes_and_pathways(
                raw_data['genes_and_pathways'], org=org)
        genes_and_pathways_descriptions = pathways_ids2pathways_descriptions(
                pathways_descriptions, genes_and_pathways_ids)
        genes_names = parse_genes_names(raw_data['genes'], org=org)
        complete_dict =  {'gene2pathway_id' : genes_and_pathways_ids, 
                        'gene2pathway_desc' : genes_and_pathways_descriptions, 
                                   'format' : format_map,
                              'genes_names' : genes_names,
                    'pathways_descriptions' : pathways_descriptions}
        parsed_data = json.dumps(complete_dict)
        with open(storage_file, 'a') as storage:
            storage.write(parsed_data)
    elif get_from == 'file':
        date = time.ctime(creation_date(storage_file))
#        date = time.ctime(os.path.getctime(storage_file))
        info_string = "Getting data from file..." 
        info_string += "\nLast modification:\n{}".format(date)
        print(info_string)
        with open(storage_file, 'r') as storage:
            storage_data = storage.read()
        complete_dict = json.loads(storage_data)
    return complete_dict


def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out

def parse_kegg_pathways_hierarchy(kegg_json_file, org='hsa'):
    '''
    Return kegg hierarchy dataframe
    '''
    with open(kegg_json_file, 'r') as input_file:
        data = input_file.read()
    kegg_dict = json.loads(data)
    flattened = flatten_json(kegg_dict)
    df = pd.DataFrame.from_dict(flattened, orient='index').reset_index()
    df.columns = ['children','pathway']
    
    # extract levels    
    df['children'] = df['children'].str.replace('children_','')
    df['children'] = df['children'].str.replace('_name','')
    pattern3 = r'(?P<l1>[\d]+)\_(?P<l2>[\d]+)?\_(?P<l3>[\d]+)'
    df[['l1','l2','l3']] = df['children'].str.extract(pattern3, expand=True)
    pattern2 =  r'(?P<l1>[\d]+)\_(?P<l2>[\d]+)'
    df[['l1','l2']] = df['children'].str.extract(pattern2, expand=True)
    pattern1 =  r'(?P<l1>[\d]+)'
    df[['l1']] = df['children'].str.extract(pattern1, expand=True)
    pattern0 = r'(?P<pathway_id>[\d]+)\s\s(?P<pathway_name>[\S\s]+)'
    df[['pathway_id','pathway_name']] = df['pathway'].str.extract(pattern0, expand=True)

    df['pathway_name'] = df['pathway_name'].fillna(df['pathway'])
    df.drop(columns=['children','pathway'], inplace=True) 
    df = df.astype({"l1": float, "l2": float, "l3": float})
    
    # change from 0-based to 1-based indexation
    df['l3'] = df['l3'] + 1
    df['l2'] = df['l2'] + 1
    df['l1'] = df['l1'] + 1
    df[['l1','l2','l3']] = df[['l1','l2','l3']].fillna(0.0)
    
    # merge l1 and l2
    df['l1_l2'] = df[['l1', 'l2']].apply(lambda x: '-'.join(x.map(int).map(str)), axis=1)
    
    # eliminate first row with all levels equal to 0
    df = df.loc[~((df['l1'] == 0.0) & (df['l2'] == 0.0) & (df['l3'] == 0.0)),:]
    
    # extract sub_dfs containing hierachic levels
    sub_df_l1 = df.loc[((df['l1'] != 0.0) & df['l2'] == 0.0) & (df['l3'] == 0.0),
                       ['l1','pathway_name']]
    sub_df_l1['l1'] = sub_df_l1['l1'].map(int).map(str)
    sub_df_l1 = sub_df_l1.set_index('l1')
    map_l1 = sub_df_l1['pathway_name'].to_dict()
    
    sub_df_l2 = df.loc[((df['l1'] != 0.0) & df['l2'] != 0.0) & (df['l3'] == 0.0),
                       ['l1','l2','pathway_name']]
    sub_df_l2['merged'] = sub_df_l2[['l1', 'l2']].apply(lambda x: '-'.join(x.map(int).map(str)), axis=1)
    sub_df_l2 = sub_df_l2[['merged', 'pathway_name']]
    sub_df_l2 = sub_df_l2.set_index('merged')
    map_l2 = sub_df_l2['pathway_name'].to_dict()
    
    # reorganize columns
    df = df.drop(columns=['l2','l3'])
    cols = ['l1', 'l1_l2', 'pathway_id', 'pathway_name']
    df = df[cols] 
    df['l1'] = df['l1'].map(int).map(str).map(map_l1)
    df['l1_l2'] = df['l1_l2'].map(map_l2)
    df.dropna(subset=['pathway_id'], inplace=True)
    df = df.rename(columns={'l1_l2':'l2'})
    
    # put hsa aprefix if human
    if org:
        df['pathway_id'] = org + df['pathway_id'].astype(str)
    
    return df

def map_kegg_pathways_hierarchy(df, group_by_level=['l1','l2'][0], 
                                  use=['pathway_id', 'pathway_name'][0]):
    '''
    Return mapping dictionary.
    '''
    df = df.copy()
    df = df.set_index(use)
    dictionary = df[group_by_level].to_dict()
    return dictionary

def kegg_pathway_hierarchy(kegg_json_file, org='hsa'):
    '''
    '''
    
    pathway_hierarchy = parse_kegg_pathways_hierarchy(kegg_json_file, org=org)
    hierarchy_level1 = map_kegg_pathways_hierarchy(pathway_hierarchy,
                                             group_by_level='l1',
                                             use='pathway_id')
    hierarchy_level2 = map_kegg_pathways_hierarchy(pathway_hierarchy,
                                             group_by_level='l2',
                                             use='pathway_id')
    hierarchy_level3 = map_kegg_pathways_hierarchy(pathway_hierarchy,
                                             group_by_level='pathway_name',
                                             use='pathway_id')
    return {1: hierarchy_level1, 2: hierarchy_level2, 3: hierarchy_level3}
