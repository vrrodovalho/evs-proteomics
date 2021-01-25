# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 11:55:19 2019

@author: vrrod

Set of useful modules
"""

import os 
import re
import sys 
import urllib
import platform
import pandas as pd
import numpy as np
from itertools import islice




def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def get_dir_above(path):
    return os.path.dirname(path)

def invert(dictionary):
    """
    Invert keys and values.
    """
    return {v: k for k, v in dictionary.items()}


def fasta_parser(fasta_file, simplify_id=True, how=('|',1), verbose=True):
    """
    
    """
    with open(fasta_file, 'r') as input_file:
        data = input_file.readlines()
    
    sequences = {}
    duplicated_ids = {}
    for line in data:
        if line.startswith(">") and line != "" and line != " ":
            fasta_id = line[1:]
            fasta_id = fasta_id.strip()
            if simplify_id:
                fasta_id = fasta_id.split(how[0])[how[1]]
            if fasta_id in sequences:
                duplicated_ids[fasta_id] = ""
                duplicated = True
            else:
                sequences[fasta_id] = ""
                duplicated = False
        else:
            if duplicated:
                duplicated_ids[fasta_id] += line.strip()
            else:
                sequences[fasta_id] += line.strip()
    if duplicated_ids:
        final = (sequences, duplicated_ids)
        if verbose:
            print("\nGot {} unique fasta ids!".format(len(sequences))) 
            print("Got {} duplicated fasta ids!".format(len(duplicated_ids)))
            print("Returned as tuple (sequences, duplicated_ids)")
    else:
        final = sequences
        if verbose:
            print("\nGot {} unique fasta ids!".format(len(sequences))) 
            print("Returned a sequences dictionary")    
    return final

def filter_sequences(sequences, min_seq_len=35, forbidden=['B','J','O','U','X','Z']):
    ""
    ""
    filtered_sequences = {}
    for fasta_id in sequences:
        seq = sequences[fasta_id].upper()
        seq_len = len(seq)
        if seq_len > min_seq_len and all(char not in seq for char in forbidden):
            filtered_sequences[fasta_id] = sequences[fasta_id]
    print("\nFiltered {} out of {} sequences.".format(
            len(sequences) - len(filtered_sequences), len(sequences)))
    return filtered_sequences

def fasta_list2dict(fasta_list):
    """
    """
    fasta_dict = {}
    fasta_id = ''
    for line in fasta_list:
        if line not in ['', ' ', '\n']:
            if line.startswith(">"):
                fasta_id = line[1:]
                fasta_dict[fasta_id] = ''
            else:
                fasta_dict[fasta_id] += line.strip()
    return fasta_dict
        

def uniprot(query='P31946', file_format='fasta', include='no', input_id="ID",
            output_id="ID", taxid='9606'):
    """
    """
    url = 'https://www.uniprot.org/uniprot/'
    params = {'query': query,
              'from': input_id, 
              'to': output_id, 
              'organism': taxid,
              'format': file_format, 
              'include': include}

    data = urllib.parse.urlencode(params).encode("utf-8")
    requested = urllib.request.Request(url, data)
    # Please set a contact email address here to help us debug in case of 
    # problems (see https://www.uniprot.org/help/privacy).
    #contact = "rodovalhovr@yahoo.com" 
    #requested.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(requested)
    page = response.read(200000)
    return str(page.decode("utf-8"))

def retrieve_data(ids_list, file_format=['fasta','tab'][0]):
    """
    Retrive fasta or tab strings from Uniprot IDs
    """
    data = []
    for uniprot_id in ids_list:
        this_query = uniprot(query=uniprot_id, file_format=file_format)
        data.append(this_query)
    return data

def parse_uniprot(tab_data, fasta_data):
    """
    Parse tabular and fasta data from Uniprot to df
    """
    tab_data_list1D = ''.join(tab_data).strip().split('\n')
    #print("tab_data_list1D",tab_data_list1D)
    tab_data_list2D = [string.split('\t') for string in tab_data_list1D]
    #print("tab_data_list2D",tab_data_list2D)
    df = pd.DataFrame(tab_data_list2D, columns=tab_data_list2D[0])
    mask = df['Entry'].isin(['Entry'])
    df = df[~mask]
    #fasta_string = ''.join(fasta_data).strip()
    #fasta_lines = fasta_string.split('\n')
    #fasta_dict = fasta_list2dict(fasta_lines)
    #new_fasta_dict = {}
    #for fasta_id in fasta_dict:
    #    old_key, new_key = fasta_id, fasta_id.split("|")[1]
    #    new_fasta_dict[new_key] = fasta_dict[old_key]
    #df["sequence"] = df["Entry"].map(new_fasta_dict)
    return df


def export(df, column, output):
    """
    """
    np.savetxt(output, df[column].values, fmt='%s')
    return list(df[column].values)

def insert_newlines(string, every=64):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

def df2fasta_dict(df, output_dir, output_file_name="",
                 id_columns=['gene_from_uniprot_splitted'], 
                 seq_column='sequence', filter_column="medium", 
                 filter_scope=["U", "Y", "2"]):
    """
    Construct fasta file from dataframe
    """
    if len(filter_column) > 0:
        df = df[df[filter_column].isin(filter_scope)]
    df2 = pd.DataFrame(columns=['fasta_id', 'sequence'])
    
    id_columns_numbers = [df.columns.get_loc(x) for x in id_columns]
    seq_column_number = df.columns.get_loc(seq_column) 
    
    df2["fasta_id"] = ( df.iloc[:,id_columns_numbers]
                       .apply(lambda x: ' '.join(x), axis=1)
                     )
    df2["sequence"] = df.iloc[:,seq_column_number]
    fasta_dict = dict(zip(df2.fasta_id, df2.sequence))

    final_string = ""    
    for fasta_id in fasta_dict:
        sequence = fasta_dict[fasta_id] + "\n"
        fasta_id = ">" + fasta_id + "\n"
        final_string += fasta_id
        final_string += insert_newlines(sequence, every=64) 
    
    output_file_name += ("EVs_Pf129_" + "".join(filter_scope) + ".fasta")
    output_file = output_dir / output_file_name

    with open(output_file, 'w') as out:
        out.write(final_string )
    
    print("Exported %d sequences to %s." % (len(fasta_dict), output_file))
    
    return fasta_dict

def search_empty_row(df, show_id = 'gi', search_column = 'COG'):
    print(df[show_id][df[search_column].isnull()])

def check_duplicates(df, column="gene_from_uniprot"):
    """
    """
    dup = ( df[
            df[column]
            .str.split(';')
            .str[0].duplicated()
                    ][column] )
    if not dup.empty:
        print("In column {} got {} duplicated rows: {}".format(
                column, dup.shape[0], dup.index))
    else:
        print("In column {}, got no duplicated rows".format(column))
        
def fasta_dict2file(fasta_dict, output):
    """
    """
    print("Exporting sequences to file {}".format(output))
    final_string = ''
    for fasta_id in fasta_dict:
        sequence = insert_newlines(fasta_dict[fasta_id], every=64)
        this_entry = ">{}\n{}\n".format(fasta_id, sequence)
        final_string += this_entry
    with open(output, 'w') as out:
        out.write(final_string)
    return None

def change_fasta_id(fasta_dict, pattern='\w+\|(\S+)\|([\S\s]+)OS=.?'):
    new_dict = {}
    for fasta_id in fasta_dict:
        old_value = fasta_dict[fasta_id]
        matches = re.match(pattern, fasta_id)
        if matches:
            new_key = '>{} {}'.format(matches.group(1), matches.group(2))
        else:
            new_key = '####{}'.format(fasta_id)
        new_dict[new_key] = old_value
    return new_dict

def export_fasta(df, output_dir, output_file_name="X.fasta", filters = {},
                 id_column='gi', seq_column='Sequence'):
    """
    Construct fasta file from dataframe
    INPUT:      dataframe, output directory and filename,
                id and seq columns, filter dict
    OUTPUT:     fasta_dict, fasta file
    """
    if len(filters) > 0:
        for filter_column in filters:
            filter_rule = filters[filter_column]
            df = df[df[filter_column].isin(filter_rule)]
    df = df.copy()
    cols = [id_column, seq_column]
    sub_df = df.loc[:,cols]
    fasta_dict = dict(zip(sub_df[id_column], sub_df[seq_column]))

    final_string = ""    
    for fasta_id in fasta_dict:
        sequence = insert_newlines(fasta_dict[fasta_id], every=64)
        this_string = ">{}\n{}\n".format(fasta_id, sequence)
        final_string += this_string
    
    output_file = output_dir / output_file_name

    with open(output_file, 'w') as out:
        out.write(final_string.strip())
    
    print("Exported %d sequences to %s." % (len(fasta_dict), output_file))
    #print(final_string)
    return 

def replace(string, char):
    """
    Function to replace multiple occurrences   
    of a character by a single character 
    """
    pattern = char + '{2,}'
    string = re.sub(pattern, char, string) 
    return string

def print_basic_info(df, columns=True, size=True, empties=False,
                     duplicates=True, 
                     check_columns=["gene_from_uniprot", 
                                    "gene_from_uniprot_splitted"]):
    if columns:
        print("\nColumns:")
        print(df.columns)
    if size:
        print("\nSize (rowns, columns):")
        print(df.shape)
    if empties:
        print("\nNumber of empty rows per column:")
        print(df.isna().sum())
    if duplicates:
        for col in check_columns:
            check_duplicates(df, column=col)

def filter_dict(freqs, lower_limit=0.01):
    """
    Filter dictionary based on value
    (Filtering by % in relation to total sum)
    """
    new_dict = {}
    total = sum(freqs.values())
    
    for key in freqs:
        porcent = freqs[key] / float(total)
        if porcent > lower_limit:
            new_dict[key] = freqs[key]
        else:
            if '+' not in new_dict:
                new_dict['+'] = freqs[key]
            else:
                new_dict['+'] += freqs[key]
    return new_dict

def separate_pairs(dictionary):
    """
    """
    new_dictionary = {}
    for key in dictionary:
        if len(key) > 1:
            #keys = key.split(sep)
            keys = list(key)
            value = dictionary[key]
            for key in keys:
                if key in new_dictionary:
                    new_dictionary[key] += value
                else:
                    new_dictionary[key] = value
        else:
            new_dictionary[key] = dictionary[key]
    return new_dictionary

def export_sequences(df, output_dir, output_file_name, list_of_proteins, 
                        mapping_column):
    new_df = df[df[mapping_column].str.contains('|'.join(list_of_proteins))]
    export_fasta(new_df, output_dir, output_file_name, id_columns=['query_name',
                 'description', 'method', 'medium_UC','medium_CR'], 
                 seq_column='Sequence', filters = {})
    return dict(zip(new_df.query_name, new_df.Sequence))

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
        
def chunks(data, SIZE=10000):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}

def split_fasta_file(fasta_file, output_dir, par_type=['n_files','n_seqs'][0], par_value=2):
    '''
    
    '''
    
    fastas = fasta_parser(fasta_file, simplify_id=False)
    n_seqs_total = len(fastas)
    if par_type == 'n_files':
        n_seqs = round(n_seqs_total / par_value)
        n_files = par_value
    elif par_type == 'n_seqs':
        n_files = round(n_seqs_total / par_value)
        n_seqs = par_value
    print("n_total = {}, n_seqs = {}, n_files = {}".format(n_seqs_total, n_seqs, n_files))
    
    splitted_dictionaries = []
    for item in chunks(fastas, n_seqs):
        splitted_dictionaries.append(item)
    
    splitted_dir = output_dir / 'splited_fastas'
    for i in range(n_files):
        this_file = str(i) + '.fasta'
        file_path = splitted_dir / this_file
        filtered_dict = filter_sequences(splitted_dictionaries[i], min_seq_len=40)
        fasta_dict2file(filtered_dict, file_path)
    return splitted_dictionaries

def join_csv(csv_list, header=True):
    '''
    '''
    all_lines = []
    for file_path in csv_list:
        with open(file_path, 'r') as input_data:
            data = input_data.readlines()
            if header:
                header_content = data[0]
                file_content = data[1:]
            else:
                header_content = None
                file_content = data
        
        if not all_lines:
            all_lines.append(header_content)
        all_lines.append(file_content)
    
    return all_lines


def simplify_fasta_keys(dictionary, how=('|',1)):
    """
    In the fasta dictionary, change keys with split / selection
    """
    new_dict = {}
    for key in dictionary:
        value = dictionary[key]
        simpler_key = key.split(how[0])[how[1]]
        new_dict[simpler_key] = value
    return new_dict

def change_dict_values(dictionary, join_char='/'):
    new_dict = {}
    for key in dictionary:
        value_list = dictionary[key]
        value_stripped = [i.strip() for i in value_list]
        value_string = join_char.join(value_stripped)
        new_dict[key] = value_string
    return new_dict


        