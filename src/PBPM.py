# -*- coding: utf-8 -*-

"""
STALICLA
PBPM Module (1)
This module contains the functions to create a general protocol to generate patient by disrupted pathways matrices
(PBPMs). It takes two main file types to generate a PBPM matrix according to user parameters passed via command line.
Author: Gabriela Martinez - airamgabriela17@gmail.com
"""

import pandas as pd
from time import gmtime, strftime
import argparse
import os
from pathlib import Path
import logging
from datetime import datetime
import getpass


# Auxiliary functions

def create_directory(folderpath):

    """ Creates a folder.

    Parameters:
        folderpath (str):The path to the folder that will be create.

    Returns:
        create_directory(folderpath):The created folder in the specified path.

    """

    if not os.path.exists(folderpath):
        os.makedirs(folderpath)


def drop_na(df):

    """ Drops records with null or 'NA' genes given by the HGNC_symbol.

    Parameters:
        df (pandas dataframe):Pandas dataframe with the column HGNC_symbol.

    Returns:
        drop_na(df):Pandas dataframe without null or invalid genes

    """

    df = df.dropna(subset=['HGNC_symbol'])
    df = df[pd.notnull(df['HGNC_symbol'])]  # Make sure that null columns are deleted
    df = df[~df['HGNC_symbol'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
    return df


def return_file(filepath):

    """ Returns a file read line by line.

    Parameters:
        filepath (str): folder where the file is located.

    Returns:
        return_file(pathfile): returns readlines() object.

    """

    if filepath != '':  # If I provided a path
        file_name = Path(filepath)
        if os.path.isfile(file_name):
            with open(file_name, 'r') as file:
                lines = file.readlines()
                return lines
    elif filepath == '':
        return ''


def edit_lines(filepath, lines, folderpath):

    """ Returns a file with its lines starting with '#" at the beginning and a single table for the remaining data.

    Parameters:
        filepath (str): folder where the file is located.
        lines (readlines() object): line by line typed file.
        folderpath (str): destination to save the generated output


    Returns:
        edit_lines(filepath, lines, folderpath): exports input as csv file and returns nothing.

    """

    comments = []
    header = ''
    data = []

    if filepath != '':  # If I provided a path
        file_name = Path(filepath)
        if os.path.isfile(file_name):
            for i in range(len(lines)):
                if not lines[i].startswith('#') and not lines[i - 1].startswith('#'):
                    data.append(lines[i])
                elif lines[i].startswith('#'):
                    comments.append(lines[i])
                elif lines[i - 1].startswith('#'):
                    header = lines[i]
            temp = os.path.join(folderpath, '_temp.csv')
            with open(temp, 'w') as f:
                f.writelines(comments)
                f.write(header)
                f.writelines(data)
            os.rename(temp, file_name)
    elif filepath == '':
        pass


def none_to_str(s):

    """ Function to convert a None object to string.

    Parameters:
        s (None): a None Python object.

    Returns:
        none_to_str(s): a string Python object.

    """

    convert = lambda i: i or ''
    if s is None:
        s = convert(s)
        return s
    if s is not None:
        return s


# Core functions


def create_arg_parser():

    """ Creates a parser to read arguments via the command line.

    Parameters:
        parser object from command line.

    Returns:
        create_arg_parser():Args dictionary mapping user inputs.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', help='(Mandatory): Path to the input file containing patient mutations file')
    parser.add_argument('pathwaysDirectory', help='(Mandatory ):Path to the directory that contains the pathway files '
                                                  'with annotated genes')
    parser.add_argument('-pathway', help='Pathway(s) to be extracted. If not given, all by default. '
                                         'A single pathway, a subset of pathways separated '
                                         'by comma without spaces or a filepath (tab text with no headers '
                                         'and a single column) can be provided. '
                                         'Examples: -pathway R-HSA-69620 | -pathway R-HSA-69620,0051705 | '
                                         '-pathway path-to-my-file.txt', default=None)
    parser.add_argument('-gene', help='Gene(s) to be extracted. '
                                      'If not given, all by default. '
                                      'A single gene, a subset of genes separated by comma without spaces or a '
                                      'filepath (tab text with no headers and a single column) can be '
                                      'provided. Examples: -gene CTR9 | -gene CTR9,NOCL2 | '
                                      '-gene path-to-my-file.txt', default=None)
    parser.add_argument('-patient', help='Patient(s) id(s) to be extracted. '
                                         'If not given, all by default. '
                                         'A single patient id, a subset of patients separated '
                                         'by comma without spaces or a filepath (tab text with no headers and a '
                                         'single column) can be provided. '
                                         'Examples: -patient Patient_X | -patient Patient_X,Patient_Y | '
                                         '-patient path-to-my-file.txt', default=None)
    parser.add_argument('-mutation', help='Mutation(s) to be extracted. '
                                          'If not given, all by default. '
                                          'A single mutation, a subset of mutations separated '
                                          'by comma without spaces or a filepath (tab text with no headers and a '
                                          'single column) can be provided. '
                                          'Examples: -mutation missense_variant | -mutation missense_variant,Intron | '
                                          '-mutation path-to-my-file.txt', default=None)
    parser.add_argument('-pli_gt', help='Filters records with values greater than or equal to a pLI threshold.'
                                        'Example: -pli_gt 0.6', default=None)
    parser.add_argument('-af_lt', help='Filters records with values less than a max_control_AF threshold.'
                                       'Example: -af_lt 0.3', default=None)
    parser.add_argument('-pph2', help='Filters records for qualifiers of pph2 predictions.'
                                      'Available options: "benign", "possibly damaging", "probably damaging".'
                                      'Qualifiers must be enclosed with quotes as they hold blank spaces.'
                                      'More than one qualifier is accepted separated by a comma'
                                      'Examples: -pph2 "possibly damaging" | '
                                      '-pph2 "probably damaging,possibly damaging"', default=None)
    parser.add_argument('-mpc_gt', help='Filters records with values greater than or equal to an MPC threshold.'
                                        'Example: -mpc_gt 0.7', default=None)
    parser.add_argument('-adj', help='Filters records by specific value(s) of adjusted consequence.'
                                     'Available options: PTV, Missense3, Missense, etc.'
                                     'More than one adjusted consequence is accepted separated by a comma'
                                     'Example: -adj PTV | -adj Missense3,Missense', default=None)
    parser.add_argument('-matrix', help='Type of matrix to be extracted.'
                                        'Available options: binary(b), numerical(n), normalized(nn).'
                                        'Example: -matrix b | -matrix n | -matrix nn.'
                                        'It is only possible to generate one matrix at a time.'
                                        'By default, the binary matrix is created', default='b')
    parser.add_argument('-r', help='Set this parameter to yes (y) to download a copy of the base matrix (raw data)'
                                   ' generated to build the PBP matrix. '
                                   'By default, this parameter is set to no (n).', default='n')
    parser.add_argument('-append', help='Set this parameter for adding a txt or csv filepath to which base matrix '
                                        'results want to be appended.', default='')
    args = parser.parse_args()
    args_dict = vars(args)
    return args_dict


def read_mutations(filepath):

    """ Reads the mutations file with gene annotations.
        Column names must be respected.

    Parameters:
        filepath (str):The path to the mutations file.

    Returns:
        read_mutations(filepath):Pandas dataframe with mutations file information.

    """

    if filepath:
        input_file = pd.read_csv(filepath, header=0, sep="\t")
        # input_file['MPC'] = np.random.randint(0, 5, size=(len(input_file), 1)) # VALIDAR CON LAURA SI ES CALCULADO
        # input_file['adjusted_consequence'] = input_file.apply(lambda row: adj_consequence(row), axis=1) # VALIDAR
        # input_file.to_csv('mutations-last.csv', index=False)
        return input_file


def extract_pathways(folderpath, pathway):

    """ Extracts directory locations for a subset of pathways to be analyzed.
        The pathway parameter should be pass like: -pathway yourdesiredpathways.
        This function process the pathway argument as follows:

        * if you type '-pathway PathwayName' only information for this pathway will be computed
        * if you type '-pathway PathwayName1,PathwayName2' pathways with those two names will be computed
        * if you type nothing, by default, all the pathways will be computed

    Parameters:
        folderpath (str):The path to the folder containing the pathways files.
        pathway (str): single string or collection of paths to be analyzed.

    Returns:
        read_mutations(folderpath, pathway):List or string containing the pathways files' paths to process

    """

    if folderpath and pathway:
        if os.path.isfile(pathway):
            file_name = Path(pathway)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            files = []
            for i in file_data:
                file_name = '{}_with_gene_annotations'.format(i)
                path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
                files.append(path_to_find)
            return files
        elif ',' not in pathway:
            files = []
            file_name = '{}_with_gene_annotations'.format(pathway)
            path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
            files.append(path_to_find)
            return files
        elif ',' in pathway:
            paths = pathway.split(',')
            files = []
            for i in paths:
                file_name = '{}_with_gene_annotations'.format(i)
                path_to_find = os.path.join(folderpath, file_name + '.' + 'txt')
                files.append(path_to_find)
            return files
    elif folderpath and pathway is None:
        files = []
        for file in os.listdir(folderpath):
            if file.endswith('_with_gene_annotations.txt'):
                files.append(os.path.join(folderpath, file))
        return files


def read_pathways(folderpath, pathway):

    """ Reads a subset of pathways files.

    Parameters:
        folderpath (str):The path to the folder containing the pathways files.
        pathway (str): single string or collection of paths to be analyzed.

    Returns:
        read_pathways(folderpath, pathway):Pandas dataframe containing the pathways files' joint data

    """

    frames = []
    paths_to_find = extract_pathways(folderpath, pathway)
    for path in paths_to_find:
        g = pd.read_csv(path, header=0, sep='\t', engine='python')
        g['Pathway'] = '{}'.format(path.split('_with_gene_annotations')[0].split('/')[-1])
        g.columns = ['Gene', 'HGNC_symbol', 'HGNC_mapping', '%RVIS_ESP_0.1%', '%RVIS_ExAC_0.01%',
                     '%RVIS_ExAC_0.05%popn', '%RVIS_ExAC_0.1%popn', 'constraint_score', 'pLI', 'pRecessive', 'pNull',
                     'pHI', 'chrN:N-N', 'Pathway']
        g = g[['Pathway', 'Gene', 'HGNC_symbol', 'HGNC_mapping', '%RVIS_ESP_0.1%', '%RVIS_ExAC_0.01%',
               '%RVIS_ExAC_0.05%popn', '%RVIS_ExAC_0.1%popn', 'constraint_score', 'pLI', 'pRecessive', 'pNull', 'pHI',
               'chrN:N-N']]
        frames.append(g)
    data = pd.concat(frames, ignore_index=True)
    data = data.sort_values(by='Pathway', ascending=True)
    return data


def genes_to_normalize_pathways(df):

    """ Extracts the count of unique genes per pathway. This will be used later to normalize the numerical matrix

    Parameters:
        df (pandas dataframe):Pandas dataframe containing joint information of pathways.

    Returns:
        genes_to_normalize_pathways(df):Pandas dataframe relating pathways and unique count of genes

    """

    pivot_df = pd.pivot_table(df, index='Pathway', values='HGNC_symbol', aggfunc=lambda x: len(x.unique()))
    flattened = pd.DataFrame(pivot_df.to_records())
    #flattened.to_csv('all_genes_to_norm.csv', index=False)
    return flattened


def join_input_data(df1, df2):

    """ Joins data coming from mutations file and pathways, on HGNC_symbol column.

    Parameters:
        df1 (pandas dataframe):Pandas dataframe with mutations and patients data.
        df2 (pandas dataframe):Pandas dataframe with the pathways data.

    Returns:
        join_input_data(df1, df2):joint pandas dataframe.

    """

    pathways = df2[['Pathway', 'HGNC_symbol']]
    joint_data = pd.merge(df1, pathways, on=['HGNC_symbol'])
    return joint_data


def create_base_matrix(df):

    """ Creates base data matrix containing summarized information per patient and mutation (given by
    chromosome, position, reference and alteration).

    Parameters:
        df (pandas dataframe):Pandas dataframe with mutations, patients and pathways data.

    Returns:
        create_base_matrix(df):joint pandas dataframe.

    """

    pivot_table = pd.pivot_table(df,
                                 index=['child_id', 'Chr', 'Position', 'Ref', 'Alt', 'consequence'],
                                 values=['HGNC_symbol', 'Pathway'],
                                 aggfunc=lambda x: ', '.join(set(x)))
    flattened = pd.DataFrame(pivot_table.to_records())
    return flattened


def download_matrix(setting, df, info, csv_name):

    """ Downloads a pandas dataframe into a csv file.

    Parameters:
        setting (str): if set to 'y' this function will be executed. If set 'n', the default value, it will not.
        df (pandas dataframe):Pandas dataframe to be downloaded.
        info (str): the command line that triggered the creation of the matrix to download.
        csv_name (str): destination where the dataframe will be stored.

    Returns:
        download_base_matrix(setting, df, info, folderpath)): resulting csv file.

    """

    if setting:
        if setting == 'y':
            file_name = csv_name
            with open(file_name, 'w') as file:
                file.write(info + '\n')
                df.to_csv(file, index=False)
                return df
        elif setting == 'n':
            return df
    elif setting is None:
        return df


def append_base_matrices(df, info, filepath):

    """ Appends a pandas dataframe to a new or existing csv file containing same type of data.

    Parameters:
        df (pandas dataframe):Pandas dataframe to be appended.
        info (str): the command line that triggered the creation of the pandas dataframe to be appended.
        filepath (str): destination file where the dataframe will be appended/stored.

    Returns:
        append_base_matrices(df, info, pathfile)): resulting csv file with appended data.

    """

    if filepath != '':  # If I provided a path
        file_name = Path(filepath)
        if os.path.isfile(file_name):
            with open(filepath, 'a') as file:
                file.write(info + '\n')
                df.to_csv(file, index=False)
                return file_name
        elif not os.path.isfile(file_name):
            with open(filepath, 'w') as file:
                file.write(info + '\n')
                df.to_csv(file, index=False)
                return file_name
    elif filepath == '':
        return ''


def analyze_patients(df, patient_id):

    """ Filters a matrix for a set of patients.
    The patient parameter should be pass like: -patient Patient_XXXX
    In case more than one patient wants to be included, type: -patient Patient_XXXX, Patient_YYYY, Patient_ZZZZ
    This function also accepts a path to a file of patients (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        patient_id (str): single id or subset of ids to be filtered.

    Returns:
        analyze_patients(df, patient_id): pandas dataframe object with filter applied.

    """

    df = df[pd.notnull(df['child_id'])]  # Make sure that null columns are deleted
    df = df[~df['child_id'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python

    if patient_id:
        if os.path.isfile(patient_id):
            file_name = Path(patient_id)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['child_id'].isin(file_data)]
            return df
        else:
            if ',' in patient_id:
                patients = patient_id.split(',')
                df = df[df['child_id'].isin(patients)]
                return df
            elif ',' not in patient_id:
                df = df[df.child_id == patient_id]
                return df
    elif patient_id is None:
        return df


def analyze_genes(df, gene_id):

    """ Filters a matrix for a set of genes.
    The gene parameter should be pass like: -gene XXXXX
    In case more than one gene wants to be included, type: -gene XXXXX,YYYYY,ZZZZZ
    This function also accepts a path to a file of genes (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        gene_id (str): single id or subset of ids to be filtered.

    Returns:
        analyze_genes(df, gene_id): pandas dataframe object with filter applied.

    """

    if gene_id:
        if os.path.isfile(gene_id):
            file_name = Path(gene_id)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['HGNC_symbol'].isin(file_data)]
            return df
        else:
            if ',' in gene_id:
                patients = gene_id.split(',')
                df = df[df['HGNC_symbol'].isin(patients)]
                return df
            elif ',' not in gene_id:
                df = df[df.HGNC_symbol == gene_id]
                return df
    elif gene_id is None:
        return df


def analyze_mutations(df, mutation):

    """ Filters a matrix for a set of mutations.
    The mutation parameter should be pass like: -mutation desiredmutation
    In case more than one mutation wants to be included, type: -mutation desiredmutation1,desiredmutation2
    This function also accepts a path to a file of mutations (txt tab delimited with no headers)

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        mutation (str): subset of mutations to be filtered.

    Returns:
        analyze_mutations(df, mutation): pandas dataframe object with filter applied.

    """

    if mutation:
        if os.path.isfile(mutation):
            file_name = Path(mutation)
            file_data = pd.read_csv(file_name, header=None, sep="\t")
            file_data = file_data.iloc[:, 0].tolist()
            df = df[df['consequence'].isin(file_data)]
            return df
        else:
            if ',' in mutation:
                patients = mutation.split(',')
                df = df[df['consequence'].isin(patients)]
                return df
            elif ',' not in mutation:
                df = df[df.consequence == mutation]
                return df
    elif mutation is None:
        return df


def analyze_pli(df, pli):

    """ Filters a matrix by pLI level.
    The parameter should be pass like: pLI desiredpLI
    Results will be shown for mutations with pLI greater or equal than the specified threshold.

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        pli (float): pLI threshold.

    Returns:
        analyze_pli(df, pli_gt): pandas dataframe object with filter applied.

    """

    if pli:
        pli = float(pli)
        df = df[pd.notnull(df['pLI'])]  # Make sure that null columns are deleted
        df = df[~df['pLI'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        df = df[df.pLI >= pli]
        return df
    elif pli is None:
        return df


def analyze_max_control_af(df, max_control_af):

    """ Filters a matrix by max_control_af frequency.
    The parameter should be pass like: -af desiredfrequency
    Results will be shown for mutations with max_control_AF less than the specified threshold.

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        max_control_af (float): max_control_af threshold.

    Returns:
        analyze_max_control_af(df, max_control_af): pandas dataframe object with filter applied.

    """

    if max_control_af:
        paf = float(max_control_af)
        df = df[pd.notnull(df['max_control_AF'])]  # Make sure that null columns are deleted
        df = df[~df['max_control_AF'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        df = df[df.max_control_AF < paf]
        return df
    elif max_control_af is None:
        return df


def analyze_pph2_prediction(df, pph2):

    """ Filters a matrix by a subset of pph2 labels.
    The parameter should be pass like: -pph2 qualifier
    Example: -pph2 possibly damaging
             -pph2 "possibly damaging, probably damaging"

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        pph2 (str): subset of labels to be filtered.

    Returns:
        analyze_pph2_prediction(df, pph2): pandas dataframe object with filter applied.

    """

    if pph2:
        df = df[pd.notnull(df['pph2_prediction'])]  # Make sure that null columns are deleted
        df = df[~df['pph2_prediction'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        if ',' in pph2:
            pph2_labels = pph2.split(',')
            df = df[df['pph2_prediction'].isin(pph2_labels)]
            return df
        elif ',' not in pph2:
            df = df[df.pph2_prediction == pph2]
            return df
    elif pph2 is None:
        return df


def analyze_mpc(df, mpc):

    """ Filters a matrix by an mpc value.
    The parameter should be pass like: -mpc desiredmpc
    Results will be shown for mutations with mpc greater or equal to the specified threshold.

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        mpc (float): mpc threshold.

    Returns:
        analyze_mpc(df, mpc): pandas dataframe object with filter applied.

    """

    if mpc:
        mpc = float(mpc)
        df = df[pd.notnull(df['MPC'])]  # Make sure that null columns are deleted
        df = df[~df['MPC'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        df = df[df.MPC >= mpc]
        return df
    elif mpc is None:
        return df


def analyze_adj_consequence(df, adj_cons):

    """ Filters a matrix by a subset of adjusted consequence labels.
    The parameter should be pass like: -adj desiredcons
    Example: -adj Missense3
             -adj Missense3, PTV

    Parameters:
        df (pandas dataframe): matrix to be filtered.
        adj_cons (str): subset of labels to be filtered.

    Returns:
        analyze_adj_consequence(df, adj_cons): pandas dataframe object with filter applied.

    """

    if adj_cons:
        df = df[pd.notnull(df['adjusted_consequence'])]  # Make sure that null columns are deleted
        df = df[~df['adjusted_consequence'].astype(str).str.contains('NA')]  # Delete NAs not recognized by Python
        if ',' in adj_cons:
            cons_labels = adj_cons.split(',')
            df = df[df['adjusted_consequence'].isin(cons_labels)]
            return df
        elif ',' not in adj_cons:
            df = df[df.adjusted_consequence == adj_cons]
            return df
    elif adj_cons is None:
        return df


def create_pbpm(setting, df, df_to_norm):

    """ Extracts PBP matrices of three types: binary, numerical and normalized.
        -matrix 'b': binary matrix
        -matrix 'n': numerical matrix
        -matrix 'nn': normalized matrix
    If the parameter -matrix is not set, it takes the argument 'b' by default.

    available_mutations = ['3UTR', '5UTR', 'frameshift_variant', 'InFrame', 'Intron', 'missense_variant',
                           'Near_3splice', 'Near_5splice', 'Promoter', 'Silent', 'splice_acceptor_variant',
                           'splice_donor_variant', 'StartCodonDisrupted', 'stop_gained']

    Parameters:
        setting (str): specifies the type of matrix to extract.
        df (pandas dataframe): input data for the matrix.
        df_to_norm (pandas dataframe): input data for the normalization.

    Returns:
        create_pbpm(setting, df, df_to_norm): pandas dataframe.

    """

    if setting == 'b':
        pivot = pd.concat([df.drop('Pathway', 1), pd.get_dummies(df.Pathway).mul(1)], axis=1)
        flattened = pd.DataFrame(pivot.to_records())
        filter_cols = [col for col in flattened if (col.startswith('path') or col.startswith('child'))]
        filtered = flattened[filter_cols]
        final_matrix = filtered.sort_values(by='child_id', ascending=True)
        final_matrix = final_matrix.groupby('child_id').max()
        final_matrix = pd.DataFrame(final_matrix.to_records())
        return final_matrix
    elif setting == 'n':
        pivot = df.pivot_table(values='consequence', index='child_id', columns='Pathway', aggfunc='count', fill_value=0)
        flattened = pd.DataFrame(pivot.to_records())
        final_matrix = flattened.sort_values(by='child_id', ascending=True)
        return final_matrix
    elif setting == 'nn':
        pivot = df.pivot_table(values='consequence', index='child_id', columns='Pathway', aggfunc='count', fill_value=0)
        flattened = pd.DataFrame(pivot.to_records())
        p = list(flattened.columns.values)[1:]
        df_to_norm = df_to_norm[df_to_norm['Pathway'].isin(p)]
        final_matrix = flattened.set_index('child_id').div(df_to_norm.set_index('Pathway')['HGNC_symbol']).reset_index()
        return final_matrix


def download_pbpm(setting, df, info, b_csv_name, n_csv_name, nn_csv_name):

    """ Downloads a PBP matrix.

    Parameters:
        setting (str): specifies the type of matrix to extract.
        df (pandas dataframe): data to be downloaded.
        info (str): the command line that triggered the creation of the dataframe to download.

    Returns:
        download_pbpm(setting, df, info, folderpath):  csv file.

    """

    if setting == 'b':
        with open(b_csv_name, 'w') as file:
            file.write(info + '\n')
            df.to_csv(file, index=False)
    elif setting == 'n':
        with open(n_csv_name, 'a') as file:
            file.write(info + '\n')
            df.to_csv(file, index=False)
    elif setting == 'nn':
        with open(nn_csv_name, 'a') as file:
            file.write(info + '\n')
            df.to_csv(file, index=False)


def logging_info(args, folderpath, folderpath2):

    """ Records an INFO type log.

    Parameters:
        args (dictionary): parsed arguments from command line.
        folderpath (str): destination of the info log.
        folderpath2 (str): folder where the final PBP matrix (and base matrix if applies) are downloaded.

    Returns:
        logging_info(args, folderpath): info log.

    """

    is_base = "Base matrix downloaded at: {}".format(folderpath2) if args['r'] == 'y' else ''
    is_append = "Append filepath for base matrix (if provided): {}".format(args['append']) if args['append'] else ''

    log_filename = os.path.join(folderpath, 'logINFO_{}.log'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    logging.basicConfig(filename=log_filename, level=logging.INFO)

    logging.info('PBPM protocol started generation on: ' + datetime.now().strftime('%d-%m-%Y %H:%M:%S'))
    logging.info('Process initiated by the user ' + getpass.getuser())
    logging.info('\n' + '\n' +
                 "Script generated with the following parameters --- " + '\n' +
                 " inputFile: " + args['inputFile'] + '\n' +
                 " pathwaysDirectory: " + args['pathwaysDirectory'] + '\n' +
                 "-pathway: " + none_to_str(args['pathway']) + '\n' +
                 "-gene: " + none_to_str(args['gene']) + '\n' +
                 "-patient: " + none_to_str(args['patient']) + '\n' +
                 "-mutation: " + none_to_str(args['mutation']) + '\n' +
                 "-pli (greater/equal than): " + none_to_str(args['pli_gt']) + '\n' +
                 "-af (less than): " + none_to_str(args['af_lt']) + '\n' +
                 "-pph2: " + none_to_str(args['pph2']) + '\n' +
                 "-mpc (greater/equal than): " + none_to_str(args['mpc_gt']) + '\n' + '\n' +
                 "Type of PBP matrix generated: " + none_to_str(args['matrix']) + '\n' +
                 "PBP matrix downloaded at: " + folderpath2 + '\n' + '\n' +
                 "Base matrix downloaded? " + none_to_str(args['r']) + '\n' +
                 is_base + is_append)

    print("---")
    print("Protocol successfully generated!")
