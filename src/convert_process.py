#'/usr/bin/env python3.10.5
"""
converts ANNOVAR files to UMAP input format.
"""
import numpy as np
import pandas as pd
import os
import zipfile
import glob
import re
import json
from functools import reduce

def convert(input_path, patient_info_file, output_csv):
    """
    Extracts ANNOVAR files, processes for each patient, and outputs a CSV for UMAP.
    """
    if zipfile.is_zipfile(input_path):
        zip_file_dir = input_path.rsplit(".", 1)[0]
        if not os.path.exists(zip_file_dir):
            os.makedirs(zip_file_dir)
        with zipfile.ZipFile(input_path, 'r') as zip_ref:
            zip_ref.extractall(zip_file_dir)
        files_path = zip_file_dir
    else:
        files_path = input_path

    # Find all relevant files
    files = sorted(glob.glob(f"{files_path}/*"), key=str.lower)
    
    # Load patients information
    patients_info = pd.read_csv(patient_info_file)
    
    # Process files for each patient
    stats_all = []
    for patient_name in sorted(list(patients_info[patients_info.columns[0]])):
        patient_files = [file for file in files if re.search(str(patient_name), file)]
        gene_count = [process_patient_file(patient_file) for patient_file in patient_files]
        
        # Merge counts for all files of a patient
        annovar_stat = reduce(lambda left, right: pd.merge(left, right, on='Gene.refGene', how='outer'), gene_count).fillna(0)
        stats_all.append(annovar_stat)
    
    # Combine all patient stats
    annovar_stat_all = reduce(lambda left, right: pd.merge(left, right, on='Gene.refGene', how='outer'), stats_all).fillna(0)
    
    # Add more info to sample names
    annovar_stat_all, categorical_columns = enhance_sample_names(annovar_stat_all, patients_info)
    
    # Convert float values to int where appropriate
    annovar_stat_all = convert_float_columns(annovar_stat_all)
    
    # Save to CSV
    annovar_stat_all.to_csv(output_csv, index=False)
    
    with open('categorical_columns.json', 'w') as f:
    	json.dump('_'.join(categorical_columns), f)

def process_patient_file(patient_file):
    """
    Processes a single patient file for nonsynonymous SNV counts.
    """
    df = pd.read_csv(patient_file, sep="\t", header=0, quotechar='"', low_memory=False)
    df = df[df['ExonicFunc.refGene'] == "nonsynonymous SNV"]
    return df.groupby('Gene.refGene').size().reset_index(name=os.path.basename(patient_file).split('.')[0])

def enhance_sample_names(annovar_stat_all, patients_info):
    """
    Replaces sample names with more descriptive identifiers.
    """
    ##record the names of all categorical columns in patients information
    patient_id_col = patients_info.columns[0]
    non_numerical_columns = patients_info.select_dtypes(exclude=['int64', 'float64']).columns.tolist()
    categorical_columns = [patient_id_col] + [col for col in non_numerical_columns if col != patient_id_col]
    categorical_df = patients_info[categorical_columns]
    # Concatenate all categorical columns to form a new 'info' column
    patients_info['info'] = categorical_df.apply(lambda x: '_'.join(x.astype(str)), axis=1) + "_"
    
    for i, patient_id in enumerate(patients_info[patient_id_col]):
        info = patients_info.loc[i, 'info']
        # Replace annovar_stat_all column names that contain this patient_id
        annovar_stat_all.columns = [col.replace(str(patient_id), info) if str(patient_id) in col else col for col in annovar_stat_all.columns]
    
    return annovar_stat_all, categorical_columns

def convert_float_columns(df):
    """
    Converts float columns to int where all values are integers.
    """
    for col in df.select_dtypes(include=['float']).columns:
        if all(df[col].dropna().apply(float.is_integer)):
            df[col] = df[col].astype(int)
    return df
