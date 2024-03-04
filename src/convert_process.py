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
    
    # Load patient stage information
    patients_stage = pd.read_csv(patient_info_file)
    
    # Process files for each patient
    stats_all = []
    for patient_name in sorted(list(patients_stage["Patient"])):
        patient_files = [file for file in files if re.search(str(patient_name), file)]
        gene_count = [process_patient_file(patient_file) for patient_file in patient_files]
        
        # Merge counts for all files of a patient
        annovar_stat = reduce(lambda left, right: pd.merge(left, right, on='Gene.refGene', how='outer'), gene_count).fillna(0)
        stats_all.append(annovar_stat)
    
    # Combine all patient stats
    annovar_stat_all = reduce(lambda left, right: pd.merge(left, right, on='Gene.refGene', how='outer'), stats_all).fillna(0)
    
    # Add more info to sample names
    annovar_stat_all = enhance_sample_names(annovar_stat_all, patients_stage)
    
    # Convert float values to int where appropriate
    annovar_stat_all = convert_float_columns(annovar_stat_all)
    
    # Save to CSV
    annovar_stat_all.to_csv(output_csv, index=False)

def process_patient_file(patient_file):
    """
    Processes a single patient file for nonsynonymous SNV counts.
    """
    df = pd.read_csv(patient_file, sep="\t", header=0, quotechar='"', low_memory=False)
    df = df[df['ExonicFunc.refGene'] == "nonsynonymous SNV"]
    return df.groupby('Gene.refGene').size().reset_index(name=os.path.basename(patient_file).split('.')[0])

def enhance_sample_names(annovar_stat_all, patients_stage):
    """
    Replaces sample names with more descriptive identifiers.
    """
    patients_stage['info'] = patients_stage['Patient'].astype(str) + "_" + patients_stage['stage'].astype(str) + "_" + patients_stage['type'].astype(str) + "_" + patients_stage['status'].str.replace(r"\(.*\)", "", regex=True).astype(str) + "_"
    
    for i, col_name in enumerate(annovar_stat_all.columns):
        if i == 0: continue  # Skip the gene column
        pattern = col_name[:5]
        matches = patients_stage['info'][patients_stage['info'].str.contains(pattern)]
        if not matches.empty:
            new_col_name = col_name.replace(pattern, matches.iloc[0], 1)
            annovar_stat_all.rename(columns={col_name: new_col_name}, inplace=True)
    return annovar_stat_all

def convert_float_columns(df):
    """
    Converts float columns to int where all values are integers.
    """
    for col in df.select_dtypes(include=['float']).columns:
        if all(df[col].dropna().apply(float.is_integer)):
            df[col] = df[col].astype(int)
    return df
