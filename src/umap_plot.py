#'/usr/bin/env python3.10.5
"""
Plot UMAP and output related info
"""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import shutil
import warnings
from venny4py.venny4py import venny4py
import matplotlib
matplotlib.use('Agg')

output_dir = "UMAP_outputs"
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir)
#suppress warnings
warnings.filterwarnings("ignore")

def umap_plot(umap_input_file, category_info, min_cells=3, min_genes=30, n_top_genes=3000, n_neighbors=60, n_pcs=40, leiden_resolution=1.5, plot_venn=True, venn='type'):
    # Load data and transpose to have genes as columns
    adata = sc.read_csv(umap_input_file).T
    sc.settings.verbosity = 'error'
    # Preprocessing
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adata, min_cells)
    sc.pp.filter_cells(adata, min_genes)
    
    metadata_fields = category_info.split('_') 
    for i, field in enumerate(metadata_fields):
        adata.obs[field] = adata.obs.index.str.split('_').map(lambda x: x[i] if len(x) > i else 'unknown')
    
    # Output adata.obs details
    adata.obs.to_csv(f"{output_dir}/adata_obs_details.txt", sep='\t')
    
    # Patient summary
    patient_summary(adata)
    
    # Violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save='_qc_metrics.png')
    
    # Filter outliers based on QC
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts.values < upper_lim]
    
    # Gene overlap analysis and plotting
    if plot_venn and venn:
    	gene_overlap = calculate_gene_overlap(adata, venn)
    	if len(gene_overlap) > 4:
        	raise ValueError("Error: Cannot plot a Venn diagram for more than 4 groups. Exiting.")
    	else:
    		plot_venn_diagram(gene_overlap)    
    
    # Normalization and log transformation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var.highly_variable]
    
    # Regression and scaling
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
    
    # PCA and neighbors
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='_pca_variance_ratio.png')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Clustering and UMAP
    sc.tl.leiden(adata, resolution=leiden_resolution)
    sc.tl.umap(adata)
    
    # UMAP plots
    sc.pl.umap(adata, color='leiden', save='_leiden.pdf')    
    for category in metadata_fields:
    	if category == "Patient":
    		sc.pl.umap(adata, color=category, legend_loc='on data', save=f"_{category}.pdf")
    	else:
    		sc.pl.umap(adata, color=category, save=f"_{category}.pdf")
    
def patient_summary(adata):
    summary_path = f"{output_dir}/patient_total_count_mutated_genes.txt"
    with open(summary_path, 'w') as summary_file:
        for i in adata.obs.Patient.unique():
            patient_data = adata.obs[adata.obs.Patient == i]
            mean_count = patient_data.total_counts.mean()
            min_count = int(patient_data.total_counts.min())
            max_count = int(patient_data.total_counts.max())
            summary_line = f"The average total count of mutated genes for Patient {i} is {round(mean_count, 2)}, with a minimum of {min_count} and a maximum of {max_count}.\n"
            summary_file.write(summary_line)
            
def calculate_gene_overlap(adata, venn):
    gene_overlap = {}
    overlap_info = f"{output_dir}/mutated_gene_info.txt"
    with open(overlap_info, "w") as overlap_file:
        unique_categories = adata.obs[venn].unique()
        
        for category in unique_categories:
            adata_subset = adata[adata.obs[venn] == category]
            gene_overlap[category] = set()
            for patient_id in adata_subset.obs.Patient.unique():
                adata_patient = adata_subset[adata_subset.obs['Patient'] == patient_id]
                gene_counts = adata_patient.X > 0
                mutated_genes = adata_patient.var_names[gene_counts.sum(axis=0) > 1]
                overlap_file.write(f"Patient {patient_id} has {len(mutated_genes)} mutated genes.\n")
                gene_overlap[category].update(set(mutated_genes))
            overlap_file.write(f"Patients with the same {venn} of {category} have {len(gene_overlap[category])} overlap mutated genes.\n")
    return gene_overlap

def plot_venn_diagram(gene_overlap):
    # Assuming venny4py is already correctly set up for your environment
    # Prepare sets for the Venn diagram
    sets = {key: gene_overlap[key] for key in sorted(gene_overlap.keys())[:len(gene_overlap)]}  # Adjust based on actual data
    venny4py(sets=sets)
            