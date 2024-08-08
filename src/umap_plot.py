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
import matplotlib.pyplot as plt
matplotlib.use('Agg')

output_dir = "UMAP_outputs"
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir)
#suppress warnings
warnings.filterwarnings("ignore")

def umap_plot(umap_input_file, category_info, min_cells=3, min_genes=30, n_top_genes=3000, n_neighbors=60, n_pcs=40, leiden_resolution=1.5, louvain_resolution=None, no_plot_venn=False, venn='type'):
    # Load data and transpose to have genes as columns
    adata = sc.read_csv(umap_input_file).T
    sc.settings.verbosity = 'error'
    # Preprocessing
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # Plot distribution of number of mutated cells per gene
    mut_cells_per_gene = np.sum(adata.X > 0, axis=0)
    plt.figure(figsize=(10, 6))
    plt.hist(mut_cells_per_gene, bins=50, color='skyblue', range=(0, np.quantile(mut_cells_per_gene, 0.999)))
    plt.axvline(x=min_cells, color='red', linestyle='dashed', label=f'Cutoff at {min_cells}')
    plt.xlabel('Number of mutated cells per gene')
    plt.ylabel('Frequency')
    plt.title('Distribution of mutated cells per gene')
    plt.legend()
    plt.savefig(f"{output_dir}/mut_cells_per_gene_distribution.png")
    plt.close()
    # Plot distribution of number of mutated genes per cell
    mut_genes_per_cell = np.sum(adata.X > 0, axis=1)
    plt.figure(figsize=(10, 6))
    plt.hist(mut_genes_per_cell, bins=50, color='lightgreen', range=(0, np.quantile(mut_genes_per_cell, 0.999)))
    plt.axvline(x=min_genes, color='red', linestyle='dashed', label=f'Cutoff at {min_genes}')
    plt.xlabel('Number of mutated genes per cell')
    plt.ylabel('Frequency')
    plt.title('Distribution of mutated genes per cell')
    plt.legend()
    plt.savefig(f"{output_dir}/mut_genes_per_cell_distribution.png")
    plt.close()
    
    sc.pp.filter_genes(adata, min_cells)
    sc.pp.filter_cells(adata, min_genes)
    
    metadata_fields = category_info.split('_') 
    for i, field in enumerate(metadata_fields):
        adata.obs[field] = adata.obs.index.str.split('_').map(lambda x: x[i] if len(x) > i else 'unknown')
    # Output adata.obs details
    adata.obs.to_csv(f"{output_dir}/adata_obs_details.txt", sep='\t')
    
    # Patient summary
    patient_summary(adata, metadata_fields)
    
    # Violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, save='_qc_metrics.png')
    
    # Filter outliers based on QC
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    
    # Violin plot for number of genes with at least 1 count in a cell 
    plt.figure(figsize=(8, 10))
    sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, show=False)
    plt.title('Violin plot of the number of genes with at least 1 count in a cell ')
    plt.axhline(y=upper_lim, color='r', linestyle='--', label=f'Cutoff at {round(upper_lim,2)}')
    plt.legend(loc='upper right')
    plt.savefig(f"{output_dir}/violin_n_genes_by_counts_per_cell.png")
    plt.close()
    
    adata = adata[adata.obs.n_genes_by_counts.values < upper_lim]
    adata.obs.columns = adata.obs.columns.str.lower()    

    # Check if the 'venn' field is a valid metadata field
    if venn not in [field_info.lower() for field_info in metadata_fields]:
        raise ValueError(f"The specified Venn category '{venn}' is not valid. Please select one of the following category: {', '.join(metadata_fields)}.")
    
    # Gene overlap analysis and plotting
    if no_plot_venn==False and venn:
    	gene_overlap = calculate_gene_overlap(adata, metadata_fields, venn)
    	if len(gene_overlap) > 4:
    		sorted_keys = sorted(gene_overlap, key=lambda k: len(gene_overlap[k]), reverse=True)
    		sets = {key: gene_overlap[key] for key in sorted_keys[:4]}
    		venny4py(sets=sets)
    		print("Warning: more than 4 groups detected. Only the top 4 groups will be plotted in the Venn diagram.")
    	else:
    		plot_venn_diagram(gene_overlap)    
    
    # Normalization and log transformation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata, show=False)
    figures = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figures:
        for ax in fig.axes:
            ax.set_xlabel('mean mutations of genes')  # Set new x-axis label for each subplot
    plt.savefig(f"{output_dir}/highly_variable_genes.png")
    plt.close()
    
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
    sc.tl.louvain(adata, resolution=louvain_resolution)
    sc.tl.umap(adata)
    
    # UMAP plots
    sc.pl.umap(adata, color='leiden', save='_leiden.pdf')
    sc.pl.umap(adata, color='louvain', save='_louvain.pdf')     
    for category in metadata_fields:
    	if category == metadata_fields[0]:
    		sc.pl.umap(adata, color=category.lower(), legend_loc='on data', save=f"_{category}.pdf")
    	else:
    		sc.pl.umap(adata, color=category.lower(), save=f"_{category}.pdf")
    
def patient_summary(adata, metadata_fields):
    summary_path = f"{output_dir}/total_count_mutated_genes.txt"
    with open(summary_path, 'w') as summary_file:
        for i in adata.obs[metadata_fields[0]].unique():
            patient_data = adata.obs[adata.obs[metadata_fields[0]] == i]
            mean_count = patient_data.total_counts.mean()
            min_count = int(patient_data.total_counts.min())
            max_count = int(patient_data.total_counts.max())
            summary_line = f"The average total count of mutated genes for {metadata_fields[0]} {i} is {round(mean_count, 2)}, with a minimum of {min_count} and a maximum of {max_count}.\n"
            summary_file.write(summary_line)
            
def calculate_gene_overlap(adata, metadata_fields, venn):
    gene_overlap = {}
    overlap_info = f"{output_dir}/mutated_gene_info.txt"
    with open(overlap_info, "w") as overlap_file:
        unique_categories = adata.obs[venn].unique()
        # Only retain the genes mutated in at least two cells for each patient and each cell type
        for category in unique_categories:
            adata_subset = adata[adata.obs[venn] == category]
            gene_overlap[category] = set()
            for patient_id in adata_subset.obs[metadata_fields[0]].unique():
                adata_patient = adata_subset[adata_subset.obs[metadata_fields[0]] == patient_id]
                gene_counts = adata_patient.X > 0
                mutated_genes = adata_patient.var_names[gene_counts.sum(axis=0) > 1]
                overlap_file.write(f"{metadata_fields[0]} {patient_id} has {len(mutated_genes)} mutated genes in at least two cells.\n")
                gene_overlap[category].update(set(mutated_genes))
            overlap_file.write(f"{metadata_fields[0]} with the same {venn} of {category} have {len(gene_overlap[category])} mutated genes.\n")
    return gene_overlap

def plot_venn_diagram(gene_overlap):
    # Assuming venny4py is already correctly set up for your environment
    # Prepare sets for the Venn diagram
    sets = {key: gene_overlap[key] for key in sorted(gene_overlap.keys())[:len(gene_overlap)]}  # Adjust based on actual data
    venny4py(sets=sets)
            
