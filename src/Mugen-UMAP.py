#'/usr/bin/env python3.10.5
import os
import argparse
from convert_process import convert
from umap_plot import umap_plot

def convert_annovar(input_path, patient_info_file, output_csv):
    """
    Converts ANNOVAR files to UMAP input format.
    """
    convert(input_path, patient_info_file, output_csv)

def umap(umap_input_file, min_cells, min_genes, n_top_genes, n_neighbors, n_pcs, leiden_resolution, plot_venn):
    """
    Plots UMAP figure and generates summary and figures.
    """
    umap_plot(umap_input_file, min_cells, min_genes, n_top_genes, n_neighbors, n_pcs, leiden_resolution, plot_venn)

def all_in_one(input_path, patient_info_file, min_cells, min_genes, n_top_genes, n_neighbors, n_pcs, leiden_resolution, plot_venn):
    """
    Combines ANNOVAR file conversion and UMAP plotting in one step.
    """
    output_csv = os.path.join(os.path.dirname(input_path), "umap_format.csv")
    convert_annovar(input_path, patient_info_file, output_csv)
    umap_plot(output_csv, min_cells, min_genes, n_top_genes, n_neighbors, n_pcs, leiden_resolution, plot_venn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ANNOVAR files and generate UMAP plots.")
    subparsers = parser.add_subparsers(dest="command")

    # Convert command
    parser_convert = subparsers.add_parser('convert', help='Convert ANNOVAR files to UMAP input format.')
    parser_convert.add_argument("-i", "--input", type=str, help="Input ANNOVAR zip file or ANNOVAR directory directly")
    parser_convert.add_argument("-p", "--patient", type=str, help="Input patient info file, including at least the info of Patient, stage, status, and type (e.g., histology type).")
    parser_convert.add_argument("-o", "--output", type=str, help="Output UMAP CSV format")

    # Umap_plot command
    parser_umap = subparsers.add_parser('umap', help='Plot UMAP figure from preprocessed data.')
    parser_umap.add_argument("-i",'--input_umap', type=str, help='Path to UMAP input CSV file.')
    parser_umap.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells for filtering (default=3).')
    parser_umap.add_argument('--min_genes', type=int, default=30, help='Minimum number of genes for filtering (default=30).')
    parser_umap.add_argument('--n_top_genes', type=int, default=3000, help='Number of top genes (default=3000).')
    parser_umap.add_argument('--n_neighbors', type=int, default=60, help='Number of neighbors for UMAP (default=60).')
    parser_umap.add_argument('--n_pcs', type=int, default=40, help='Number of principal components (default=40).')
    parser_umap.add_argument('--leiden_resolution', type=float, default=1.5, help='Leiden algorithm resolution (default=1.5).')
    parser_umap.add_argument('--plot_venn', action='store_false', help='Whether to plot Venn diagram (default=True).')

    # All command
    parser_all = subparsers.add_parser('all', help='Execute full pipeline from ANNOVAR files to UMAP plotting.')
    parser_all.add_argument("-i", "--input", type=str, help="Input ANNOVAR zip file or ANNOVAR directory directly")
    parser_all.add_argument("-p", "--patient", type=str, help="Input patient info file, including at least the info of Patient, stage, status, and type (e.g., histology type).")
    parser_all.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells for filtering (default=3).')
    parser_all.add_argument('--min_genes', type=int, default=30, help='Minimum number of genes for filtering (default=30).')
    parser_all.add_argument('--n_top_genes', type=int, default=3000, help='Number of top genes (default=3000).')
    parser_all.add_argument('--n_neighbors', type=int, default=60, help='Number of neighbors for UMAP (default=60).')
    parser_all.add_argument('--n_pcs', type=int, default=40, help='Number of principal components (default=40).')
    parser_all.add_argument('--leiden_resolution', type=float, default=1.5, help='Leiden algorithm resolution (default=1.5).')
    parser_all.add_argument('--plot_venn', action='store_false', help='Whether to plot Venn diagram (default=True).')

    args = parser.parse_args()

    if args.command == 'convert':
        convert_annovar(args.input, args.patient, args.output)
    elif args.command == 'umap':
        umap(args.input_umap, args.min_cells, args.min_genes, args.n_top_genes, args.n_neighbors, args.n_pcs, args.leiden_resolution, args.plot_venn)
    elif args.command == 'all':
        all_in_one(args.input, args.patient, args.min_cells, args.min_genes, args.n_top_genes, args.n_neighbors, args.n_pcs, args.leiden_resolution, args.plot_venn)
