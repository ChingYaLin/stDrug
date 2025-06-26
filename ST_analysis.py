"""
Author: Ching-Ya Lin
Last update: 20250626
Description: spatial transcriptome analysis pipeline.
Usage:
python3 ST_analysis.py -i /media/data1/PhD/ChingYaLin/VisiumHD_Colon/binned_outputs/square_008um \
                -o /media/data1/PhD/ChingYaLin/writing_reso1 \
                -p /media/data1/PhD/ChingYaLin/VisiumHD_Colon/8um_squares_annotation_pathologist.csv
"""

### Package import
# Single cell and Spatial transcript analysis
import scanpy as sc
import squidpy as sq
import scipy.sparse 
from scipy.sparse import csc_matrix # For sparse matrix

# Data processing
import numpy as np
import pandas as pd

# Data visualization and export to PDF
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

# Other
from pathlib import Path # for parquet file convert 
from sklearn.metrics import silhouette_score # for esrimate the cluster quality
import argparse, os
import time

def readFile(input_path, pathologist_path):
    """
    input_path = The path of spatial data
    """
    # Output: adata with raw count and metadata
    print('Reading files...')
    ## Transfer the parquet table into csv that squidpy needed. (Only Visium HD use parquet file)
    bin_dir = Path(input_path)
    df = pd.read_parquet(bin_dir / "spatial/tissue_positions.parquet")
    csv_path = bin_dir / "spatial/tissue_positions.csv"
    df.to_csv(csv_path, header=False, index=False) # Squidpy don't want header
    print("Wrote:", csv_path)

    ## Read visiumHD data
    adata = sq.read.visium(
        path = input_path,
        counts_file="filtered_feature_bc_matrix.h5",   # 指定 filtered
        load_images=True
    )

    # Add the pathologist annotation to adata table
    anno_df = pd.read_csv('/media/data1/PhD/ChingYaLin/VisiumHD_Colon/8um_squares_annotation_pathologist.csv', header=None, sep='\t', index_col=0)
    anno_df.columns = ['tissue_annotation_pathologist']  # 幫它加一個欄名
    adata.obs['tissue_annotation_pathologist'] = anno_df['tissue_annotation_pathologist']
    adata.obs['tissue_annotation_pathologist'] = adata.obs['tissue_annotation_pathologist'].astype('category')
    return adata

# process arguments
parser = argparse.ArgumentParser(description="Spatial transcriptomic data analysis")
parser.add_argument("-i", "--input", required=True, help="path to input Visium HD data")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-p", "--pathologist", required=True, help="csv file path to pathologist annotation.")


### Read input
global args
args = parser.parse_args()


# read input files
global pdf_save, results_file
pdf_save = PdfPages(f'{args.output}/results.pdf') 

### Main
adata = readFile(args.input, args.pathologist)







# close PDF file
pdf_save.close()
