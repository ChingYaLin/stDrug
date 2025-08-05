"""
Author: Ching-Ya Lin
Last update: 20250626
Description: spatial transcriptome analysis pipeline.
Usage:
python3 ST_analysis.py -i /media/data1/PhD/ChingYaLin/VisiumHD_Colon/binned_outputs/square_008um \
                -o /media/data1/PhD/ChingYaLin/ST_analysis_test/write_first \
                -p /media/data1/PhD/ChingYaLin/VisiumHD_Colon/8um_squares_annotation_pathologist.csv \
                -t Neoplasm \
                -r 0.6 \
                --cpus 4
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
import argparse, os, sys
import time
import warnings


### Function part
def readFile(input_path, pathologist_path):
    """
    The function to read spatial transcriptome data to adata structure.
    # Input
    1. input_path:The path of spatial data.
    2. pathologist_path: The file path of pathologist annotation of each spot or pixel.
    # Output: original adata (un-preprocessing)
    """
    print('Reading files...')
    start = time.time()
    if args.format == 'VisiumHD':
        ## Transfer the parquet table into csv that squidpy needed. (Only Visium HD use parquet file)
        bin_dir = Path(input_path)
        df = pd.read_parquet(bin_dir / "spatial/tissue_positions.parquet")
        csv_path = bin_dir / "spatial/tissue_positions.csv"
        df.to_csv(csv_path, header=False, index=False) # Squidpy don't want header
        print("Create:", csv_path)

        ## Read visiumHD data
        adata = sq.read.visium(
            path = input_path,
            counts_file="filtered_feature_bc_matrix.h5",   # 指定 filtered
            load_images=True
        )
        adata.var_names_make_unique() # unique duplicate gene name

        # Add the pathologist annotation to adata table
        if pathologist_path != "none":
            anno_df = pd.read_csv(pathologist_path, header=None, sep='\t', index_col=0)
            anno_df.columns = ['tissue_annotation_pathologist']  # 幫它加一個欄名
            adata.obs['tissue_annotation_pathologist'] = anno_df['tissue_annotation_pathologist']
            adata.obs['tissue_annotation_pathologist'] = adata.obs['tissue_annotation_pathologist'].astype('category')
        print('Reading files...Done')
        print(f'n_obs x n_vars = {adata.n_obs} x {adata.n_vars}')

    # read main file
    if args.format == 'h5ad':
        adata = sc.read(args.input)
        if not adata.raw is None:
            adata = adata.raw.to_adata()
        if adata.X.max() < 50:
            print("Warning: The input h5ad doesn't contain raw count information.")
            adata.X = np.expm1(adata.X)
    
    # Filter clustering to post-analysis
    if args.clusters:
        clusters = [x.strip() for x in args.clusters.split(',')]
        if args.cname in adata.obs:
            adata = adata[adata.obs[args.cname].isin(clusters)]
        else:
            sys.exit(f"{args.cname} cannot be found in data.")
    end = time.time()
    print('time: {}', end-start)
    return adata

def preprocessing(adata, impute = False):
    """
    The function to preprocess spatial transcriptome data.
    """
    print('Preprocessing...')
    print(f'adata before filtering: n_obs x n_vars = {adata.n_obs} x {adata.n_vars}')
    # Filter the region that don't contain the tissue
    print("Filter the spot don't contain tissue")
    adata = adata[adata.obs["in_tissue"] == 1].copy()

    ## Remove mitochondrion
    print("Filter high mitochondrion spots")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # After calculate qc, we got those clumns can plot violin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax = sc.pl.violin(adata[adata.obs.sample(5000).index].copy(), # random choose 5000 spot to plot
                ["total_counts", "n_genes_by_counts", "pct_counts_mt"], jitter=0.4, palette="pastel", show=False)
        fig = ax[0].get_figure() if isinstance(ax, list) else ax.get_figure()
        fig.suptitle("QC Metrics Violin Plot")  # setting the title
        fig.savefig(os.path.join(figure_save, 'QC_metrix_violin_plots.png'), dpi=300, bbox_inches='tight')

    if not (adata.obs.pct_counts_mt == 0).all():
        adata = adata[adata.obs.pct_counts_mt < 30, :]

    ## Filter low express genes and spot
    print("Filter low-exp genes and spots")
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f'adata after filtering: n_obs x n_vars = {adata.n_obs} x {adata.n_vars}')

    
    
    ## Normalization and log-transformation
    print("Normalize and log-transform")
    
    sc.pp.normalize_total(adata, target_sum=1e4)   # CPM normalization
    adata.raw=adata.copy()
    sc.pp.log1p(adata)

    if impute:
        print('Doing imputation...')
        start = time.time()
        bk_adata_raw = adata.raw.to_adata()
        sc.external.pp.magic(adata)
        adata.raw = bk_adata_raw
        end = time.time()
        print('time: {}', end-start)

    print("Find highly variable genes")
    sc.pp.highly_variable_genes(adata)
    if 'features' in adata.raw.var.columns:
        adata.raw.var.index = adata.raw.var['features']
    adata = adata[:, adata.var.highly_variable].copy()
    # Remove comfunding variables
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # PCA dimontion reduction
    print("Scale and PCA")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    print('Preprocessing...Done')
    return adata

def autoResolution(adata):
    print("Automatically determine clustering resolution...")   
    start = time.time()    
    rep_n = 3
    sample_n = len(adata.obs)
    subsample_n = min(5000,int(sample_n))
    resolutions = np.linspace(0.1, 1.5, 15)
    silhouette_avg = {}
    np.random.seed(1)
    best_resolution = 0
    highest_sil = 0
    for r in resolutions:
        r = np.round(r, 1)
        print("Clustering test: resolution = ", r)
        sub_start = time.time()
        score_list = []
        for i in range(rep_n):
            random_indices = np.random.choice(adata.obs_names, size=subsample_n, replace=False)
            subadata = adata[random_indices].copy()
            sc.pp.neighbors(subadata, n_neighbors=15, n_pcs=20)
            sc.tl.leiden(subadata, resolution=r, flavor="igraph")
            X = subadata.obsm["X_pca"][:, :20]
            labels = subadata.obs['leiden']
            if len(set(labels)) <= 1:
                print(f"Warning: Only one cluster found at resolution={r}, replicate={i}. Skipping silhouette.")
                continue
            score = silhouette_score(X, labels)
            score_list.append(score)
        
        silhouette_avg[str(r)] = np.mean(score_list)
        if silhouette_avg[str(r)] > highest_sil:
            highest_sil = silhouette_avg[str(r)]
            best_resolution = r
        print("robustness score = ", silhouette_avg[str(r)])
        sub_end = time.time()
        print(f'time: {sub_end - sub_start}')
        print()

    print("resolution with highest score: ", best_resolution)
    res = best_resolution
    sc.tl.leiden(adata, resolution=best_resolution, flavor="igraph")
    # write silhouette record to uns and remove the clustering results except for the one with the best resolution
    adata.uns['sihouette score'] = silhouette_avg
    # draw lineplot
    df_sil = pd.DataFrame(list(silhouette_avg.values()), columns=['silhouette score'], index=[float(x) for x in silhouette_avg.keys()])
    df_sil.plot.line(style='.-', color='green', title='Auto Resolution', xticks=resolutions, xlabel='resolution', ylabel='silhouette score', legend=False)
    plt.xticks(resolutions)
    plt.tight_layout()
    plt.savefig(os.path.join(figure_save, 'Auto_resolution.png'), dpi=300, bbox_inches='tight')
    plt.close()
    end = time.time()
    print(f'time: {end-start}')
    return adata, res



def clustering(adata, resolution):
    print("Clustering...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata)

    print(f"Clustering with resolution = {resolution}")
    # Leiden cluster
    sc.tl.leiden(adata,resolution=0.8, flavor="igraph") # if set random_state=0 can reproduce the clustering results.

    '''
    if args.auto_resolution:
        adata, res = autoResolution(adata)
    else:
        print("Clustering with resolution = ", args.resolution)
        sc.tl.louvain(adata, resolution=args.resolution)
        res = args.resolution
    '''
    
    print('Exporting UMAP...')
    fig, ax = plt.subplots(dpi=300, figsize=(12,12))
    sc.pl.umap(adata, color='leiden', ax=ax, legend_loc='on data',
           frameon=False, legend_fontsize='small', legend_fontoutline=2, legend_fontweight='normal',
           use_raw=False, show=False, title='leiden, resolution='+str(resolution))
    
    fig.savefig(os.path.join(figure_save, '1.uamp_with_cluster.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Exporting spatial plot...")
    # Cluster result on tissue location
    fig, ax = plt.subplots(dpi=300)
    sq.pl.spatial_scatter(adata, color="leiden",fig=fig, ax=ax, img= False, size=2)
    fig.savefig(os.path.join(figure_save, '2.tissue_with_cluster.png'), dpi=300, bbox_inches='tight', transparent=False)
    plt.close()
    print("Clustering...Done")

    return adata

def annotation(adata, groups, species, output, cpus):
    if not adata.raw:
        print('Skip annotation since the process needs the expression data for all genes.')
    else:
        print('Cell type annotation...')
        
        # Export csv used by scMatch
        mat = np.zeros((len(adata.raw.var.index), len(groups)), dtype=float)
        for group in groups:
            mask = (adata.obs['leiden'] == group).values
            mat[:, int(group)] = adata.raw.X[mask, :].mean(axis=0).A1
        dat = pd.DataFrame(mat, index = adata.raw.var.index, columns = groups)
        dat.to_csv(os.path.join(output, 'cluster_mean_exp.csv'))
        
        os.system('python /media/data1/PhD/ChingYaLin/scMatch/scMatch.py --refDS /media/data1/PhD/ChingYaLin/scMatch/refDB/FANTOM5 \
                --dFormat csv --refType {} --testType {} --testDS {} --coreNum {}'.format(
                species, species, os.path.join( output, 'cluster_mean_exp.csv'), cpus))
        
        # Cell annotation result
        scMatch_cluster_df = pd.read_csv(os.path.join(output, 'cluster_mean_exp') + '/annotation_result_keep_all_genes/{}_Spearman_top_ann.csv'.format(species))
        scMatch_cluster_names = [group + " " + scMatch_cluster_df.loc[scMatch_cluster_df['cell']==int(group)]\
                                ['cell type'].tolist()[0] for group in groups]
        adata.obs['cell_type'] = adata.obs['leiden'].cat.rename_categories(scMatch_cluster_names)
        scMatch_candidate_df = pd.read_excel(os.path.join(output, 'cluster_mean_exp') + '/annotation_result_keep_all_genes/{}_Spearman.xlsx'.format(species), skiprows=4, header=None, index_col=0)
        for i in range(len(scMatch_candidate_df.columns)):
            if i%2 == 0:
                scMatch_candidate_df.iloc[:, i] = [x.split(',',1)[0].split(':',1)[0] for x in scMatch_candidate_df.iloc[:, i]]
        dict_candidates = {}
        for i in range(int(len(scMatch_candidate_df.columns)/2)):
            candidates = list(dict.fromkeys(scMatch_candidate_df.iloc[:5, i*2]))
            idx = 5
            while len(candidates) < 5:
                cell = scMatch_candidate_df.iloc[idx, i*2]
                if not cell in candidates:
                    candidates.append(cell)
                idx += 1
            dict_candidates[str(i)] = candidates
        df_candidate = pd.DataFrame(dict_candidates).T.reset_index().rename(columns={'index':'cluster'})
        del scMatch_candidate_df

        # UMAP with custom legend
        fig, ax = plt.subplots(dpi=300, figsize=(12, 12))
        sc.pl.umap(
            adata, color='leiden', ax=ax, legend_loc='on data', frameon=False, size=5,
            legend_fontsize='small', legend_fontoutline=2, legend_fontweight='normal')
        labels = sorted(adata.obs['cell_type'].unique(), key=lambda y: int(y.split(' ', 1)[0]))
        for label in labels:
            ind = int(label.split(' ', 1)[0])
            ax.scatter([], [], c=adata.uns['leiden_colors'][ind], label=label)
        ax.legend(
            frameon=False, loc='center left', bbox_to_anchor=(1, 0.5),
            ncol=1, fontsize=16
        )
        fig.subplots_adjust(right=0.85)  # 給 legend 多一點空間
        fig.savefig(os.path.join(figure_save, '3.umap_with_cluster(custom).png'), dpi=300, bbox_inches='tight')
        plt.close()

        # Table: top 5 annotation
        fig, ax = plt.subplots(dpi=300)
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        table = ax.table(cellText=df_candidate.values, colLabels=df_candidate.columns, loc='center')
        table.auto_set_column_width(col=list(range(len(df_candidate.columns))))
        table.scale(1, 1.8)

        for cell in table._cells:
            if cell[0] == 0:
                table._cells[cell].set_color('lightblue')
                table._cells[cell].set_height(0.05)

        ax.set_title('Top 5 annotation')
        fig.subplots_adjust(left=0.1, right=0.9)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(os.path.join(figure_save, '4.top_5_annotation.png'), dpi=300, bbox_inches='tight')
        plt.close()

        # Spatial plot with annotation
        fig, ax = plt.subplots(dpi=300, figsize=(5, 5))
        sq.pl.spatial_scatter(adata, ax=ax, color="cell_type",img= False, size=2)
        fig.savefig(os.path.join(figure_save, '5.tissue_with_annotation.png'), dpi=300, bbox_inches='tight', transparent=False)
        plt.close()

    return adata

def findDEG(adata, groups, output):
    # Finding differentially expressed genes
    print('Finding Differentially Expressed Genes...')
    method = "t-test"
    if adata.raw:
        if adata.raw.X.max() >= 50:
            # haven't been log-transform
            adata_raw = adata.raw.to_adata()
            adata_raw_bk = adata_raw.copy()
            adata_raw.X = np.log1p(adata_raw.X)
            adata.raw = adata_raw
            sc.tl.rank_genes_groups(adata, 'leiden', method=method, pts=True, use_raw=True)
            adata.raw = adata_raw_bk
    else:
        sc.tl.rank_genes_groups(adata, 'leiden', method=method, pts=True)


    # cluster DEGs
    result = adata.uns['rank_genes_groups']
    dat = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    dat.to_csv(os.path.join(output, 'cluster_DEGs.csv'))
    print('Finding Differentially Expressed Genes...Done')

    return adata

def pathologist_plot(adata, ratio_taget, output):
    fig, axs = plt.subplots(1, 2, figsize=(10, 5), dpi=300)

    # plot umap
    sc.pl.umap(adata, color='tissue_annotation_pathologist', show=False, ax=axs[0], legend_loc='none')
    axs[0].set_title("UMAP")

    # plot Spatial scatter
    sq.pl.spatial_scatter(adata, ax=axs[1], color="tissue_annotation_pathologist", img=False, size=2)
    axs[1].set_title("Spatial")

    plt.tight_layout()
    fig.savefig(os.path.join(figure_save, 'combined_umap_spatial.png'), dpi=300, bbox_inches='tight', transparent=False)
    plt.close()

    # Calculate the cluster's tumor ratio with pathologist annotation
    leiden = adata.obs['leiden']
    annotation = adata.obs['tissue_annotation_pathologist']
    df = pd.DataFrame({'leiden': leiden, 'annotation': annotation})
    total_per_cluster = df['leiden'].value_counts().sort_index()
    neoplasm_per_cluster = df[df['annotation'] == ratio_taget]['leiden'].value_counts().sort_index()
    ratio_table = pd.DataFrame({'total_spots': total_per_cluster,'pathologist_cancer_spots': neoplasm_per_cluster})
    ratio_table = ratio_table.fillna(0).astype(int)
    ratio_table['cancer_spot_ratio'] = ratio_table['pathologist_cancer_spots'] / ratio_table['total_spots']
    ratio_table.to_csv(os.path.join(output, 'pathologist_annotation_inCluster_ratio.csv'))

    fig, ax = plt.subplots(dpi=300, figsize=(5, 5))
    ratio_table['cancer_spot_ratio'].plot(kind='bar')
    plt.axhline(y=0.8, color='red', linestyle='--', linewidth=1)  # 紅色虛線
    plt.ylabel('Proportion of Cancer spots')
    plt.xlabel('Leiden Cluster')
    plt.title('Proportion of cancer spots in Each Leiden Cluster')
    plt.tight_layout()
    fig.savefig(os.path.join(figure_save, 'pathologist_annotation_ratio.png'), dpi=300, bbox_inches='tight', transparent=False)
    plt.close()
    




### process arguments
parser = argparse.ArgumentParser(description="Spatial transcriptomic data analysis")
parser.add_argument("-i", "--input", required=True, help="path to input Visium HD data")
parser.add_argument("-f", "--format", default='VisiumHD', help="input format, VisiumHD (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-p", "--pathologist", default='none', help="csv file path to pathologist annotation.")
parser.add_argument("-t", "--pathologist_target", default="Neoplasm", help = "the target of pathologist annotation that indicate tumor, default='Neoplasm'")
parser.add_argument("--impute", action="store_true", help="do imputation. default: no")
parser.add_argument("-r", "--resolution", type=float, default=0.8, help="resolution for clustering, default=0.8")
parser.add_argument("--species", default="human", help="sample species. Options: human (default) | mouse")
parser.add_argument("--cpus", default=1, type=int, help="number of CPU used for auto-resolution and annotation, default=1")
parser.add_argument("-c", "--clusters", default=None, help="perform single cell analysis only on specified clusters, e.g. '1,3,8,9'")
parser.add_argument("--cname", default='leiden', help="which variable should be used when selecting clusters; required when clusters are provided. Default: 'leiden'")



### Read input
global args
args = parser.parse_args()


# read input files
global figure_save , results_file
figure_save = os.path.join(args.output, 'figures')
os.makedirs(figure_save, exist_ok=True) 
results_file = os.path.join(args.output, 'scanpyobj.h5ad')

### Main
adata = readFile(args.input, args.pathologist)
adata = preprocessing(adata, impute = args.impute)
adata = clustering(adata, resolution=args.resolution)

groups = sorted(adata.obs['leiden'].unique(), key=int)
adata = annotation(adata, groups, args.species, args.output, args.cpus)
adata = findDEG(adata, groups, args.output)
if args.pathologist != "none":
    pathologist_plot(adata, args.pathologist_target, args.output)

# plot tissue and annotation spatial plot (save in pdf)
#sc.pl.spatial(adata,color="leiden",size=2,alpha=1,alpha_img=0.3)




if args.clusters:
    results_file = '{}.sub.h5ad'.format(results_file.rsplit('.',1)[0])
# save h5ad file
adata.write(results_file)
