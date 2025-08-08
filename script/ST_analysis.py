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
import pickle
import math
from os import listdir
from os.path import isfile, join
import collections
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter, CoxPHFitter


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



def clustering(adata, resolution, auto_reso):
    print("Clustering...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata)

    # Whether use auto-resolution selection
    if auto_reso:
        adata, res = autoResolution(adata)
    else:
        print(f"Clustering with resolution = {resolution}")
        sc.tl.leiden(adata,resolution, flavor="igraph") # if set random_state=0 can reproduce the clustering results.
        res = resolution
    
    print('Exporting UMAP...')
    fig, ax = plt.subplots(dpi=300, figsize=(12,12))
    sc.pl.umap(adata, color='leiden', ax=ax, legend_loc='on data',
           frameon=False, legend_fontsize='small', legend_fontoutline=2, legend_fontweight='normal',
           use_raw=False, show=False, title='leiden, resolution='+str(res))
    
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
    
def runGSEAPY(adata, output, group_by='leiden', cutoff=0.05, logfc_threshold=2):
    import gseapy as gp

    print('Running GSEAPY...')
    start = time.time()
    
    with open('./data/GO_Biological_Process_2021.pkl', 'rb') as handle:
        gene_sets = pickle.load(handle)

    df_list = []
    cluster_list = []
    celltypes = sorted(adata.obs[group_by].unique())

    for celltype in celltypes:
        degs = sc.get.rank_genes_groups_df(adata, group=celltype, key='rank_genes_groups', log2fc_min=logfc_threshold, 
                                    pval_cutoff=cutoff)['names'].squeeze()
        if isinstance(degs, str):
            degs = [degs.strip()]
        else:
            degs = degs.str.strip().tolist()
        
        if not degs:
            continue

        enr = gp.enrichr(gene_list=degs,
                gene_sets=gene_sets,
                no_plot=True
                )
        if (enr is not None) and hasattr(enr, 'res2d') and (enr.res2d.shape[0] > 0):
            df_list.append(enr.res2d)
            cluster_list.append(celltype)

    columns = ['Cluster', 'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Genes']

    df = pd.DataFrame(columns = columns)
    for cluster_ind, df_ in zip(cluster_list, df_list):
        df_ = df_[df_['Adjusted P-value'] <= cutoff]
        df_ = df_.assign(Cluster = cluster_ind)
        if (df_.shape[0] > 0):
            df = pd.concat([df, df_[columns]], sort=False)
            df_tmp = df_.loc[:, ['Term', 'Adjusted P-value']][:min(10, df_.shape[0])]
            df_tmp['Term'] = [x.split('(',1)[0] for x in df_tmp['Term']]
            df_tmp['-log_adj_p'] = - np.log10(df_tmp['Adjusted P-value'])
            df_tmp = df_tmp.sort_values(by='-log_adj_p', ascending=True)
            ax = df_tmp.plot.barh(y='-log_adj_p', x='Term', legend=False, grid=False, figsize=(12,4))
            ax.set_title('Cluster {}'.format(cluster_ind))
            ax.set_ylabel('')
            ax.set_xlabel('-log(Adjusted P-value)')
            fig = ax.figure
            fig.savefig(os.path.join(figure_save, f'GSEA_cluster{cluster_ind}.png'), dpi=300, bbox_inches='tight', transparent=False)
            plt.close()
        else:
            print('No pathway with an adjusted P-value less than the cutoff (={}) for cluster {}'.format(cutoff, cluster_ind))
    df.to_csv(os.path.join(output, 'GSEA_results.csv'))

    end = time.time()
    print('time: {}', end-start)

### Survival analysis
def getBulkProfile(bulkpath, gencode_table):
    bk_gep = pd.read_csv(bulkpath, sep=',', compression='gzip')
    bk_gep_name = bk_gep.merge(gencode_table, how='left', on='id')
    bk_gep_name.drop_duplicates(subset='name', inplace=True)
    bk_gep_name = bk_gep_name.set_index('name')
    bk_gep_name = bk_gep_name.drop(columns=['id'])
    return bk_gep_name

def getSpecCellDict(bk_gep_name, dict_deg):
    col_dict = {}
    for gene in bk_gep_name.index:
        median = bk_gep_name.loc[gene, :].median()
        col_dict[gene] = bk_gep_name.loc[gene, :] >= median
    score_table = pd.DataFrame(col_dict).T

    values_list = []
    for sample in score_table.columns:
        values = []
        for _, genes in dict_deg.items():
            values.append(score_table.loc[genes, sample].sum())
        values_list.append(values)

    spec_score_table = pd.DataFrame(values_list, index=score_table.columns, columns=dict_deg.keys())

    dict_low = spec_score_table.quantile(q=0.25)
    dict_high = spec_score_table.quantile(q=0.75)

    dict_celltype = collections.defaultdict(lambda: collections.defaultdict(dict))
    for c in spec_score_table.columns:
        dict_celltype[c]['high'] = spec_score_table.index[spec_score_table[c] >= dict_high[c]].tolist()
        dict_celltype[c]['low'] = spec_score_table.index[spec_score_table[c] <= dict_low[c]].tolist()

    return dict_celltype

def drawSurvivalPlot(dict_celltype, clinical_df, project_id):
    n_types = len(dict_celltype.keys())
    n_cols = 3
    n_rows = math.ceil(n_types / n_cols)
    dict_group = {'high': 1, 'low': 0}
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(2.5 * n_cols, 2.5 * n_rows))
    axs = axs.flatten() if n_rows * n_cols > 1 else [axs]  # 保證 axs 是 1D array
    i_ax = 0

    alpha_val = 0.05
    cph_selected_cols = ['project', 'cell type', 'coef', 'exp(coef)', 'se(coef)', 'coef lower 95%', 'coef upper 95%',
                         'exp(coef) lower 95%', 'exp(coef) upper 95%', 'z', 'p', '-log2(p)']
    df_hazard = pd.DataFrame(columns=cph_selected_cols)

    for c in sorted(dict_celltype.keys(), key=int):
        if i_ax >= len(axs):
            break
        ax = axs[i_ax]
        i_ax += 1

        dict_km = {'time': [], 'died': [], 'group': [], 'abundance': []}
        for g in ['high', 'low']:
            patients = [x for x in dict_celltype[c][g] if x in clinical_df.index]

            death_time = pd.to_numeric(
                clinical_df.loc[patients, '_event_time_'],
                errors='coerce'
            ).fillna(0).astype(int).tolist()

            died_series = clinical_df.loc[patients, 'days_to_death']
            died_flags = died_series.notna() & (died_series != '')  # 判斷是否為死亡事件
            dict_km['died'].extend(died_flags.astype(int).tolist())

            dict_km['time'].extend(death_time)
            dict_km['group'].extend([g] * len(death_time))
            dict_km['abundance'].extend([dict_group[g]] * len(death_time))

        df = pd.DataFrame(dict_km)

        # 確保欄位皆為數值型態
        df = df.astype({'time': int, 'died': int, 'abundance': int})

        T = df['time']
        E = df['died']
        dem = (df['group'] == 'high')

        # Kaplan-Meier plot
        kmf = KaplanMeierFitter(alpha=alpha_val)
        kmf.fit(T[dem], event_observed=E[dem], label='High')
        kmf.plot_survival_function(ax=ax, ci_show=False)
        kmf.fit(T[~dem], event_observed=E[~dem], label='Low')
        kmf.plot_survival_function(ax=ax, ci_show=False)

        out_lr = logrank_test(T[dem], T[~dem], E[dem], E[~dem], alpha=1 - alpha_val)
        p_value = out_lr.p_value

        # Cox regression
        df.drop(['group'], axis=1, inplace=True)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='time', event_col='died', show_progress=False)
        hr = cph.hazard_ratios_.values[0]
        dict_tmp = cph.summary.to_dict(orient='records')[0]
        dict_tmp['project'] = project_id
        dict_tmp['cell type'] = c
        df_new_row = pd.DataFrame([dict_tmp]).dropna(axis=1, how='all')
        df_hazard = pd.concat([df_hazard, df_new_row], ignore_index=True)  

        ax.set_title(f'C{c} (Log-rank P={p_value:.2f}, HR={hr:.2f})', fontsize=8, y=1.08)

    # 刪除空白 subplot
    for j in range(i_ax, len(axs)):
        fig.delaxes(axs[j])

    fig.suptitle('Survival Analysis: ' + project_id.rsplit('-', 1)[1], fontsize=12, y=1)
    fig.tight_layout()
    plot_name = 'survival_analysis_'+project_id.rsplit('-', 1)[1]+'.png'
    fig.savefig(os.path.join(figure_save, plot_name), dpi=300, bbox_inches='tight')
    plt.close('all')

    return df_hazard

def survivalAnalysis(adata, clinicalpath, gencode, tcga, no_treat, id):
    print('Survival Analysis...')

    # 判斷是否要排除已接受治療的樣本（由 args.not_treated 控制）
    if no_treat:
        treatment_status = 'no'
    else:
        treatment_status = ''

    # 取得所有可用的 TCGA project 檔名（不含 .csv.gz）
    project_ids = [f.rsplit('.csv',1)[0] for f in listdir(tcga) if isfile(join(tcga, f)) and f.startswith('TCGA')]

    # 若 user 有指定特定 project id，確認是否存在
    if id:
        if not id in project_ids:
            print(f"Cannot find the TCGA file for the specified project_id: {id}\nCandidates:{project_ids}")
            return adata

    # 讀取 GENCODE gene id → gene name 對照表
    gencode_table = pd.read_csv(gencode, names=['id','name'])

    # 讀取 TCGA 臨床生存資料（包含 death time、follow-up 等）
    clinical_df_all = pd.read_csv(clinicalpath, sep='\t', index_col=0)

    # 欲保留的 Cox regression hazard 結果欄位
    cph_selected_cols = ['project', 'cell type', 'coef', 'exp(coef)', 'se(coef)', 
                         'coef lower 95%', 'coef upper 95%', 'exp(coef) lower 95%', 
                         'exp(coef) upper 95%', 'z', 'p', '-log2(p)']
    df_hazard = pd.DataFrame(columns=cph_selected_cols)

    # 判斷是全部 project 還是 user 指定某一個
    run_project_ids = project_ids if not id else [id]

    for project_id in run_project_ids:
        print(project_id)

        # 讀取 TCGA bulk GEP，並將 Ensembl ID 映射成 gene name
        bk_gep_name = getBulkProfile(f'{tcga}/{project_id}.csv.gz', gencode_table)
        

        # 針對每個 cluster 抽出 top 20 DEGs（排除 MT 基因與不在 TCGA GEP 中的基因）
        dict_deg = {}
        all_degs = []
        for celltype in adata.obs['leiden'].unique().tolist():
            degs = [x for x in adata.uns['rank_genes_groups']['names'][celltype] 
                    if not (x.startswith('MT-') or not x in bk_gep_name.index)][:20]
            dict_deg[celltype] = degs
            all_degs.extend(degs)
        all_degs = list(set(all_degs))  # 移除重複基因

        # 將 TCGA bulk 表現矩陣只保留有用的 DEG
        bk_gep_name = bk_gep_name.loc[all_degs, :]

        # 根據是否排除 treated，過濾對應的臨床資料
        if treatment_status == 'no':
            clinical_df = clinical_df_all[(clinical_df_all['project_id'] == project_id) & 
                                          (clinical_df_all['isTreated'] == treatment_status)].copy()
        else:
            clinical_df = clinical_df_all[clinical_df_all['project_id'] == project_id].copy()

        # 過濾 TCGA bulk GEP 欄位，只保留 clinical info 中有的 sample
        bk_gep_name = bk_gep_name.loc[:, 
                        [x for x in clinical_df['case_submitter_id'].tolist() if x in bk_gep_name.columns]]

        # 若病人數不足，跳過該 project
        if bk_gep_name.shape[1] < 10:
            print(f'Skipped: the number of {project_id} patients (treated: {treatment_status}) is less than 10.')
            continue

        # 補上生存時間欄位：優先用 days_to_death，若沒有則用 follow-up
        clinical_df['_event_time_'] = clinical_df['days_to_death']
        for ind in clinical_df.index:
            if not clinical_df.loc[ind, '_event_time_'].isnumeric():
                clinical_df.loc[ind, '_event_time_'] = clinical_df.loc[ind, 'days_to_last_follow_up']

        # 根據 DEGs 與 TCGA bulk expression，定義 high/low 表現的病人群體
        dict_celltype = getSpecCellDict(bk_gep_name, dict_deg)

        # 計算 hazard ratio、畫生存圖
        df_new_hazard = pd.DataFrame()
        try:
            df_new_hazard = drawSurvivalPlot(dict_celltype, clinical_df, project_id)
            df_new_hazard.iloc[:,2:] = df_new_hazard.iloc[:,2:].round(2)

                # 畫表格圖，並儲存成 PDF
            fig, ax = plt.subplots()
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')
            table = ax.table(cellText=df_new_hazard.iloc[:,1:].values,
                                colLabels=df_new_hazard.columns[1:], loc='center')
            table.auto_set_font_size(False)
            table.auto_set_column_width(col=list(range(len(df_new_hazard.columns[1:]))))
            table.set_fontsize(8)
            for cell in table._cells:
                if cell[0] == 0:
                    table._cells[cell].set_fontsize(6)
                    table._cells[cell].set_color('lightblue')
                    table._cells[cell].set_height(.05)

            ax.set_title(project_id, y=1.08)
            fig.tight_layout()
            filename = 'survival_table_'+project_id
            fig.savefig(os.path.join(figure_save, filename), dpi=300, bbox_inches='tight')
            plt.close('all')
        except Exception as e:
            print(e)
            print(f'Survival analysis failed for {project_id}. The reason might be that the abundances, when conditioned on the presence or absence of death events, have very low variance and thus fail to converge.')
            continue
        else:
            df_hazard = pd.concat([df_hazard, df_new_hazard], ignore_index=True, join="inner")

    # 儲存全部 survival analysis 的結果
    adata.uns['survival_analysis'] = df_hazard
    df_hazard.to_csv(f'{output}/HR_{treatment_status}.csv')

    return adata




### process arguments
parser = argparse.ArgumentParser(description="Spatial transcriptomic data analysis")
parser.add_argument("-i", "--input", required=True, help="path to input Visium HD data")
parser.add_argument("-f", "--format", default='VisiumHD', help="input format, VisiumHD (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)")
parser.add_argument("-o", "--output", default='./', help="path to output directory, default='./'")
parser.add_argument("-p", "--pathologist", default='none', help="csv file path to pathologist annotation.")
parser.add_argument("-t", "--pathologist_target", default="Neoplasm", help = "the target of pathologist annotation that indicate tumor, default='Neoplasm'")
parser.add_argument("--impute", action="store_true", help="do imputation. default: no")
parser.add_argument("-r", "--resolution", type=float, default=0.8, help="resolution for clustering, default=0.8")
parser.add_argument("--auto-resolution", action="store_true", help="automatically determine resolution for clustering")
parser.add_argument("--species", default="human", help="sample species. Options: human (default) | mouse")
parser.add_argument("--cpus", default=1, type=int, help="number of CPU used for auto-resolution and annotation, default=1")
parser.add_argument("-c", "--clusters", default=None, help="perform single cell analysis only on specified clusters, e.g. '1,3,8,9'")
parser.add_argument("--cname", default='leiden', help="which variable should be used when selecting clusters; required when clusters are provided. Default: 'leiden'")
parser.add_argument("--gsea", action="store_true", help="perform gene set enrichment analysis (GSEA)")
parser.add_argument("--survival", action="store_true", help="perform survival analysis")
parser.add_argument("--tcga", default='/stDrug/data/TCGA/', help="path to TCGA data")
parser.add_argument("--id", default=None, help='Specify TCGA project id in the format "TCGA-xxxx", e.g., "TCGA-LIHC"')
parser.add_argument("--not_treated", action="store_true", help='only consider untreated samples from TCGA for survival analysis.')



### Read input
global args
args = parser.parse_args()

# check option
if args.survival:
    if not os.path.isdir(args.tcga):
            sys.exit("The path to TCGA files does not exist.")
    if not os.path.isfile(f'{args.tcga}/clinical.tsv'):
        sys.exit("The TCGA clinical file does not exist.")
    else:
        clinicalpath = f'{args.tcga}/clinical.tsv'
    if not os.path.isfile(f'{args.tcga}/gencode.v22.annotation.id.name.gtf'):
        sys.exit("The gencode file for id conversion does not exist.")
    else:
        gencode = f'{args.tcga}/gencode.v22.annotation.id.name.gtf'


# read input files
global figure_save , results_file
figure_save = os.path.join(args.output, 'figures')
os.makedirs(figure_save, exist_ok=True) 
results_file = os.path.join(args.output, 'scanpyobj.h5ad')

### Main
adata = readFile(args.input, args.pathologist)
adata = preprocessing(adata, impute = args.impute)
adata = clustering(adata, args.resolution, args.auto_resolution)

groups = sorted(adata.obs['leiden'].unique(), key=int)
adata = annotation(adata, groups, args.species, args.output, args.cpus)
adata = findDEG(adata, groups, args.output)
if args.pathologist != "none":
    pathologist_plot(adata, args.pathologist_target, args.output)
if args.gsea:
    runGSEAPY(adata, args.output)
if args.survival:
    adata = survivalAnalysis(adata, clinicalpath, gencode, args.tcga, args.not_treated, args.id)

# plot tissue and annotation spatial plot (save in pdf)
#sc.pl.spatial(adata,color="leiden",size=2,alpha=1,alpha_img=0.3)




if args.clusters:
    results_file = '{}.sub.h5ad'.format(results_file.rsplit('.',1)[0])
# save h5ad file
adata.write(results_file)
