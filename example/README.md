# Script usage (detail)
## ST_analysis.py

```
python3 ST_analysis.py -i /media/data1/PhD/ChingYaLin/VisiumHD_Colon/binned_outputs/square_008um \
                -o /media/data1/PhD/ChingYaLin/ST_analysis_test/write_first \
                -p /media/data1/PhD/ChingYaLin/VisiumHD_Colon/8um_squares_annotation_pathologist.csv \
                -t Neoplasm \
                --auto-resolution \
                --cpus 4 \
                --annotation \
                --gsea \
                --survival \
                --tcga /media/data1/PhD/ChingYaLin/TCGA_survival_data/TCGA \
                --id TCGA-COAD \
                --not_treated
```
## Sub-cluster analysis

```
python3 ST_analysis.py -i /media/data1/PhD/ChingYaLin/ST_analysis_test/write_first/scanpyobj.h5ad \
                -f h5ad \
                -o /media/data1/PhD/ChingYaLin/ST_analysis_test/write_sub \
                -r 0.6 \
                --cpus 4 \
                -c '0,2' \
                --annotation \
                --gsea \
                --survival \
                --tcga /media/data1/PhD/ChingYaLin/TCGA_survival_data/TCGA \
                --id TCGA-COAD \
                --not_treated
```